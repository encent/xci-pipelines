"""Call H3K27me3 valleys with a 2-state Gaussian HMM.

This is the pipeline's primary exact-equality test target. Every constant below
is carried verbatim from the verified original
(``TEST_03_00_H3K27me3_boundary_calling.ipynb``, frozen ground-truth-config
snapshot in ``XCI_valleys_check/``; algorithm transcribed line by line in
``the archaeology notes`` section 2.1).

The algorithm, in order
-----------------------
1. ``train = bioframe.subtract(chromsizes[chrom], blacklist[chrom])`` keeping
   intervals of ``size >= region_size_threshold`` (500 kb).
2. ``bioframe.read_bigwig`` -> ``bioframe.binnify(resolution)`` ->
   ``bioframe.overlap`` -> ``fillna(0.0)``.
3. One ``GaussianHMM`` per bigWig, fitted on the *concatenated training chunks*
   with ``lengths=`` -- never on the whole chromosome, and never pooled across
   samples.
4. ``valley_state = int(np.argmin(hmm.means_))``.
5. Viterbi is **restarted inside each training interval**; bins outside every
   training interval keep state ``-1`` and can never be valleys.
6. ``bioframe.merge`` on the valley bins -- adjacency only (``min_dist=0``).

There is NO smoothing, NO signal transform, NO threshold, NO minimum valley
length. A single 5 kb bin is a legal valley. Do not "improve" any of this.

Two things that look like bugs and are not
------------------------------------------
*The leading ``chrX:0-3,285,000`` call.* mm10 begins with an assembly gap. It
has no reads, ``fillna(0.0)`` makes it indistinguishable from a deep valley, and
it sits just OUTSIDE the first ENCODE blacklist interval (which starts at
3,286,700), so masking does not remove it. It is a genuine artefact. We
reproduce it anyway (ruling O-7): it survives the paper's own 3-gene filter and
is counted in the published Fig 2g bar totals, so dropping it would put every
clone one valley below the paper. ``valleys.drop_leading_gap: true`` removes it
for future data.

*Decision D-08 -- 18 files where ground truth has 17.* The original wrapped
``bioframe.merge`` in a bare ``try/except`` that printed "No valleys found!" and
``continue``d, which is why ``H3K27me3_F3_CTCF-NodTAG_Gall_valleys.{bed,bw}`` is
absent from the ground truth (35 files instead of 36; 34 locus figures instead
of 36). We do NOT swallow it: an empty valley set produces a valid empty BED, an
empty bigWig and a loud warning in the log. A correct run therefore emits **18**
degron valley BEDs where the ground truth has 17, and the visualisation stage
emits **36** locus figures where the ground truth has 34. That is expected, not
a failure.

Acceptance numbers (raw Xi valleys, ``the ground-truth inventory`` 4.3)
-----------------------------------------------------------------------------
WT      B1 261, C5 304, CL30 324, E6 377, JTG 302
degron  B1621-NodTAG 284, B1621-dTAG 281, E6A7-NodTAG 321, E6A7-dTAG 281,
        F3-NodTAG 351, F3-dTAG 205

Emits
-----
``{track}_valleys.bed``    BED3, tab-separated, no header, ONE chromosome.
                           This is the byte-comparison target -- no name column,
                           no score, no strand. (``valley_boundaries`` writes
                           BED4; the two formats differ deliberately.)
``{track}_valleys.bw``     the same intervals with constant value 1.0.
``{track}_states.tsv.gz``  per-bin chrom/start/end/value/state, so a moved
                           valley edge is localisable to a bin instead of
                           diffed as two BED files.
"""

import os
import os as _os
import shutil
import sys

import bioframe
import threadpoolctl
import numpy as np
import pandas as pd
import pyBigWig
from hmmlearn.hmm import GaussianHMM

chrom = snakemake.wildcards.chrom
track = snakemake.wildcards.track
resolution = int(snakemake.params.resolution)
hmm_cfg = dict(snakemake.params.hmm)
merge_gap = int(snakemake.params.merge_gap)
drop_leading_gap = bool(snakemake.params.drop_leading_gap)
blas_threads = int(snakemake.params.blas_threads)
size_threshold = int(snakemake.params.region_size_threshold)

log_path = snakemake.log[0]
os.makedirs(os.path.dirname(log_path), exist_ok=True)
for out in (snakemake.output.bed, snakemake.output.bw, snakemake.output.states):
    os.makedirs(os.path.dirname(out), exist_ok=True)


with open(log_path, "w") as log:

    def say(msg: str) -> None:
        log.write(msg + "\n")
        log.flush()
        print(msg, file=sys.stderr)

    # -- training regions: chromosome minus blacklist, >= 500 kb ---------
    chromsizes_bed = bioframe.read_chromsizes(snakemake.input.chrom_sizes, as_bed=True)
    chromsizes_bed = chromsizes_bed[chromsizes_bed["chrom"] == chrom]
    if chromsizes_bed.empty:
        raise SystemExit(
            f"{chrom} is not in {snakemake.input.chrom_sizes}.\n"
            "Check config.chromosomes.downstream / valley_control."
        )

    blacklist = pd.read_csv(
        snakemake.input.blacklist,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "name"],
    )
    blacklist = blacklist[blacklist["chrom"] == chrom]

    train = bioframe.subtract(chromsizes_bed, blacklist)
    train["size"] = train["end"] - train["start"]
    train = train[train["size"] >= size_threshold].copy()
    train.reset_index(drop=True, inplace=True)
    say(
        f"{track} {chrom}: {len(blacklist)} blacklist intervals -> "
        f"{len(train)} training intervals >= {size_threshold} bp "
        f"({int(train['size'].sum())} bp total)"
    )
    if train.empty:
        raise SystemExit(
            f"No training interval on {chrom} survives the "
            f"{size_threshold} bp size threshold. Either the blacklist is wrong "
            "or valleys.region_size_threshold is too large."
        )

    # -- signal: bigWig -> 5 kb bins, missing bins are 0.0 ---------------
    chromsizes = bioframe.read_chromsizes(
        snakemake.input.chrom_sizes, chrom_patterns=[f"^{chrom}$"]
    )
    chrom_len = int(chromsizes[chrom])

    signal = bioframe.read_bigwig(snakemake.input.bw, chrom, start=0, end=chrom_len)
    bins = bioframe.binnify(chromsizes, resolution)
    binned = bioframe.overlap(bins, signal)[["chrom", "start", "end", "value_"]]
    binned.rename(columns={"value_": "value"}, inplace=True)
    n_missing = int(binned["value"].isna().sum())
    binned["value"] = binned["value"].fillna(0.0)
    binned.reset_index(drop=True, inplace=True)

    if len(binned) != len(bins):
        say(
            f"WARNING: {len(binned)} rows from {len(bins)} bins -- some bin "
            "overlaps more than one bigWig interval. The original had the same "
            "behaviour; downstream indexing follows the overlap frame."
        )
    say(
        f"signal: {len(binned)} bins @ {resolution} bp, {n_missing} with no "
        f"bigWig coverage -> 0.0"
    )

    # -- fit: one HMM per bigWig, on the concatenated training chunks ----
    sequences, lengths = [], []
    for _, chunk in train.iterrows():
        inside = (binned["start"] >= chunk["start"]) & (binned["end"] <= chunk["end"])
        chunk_signal = binned.loc[inside, "value"].values.reshape(-1, 1)
        if len(chunk_signal) == 0:
            say(
                f"WARNING: training interval {chrom}:{chunk['start']}-{chunk['end']} "
                "contains no whole bin; skipped."
            )
            continue
        sequences.append(chunk_signal)
        lengths.append(len(chunk_signal))

    observations = np.vstack(sequences)
    model = GaussianHMM(
        n_components=int(hmm_cfg["n_components"]),
        covariance_type=str(hmm_cfg["covariance_type"]),
        n_iter=int(hmm_cfg["n_iter"]),
        random_state=int(hmm_cfg["random_state"]),
        min_covar=float(hmm_cfg["min_covar"]),
        tol=float(hmm_cfg["tol"]),
    )
    # ---- V-C: the fit MUST run under a pinned BLAS thread count ---------
    # hmmlearn 0.3.3 initialises GaussianHMM emission means with
    # sklearn.cluster.KMeans (hmm.py:311-315). KMeans' chunked parallel
    # reduction sums partial results in a THREAD-COUNT-DEPENDENT order, so a
    # different ambient thread count produces a different initialisation, EM
    # settles in a different local optimum, and a 2-state HMM turns that into
    # whole valleys appearing or vanishing. Measured on H3K27me3_E6_WT_Xi:
    #
    #   threads=1  378 valleys  md5 954d3093...
    #   threads=2  377 valleys  md5 5336c35f...  == ground truth, byte-identical
    #   threads=4  377 valleys  md5 5336c35f...  == ground truth, byte-identical
    #   threads=8  320 valleys  md5 b100d151...  means in the OPPOSITE order,
    #                                            so argmin(means) itself flips
    #
    # This is NOT the n_iter=15 cap: all 51 fits report converged=True.
    #
    # `threads:` on the Snakemake rule does NOT help -- it only tells the
    # scheduler how to allocate slots and never touches OpenBLAS. Nor do env
    # vars set at import time, which are global and cannot be scoped. So the
    # limit is applied here, around the fit, with threadpoolctl.
    #
    # Left unpinned this is a SHIPPING defect, not a ground-truth question: this
    # box has 48 cores and OpenBLAS takes all of them, so a reviewer on a laptop
    # would get different valley calls from us and neither would be
    # reproducible. That is fatal for a pipeline whose purpose is reproduction.
    with threadpoolctl.threadpool_limits(limits=blas_threads):
        observed = {
            f"{d['user_api']}/{d['internal_api']}": d["num_threads"]
            for d in threadpoolctl.threadpool_info()
        }
        libs = {
            d["internal_api"]: _os.path.basename(d["filepath"] or "?")
            for d in threadpoolctl.threadpool_info()
        }
        say(
            f"BLAS: requested={blas_threads} observed={observed} libs={libs}"
        )
        for api, n in observed.items():
            if n != blas_threads:
                say(
                    f"WARNING: {api} reports {n} threads, not the requested "
                    f"{blas_threads}. The valley calls will NOT be comparable "
                    "to a run at the requested count. Investigate before "
                    "trusting this output."
                )
        model.fit(observations, lengths=lengths)


    means = model.means_.flatten()
    valley_state = int(np.argmin(means))
    hill_state = int(np.argmax(means))
    if valley_state == hill_state:
        # Degenerate fit: the two components collapsed onto the same mean, so
        # argmin == argmax and the original's relabelling dict would silently
        # map one state to None. Say so, loudly, and keep the run deterministic.
        hill_state = 1 - valley_state
        say(
            f"WARNING: the two HMM components have identical means "
            f"({means.tolist()}). Valley/hill assignment is arbitrary; treat "
            "this sample's valleys as meaningless."
        )
    say(
        f"HMM: converged={model.monitor_.converged} "
        f"iterations={len(model.monitor_.history)}/{int(hmm_cfg['n_iter'])} "
        f"means={np.round(means, 6).tolist()} valley_state={valley_state}"
    )
    if not model.monitor_.converged:
        say(
            "NOTE: the EM fit hit the n_iter cap WITHOUT converging. That is the "
            "original's behaviour (n_iter=15, tol=0.01) and is reproduced on "
            "purpose -- but an unconverged fit can settle on a different local "
            "optimum under a different BLAS build or thread count, which is the "
            "leading explanation for the published valley BEDs not being "
            "regenerable. BLAS threading is pinned to 1 at the top of this "
            "script to make OUR output machine-independent. If this track's "
            "count disagrees with the published one, check whether the tracks "
            "that DO reproduce are the converged ones."
        )

    # -- predict: Viterbi restarted inside each training interval --------
    # -1 = never evaluated (blacklist, or a non-blacklisted stretch < 500 kb,
    # or a bin straddling a training-interval edge). Those are never valleys.
    predicted = np.full(binned.shape[0], -1)
    relabel = {valley_state: 0, hill_state: 1}
    for _, chunk in train.iterrows():
        inside = (binned["start"] >= chunk["start"]) & (binned["end"] <= chunk["end"])
        chunk_signal = binned.loc[inside, "value"].values.reshape(-1, 1)
        if len(chunk_signal) == 0:
            continue
        raw_states = model.predict(chunk_signal)
        predicted[inside.values] = np.vectorize(relabel.get)(raw_states)

    valley_mask = predicted == 0
    n_uncalled = int((predicted == -1).sum())
    say(
        f"states: {int(valley_mask.sum())} valley bins, "
        f"{int((predicted == 1).sum())} hill bins, {n_uncalled} not evaluated"
    )

    # -- merge adjacent valley bins --------------------------------------
    if valley_mask.any():
        merged = bioframe.merge(binned.loc[valley_mask], min_dist=merge_gap)
        merged = merged[["chrom", "start", "end"]].copy()
        merged.reset_index(drop=True, inplace=True)
    else:
        merged = pd.DataFrame(columns=["chrom", "start", "end"])
        say(
            "WARNING: no valleys found. The original printed 'No valleys "
            "found!' and wrote nothing; we write a valid empty BED and an "
            "empty bigWig (decision D-08). Expect one MORE file here than in "
            "the ground truth."
        )

    if drop_leading_gap and len(merged) and int(merged.iloc[0]["start"]) == 0:
        dropped_end = int(merged.iloc[0]["end"])
        merged = merged.iloc[1:].reset_index(drop=True)
        say(
            f"drop_leading_gap: removed {chrom}:0-{dropped_end} (mm10 assembly "
            "gap). NOTE: this puts the valley count one BELOW the published "
            "Fig 2g bars, which include it."
        )
    elif len(merged) and int(merged.iloc[0]["start"]) == 0:
        say(
            f"NOTE: leading call {chrom}:0-{int(merged.iloc[0]['end'])} is the "
            "mm10 assembly gap. Kept on purpose -- it survives the 3-gene "
            "filter and is counted in the published Fig 2g bars (ruling O-7)."
        )

    say(f"called {len(merged)} valleys")

    # -- outputs ---------------------------------------------------------
    # BED3. No name column: this file is byte-compared against ground truth.
    merged[["chrom", "start", "end"]].to_csv(
        snakemake.output.bed, sep="\t", index=False, header=False
    )

    if len(merged):
        binary = shutil.which("bedGraphToBigWig")
        if binary is None:
            raise SystemExit(
                "bedGraphToBigWig is not on PATH.\n"
                "It ships with the `xci-pipeline` env (ucsc-bedgraphtobigwig); "
                "activate the env, or install it with\n"
                "    conda install -c bioconda ucsc-bedgraphtobigwig"
            )
        as_bigwig = merged.copy()
        as_bigwig["value"] = 1.0
        bioframe.to_bigwig(
            as_bigwig, chromsizes, snakemake.output.bw, path_to_binary=binary
        )
    else:
        # bedGraphToBigWig refuses an empty bedGraph; write the header only.
        empty = pyBigWig.open(snakemake.output.bw, "w")
        empty.addHeader([(chrom, chrom_len)])
        empty.close()

    states_out = binned[["chrom", "start", "end", "value"]].copy()
    states_out["state"] = predicted
    states_out.to_csv(snakemake.output.states, sep="\t", index=False, compression="gzip")

    say(f"wrote {snakemake.output.bed}")
