"""Per-loop observed/expected contact frequency, and the paired Wilcoxon test.

Reimplements the data half of ``01_03_loops_chromosight_stats.py``. The
original computed exactly this and then threw it away -- the numbers only ever
appeared inside violin-plot titles, and its own `loop_stats*.tsv` writes are
commented out. Here the per-loop table is a first-class output.

The arithmetic, verbatim:

1. ``expected_cis(clr, nproc, ignore_diags=0, intra_only=True,
   view_df=<balanced span>, chunksize=1e6)`` -- the view is the first to last
   bin carrying a non-null balancing weight, NOT the ROI.
2. divide each pixel's ``balanced`` by ``balanced.avg`` at its genomic
   distance; +/-inf becomes NaN.
3. per loop, ``np.nanmean`` over the 3x3 block of O/E pixels centred on
   (bin1, bin2). Pixels absent from the sparse matrix are simply not in the
   mean (the original's per-neighbour try/except).
4. Wilcoxon signed-rank between the two groups of an experiment, paired by
   loop index: two-sided plus both one-sided alternatives.

The ROI selects loops with the original's asymmetric predicate --
``start1 > start and end1 > start and start2 < end and end2 < end``.
"""
import os
import sys

import cooler
import numpy as np
import pandas as pd
from cooltools import expected_cis
from scipy.stats import wilcoxon

CENTER_OFFSETS = [-1, 0, 1]


def balanced_span(clr):
    """The view the original used: first to last bin with a balancing weight."""
    bins = clr.bins()[:]
    kept = bins[bins["weight"].notnull()]
    return pd.DataFrame([{
        "chrom": kept.iloc[0]["chrom"],
        "start": int(kept.iloc[0]["start"]),
        "end": int(kept.iloc[-1]["end"]),
        "name": kept.iloc[0]["chrom"],
    }])


def oe_matrix(clr, nproc, ignore_diags, intra_only):
    """Sparse pixel table with `balanced` replaced by observed/expected."""
    view = balanced_span(clr)
    cvd = expected_cis(clr, nproc=nproc, ignore_diags=ignore_diags,
                       intra_only=intra_only, view_df=view, chunksize=1000000)
    pix = clr.matrix(as_pixels=True)[:].copy()
    pix["dist"] = (pix["bin2_id"] - pix["bin1_id"]).abs()
    merged = pix.merge(cvd[["dist", "balanced.avg"]], on="dist", how="left")
    merged["balanced"] = merged["balanced"] / merged["balanced.avg"]
    merged["balanced"] = merged["balanced"].replace([np.inf, -np.inf], np.nan)
    return {
        (int(b1), int(b2)): v
        for b1, b2, v in zip(merged["bin1_id"], merged["bin2_id"], merged["balanced"])
    }


def centre_mean(oe, bin1, bin2, center):
    """`np.nanmean` over the center x center block, missing pixels excluded."""
    half = center // 2
    offsets = range(-half, half + 1)
    vals = [oe[(bin1 + i, bin2 + j)] for i in offsets for j in offsets
            if (bin1 + i, bin2 + j) in oe]
    return np.nanmean(vals) if vals else np.nan


def select_loops(loops, start, end):
    """The original's asymmetric ROI predicate, reproduced exactly."""
    m = ((loops["start1"] > start) & (loops["end1"] > start)
         & (loops["start2"] < end) & (loops["end2"] < end))
    return loops[m].copy().reset_index(drop=True)


def main():
    smk = snakemake  # noqa: F821
    comparison = smk.wildcards.comparison
    roi = smk.wildcards.roi
    suffix = smk.wildcards.suffix
    called_on = suffix.lstrip("_") or None
    experiments = dict(smk.params.experiments)
    center = int(smk.params.center)
    exp_kw = dict(smk.params.expected)

    for path in (smk.output.stats, smk.log[0]):
        os.makedirs(os.path.dirname(path), exist_ok=True)

    rois = pd.read_csv(smk.input.roi_table, sep="\t")
    loops_of = {os.path.basename(p).rsplit(".", 1)[0]: p for p in smk.input.loops}
    mcool_of = {os.path.basename(p)[: -len(".mcool")]: p for p in smk.input.mcools}

    rows, per_loop, tests = [], [], []
    with open(smk.log[0], "w") as fh:
        for exp, members in experiments.items():
            locus = exp.split("_")[0]
            r = rois[(rois["locus"] == locus) & (rois["roi"] == roi)].iloc[0]
            start, end = int(r["start"]), int(r["end"])

            group_freqs = {}
            for cool_name in members:
                # Which loop set: the allele under test for Xa_vs_Xi, the
                # NodTAG partner's allele for dTAG_vs_NodTAG.
                allele = called_on or ("Xa" if "_Xa" in cool_name else "Xi")
                loop_set = f"{locus}_{allele}"
                if loop_set not in loops_of:
                    fh.write(f"SKIP {exp}/{cool_name}: no loop set {loop_set}\n")
                    continue
                loops = select_loops(
                    pd.read_csv(loops_of[loop_set], sep="\t"), start, end)
                if loops.empty:
                    fh.write(f"EMPTY {exp}/{cool_name}: 0 loops in {locus}/{roi}\n")
                    group_freqs[cool_name] = pd.Series(dtype=float)
                    rows.append({"file": cool_name, "experiment": exp,
                                 "loop_set": loop_set, "n_loops": 0,
                                 "median": np.nan, "mean": np.nan, "std": np.nan})
                    continue

                clr = cooler.Cooler(
                    f"{mcool_of[cool_name]}::resolutions/{smk.params.resolution}")
                oe = oe_matrix(clr, smk.threads, exp_kw["ignore_diags"],
                               exp_kw["intra_only"])
                freqs = np.array([
                    centre_mean(oe, int(b1), int(b2), center)
                    for b1, b2 in zip(loops["bin1"], loops["bin2"])
                ])
                group_freqs[cool_name] = pd.Series(freqs, index=loops.index)
                rows.append({
                    "file": cool_name, "experiment": exp, "loop_set": loop_set,
                    "n_loops": int(len(freqs)),
                    "median": np.nanmedian(freqs), "mean": np.nanmean(freqs),
                    "std": np.nanstd(freqs),
                })
                for i, v in enumerate(freqs):
                    per_loop.append({
                        "experiment": exp, "group": cool_name, "loop_set": loop_set,
                        "loop_index": i, "loop_strength": v,
                        "present": cool_name.endswith(f"_{allele}"),
                    })
                fh.write(f"OK   {exp}/{cool_name}: {len(freqs)} loops, "
                         f"mean {np.nanmean(freqs):.4f}\n")

            names = [n for n in members if n in group_freqs]
            if len(names) == 2:
                a, b = group_freqs[names[0]], group_freqs[names[1]]
                paired = pd.concat([a, b], axis=1).dropna()
                if len(paired) >= 1:
                    x, y = paired.iloc[:, 0], paired.iloc[:, 1]
                    tests.append({
                        "experiment": exp, "group1": names[0], "group2": names[1],
                        "n": int(len(paired)),
                        "wilcoxon_two_sided": wilcoxon(x, y, alternative="two-sided").pvalue,
                        "wilcoxon_greater": wilcoxon(x, y, alternative="greater").pvalue,
                        "wilcoxon_less": wilcoxon(x, y, alternative="less").pvalue,
                    })

        fh.write(f"\n{comparison}/{roi}{suffix}: {len(rows)} samples, "
                 f"{len(tests)} paired tests\n")

    pd.DataFrame(rows, columns=["file", "experiment", "loop_set", "n_loops",
                                "median", "mean", "std"]).to_csv(
        smk.output.stats, sep="\t", index=False)
    pd.DataFrame(per_loop, columns=["experiment", "group", "loop_set",
                                    "loop_index", "loop_strength", "present"]).to_csv(
        smk.output.per_loop, sep="\t", index=False)
    pd.DataFrame(tests, columns=["experiment", "group1", "group2", "n",
                                 "wilcoxon_two_sided", "wilcoxon_greater",
                                 "wilcoxon_less"]).to_csv(
        smk.output.wilcoxon, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
