"""coolpuppy pile-up of one cooler over one refined loop set.

Reimplements the pile-up half of ``01_04_loops_pileups_and_stackups.py``,
emitting the aggregate matrix as `.npz` and its scores as `.tsv` instead of
drawing a PNG/SVG. `panel_pileup` and `panel_loops_comps_grid` do the drawing,
with the original's arguments:

    plotpup.plot(pup, score=True, cmap='coolwarm', scale='log', sym=True,
                 vmax=2, vmin=None, height=4, plot_ticks=True, center=3)

`view_df` is the span of BALANCED bins, not the ROI -- the ROI only filters
which loops enter the pile-up. This matters: passing the ROI as the view would
silently change the expected-value normalisation.

An empty loop set is not an error. The original hit `continue` and produced no
file; a Snakemake output has to exist, so we write an all-NaN matrix and a
score row of NaN with `n_loops = 0`, and say so in the log.
"""
import os
import sys

import cooler
import numpy as np
import pandas as pd
from cooltools import expected_cis
from coolpuppy import coolpup


def balanced_span(clr):
    bins = clr.bins()[:]
    kept = bins[bins["weight"].notnull()]
    return pd.DataFrame([{
        "chrom": kept.iloc[0]["chrom"],
        "start": int(kept.iloc[0]["start"]),
        "end": int(kept.iloc[-1]["end"]),
        "name": kept.iloc[0]["chrom"],
    }])


def enrichment(amap, n, statistic):
    """`get_score` from the original, verbatim: the central n x n block."""
    c = amap.shape[0] // 2
    if c < n:
        raise ValueError(f"central pixel value {n} too large, max {c}")
    lo, hi = c - n // 2, c + n // 2 + 1
    block = amap[lo:hi, lo:hi]
    return {"mean": np.nanmean, "median": np.nanmedian, "std": np.nanstd}[statistic](block)


def main():
    smk = snakemake  # noqa: F821
    center = int(smk.params.score["center"])
    exp_kw = dict(smk.params.expected)
    sample = smk.wildcards.sample

    for path in (smk.output.npz, smk.log[0]):
        os.makedirs(os.path.dirname(path), exist_ok=True)

    rois = pd.read_csv(smk.input.roi_table, sep="\t")
    r = rois[(rois["locus"] == smk.params.locus)
             & (rois["roi"] == smk.wildcards.roi)].iloc[0]
    start, end = int(r["start"]), int(r["end"])

    loops = pd.read_csv(smk.input.loops, sep="\t")
    loops = loops[(loops["start1"] > start) & (loops["end1"] > start)
                  & (loops["start2"] < end) & (loops["end2"] < end)]
    loops = loops.copy().reset_index(drop=True)

    clr = cooler.Cooler(f"{smk.input.mcool}::resolutions/{smk.params.resolution}")
    view = balanced_span(clr)

    with open(smk.log[0], "w") as fh:
        fh.write(f"{sample}: {len(loops)} loops in "
                 f"{smk.params.locus}/{smk.wildcards.roi} "
                 f"({r['chrom']}:{start}-{end})\n")
        fh.write(f"view_df (balanced span): {view.iloc[0].to_dict()}\n")

        if loops.empty:
            fh.write("no loops in this ROI -- writing an empty pile-up "
                     "(the original would have skipped the file entirely).\n")
            side = 2 * int(smk.params.flank) // clr.binsize + 1
            np.savez_compressed(smk.output.npz,
                                data=np.full((side, side), np.nan),
                                n_loops=0, flank=int(smk.params.flank))
            pd.DataFrame([{"file": sample, "cooler": smk.params.cooler,
                           "n_loops": 0, "median": np.nan, "mean": np.nan,
                           "std": np.nan}]).to_csv(
                smk.output.score, sep="\t", index=False)
            return 0

        cvd = expected_cis(clr, nproc=smk.threads,
                           ignore_diags=exp_kw["ignore_diags"],
                           intra_only=exp_kw["intra_only"],
                           view_df=view, chunksize=1000000)
        pup = coolpup.pileup(
            clr, loops,
            features_format="bedpe",
            view_df=view,
            flank=int(smk.params.flank),
            min_diag=int(smk.params.min_diag),
            nproc=smk.threads,
            seed=int(smk.params.seed),
            store_stripes=False,
            expected_df=cvd,
            nshifts=int(smk.params.nshifts),
            minshift=int(smk.params.minshift),
            maxshift=int(smk.params.maxshift),
            clr_weight_name="weight",
            local=False,
            trans=False,
            by_strand=False,
            by_distance=False,
            mindist="auto",
            maxdist=None,
            flip_negative_strand=False,
        )
        amap = np.asarray(pup["data"].values[0], dtype=float)
        scores = {s: float(enrichment(amap, center, s))
                  for s in ("mean", "median", "std")}
        fh.write(f"pile-up {amap.shape}, central {center}x{center} "
                 f"mean {scores['mean']:.4f}\n")

        np.savez_compressed(
            smk.output.npz,
            data=amap,
            n_loops=int(len(loops)),
            flank=int(smk.params.flank),
            n=int(pup["n"].values[0]) if "n" in pup else int(len(loops)),
        )
        # Columns are the ground truth's loop_stats*.tsv columns.
        pd.DataFrame([{
            "file": sample, "cooler": smk.params.cooler, "n_loops": int(len(loops)),
            "median": scores["median"], "mean": scores["mean"], "std": scores["std"],
        }]).to_csv(smk.output.score, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
