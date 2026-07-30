"""Saddle matrices and compartment strengths for all three eigenvectors.

Reimplements the saddle half of ``01_05_compartments_cooltools.py``:

    cvd = cooltools.expected_cis(clr, view_df=view_df)          # defaults
    interaction_sum, interaction_count = cooltools.saddle(
        clr, cvd, eig_track[['chrom','start','end',lbl]], 'cis',
        n_bins=38, qrange=(0.025, 0.975), view_df=view_df)
    saddle_freqs = interaction_sum / interaction_count

TWO strength numbers, and they are NOT the same quantity:

`Saddle_values_{name}.tsv` -- one row, columns E1/E2/E3, holding the original's
own score::

    (mean(S[:8, :8]) + mean(S[32:, 32:])) / (mean(S[:8, 32:]) + mean(S[32:, :8]))

This is what `01_06` selects from, and therefore the number behind EFig 3b and
Fig 3c. Note the hard-coded 8 and 32: with `n_bins = 38` the matrix is 40x40
(two flanking outlier bins), so `:8` is the extreme-B corner and `32:` the
extreme-A one -- `strength_extent` in config.

`saddle_strength_{E}.tsv` -- `cooltools.api.saddle.saddle_strength(sum, count)`,
the full profile over extents, plus the value at index 8 that the original
computed and then commented out of its output.

NO PLOT: the `.npz` carries `interaction_sum`, `interaction_count` and the
ratio, which is everything `panel_saddle` needs.
"""
import os
import sys

import cooler
import cooltools
import numpy as np
import pandas as pd
from cooltools.api.saddle import saddle_strength


def corner_score(S, extent):
    """The original's (AA + BB) / (AB + BA), with its double np.nanmean."""
    top_left = np.nanmean(np.nanmean(S[:extent, :extent]))
    bottom_right = np.nanmean(np.nanmean(S[-extent:, -extent:]))
    top_right = np.nanmean(np.nanmean(S[:extent, -extent:]))
    bottom_left = np.nanmean(np.nanmean(S[-extent:, :extent]))
    return {
        "top_left": float(top_left), "bottom_right": float(bottom_right),
        "top_right": float(top_right), "bottom_left": float(bottom_left),
        "score": float((top_left + bottom_right) / (top_right + bottom_left)),
    }


def main():
    smk = snakemake  # noqa: F821
    name = smk.wildcards.name
    roi = smk.wildcards.roi
    eig_labels = list(smk.params.eigs)
    n_bins = int(smk.params.n_bins)
    q_lo, q_hi = [float(q) for q in smk.params.qrange]
    extent = int(smk.params.extent)

    os.makedirs(os.path.dirname(smk.output.values), exist_ok=True)
    os.makedirs(os.path.dirname(smk.log[0]), exist_ok=True)

    rois = pd.read_csv(smk.input.roi_table, sep="\t")
    r = rois[(rois["locus"] == smk.params.locus) & (rois["roi"] == roi)].iloc[0]
    chrom, start, end = r["chrom"], int(r["start"]), int(r["end"])
    view_df = pd.DataFrame({"chrom": [chrom], "start": [start],
                            "end": [end], "name": [chrom]})

    clr = cooler.Cooler(f"{smk.input.mcool}::resolutions/{smk.params.resolution}")
    track = pd.read_csv(smk.input.eigs, sep="\t")
    cvd = cooltools.expected_cis(clr=clr, view_df=view_df)

    scores, strengths_at_extent = {}, {}
    with open(smk.log[0], "w") as fh:
        fh.write(f"{name} / {roi}: n_bins={n_bins} qrange=({q_lo}, {q_hi}) "
                 f"extent={extent}\n")
        for lbl, npz_path, strength_path in zip(
                eig_labels, smk.output.npz, smk.output.strength):
            try:
                isum, icount = cooltools.saddle(
                    clr, cvd, track[["chrom", "start", "end", lbl]], "cis",
                    n_bins=n_bins, qrange=(q_lo, q_hi), view_df=view_df)
            except Exception as exc:                      # noqa: BLE001
                # The original wrapped all three in one try and `continue`d out
                # of the whole cooler; we fail only the eigenvector that failed.
                fh.write(f"  {lbl}: FAILED ({exc}); writing NaN\n")
                side = n_bins + 2
                nan = np.full((side, side), np.nan)
                np.savez_compressed(npz_path, interaction_sum=nan,
                                    interaction_count=nan, saddle=nan,
                                    n_bins=n_bins, qrange=(q_lo, q_hi))
                pd.DataFrame([{"eigenvector": lbl, "extent": extent,
                               "saddle_strength": np.nan}]).to_csv(
                    strength_path, sep="\t", index=False)
                scores[lbl] = np.nan
                strengths_at_extent[lbl] = np.nan
                continue

            S = isum / icount
            c = corner_score(S, extent)
            scores[lbl] = c["score"]

            profile = saddle_strength(isum, icount)
            strengths_at_extent[lbl] = float(profile[extent])

            np.savez_compressed(
                npz_path,
                interaction_sum=isum, interaction_count=icount, saddle=S,
                n_bins=n_bins, qrange=(q_lo, q_hi),
                saddle_strength_profile=profile,
                **{k: v for k, v in c.items()},
            )
            pd.DataFrame({
                "eigenvector": lbl,
                "extent": np.arange(len(profile)),
                "saddle_strength": profile,
            }).to_csv(strength_path, sep="\t", index=False)
            fh.write(f"  {lbl}: (AA+BB)/(AB+BA) = {c['score']:.4f}   "
                     f"saddle_strength[{extent}] = {profile[extent]:.4f}\n")

        # The original's Saddle_values TSV: one row, the CORNER score per E.
        pd.DataFrame([dict(file=name, **{e: scores[e] for e in eig_labels})]).to_csv(
            smk.output.values, sep="\t", index=False)
        fh.write(f"wrote {smk.output.values} (corner scores; "
                 f"saddle_strength[{extent}] is in the per-E tables)\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
