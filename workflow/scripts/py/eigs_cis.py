"""Cis eigenvectors over one ROI, with a RECORDED sign convention.

Reimplements the eigenvector half of ``01_05_compartments_cooltools.py``:

    view_df = single row  chrX:<roi.start>-<roi.end>  named 'chrX'
    cooltools.eigs_cis(clr, gc_cov, view_df=view_df, n_eigs=3)
    -> Comp_{name}.tsv  (chrom, start, end, E1, E2, E3)
    -> Comp_E{1,2,3}_{name}.bw  via bioframe.to_bigwig

PLAN 9.5 -- THE SIGN.
`eigs_cis` returns eigenvectors whose global sign is arbitrary: the same data
can yield E and -E on two runs. The curation fixture `compartments*.json`
records only WHICH eigenvector was chosen, never its ORIENTATION. If our
orientation differs from the ground truth's, the eigenvector track is inverted
and every saddle plot built on it is mirrored -- a wrong figure that looks
entirely plausible. That is risk R-3, and it reaches 15 of the 25 panels.

So: each eigenvector is oriented such that Spearman(E, GC) > 0, i.e. GC-rich is
positive (the A compartment, by the usual convention), and the flip -- applied
or not -- is written to `orientation.tsv` next to the track. The three
correlations are recorded either way, because `eigenvector_selection: auto`
selects on max |Spearman(E, GC)| and a reader needs to see how close the
runners-up were.

`canonicalise_sign: false` leaves cooltools' raw output alone and still records
the correlations, for a bit-for-bit diff against the ground truth.
"""
import os
import sys

import bioframe
import cooler
import cooltools
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def main():
    smk = snakemake  # noqa: F821
    name = smk.wildcards.name
    roi = smk.wildcards.roi
    eig_labels = list(smk.params.eigs)

    os.makedirs(os.path.dirname(smk.output.tsv), exist_ok=True)
    os.makedirs(os.path.dirname(smk.log[0]), exist_ok=True)

    rois = pd.read_csv(smk.input.roi_table, sep="\t")
    r = rois[(rois["locus"] == smk.params.locus) & (rois["roi"] == roi)].iloc[0]
    chrom, start, end = r["chrom"], int(r["start"]), int(r["end"])
    view_df = pd.DataFrame({"chrom": [chrom], "start": [start],
                            "end": [end], "name": [chrom]})

    clr = cooler.Cooler(f"{smk.input.mcool}::resolutions/{smk.params.resolution}")
    gc_cov = pd.read_csv(smk.input.gc, sep="\t")
    chromsizes = bioframe.read_chromsizes(smk.input.chrom_sizes)

    with open(smk.log[0], "w") as fh:
        fh.write(f"{name} / {roi}: view {chrom}:{start}-{end} (single row, "
                 f"named {chrom!r} -- no arms file)\n")

        cis_eigs = cooltools.eigs_cis(clr, gc_cov, view_df=view_df,
                                      n_eigs=int(smk.params.n_eigs))
        track = cis_eigs[1][["chrom", "start", "end"] + eig_labels].copy()

        # --- sign convention (plan 9.5) --------------------------------
        gc = gc_cov.set_index(["chrom", "start", "end"])["GC"]
        joined = track.set_index(["chrom", "start", "end"]).join(gc, how="left")
        rows = []
        for lbl in eig_labels:
            ok = joined[lbl].notnull() & joined["GC"].notnull()
            if int(ok.sum()) >= 3:
                rho = float(spearmanr(joined.loc[ok, lbl], joined.loc[ok, "GC"]).correlation)
            else:
                rho = float("nan")
            flip = bool(smk.params.canonicalise) and np.isfinite(rho) and rho < 0
            if flip:
                track[lbl] = -track[lbl]
            rows.append({
                "roi": roi, "cooler": name, "eigenvector": lbl,
                "spearman_gc_raw": rho,
                "spearman_gc_final": -rho if flip else rho,
                "n_bins": int(ok.sum()),
                "sign_flipped": flip,
                "canonicalise_sign": bool(smk.params.canonicalise),
            })
            fh.write(f"  {lbl}: Spearman(E, GC) = {rho:+.4f}"
                     f"{'   FLIPPED to GC-rich-positive' if flip else ''}\n")

        if any(r_["sign_flipped"] for r_ in rows):
            fh.write("NOTE: a sign flip was applied. An UNRECORDED flip mirrors "
                     "every saddle plot built on this track (plan 9.5, risk "
                     "R-3); this one is recorded in orientation.tsv.\n")

        pd.DataFrame(rows).to_csv(smk.output.orientation, sep="\t", index=False)
        track.to_csv(smk.output.tsv, index=False, sep="\t")
        for lbl, path in zip(eig_labels, smk.output.bw):
            bioframe.to_bigwig(track, chromsizes, path, value_field=lbl)
            fh.write(f"  wrote {path}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
