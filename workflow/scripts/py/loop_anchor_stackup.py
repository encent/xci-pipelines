"""Signal around both feet of every loop.

Reimplements the stackup half of ``01_04_loops_pileups_and_stackups.py``:

    sub = concat(loops[chrom1, start1, end1], loops[chrom2, start2, end2])
    mid = (start + end) // 2
    bbi.stackup(bw, sub.chrom, mid - 100_000, mid + 100_000, bins=100)
    profile = np.nanmean(stackup, axis=0)

DECLARED DEVIATION D-9. The original read `data/bigwigs/merged_bw/` -- deepTools
means of the bigWigs as originally delivered, which were never csaw/TMM
normalised. This reads `merged20`, the csaw-normalised replicate means, on the
archaeologist's recommendation (the archaeology notes section J.5).

Scope: this is the ONLY place the switch happens. `loop_pileup` reads coolers,
so every paper pile-up is unaffected, and these stackups never reached the
paper. Expect a structural test, not a ground-truth diff.
"""
import os
import sys

import bbi
import numpy as np
import pandas as pd


def anchors(loops):
    """Both feet of every loop, concatenated -- foot 1 of all loops, then foot 2."""
    a = loops[["chrom1", "start1", "end1"]].copy()
    a.columns = ["chrom", "start", "end"]
    b = loops[["chrom2", "start2", "end2"]].copy()
    b.columns = ["chrom", "start", "end"]
    out = pd.concat([a, b], ignore_index=True)
    out["midpoint"] = (out["start"] + out["end"]) // 2
    return out


def main():
    smk = snakemake  # noqa: F821
    flank = int(smk.params.flank)
    nbins = int(smk.params.nbins)
    tracks = list(smk.params.tracks)
    bigwigs = [smk.input.bw1, smk.input.bw2]

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

    with open(smk.log[0], "w") as fh:
        fh.write(f"signal={smk.wildcards.signal} pair={smk.wildcards.pair}\n")
        fh.write(f"tracks (D-9: merged20, csaw-normalised): {tracks}\n")
        fh.write(f"{len(loops)} loops -> {2 * len(loops)} anchors in "
                 f"{smk.params.locus}/{smk.wildcards.roi}\n")

        if loops.empty:
            fh.write("no loops in this ROI -- writing an empty stackup.\n")
            np.savez_compressed(
                smk.output.npz,
                matrices=np.full((2, 0, nbins), np.nan),
                profiles=np.full((2, nbins), np.nan),
                tracks=np.array(tracks), n_anchors=0, flank=flank, nbins=nbins)
            return 0

        sub = anchors(loops)
        mats, profs = [], []
        for track, bw in zip(tracks, bigwigs):
            m = bbi.stackup(bw, sub["chrom"], sub["midpoint"] - flank,
                            sub["midpoint"] + flank, bins=nbins)
            mats.append(np.asarray(m, dtype=float))
            profs.append(np.nanmean(m, axis=0))
            fh.write(f"  {track}: stackup {np.shape(m)}\n")

        np.savez_compressed(
            smk.output.npz,
            matrices=np.stack(mats),
            profiles=np.stack(profs),
            tracks=np.array(tracks),
            n_anchors=int(len(sub)),
            n_loops=int(len(loops)),
            flank=flank,
            nbins=nbins,
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
