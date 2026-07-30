"""Signal enrichment over an anchor set -> .npz (feeds Fig 1e, 2a, 2d).

Replaces 08_plot_CTCF_enrichment.py / 10_ / 11_, which drew heatmaps directly.
bbi.stackup(bw, chrom, centre +/- flank, bins=nbins); mean profile per track.
"""
import os

import bbi
import numpy as np
import pandas as pd

flank = int(snakemake.params.flank)
nbins = int(snakemake.params.nbins)

anchors = pd.read_csv(
    snakemake.input.anchors, sep="\t", header=None, usecols=[0, 1, 2],
    names=["chrom", "start", "end"], comment="#",
)
chromsizes = dict(
    pd.read_csv(snakemake.input.chrom_sizes, sep="\t", header=None,
                names=["chrom", "length"]).itertuples(index=False, name=None)
)
anchors = anchors[anchors["chrom"].isin(chromsizes)]
centre = ((anchors["start"] + anchors["end"]) // 2).to_numpy()
lo, hi = centre - flank, centre + flank
keep = (lo >= 0) & np.array(
    [hi[i] <= chromsizes[c] for i, c in enumerate(anchors["chrom"])]
)
anchors, lo, hi = anchors[keep], lo[keep], hi[keep]

profiles, labels = {}, []
for bw in snakemake.input.bws:
    track = os.path.basename(bw)[: -len(".bw")]
    stack = bbi.stackup(
        bw, anchors["chrom"].tolist(), lo.tolist(), hi.tolist(), bins=nbins
    )
    profiles[track] = np.nanmean(stack, axis=0)
    labels.append(track)

np.savez_compressed(
    snakemake.output.npz,
    tracks=np.array(labels),
    profiles=np.vstack([profiles[t] for t in labels]),
    flank=flank, nbins=nbins, n_anchors=len(anchors),
)
