"""Binned H3K27me3 coverage, NodTAG vs dTAG, inside or outside valleys -- EFig 6a.

From ``xci_density_plots.ipynb``. The question is whether degron depletion
changes H3K27me3 globally or only at valley edges, so the same two bigWigs are
compared over three different bin sets:

    valley       tiles of the NodTAG valley intervals
    antivalley   tiles of the complement of those intervals
    allcoverage  tiles of the whole chromosome

Bins are always defined by the NodTAG valley call -- the point is to compare two
conditions over one fixed partition. Tiling starts at each interval's own start
and the final tile of an interval is truncated to the interval end, so bins are
not all ``win`` wide; that is the original's ``tile_intervals``.

``snakemake.input.valleys`` is a LIST: the gene-filtered chrX set plus the raw
set for each autosome control chromosome. Every panel in the original has three
columns -- chr7, chrX, and the two pooled -- so both chromosomes have to be in
the table. Rows carry their ``chrom`` and the three statistics are computed per
chromosome; pooling is the panel's job.

Signal is ``pyBigWig.stats(type="mean")`` per bin, with no-data bins becoming
NaN (NOT zero -- these are excluded from the statistics rather than dragged to
the floor, which matters because the antivalley set covers the blacklist).

Three summary statistics per chromosome, carried as constant columns:

    pearson_r    scipy.stats.pearsonr over finite pairs
    ccc          Lin's concordance correlation coefficient,
                 2*cov / (var_x + var_y + (mean_x - mean_y)^2), ddof = 0
    median_lfc   median of log2(dtag / nodtag) over bins where both are > 0

Not a figure. The MA plot and the density scatter are drawn later from this
table. Emits one row per bin::

    chrom start end nodtag dtag log2_ratio pearson_r ccc median_lfc
"""

import os
import sys

import numpy as np
import pandas as pd
import pyBigWig
from scipy import stats

log_path = snakemake.log[0]
os.makedirs(os.path.dirname(log_path), exist_ok=True)
os.makedirs(os.path.dirname(snakemake.output.tsv), exist_ok=True)

window = int(snakemake.params.window)
mask = str(snakemake.params.mask)
clone = snakemake.wildcards.clone
allele = snakemake.wildcards.allele


def complement(intervals, chrom, chrom_size):
    out, pos = [], 0
    for start, end in sorted(intervals):
        if start > pos:
            out.append((pos, start))
        pos = max(pos, end)
    if pos < chrom_size:
        out.append((pos, chrom_size))
    return out


def tile(intervals, bin_size):
    bins = []
    for start, end in intervals:
        pos = start
        while pos < end:
            bins.append((pos, min(pos + bin_size, end)))
            pos += bin_size
    return bins


def bw_means(path, chrom, bins):
    handle = pyBigWig.open(path)
    out = np.empty(len(bins), dtype=float)
    for i, (start, end) in enumerate(bins):
        if end <= start:
            out[i] = np.nan
            continue
        try:
            value = handle.stats(chrom, start, end, type="mean")[0]
            out[i] = np.nan if value is None else value
        except (RuntimeError, ValueError):
            out[i] = np.nan
    handle.close()
    return out


def lin_ccc(x, y):
    sxy = np.cov(x, y, ddof=0)[0, 1]
    return 2 * sxy / (x.var() + y.var() + (x.mean() - y.mean()) ** 2)


with open(log_path, "w") as log:

    def say(msg: str) -> None:
        log.write(msg + "\n")
        log.flush()
        print(msg, file=sys.stderr)

    chrom_sizes = pd.read_csv(
        snakemake.input.chrom_sizes, sep="\t", header=None, names=["chrom", "size"]
    ).set_index("chrom")["size"]

    valley_files = snakemake.input.valleys
    if isinstance(valley_files, str):
        valley_files = [valley_files]
    valleys = pd.concat(
        [
            pd.read_csv(path, sep="\t", header=None, names=["chrom", "start", "end"])
            for path in valley_files
        ],
        ignore_index=True,
    )
    say(
        f"{clone} {allele} {mask} @ {window} bp: {len(valleys)} NodTAG valley "
        f"intervals from {len(valley_files)} file(s) covering "
        f"{sorted(valleys['chrom'].unique())}"
    )

    frames = []
    for chrom in sorted(valleys["chrom"].unique()) or []:
        chrom_size = int(chrom_sizes[chrom])
        own = [
            (int(r["start"]), int(r["end"]))
            for _, r in valleys[valleys["chrom"] == chrom].iterrows()
        ]

        if mask == "valley":
            bins = tile(own, window)
        elif mask == "antivalley":
            bins = tile(complement(own, chrom, chrom_size), window)
        elif mask == "allcoverage":
            bins = tile([(0, chrom_size)], window)
        else:
            raise SystemExit(f"unknown mask {mask!r}; expected valley|antivalley|allcoverage")

        nodtag = bw_means(snakemake.input.nodtag_bw, chrom, bins)
        dtag = bw_means(snakemake.input.dtag_bw, chrom, bins)

        df = pd.DataFrame(
            {
                "chrom": chrom,
                "start": [b[0] for b in bins],
                "end": [b[1] for b in bins],
                "nodtag": nodtag,
                "dtag": dtag,
            }
        )

        both = np.isfinite(nodtag) & np.isfinite(dtag)
        positive = both & (nodtag > 0) & (dtag > 0)
        with np.errstate(divide="ignore", invalid="ignore"):
            df["log2_ratio"] = np.where(positive, np.log2(dtag / nodtag), np.nan)

        if both.sum() >= 2:
            r = stats.pearsonr(nodtag[both], dtag[both])[0]
            ccc = lin_ccc(nodtag[both], dtag[both])
        else:
            r = ccc = np.nan
        median_lfc = (
            float(np.median(df.loc[positive, "log2_ratio"])) if positive.sum() else np.nan
        )

        df["pearson_r"] = r
        df["ccc"] = ccc
        df["median_lfc"] = median_lfc
        frames.append(df)

        say(
            f"  {chrom}: {len(bins)} bins, {int(both.sum())} with signal in both, "
            f"r = {r:.4f}  CCC = {ccc:.4f}  median LFC = {median_lfc:+.4f}"
        )

    out = (
        pd.concat(frames, ignore_index=True)
        if frames
        else pd.DataFrame(
            columns=[
                "chrom", "start", "end", "nodtag", "dtag",
                "log2_ratio", "pearson_r", "ccc", "median_lfc",
            ]
        )
    )
    out.to_csv(snakemake.output.tsv, sep="\t", index=False, na_rep="NA")
    say(f"wrote {len(out)} rows to {snakemake.output.tsv}")
