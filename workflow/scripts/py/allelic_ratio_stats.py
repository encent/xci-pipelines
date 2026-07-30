"""Allelic-ratio summary tables per locus and ROI.

Reimplements the `AR_stats` half of ``01_07_get_allelic_ratio_gtf.py``. Nine
rows -- sum / mean / median crossed with three gene subsets -- one column per
clone, for the genes lying fully inside the ROI (``start >= roi.start and
end <= roi.end``, both inclusive, as the original wrote it).

The three subsets, verbatim::

    all       ratio >= 0
    non_zero  ratio >  0
    >0.1      ratio >  0.1

Written with ``sep='\\t'`` despite the ``.csv`` extension -- that is what 01_07
did and what `02_02_scatter_analysis_metaloci.ipynb` reads back. Changing it
would silently break Fig 3c.

The degron table carries no coordinates: 01_07 back-fills `start`/`end` by
joining its `X` column against `gene_id` in the GTF, keeping a coordinate only
when EXACTLY ONE gene matches. That single-match rule is preserved -- it is the
reason a handful of degron genes drop out of every ROI window.
"""
import os
import sys

import bioframe
import bioframe.sandbox.gtf_io
import pandas as pd

ROW_ORDER = ["sum_all", "mean_all", "median_all",
             "sum_non_zero", "mean_non_zero", "median_non_zero",
             "sum_>0.1", "mean_>0.1", "median_>0.1"]
SUBSETS = {"all": lambda s: s >= 0, "non_zero": lambda s: s > 0,
           ">0.1": lambda s: s > 0.1}


def gene_coordinates(gtf_path, chrom):
    """gene_id -> (start, end), only where exactly one gene carries that id."""
    genes = (bioframe.read_table(gtf_path, schema="gtf")
             .query('feature == "gene" and chrom == @chrom')
             .reset_index(drop=True))
    attrs = bioframe.sandbox.gtf_io.parse_gtf_attributes(
        genes["attributes"], kv_sep=" ", item_sep=";")
    df = pd.concat([genes[["start", "end"]], attrs[["gene_id"]]], axis=1)
    counts = df["gene_id"].value_counts()
    unique = df[df["gene_id"].isin(counts[counts == 1].index)]
    return unique.set_index("gene_id")[["start", "end"]]


def summarise(ar, columns, start, end):
    """The 9 x n_clones table, in the original's row order."""
    out = pd.DataFrame(index=ROW_ORDER, columns=columns, dtype=object)
    inside = ar[(ar["start"] >= start) & (ar["end"] <= end)]
    for col in columns:
        for label, predicate in SUBSETS.items():
            values = inside.loc[predicate(inside[col]), col]
            out.loc[f"sum_{label}", col] = values.sum()
            out.loc[f"mean_{label}", col] = values.mean()
            out.loc[f"median_{label}", col] = values.median()
    return out


def main():
    smk = snakemake  # noqa: F821
    locus, roi = smk.wildcards.locus, smk.wildcards.roi
    wt_map, degron_map = dict(smk.params.wt_map), dict(smk.params.degron_map)

    os.makedirs(os.path.dirname(smk.output.wt), exist_ok=True)
    os.makedirs(os.path.dirname(smk.log[0]), exist_ok=True)

    rois = pd.read_csv(smk.input.roi_table, sep="\t")
    r = rois[(rois["locus"] == locus) & (rois["roi"] == roi)].iloc[0]
    start, end = int(r["start"]), int(r["end"])

    with open(smk.log[0], "w") as fh:
        fh.write(f"{locus}/{roi}: {r['chrom']}:{start}-{end}\n")

        wt = pd.read_csv(smk.input.wt, sep=",", index_col=0).reset_index(drop=True)
        wt_cols = [c for c in wt_map if c in wt.columns]
        missing = [c for c in wt_map if c not in wt.columns]
        if missing:
            fh.write(f"WARNING: WT table has no column(s) {missing}\n")
        summarise(wt, wt_cols, start, end).to_csv(
            smk.output.wt, sep="\t", index=True, header=True)
        fh.write(f"WT: {len(wt)} genes, {len(wt_cols)} clones -> "
                 f"{smk.output.wt}\n")

        degron = pd.read_csv(smk.input.degron, sep=",").reset_index(drop=True)
        coords = gene_coordinates(smk.input.gtf, smk.params.chrom)
        joined = degron.join(coords, on="X")
        degron["start"] = joined["start"]
        degron["end"] = joined["end"]
        n_placed = int(degron["start"].notnull().sum())
        fh.write(f"degron: {len(degron)} genes, {n_placed} with a unique "
                 f"gene_id match in the GTF\n")
        degron_cols = [c for c in degron_map if c in degron.columns]
        missing = [c for c in degron_map if c not in degron.columns]
        if missing:
            fh.write(f"WARNING: degron table has no column(s) {missing}\n")
        summarise(degron, degron_cols, start, end).to_csv(
            smk.output.degron, sep="\t", index=True, header=True)
        fh.write(f"degron: {len(degron_cols)} clones -> {smk.output.degron}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
