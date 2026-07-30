"""Pairwise peak-overlap counts within a clone/condition (QC).

Emits the three numbers a Venn needs (a_only, shared, b_only) so the viz stage
does not have to re-derive them.
"""
import itertools
import os

import bioframe as bf
import pandas as pd

def load(path):
    df = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2],
                     names=["chrom", "start", "end"], comment="#")
    return bf.merge(df.sort_values(["chrom", "start"]))

peaks = {os.path.basename(p)[: -len("_peaks.bed")]: load(p) for p in snakemake.input.peaks}

rows = []
for a, b in itertools.combinations(sorted(peaks), 2):
    da, db = peaks[a], peaks[b]
    ov = bf.overlap(da, db, how="left", return_index=True)
    shared_a = ov.loc[ov["index_"].notna(), "index"].nunique()
    ov2 = bf.overlap(db, da, how="left", return_index=True)
    shared_b = ov2.loc[ov2["index_"].notna(), "index"].nunique()
    rows.append({"a": a, "b": b, "n_a": len(da), "n_b": len(db),
                 "a_overlapping_b": shared_a, "b_overlapping_a": shared_b,
                 "a_only": len(da) - shared_a, "b_only": len(db) - shared_b})

pd.DataFrame(rows).to_csv(snakemake.output.tsv, sep="\t", index=False)
