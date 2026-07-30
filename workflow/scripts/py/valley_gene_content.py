"""Classify each valley by the genes it contains -- the data behind Fig 2g / 4d.

Reads the GENE-FILTERED valley set (plan correction P-2). The published bars are
the filtered counts -- E6 375, C5 301, B1 258, JTG 300, CL30 321 -- not the raw
HMM calls (377 / 304 / 261 / 302 / 324). Pointing this rule at
``P.valleys()`` instead of ``P.valleys_filtered()`` makes every bar two or
three valleys too tall and nothing else looks wrong.

Classification, priority order, highest wins (verbatim from
``xci_valleys_check_LAST.ipynb`` cell 10)::

    escaping-gene-valley       overlaps >= 1 escapee gene
    silent-gene-valley         overlaps >= 1 silent gene, no escapee
    no-expressed-gene-valley   overlaps >= 1 GTF gene body, neither of the above
    no-gene                    overlaps no annotated gene body at all

"Escapee" and "silent" are not computed here. They arrive as per-clone BED6
deliveries in ``resources/fixtures/{escapees,silent}_per_clone/`` -- the
collaborator's RNA-seq allelic-ratio call, escapee := allelic ratio > 0.1 in
BOTH replicates. The alias map from our clone/condition grammar to those
filenames lives at the top of ``30_valleys.smk``; never re-derive it here.

Gene bodies come from the Ensembl GRCm38.102 GTF, ``feature == "gene"``, with
Ensembl chromosome names (``X``, not ``chrX``) and 1-based closed coordinates,
so ``start`` is decremented by one on load. The file is ~1 GB; it is read in
chunks and immediately reduced to one chromosome.

Also emits the size class used by the (non-paper) size-vs-category confusion
matrix, with the original's boundaries: <=25 kb, <=50 kb, <=100 kb, <=500 kb,
over 500 kb.

Emits one row per valley::

    chrom start end size has_escapee has_silent has_gene category size_class
"""

import os
import sys

import bioframe
import pandas as pd

log_path = snakemake.log[0]
os.makedirs(os.path.dirname(log_path), exist_ok=True)
os.makedirs(os.path.dirname(snakemake.output.tsv), exist_ok=True)

GENE_BED_COLUMNS = ["chrom", "start", "end", "gene_name", "score", "strand"]
SIZE_CLASSES = [
    (25_000, "under 25 Kbp"),
    (50_000, "25-50 Kbp"),
    (100_000, "50-100 Kbp"),
    (500_000, "100-500 Kbp"),
]


def size_class(width: int) -> str:
    for limit, label in SIZE_CLASSES:
        if width <= limit:
            return label
    return "over 500 Kbp"


def read_gene_bed(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None, names=GENE_BED_COLUMNS)
    return df[["chrom", "start", "end", "gene_name"]]


def read_gtf_genes(path: str, chrom: str) -> pd.DataFrame:
    """chrX gene bodies from an Ensembl GTF, converted to 0-based half-open."""
    ensembl = chrom[3:] if chrom.startswith("chr") else chrom
    parts = []
    for chunk in pd.read_csv(
        path,
        sep="\t",
        comment="#",
        header=None,
        usecols=[0, 2, 3, 4],
        names=["chrom", "feature", "start", "end"],
        dtype={"chrom": str, "feature": str, "start": int, "end": int},
        chunksize=200_000,
    ):
        parts.append(chunk[(chunk["chrom"] == ensembl) & (chunk["feature"] == "gene")])
    genes = pd.concat(parts, ignore_index=True)[["chrom", "start", "end"]].copy()
    genes["chrom"] = chrom
    genes["start"] -= 1
    return genes


def overlaps_any(valleys: pd.DataFrame, other: pd.DataFrame, probe: str) -> pd.Series:
    """True per valley row if it overlaps at least one interval of `other`."""
    if len(other) == 0:
        return pd.Series(False, index=valleys.index)
    hits = bioframe.overlap(valleys, other, how="left", return_index=True)
    touched = set(hits.loc[hits[probe].notna(), "index"])
    return pd.Series([i in touched for i in valleys.index], index=valleys.index)


with open(log_path, "w") as log:

    def say(msg: str) -> None:
        log.write(msg + "\n")
        log.flush()
        print(msg, file=sys.stderr)

    valleys = bioframe.read_table(snakemake.input.bed, schema="bed3")
    valleys = valleys[["chrom", "start", "end"]].reset_index(drop=True)
    chrom = valleys["chrom"].iloc[0] if len(valleys) else "chrX"
    say(f"{len(valleys)} valleys (gene-filtered set) on {chrom}")

    escapees = read_gene_bed(snakemake.input.escapees)
    silent = read_gene_bed(snakemake.input.silent)
    genes = read_gtf_genes(snakemake.input.gtf, chrom)
    say(
        f"{len(escapees)} escapee genes, {len(silent)} silent genes, "
        f"{len(genes)} annotated {chrom} gene bodies"
    )

    has_escapee = overlaps_any(valleys, escapees[["chrom", "start", "end"]], "index_")
    has_silent = overlaps_any(valleys, silent[["chrom", "start", "end"]], "index_")
    has_gene = overlaps_any(valleys, genes, "index_")

    out = valleys.copy()
    out["size"] = out["end"] - out["start"]
    out["has_escapee"] = has_escapee
    out["has_silent"] = has_silent
    out["has_gene"] = has_gene
    out["category"] = [
        "escaping-gene-valley"
        if e
        else "silent-gene-valley"
        if s
        else "no-expressed-gene-valley"
        if g
        else "no-gene"
        for e, s, g in zip(has_escapee, has_silent, has_gene)
    ]
    out["size_class"] = [size_class(int(w)) for w in out["size"]]

    for category, n in out["category"].value_counts().items():
        say(f"  {category:26s} {n}")

    out.to_csv(snakemake.output.tsv, sep="\t", index=False)
    say(f"wrote {snakemake.output.tsv}")
