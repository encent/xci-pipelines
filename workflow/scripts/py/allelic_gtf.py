"""Per-clone GTF with the allele-specific expression ratio injected.

Reimplements the GTF half of ``01_07_get_allelic_ratio_gtf.py``: read the
Ensembl GTF, keep `feature == "gene"` on chromosome X, and append
``allelic_ratio "<value>";`` to the attributes column of every gene the clone's
RNA-seq table names. Missing values become the literal string ``NA``.

This is the "D-score" source for the green colour scale in Fig 3a/b, 5b/c/g/h,
EFig 8c/d and 10c/d.

THE TWO JOINS ARE NOT SYMMETRIC, and that asymmetry is the original's:

* **degron** clones join on ``gene_id`` (the CSV's `X` column).
* **WT** clones join on ``gene_name`` **OR** an exact ``(start, end)`` match
  (the CSV's `name`, `start`, `end` columns). The `or` makes the WT join looser
  and can annotate more than one gene per row; reproduced as-is.

Output is written with ``quoting=QUOTE_NONE, escapechar='\\\\'`` and no header,
like the original, then bgzipped and tabixed for coolbox.

`GRCm38.102_NC.gtf` -- 01_07's other output -- is NOT written here;
`make_noncoding_gtf` in 00_reference.smk owns that path.
"""
import csv
import os
import subprocess
import sys

import bioframe
import bioframe.sandbox.gtf_io
import pandas as pd


def load_genes(gtf_path, chrom):
    genes = (bioframe.read_table(gtf_path, schema="gtf")
             .query('feature == "gene" and chrom == @chrom')
             .reset_index(drop=True))
    attrs = bioframe.sandbox.gtf_io.parse_gtf_attributes(
        genes["attributes"], kv_sep=" ", item_sep=";")
    return genes, pd.concat([genes, attrs], axis=1)


def as_text(value):
    return "NA" if pd.isna(value) else str(value)


def main():
    smk = snakemake  # noqa: F821
    target = smk.wildcards.ar_clone
    wt_map = dict(smk.params.wt_map)
    degron_map = dict(smk.params.degron_map)
    chrom = smk.params.chrom

    os.makedirs(os.path.dirname(smk.output.gtf), exist_ok=True)
    os.makedirs(os.path.dirname(smk.log[0]), exist_ok=True)

    inverse_wt = {v: k for k, v in wt_map.items()}
    inverse_degron = {v: k for k, v in degron_map.items()}

    with open(smk.log[0], "w") as fh:
        genes, concat = load_genes(smk.input.gtf, chrom)
        fh.write(f"{target}: {len(genes)} genes on chromosome {chrom!r} "
                 f"from {smk.input.gtf}\n")

        clone = concat.copy()
        clone["start"] = clone["start"].astype(int)
        clone["end"] = clone["end"].astype(int)
        annotated = 0

        if target in inverse_degron:
            column = inverse_degron[target]
            ar = pd.read_csv(smk.input.degron, sep=",").reset_index(drop=True)
            fh.write(f"degron table column {column!r}, {len(ar)} rows; "
                     f"join on gene_id\n")
            for _, row in ar.iterrows():
                hit = clone[clone["gene_id"] == row["X"]]
                value = as_text(row[column])
                for idx in hit.index:
                    clone.loc[idx, "attributes"] += f' allelic_ratio "{value}";'
                    annotated += 1
        elif target in inverse_wt:
            column = inverse_wt[target]
            ar = pd.read_csv(smk.input.wt, sep=",", index_col=0).reset_index(drop=True)
            fh.write(f"WT table column {column!r}, {len(ar)} rows; "
                     f"join on gene_name OR exact (start, end)\n")
            for _, row in ar.iterrows():
                hit = clone[(clone["gene_name"] == row["name"])
                            | ((clone["start"] == row["start"])
                               & (clone["end"] == row["end"]))]
                value = as_text(row[column])
                for idx in hit.index:
                    clone.loc[idx, "attributes"] += f' allelic_ratio "{value}";'
                    annotated += 1
        else:
            raise ValueError(
                f"{target!r} is in neither hic.allelic_ratio.wt_clone_map nor "
                f"degron_clone_map; known clones are "
                f"{sorted(set(wt_map.values()) | set(degron_map.values()))}"
            )

        clone[list(genes.columns)].to_csv(
            smk.output.gtf, sep="\t", header=None, index=None,
            quoting=csv.QUOTE_NONE, escapechar="\\")
        fh.write(f"annotated {annotated} gene row(s); wrote {smk.output.gtf}\n")

        with open(smk.output.bgz, "wb") as out:
            sort = subprocess.Popen(
                ["sort", "-k1,1", "-k4,4n", smk.output.gtf], stdout=subprocess.PIPE)
            subprocess.run(["bgzip", "-c"], stdin=sort.stdout, stdout=out, check=True)
            sort.wait()
        subprocess.run(["tabix", "-p", "gff", "-f", smk.output.bgz], check=True)
        fh.write(f"indexed {smk.output.bgz}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
