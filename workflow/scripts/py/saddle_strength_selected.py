"""Collect the selected saddle strengths of one scope + ROI into one table.

The tail of ``01_06_compartments_refine.py``::

    saddle_strength_selected.to_csv(
        f"{output_folder}/saddle_strength_selected_{locus_type}.tsv",
        sep="\\t", index=False)

with columns `file` and `value` and one row per cooler, in sorted-basename
order. That exact shape is preserved -- it is what `02_00_analysis_compartments
_saddle.ipynb` (EFig 3b) and `02_02_scatter_analysis_metaloci.ipynb` (Fig 3c)
read back.

The richer per-cooler selection table is written beside it, so a reader of the
strengths can see which eigenvector each came from, whether it was
sign-flipped, and how strongly it correlated with GC.
"""
import os
import sys

import pandas as pd


def main():
    smk = snakemake  # noqa: F821
    os.makedirs(os.path.dirname(smk.output.selected), exist_ok=True)
    os.makedirs(os.path.dirname(smk.log[0]), exist_ok=True)

    frames = [pd.read_csv(p, sep="\t") for p in smk.input.selection]
    table = pd.concat(frames, ignore_index=True).sort_values("file").reset_index(drop=True)

    table[["file", "value"]].to_csv(smk.output.selected, sep="\t", index=False)
    table.to_csv(smk.output.choice, sep="\t", index=False)

    with open(smk.log[0], "w") as fh:
        fh.write(f"{smk.wildcards.cscope} / {smk.wildcards.roi}: "
                 f"{len(table)} coolers\n")
        by_eig = table["eigenvector"].value_counts().to_dict()
        fh.write(f"eigenvector choices: {by_eig}\n")
        n_flipped = int(table["sign_flipped"].sum())
        fh.write(f"sign-flipped to GC-rich-positive: {n_flipped}/{len(table)}\n")
        if n_flipped:
            fh.write("Those flips are recorded, not silent. If the ground "
                     "truth's orientation differs, its saddle plots are "
                     "mirrored relative to ours (plan 9.5, risk R-3).\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
