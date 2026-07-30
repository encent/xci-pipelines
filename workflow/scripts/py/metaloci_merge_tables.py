"""Concatenate the per-run `compartmentalization_full.tsv` into one table.

Reproduces `compartmentalization_full_final.tsv`, the table
`02_01_analysis_compartments_metaloci.ipynb` and
`02_02_scatter_analysis_metaloci.ipynb` read. Column set is identical to the
per-run tables -- there is deliberately no `run` column, because the dataset
names already partition cleanly (`_WT_` = wt, `TAG` = degron, group names =
consensus) and adding one would break every existing consumer.

DEVIATION FROM THE GROUND TRUTH, STATED UP FRONT
------------------------------------------------
The ground-truth `compartmentalization_full_final.tsv` (identical copies at the
root of `METALOCI_NEW/` and `METALOCI_NEW_DEGRON/`, both 43,383 B, both
2026-03-04 20:44) contains **188 rows over 44 datasets: wt + degron only**. The
8 consensus datasets are absent, although `METALOCI_NEW_CONSENSUS/
compartmentalization_full.tsv` exists and holds 40 rows.

`the design notes` section 2.2 specifies "the three runs", so this
rule merges whatever is in `config.metaloci.runs` -- by default all three,
giving 228 rows. The extra 40 are the consensus rows the original never merged.
Testers must compare on the wt + degron subset. Setting
`metaloci.runs: [wt, degron]` reproduces the ground-truth row set exactly.

Row order follows `config.metaloci.runs`, then each run's own order. The
ground-truth order came from `os.listdir()` and is not reproducible; compare on
the `(dataset, signal, merge)` key.
"""

import os
import sys

import pandas as pd

tables = list(snakemake.input.tables)
runs = list(snakemake.params.runs)
out_tsv = snakemake.output.tsv

os.makedirs(os.path.dirname(snakemake.log[0]), exist_ok=True)
log = open(snakemake.log[0], "w")

frames = []
for run, path in zip(runs, tables):
    df = pd.read_csv(path, sep="\t")
    log.write(f"{run}\t{path}\t{len(df)} rows\n")
    frames.append(df)

if not frames:
    log.write("ERROR no input tables\n")
    log.close()
    sys.exit("metaloci_merge_tables: no input tables")

merged = pd.concat(frames, ignore_index=True)

# Guard the one thing a silent schema drift would ruin: the consumers index by
# position in places, so the column set must be identical across runs.
for run, df in zip(runs, frames):
    if list(df.columns) != list(merged.columns):
        log.write(f"ERROR column mismatch in {run}: {list(df.columns)}\n")
        log.close()
        sys.exit(f"metaloci_merge_tables: column mismatch in run {run!r}")

os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
merged.to_csv(out_tsv, index=False, sep="\t")
log.write(
    f"wrote {len(merged)} rows over {merged['dataset'].nunique()} datasets "
    f"to {out_tsv}\n"
)
log.close()
