"""Drop valleys overlapping the three excluded gene bodies.

Behind Fig 2g / 4d. The paper's bar totals are the FILTERED counts, not the raw
HMM calls:

    filtered   E6 375, C5 301, B1 258, JTG 300, CL30 321
    raw        E6 377, C5 304, B1 261, JTG 302, CL30 324

Both sets are kept on disk. Everything gene-related downstream --
``valley_gene_content``, ``valley_overlap_dtag``, ``valley_xa_xi_overlap`` --
reads THIS output (plan correction P-2). Getting that wrong is the single most
likely silent divergence in the pipeline: the numbers stay plausible and every
bar is 2-3 valleys too tall.

The excluded genes (``config.valleys.exclude_genes``) are Mid1, Tmem29 and
Firre: three loci whose gene bodies are constitutively H3K27me3-depleted for
reasons unrelated to X inactivation, so a valley called on them is an artefact
of the assay rather than a silencing feature.

Note what is NOT dropped: ``chrX:0-3,285,000``. That interval is the mm10
assembly gap, it is an artefact too -- but it does not overlap any of the three
genes, so it survives here and is counted in the published bars. See ruling
O-7 and ``call_valleys.py``.

A valley is removed if it overlaps an excluded gene by even one base
(``bioframe.overlap(..., how='inner')``), matching the original.

Emits BED3, same format as the input.
"""

import os
import sys

import bioframe
import pandas as pd

log_path = snakemake.log[0]
os.makedirs(os.path.dirname(log_path), exist_ok=True)
os.makedirs(os.path.dirname(snakemake.output.bed), exist_ok=True)

with open(log_path, "w") as log:

    def say(msg: str) -> None:
        log.write(msg + "\n")
        log.flush()
        print(msg, file=sys.stderr)

    _src = getattr(snakemake.params, "source", "computed")
    say(
        f"valley source: {_src.upper()}"
        + (
            "  <- Fig 4d fixture (resources/fixtures/valleys_fig4d/). This panel "
            "is drawn AS PUBLISHED, not from our computation. See D-18."
            if _src == "fixture"
            else ""
        )
    )

    valleys = bioframe.read_table(snakemake.input.bed, schema="bed3")
    n_before = len(valleys)

    excluded = pd.DataFrame(
        [
            {
                "chrom": str(g["chrom"]),
                "start": int(g["start"]),
                "end": int(g["end"]),
                "name": str(g.get("name", "")),
            }
            for g in snakemake.params.exclude
        ]
    )
    say(f"excluding valleys overlapping {len(excluded)} gene bodies:")
    for _, g in excluded.iterrows():
        say(f"  {g['name']:8s} {g['chrom']}:{g['start']}-{g['end']}")

    if n_before == 0:
        filtered = valleys
        say("input has no valleys; nothing to filter")
    else:
        hits = bioframe.overlap(
            valleys[["chrom", "start", "end"]],
            excluded[["chrom", "start", "end"]],
            how="inner",
        )
        # Keyed on the coordinate triple, exactly as the original did -- valley
        # intervals are unique within a file, so this is a safe key.
        bad = set(zip(hits["chrom"], hits["start"], hits["end"]))
        keep = [
            (r["chrom"], r["start"], r["end"]) not in bad
            for _, r in valleys.iterrows()
        ]
        filtered = valleys[keep].reset_index(drop=True)

    say(f"{n_before} -> {len(filtered)} valleys ({n_before - len(filtered)} removed)")

    filtered[["chrom", "start", "end"]].to_csv(
        snakemake.output.bed, sep="\t", index=False, header=False
    )
