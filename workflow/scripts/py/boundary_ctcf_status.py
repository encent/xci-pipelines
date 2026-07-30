"""CTCF-peak status at every valley boundary -- the data behind Fig 2h / EFig 6b.

Boundaries are re-derived HERE from the gene-content table, i.e. from the
GENE-FILTERED valley set carrying its gene class (plan correction P-2). They are
deliberately NOT read from ``boundaries/all/``: those files are ground truth
from the ``Stackups_*`` lineage, which predates the 3-gene filter and carries no
gene class. Both sets exist on purpose.

The proximity rule is SIDE-AWARE and asymmetric -- 50 kb inward (into the
valley), 10 kb outward. Verbatim from ``xci_valleys_check_LAST.ipynb``::

    L:  same_chrom & (cp_end   > pos - outside) & (cp_start < pos + inside)
    R:  same_chrom & (cp_start < pos + outside) & (cp_end   > pos - inside)

There is NO intra-valley CTCF filter. An earlier version restricted peaks to
those falling inside a called valley; the final convention (the one that
produces the published Fig 2h counts, and the one that reconciles with the
collaborator's independent bedtools analysis) uses every chrX peak. See
``xci_valleys_check_LAST.ipynb``: "intra-valley CTCF restriction removed".

Peaks come from ``resources/fixtures/CTCFpeak_per_clone/`` -- the collaborator's
consensus peak delivery -- NEVER from MACS2 (plan correction N-5). The alias map
from our grammar to those filenames lives at the top of ``30_valleys.smk``.

The 2x2 Fisher exact test (escaping vs other x CTCF vs no CTCF, TWO-SIDED,
``scipy.stats.fisher_exact`` default) is computed here and carried on every row
as constant columns, so the panel and ``paper_numbers.tsv`` read one file.

Acceptance: escaping-boundary counts per WT clone are 110, 96, 74, 76, 76
(Fig 2h); E6A7 110 (80 CTCF+ / 30 CTCF-) and F3 36 (23 / 13) for EFig 6b.

Emits one row per boundary::

    chrom start end side valley_class is_escaping has_ctcf
    fisher_p fisher_odds_ratio n_boundaries n_escaping
"""

import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

log_path = snakemake.log[0]
os.makedirs(os.path.dirname(log_path), exist_ok=True)
os.makedirs(os.path.dirname(snakemake.output.tsv), exist_ok=True)

inside_flank = int(snakemake.params.inside)
outside_flank = int(snakemake.params.outside)
ESCAPING = "escaping-gene-valley"

with open(log_path, "w") as log:

    def say(msg: str) -> None:
        log.write(msg + "\n")
        log.flush()
        print(msg, file=sys.stderr)

    valleys = pd.read_csv(snakemake.input.gene_content, sep="\t")

    peaks = pd.read_csv(
        snakemake.input.peaks,
        sep="\t",
        header=None,
        usecols=[0, 1, 2],
        names=["chrom", "start", "end"],
        comment="#",
    )
    say(f"{len(valleys)} valleys, {len(peaks)} CTCF consensus peaks")

    rows = []
    for _, v in valleys.iterrows():
        for side, pos in (("L", int(v["start"])), ("R", int(v["end"]))):
            rows.append(
                {
                    "chrom": v["chrom"],
                    "start": pos,
                    "end": pos + 1,
                    "side": side,
                    "valley_class": v["category"],
                }
            )
    boundaries = pd.DataFrame(
        rows, columns=["chrom", "start", "end", "side", "valley_class"]
    )

    cp_chrom = peaks["chrom"].values
    cp_start = peaks["start"].values
    cp_end = peaks["end"].values

    has_ctcf = np.zeros(len(boundaries), dtype=bool)
    for i, (_, row) in enumerate(boundaries.iterrows()):
        pos = int(row["start"])
        same = cp_chrom == row["chrom"]
        if row["side"] == "L":
            near = (
                same
                & (cp_end > pos - outside_flank)
                & (cp_start < pos + inside_flank)
            )
        else:
            near = (
                same
                & (cp_start < pos + outside_flank)
                & (cp_end > pos - inside_flank)
            )
        has_ctcf[i] = bool(near.any())

    boundaries["is_escaping"] = boundaries["valley_class"] == ESCAPING
    boundaries["has_ctcf"] = has_ctcf

    escaping = boundaries[boundaries["is_escaping"]]
    other = boundaries[~boundaries["is_escaping"]]
    a = int(escaping["has_ctcf"].sum())
    b = int(other["has_ctcf"].sum())
    c = int((~escaping["has_ctcf"]).sum())
    d = int((~other["has_ctcf"]).sum())

    if min(a + c, b + d) == 0:
        odds, pval = float("nan"), float("nan")
        say("WARNING: one margin of the 2x2 table is empty; Fisher test skipped")
    else:
        odds, pval = fisher_exact([[a, b], [c, d]])

    say(
        f"boundaries: {len(boundaries)} total, {len(escaping)} escaping "
        f"({a} CTCF+ / {c} CTCF-), {len(other)} other ({b} CTCF+ / {d} CTCF-)"
    )
    say(
        f"Fisher exact (two-sided): p = {pval:.3e}  odds ratio = {odds:.4f}  "
        f"flanks {inside_flank} inward / {outside_flank} outward"
    )

    boundaries["fisher_p"] = pval
    boundaries["fisher_odds_ratio"] = odds
    boundaries["n_boundaries"] = len(boundaries)
    boundaries["n_escaping"] = len(escaping)

    boundaries.to_csv(snakemake.output.tsv, sep="\t", index=False)
    say(f"wrote {snakemake.output.tsv}")
