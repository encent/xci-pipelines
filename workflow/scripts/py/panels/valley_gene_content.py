"""Fig 2g (WT) and Fig 4d (degron) — gene content of H3K27me3 valleys.

Stacked horizontal bars, one per clone, split into the four gene classes, with
the total printed at the end of each bar.

THE TRAP (P-2 / the archaeology notes J.3)
---------------------------
The published bars are the **gene-filtered** valley counts:

    E6 375   C5 301   B1 258   JTG 300   CL30 321

The raw HMM calls are 377 / 304 / 261 / 302 / 324. Both sets of bars look
entirely plausible; nothing about the figure says which one you are looking at.
This panel therefore reads ``P.gene_content(track)``, which the
``valley_gene_content`` rule builds from ``P.valleys_filtered`` — the output of
``filter_valleys_genes``, which removes valleys overlapping Mid1, Tmem29 and
Firre. ``paper_numbers.tsv`` asserts the five totals so a regression here fails
the run instead of quietly redrawing the paper.

Two things that look wrong and are not:

* the leading ``chrX:0-3,285,000`` valley is the mm10 assembly gap. It survives
  the 3-gene filter and IS counted in the published bars, so it stays
  (``valleys.drop_leading_gap: false``);
* row order is E6, C5, B1, JTG, CL30 — the published order, not sorted by size.
  For the degron cohort it is paired, NodTAG above dTAG, clones B1621, E6A7,
  F3 from top.

Verbatim from ``xci_valleys_check_LAST.ipynb`` cell 18: bar height 0.6, total
label at ``x = total + 3``, xlim ``0..max*1.15``, figsize ``(7, max(3, 0.7n+1.5))``.
"""

import os
import sys

sys.path.insert(0, os.path.join(snakemake.params.spec["libdir"],
                                "scripts", "py", "panels"))
sys.path.insert(0, snakemake.params.spec["libdir"])

import numpy as np                                             # noqa: E402
import _common as C                                            # noqa: E402
from lib import plotting as pl                                 # noqa: E402

spec = C.bootstrap(snakemake)
out_path = C.outputs(snakemake)

#: Published row order, top to bottom.
WT_ORDER = ["E6", "C5", "B1", "JTG", "CL30"]
DEGRON_CLONE_ORDER = ["B1621", "E6A7", "F3"]


def _rows_wt(summary):
    order = {c: i for i, c in enumerate(WT_ORDER)}
    rows = sorted(summary, key=lambda r: order.get(r["clone"], 99))
    return list(reversed(rows))          # barh: y=0 is the bottom


def _rows_degron(summary):
    """NodTAG above dTAG within each clone, clones top-to-bottom as published."""
    by_clone = {}
    for row in summary:
        by_clone.setdefault(row["clone"], []).append(row)
    ordered = []
    for clone in reversed(
        DEGRON_CLONE_ORDER + [c for c in sorted(by_clone) if c not in DEGRON_CLONE_ORDER]
    ):
        pair = by_clone.get(clone)
        if not pair:
            continue
        dtag = [r for r in pair if "-dTAG" in r["condition"]]
        nodtag = [r for r in pair if "-NodTAG" in r["condition"]]
        ordered.extend(dtag)             # lower y
        ordered.extend(nodtag)           # higher y
    return ordered


with C.logging(snakemake) as say:
    pl.style()
    import matplotlib.pyplot as plt

    cohort = spec.get("cohort", "WT")
    say("panel {} (cohort {})".format(spec["panel_id"], cohort))

    summary = []
    for path in map(str, snakemake.input):
        track = os.path.basename(path)[: -len(".tsv")]
        table = C.read_table(path)
        category = C.need(table, "category", path)
        counts = category.value_counts()
        summary.append({
            "track": track,
            "clone": C.clone_of(track),
            "condition": C.condition_of(track),
            "label": C.track_label(track),
            "no_gene": int(counts.get("no-gene", 0)),
            "no_expressed": int(counts.get("no-expressed-gene-valley", 0)),
            "silent": int(counts.get("silent-gene-valley", 0)),
            "escaping": int(counts.get("escaping-gene-valley", 0)),
            "total": int(len(table)),
        })
        say("  {:28s} total={:4d}  no-gene={:3d}  no-expressed={:3d}  "
            "silent={:3d}  escaping={:3d}".format(
                summary[-1]["label"], summary[-1]["total"],
                summary[-1]["no_gene"], summary[-1]["no_expressed"],
                summary[-1]["silent"], summary[-1]["escaping"]))

    if not summary:
        pl.empty_panel(out_path, "no gene-content tables for cohort " + cohort)
        C.write_n_items(snakemake, 0)
        raise SystemExit(0)

    rows = _rows_wt(summary) if cohort == "WT" else _rows_degron(summary)

    y = np.arange(len(rows))
    g = np.array([r["no_gene"] for r in rows], dtype=float)
    d = np.array([r["no_expressed"] for r in rows], dtype=float)
    s = np.array([r["silent"] for r in rows], dtype=float)
    e = np.array([r["escaping"] for r in rows], dtype=float)
    totals = g + d + s + e

    fig_h = max(3.0, 0.7 * len(rows) + 1.5)
    fig, ax = plt.subplots(figsize=(7, fig_h))

    ax.barh(y, g, height=0.6, label="No gene (intergenic)",
            color=pl.CAT_COLORS["no-gene"])
    ax.barh(y, d, left=g, height=0.6, label="No expressed gene",
            color=pl.CAT_COLORS["no-expressed-gene-valley"])
    ax.barh(y, s, left=g + d, height=0.6, label="Silent gene(s)",
            color=pl.CAT_COLORS["silent-gene-valley"])
    ax.barh(y, e, left=g + d + s, height=0.6, label="Escaping gene(s)",
            color=pl.CAT_COLORS["escaping-gene-valley"])

    for yi, total in zip(y, totals):
        ax.text(total + 3, yi, str(int(total)), va="center", ha="left",
                fontsize=8)

    ax.set_yticks(y)
    ax.set_yticklabels([r["label"] for r in rows], fontsize=9)
    ax.set_xlabel("Number of valleys")
    ax.set_xlim(0, (totals.max() if len(totals) else 1) * 1.15)
    ax.set_title(
        "Gene content of H3K27me3 valleys - {} (chrX)\n"
        "gene-filtered valley set".format(
            "WT clones" if cohort == "WT" else "degron clones +/-dTAG"),
        fontsize=10,
    )
    pl.despine(ax)
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0,
              frameon=False, fontsize=8)
    fig.tight_layout()

    pl.save(fig, out_path, dpi=spec.get("dpi", 300))
    C.write_n_items(snakemake, len(rows),
                    extra={r["label"]: r["total"] for r in rows})
    say("wrote {}".format(out_path))
