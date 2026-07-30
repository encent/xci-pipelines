"""EFig 3a — paired Xa-vs-Xi loop-strength violin, N = 15 clone-alleles.

Nine Mecp2 plus six Jarid clone-alleles. One dot pair per clone joined by a
grey line, because the test is paired: Wilcoxon signed-rank, two-sided plus
both one-sided alternatives, all three printed in the title exactly as
`02_00_analysis_loop_chromosight.ipynb` did.

WHERE THE NUMBERS COME FROM, AND WHY NOT FROM `loop_stats*.tsv`
--------------------------------------------------------------
Plan section 2.2 assigns `loop_stats{,_Xa,_Xi}.tsv` to `loop_strength`, the
port of `01_03`. But `01_03`'s writes of that filename are **commented out**
(lines 210, 319-320): the `loop_stats*.tsv` files that exist on disk, and that
EFig 3a was made from, were written by `01_04` and hold **pile-up scores** — a
different quantity carrying the same name.

So this panel reads `P.pileup_score(comparison, roi)`. `P.loop_stats` is still
produced and still compared against ground truth; it is simply not what this
figure plots. Reading the wrong one gives a violin of plausible numbers.

Statistics are computed HERE, with scipy, and the same numbers go into
`paper_numbers.tsv`. The figure and the table therefore cannot disagree. This
is also why we do not reproduce the original's `annotate_violin_mwu`, which
announced Mann-Whitney in the title and called `ttest_rel` (deviation D-5):
`config.legacy.stackup_violin_test` exists for the stackup violins that had
that bug; the loop violins always used Wilcoxon and always will.

Recipe, verbatim: figsize (8, 6), `inner="quart"`, `color="0.9"`, dots at
alpha 0.2 with the "called on" allele in red, ylim (0, 4.0).

A companion heatmap of the same table is emitted as `{id}__heatmap`.
"""

import os
import re
import sys

sys.path.insert(0, os.path.join(snakemake.params.spec["libdir"],
                                "scripts", "py", "panels"))
sys.path.insert(0, snakemake.params.spec["libdir"])

import numpy as np                                             # noqa: E402
import pandas as pd                                            # noqa: E402
import _common as C                                            # noqa: E402
from lib import plotting as pl                                 # noqa: E402

spec = C.bootstrap(snakemake)
out_path = C.outputs(snakemake)


def _paired_table(paths, comparison, stat, called_on):
    """Wide table: one row per clone-experiment, columns = the two groups.

    `scores.tsv` holds every pile-up of the comparison in one file. For
    `Xa_vs_Xi` the `file` column carries the allele the loops were CALLED on as
    a prefix (`Xa_Jarid_B1_WT_G1_Xa`), reproducing the original's
    `pileup_{Xa|Xi}_{cooler}`; the published EFig 3a is the `called_on = Xa`
    set. Within that set, each cooler is then paired with its partner allele.
    """
    scores = C.read_table(paths[0])
    path = paths[0]
    files = C.need(scores, "file", path).astype(str)
    values = C.need(scores, stat, path, alternatives=("mean", "median"))
    cooler = (scores["cooler"].astype(str) if "cooler" in scores.columns
              else files)

    rows = {}
    for name, cool, value in zip(files, cooler, values):
        if comparison == "Xa_vs_Xi":
            if not name.startswith(called_on + "_"):
                continue
            key = re.sub(r"_(G1|G2)_(Xa|Xi)(_rep\d)?$", "", str(cool))
            group = "Xa" if str(cool).endswith("Xa") else (
                "Xi" if str(cool).endswith("Xi") else None)
        else:
            key = re.sub(r"-(No)?dTAG", "", str(cool))
            group = "dTAG" if "-dTAG" in str(cool) else (
                "NodTAG" if "-NodTAG" in str(cool) else None)
        if group is None:
            continue
        rows.setdefault(key, {})[group] = float(value)

    g1, g2 = ("Xa", "Xi") if comparison == "Xa_vs_Xi" else ("NodTAG", "dTAG")
    paired = {k: v for k, v in rows.items() if g1 in v and g2 in v}
    wide = pd.DataFrame(paired).T
    wide = wide[[g1, g2]] if not wide.empty else wide
    var_name = "Chromosome" if comparison == "Xa_vs_Xi" else "Genotype"
    return wide, g1, g2, var_name


with C.logging(snakemake) as say:
    pl.style()
    import matplotlib.pyplot as plt
    from scipy.stats import wilcoxon

    comparison = spec.get("comparison", "Xa_vs_Xi")
    roi = spec.get("roi", "full")
    stat = spec.get("stat", "mean")
    called_on = spec.get("called_on", "Xa")
    mark = spec.get("mark")
    candidate = spec.get("candidate")
    paths = sorted(map(str, snakemake.input))
    say("panel {} ({} / {} / {})".format(
        spec["panel_id"], comparison, roi, stat))

    table, g1, g2, var_name = _paired_table(paths, comparison, stat, called_on)
    if mark:
        keep = [i for i in table.index if mark.split("-")[0] in str(i)
                or mark in str(i)]
        if keep:
            table = table.loc[keep]
    table = table.dropna() if not table.empty else table
    n = len(table)
    say("  N = {} clone-alleles".format(n))
    for idx in table.index:
        say("    {:38s} {} {:.4f}   {} {:.4f}".format(
            str(idx), g1, table.loc[idx, g1], g2, table.loc[idx, g2]))

    if n < 2:
        pl.empty_panel(out_path,
                       "only {} paired observation(s) for {}".format(n, comparison))
        C.write_n_items(snakemake, n)
        raise SystemExit(0)

    p_two = wilcoxon(table[g1], table[g2], alternative="two-sided").pvalue
    p_gt = wilcoxon(table[g1], table[g2], alternative="greater").pvalue
    p_lt = wilcoxon(table[g1], table[g2], alternative="less").pvalue
    say("  Wilcoxon two-sided p = {:.6g}   {}>{} p = {:.6g}   {}<{} p = {:.6g}"
        .format(p_two, g1, g2, p_gt, g1, g2, p_lt))

    value_label = "{} Loop Strength".format(stat)

    if candidate == "heatmap":
        fig, ax = plt.subplots(figsize=(12, 4))
        data = table[[g1, g2]].T
        im = ax.imshow(data.values, cmap="coolwarm", aspect="auto")
        ax.set_xticks(range(data.shape[1]))
        ax.set_xticklabels(data.columns, rotation=90, fontsize=6)
        ax.set_yticks(range(data.shape[0]))
        ax.set_yticklabels(data.index)
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                ax.text(j, i, "{:.2f}".format(data.values[i, j]),
                        ha="center", va="center", fontsize=5)
        fig.colorbar(im, ax=ax, label=value_label)
        ax.set_title("{} -- {} / {} / N={} clone-alleles".format(
            value_label, comparison, roi, n), fontsize=10)
    else:
        fig, ax = plt.subplots(figsize=(8, 6))
        pl.paired_violin(ax, table, g1, g2, value_label,
                         highlight=called_on if comparison == "Xa_vs_Xi"
                         else "NodTAG",
                         ylim=(0, 4.0))
        ax.set_xlabel(var_name)
        ax.set_title(
            "{} -- comparison: {} locus type: {}\n"
            "loops called on: {}  N={} clone-alleles\n"
            "Wilcoxon p-values: two-sided={:.5f}\n{}>{}={:.5f}\n{}<{}={:.5f}"
            .format(value_label, comparison, roi, called_on, n,
                    p_two, g1, g2, p_gt, g1, g2, p_lt),
            fontsize=9)

    fig.tight_layout()
    pl.save(fig, out_path, dpi=spec.get("dpi", 300))
    C.write_n_items(snakemake, n, extra={
        "wilcoxon_two_sided": p_two, "wilcoxon_greater": p_gt,
        "wilcoxon_less": p_lt, "group1": g1, "group2": g2,
    })
    say("wrote {}".format(out_path))
