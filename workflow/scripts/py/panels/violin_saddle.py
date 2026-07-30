"""EFig 3b — paired Xa-vs-Xi compartment (saddle) strength violin.

Same recipe as `violin_loops.py`; the quantity is the saddle corner score
(AA + BB) / (AB + BA) of the SELECTED eigenvector, which
`saddle_strength_selected` picks per (roi, cooler) using the
`compartments*.json` fixtures.

The input is one two-column table (`file`, `value`) covering every cooler, so
the pairing has to be done here: strip the allele token off each file name and
pair Xa with Xi (or NodTAG with dTAG). The dropped rows are logged -- an
unpaired cooler is usually a sheet problem, not a plotting problem.

Read this panel together with `eigenvector_orientation.tsv`: the saddle score
is invariant under a global sign flip only if the flip is global. A per-cooler
flip silently mirrors individual saddles and moves this violin.
"""

import os
import re
import sys

sys.path.insert(0, os.path.join(snakemake.params.spec["libdir"],
                                "scripts", "py", "panels"))
sys.path.insert(0, snakemake.params.spec["libdir"])

import pandas as pd                                            # noqa: E402
import _common as C                                            # noqa: E402
from lib import plotting as pl                                 # noqa: E402

spec = C.bootstrap(snakemake)
out_path = C.outputs(snakemake)


def _pair_key(name, comparison):
    """Strip the token that distinguishes the two members of a pair."""
    if comparison == "Xa_vs_Xi":
        key = re.sub(r"_(G1|G2)_(Xa|Xi)", "", str(name))
        group = "Xa" if str(name).endswith("Xa") else (
            "Xi" if str(name).endswith("Xi") else None)
    else:
        key = re.sub(r"-(No)?dTAG", "", str(name))
        group = "dTAG" if "-dTAG" in str(name) else (
            "NodTAG" if "-NodTAG" in str(name) else None)
    return key, group


with C.logging(snakemake) as say:
    pl.style()
    import matplotlib.pyplot as plt
    from scipy.stats import wilcoxon

    comparison = spec.get("comparison", "Xa_vs_Xi")
    roi = spec.get("roi", "full")
    mark = spec.get("mark")
    candidate = spec.get("candidate")
    say("panel {} ({} / {})".format(spec["panel_id"], comparison, roi))

    table = C.read_table(str(snakemake.input[0]))
    files = C.need(table, "file", str(snakemake.input[0]))
    values = C.need(table, "value", str(snakemake.input[0]),
                    alternatives=("saddle_strength", "score", "E1"))

    rows = {}
    for name, value in zip(files, values):
        key, group = _pair_key(name, comparison)
        if group is None:
            say("  skipping {} (no {} group token)".format(name, comparison))
            continue
        rows.setdefault(key, {})[group] = float(value)

    g1, g2 = ("Xa", "Xi") if comparison == "Xa_vs_Xi" else ("NodTAG", "dTAG")
    paired = {k: v for k, v in rows.items() if g1 in v and g2 in v}
    dropped = sorted(set(rows) - set(paired))
    for key in dropped:
        say("  unpaired, dropped: {} (has {})".format(
            key, ", ".join(sorted(rows[key]))))

    if mark:
        paired = {k: v for k, v in paired.items()
                  if mark.split("-")[0] in k or mark in k} or paired

    wide = pd.DataFrame(paired).T[[g1, g2]].dropna()
    n = len(wide)
    say("  N = {} paired coolers".format(n))
    for idx in wide.index:
        say("    {:38s} {} {:.4f}   {} {:.4f}".format(
            str(idx), g1, wide.loc[idx, g1], g2, wide.loc[idx, g2]))

    if n < 2:
        pl.empty_panel(out_path,
                       "only {} paired saddle score(s) for {}".format(n, comparison))
        C.write_n_items(snakemake, n)
        raise SystemExit(0)

    p_two = wilcoxon(wide[g1], wide[g2], alternative="two-sided").pvalue
    p_gt = wilcoxon(wide[g1], wide[g2], alternative="greater").pvalue
    p_lt = wilcoxon(wide[g1], wide[g2], alternative="less").pvalue
    say("  Wilcoxon two-sided p = {:.6g}   {}>{} p = {:.6g}   {}<{} p = {:.6g}"
        .format(p_two, g1, g2, p_gt, g1, g2, p_lt))

    value_label = "Compartment (saddle) strength"

    if candidate == "heatmap":
        fig, ax = plt.subplots(figsize=(12, 4))
        data = wide.T
        im = ax.imshow(data.values, cmap="coolwarm", aspect="auto")
        ax.set_xticks(range(data.shape[1]))
        ax.set_xticklabels(data.columns, rotation=90, fontsize=6)
        ax.set_yticks(range(data.shape[0]))
        ax.set_yticklabels(data.index)
        fig.colorbar(im, ax=ax, label=value_label)
        ax.set_title("{} -- {} / {} / N={}".format(
            value_label, comparison, roi, n), fontsize=10)
    else:
        fig, ax = plt.subplots(figsize=(8, 6))
        top = float(max(wide.max())) * 1.15
        pl.paired_violin(ax, wide, g1, g2, value_label,
                         highlight=g1 if comparison == "Xa_vs_Xi" else "NodTAG",
                         ylim=(0, max(top, 1.0)))
        ax.set_xlabel("Chromosome" if comparison == "Xa_vs_Xi" else "Genotype")
        ax.set_title(
            "{} -- comparison: {} locus type: {}  N={}\n"
            "Wilcoxon p-values: two-sided={:.5f}\n{}>{}={:.5f}\n{}<{}={:.5f}"
            .format(value_label, comparison, roi, n,
                    p_two, g1, g2, p_gt, g1, g2, p_lt),
            fontsize=9)

    fig.tight_layout()
    pl.save(fig, out_path, dpi=spec.get("dpi", 300))
    C.write_n_items(snakemake, n, extra={
        "wilcoxon_two_sided": p_two, "wilcoxon_greater": p_gt,
        "wilcoxon_less": p_lt, "n_unpaired_dropped": len(dropped),
    })
    say("wrote {}".format(out_path))
