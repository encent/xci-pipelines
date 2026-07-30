"""METALoci violins, heatmaps and scatters from the compartmentalization TSVs.

    panel_metaloci_violin   Xa vs Xi compartment-like strength, per signal
    panel_metaloci_heatmap  the same table as a clone x signal heatmap
    panel_metaloci_scatter  compartmentalization against loop/saddle strength

`panel_metaloci_violin` is not merely an extra: its PDFs are column 0 of the
Fig 6 grid, so `panel_metaloci_composite` depends on them. Only H3K27ac and
H3K27me3 reached the paper as violins; the other four signals are produced
because the original produced them.

The compartmentalization value is
`(sq1 + sq3) / (q1 + q2 + q3 + q4) * 100` -- the share of bins in the
significant LMI quadrants 1 and 3.
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


def _load(paths):
    frames = []
    for path in paths:
        if "compartmentalization" not in os.path.basename(path):
            continue
        try:
            table = C.read_table(path)
        except Exception:
            continue
        if {"dataset", "signal", "compartmentalization"} <= set(table.columns):
            frames.append(table)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def _allele(dataset):
    if str(dataset).endswith("Xa") or "_Xa" in str(dataset):
        return "Xa"
    if str(dataset).endswith("Xi") or "_Xi" in str(dataset):
        return "Xi"
    return None


def _pair_key(dataset):
    return re.sub(r"_(G1|G2)?_?(Xa|Xi)(_rep\d)?$", "", str(dataset))


with C.logging(snakemake) as say:
    pl.style()
    import matplotlib.pyplot as plt

    producer = spec.get("producer", "")
    locus = spec.get("locus", "")
    signal = spec.get("signal", "")
    say("panel {} ({}) locus={} signal={}".format(
        spec["panel_id"], producer, locus, signal))

    table = _load([str(p) for p in snakemake.input])
    if table.empty:
        pl.empty_panel(out_path, "no compartmentalization tables")
        C.write_n_items(snakemake, 0)
        raise SystemExit(0)

    if locus:
        table = table[table["dataset"].astype(str).str.startswith(locus + "_")]
    if signal and producer != "panel_metaloci_scatter":
        table = table[table["signal"] == signal]
    if table.empty:
        pl.empty_panel(out_path, "no rows for locus={} signal={}".format(
            locus, signal))
        C.write_n_items(snakemake, 0)
        raise SystemExit(0)

    table = table.copy()
    table["allele"] = table["dataset"].map(_allele)
    table["pair"] = table["dataset"].map(_pair_key)

    if producer == "panel_metaloci_heatmap":
        wide = table.pivot_table(index="dataset", columns="signal",
                                 values="compartmentalization")
        fig, ax = plt.subplots(
            figsize=(1.2 * max(3, wide.shape[1]) + 3,
                     0.28 * max(6, wide.shape[0]) + 2))
        im = ax.imshow(wide.values, cmap="viridis", aspect="auto")
        ax.set_xticks(range(wide.shape[1]))
        ax.set_xticklabels(wide.columns, rotation=90, fontsize=7)
        ax.set_yticks(range(wide.shape[0]))
        ax.set_yticklabels(wide.index, fontsize=6)
        for i in range(wide.shape[0]):
            for j in range(wide.shape[1]):
                value = wide.values[i, j]
                if np.isfinite(value):
                    ax.text(j, i, "{:.1f}".format(value), ha="center",
                            va="center", fontsize=5, color="white")
        fig.colorbar(im, ax=ax, label="compartmentalization (%)")
        ax.set_title("METALoci compartmentalization - {} {}".format(
            C.display_locus(locus, spec.get("display_names")), signal),
            fontsize=10)
        n = int(wide.shape[0])

    elif producer == "panel_metaloci_scatter":
        wide = table.pivot_table(index="dataset", columns="signal",
                                 values="compartmentalization")
        signals = list(wide.columns)
        if len(signals) < 2:
            pl.empty_panel(out_path, "need >= 2 signals to scatter")
            C.write_n_items(snakemake, 0)
            raise SystemExit(0)
        xsig, ysig = signals[0], signals[1]
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.scatter(wide[xsig], wide[ysig], s=70, edgecolor="black",
                   color="#4477aa")
        for name, row in wide.iterrows():
            ax.text(row[xsig], row[ysig], str(name), fontsize=5, ha="center",
                    va="bottom")
        ax.set_xlabel("compartmentalization, {} (%)".format(xsig))
        ax.set_ylabel("compartmentalization, {} (%)".format(ysig))
        ax.set_title("METALoci compartmentalization - {}".format(
            C.display_locus(locus, spec.get("display_names"))), fontsize=10)
        pl.despine(ax)
        n = int(wide.shape[0])

    else:
        wide = table.pivot_table(index="pair", columns="allele",
                                 values="compartmentalization")
        wide = wide.dropna(subset=[c for c in ("Xa", "Xi") if c in wide])
        if not {"Xa", "Xi"} <= set(wide.columns) or len(wide) < 2:
            pl.empty_panel(
                out_path,
                "only {} paired Xa/Xi dataset(s) for {} {}".format(
                    len(wide), locus, signal))
            C.write_n_items(snakemake, len(wide))
            raise SystemExit(0)
        from scipy.stats import wilcoxon
        p_two = wilcoxon(wide["Xa"], wide["Xi"],
                         alternative="two-sided").pvalue
        fig, ax = plt.subplots(figsize=(6, 6))
        top = float(wide.max().max()) * 1.2
        pl.paired_violin(ax, wide, "Xa", "Xi",
                         "compartmentalization (%)", highlight="Xa",
                         ylim=(0, max(top, 1.0)))
        ax.set_xlabel("Chromosome")
        ax.set_title("{} {}  N={}\nWilcoxon two-sided p = {:.5f}".format(
            C.display_locus(locus, spec.get("display_names")), signal,
            len(wide), p_two), fontsize=9)
        say("  N = {}  Wilcoxon p = {:.6g}".format(len(wide), p_two))
        n = int(len(wide))

    fig.tight_layout()
    pl.save(fig, out_path, dpi=spec.get("dpi", 300))
    C.write_n_items(snakemake, n)
    say("wrote {}".format(out_path))
