"""Fig 3c — delta-compartment against delta-loop strength. MEDIUM CONFIDENCE.

`02_02_scatter_analysis_metaloci.ipynb` writes about ten scatters into
`results/final_analyses/` and nothing marks which one became Fig 3c. The
archaeology short-listed four; all four are rendered as
`scatter_delta_loop_comp__{candidate}.svg` and a human picks against
`paper_figures/20260722_Figure_3.pdf`.

    loops     Loop strength Xi vs Loop strength Xa
    saddle    Saddle strength Xi vs Saddle strength Xa
    xi_xa     Loop strength (Xi - Xa) vs saddle strength (Xi - Xa)
    ar_loops  METALoci compartmentalization vs loop strength delta

Recipe, verbatim: seaborn scatterplot, `hue` = mean allelic ratio,
`palette="Greens"`, `s=100`, `edgecolor="black"`, a text label per point,
a manually built colorbar, `set_aspect("equal")`, grid off.

The colour normalisation is NOT fixture F-8 here: the original called
`plt.Normalize(total_stats['Allelic ratio'].min(), max())` for these scatters
and the hard-coded pair only for the Fig 3a/b circles. Reproduced as-is --
the two figures genuinely use different scales.
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

AXES = {
    "loops": ("Loop strength Xi", "Loop strength Xa"),
    "saddle": ("Saddle strength Xi", "Saddle strength Xa"),
    "xi_xa": ("Loop strength Xi - Xa", "Saddle strength Xi - Xa"),
    "ar_loops": ("Compartmentalization", "Loop strength Xi - Xa"),
}


def _pair_key(name):
    """`Mecp2_E6_WT_G1_Xi` -> (`Mecp2_E6_WT`, `Xi`)."""
    key = re.sub(r"_(G1|G2)_(Xa|Xi)$", "", str(name))
    if str(name).endswith("Xa"):
        return key, "Xa"
    if str(name).endswith("Xi"):
        return key, "Xi"
    return key, None


with C.logging(snakemake) as say:
    pl.style()
    import matplotlib.pyplot as plt

    candidate = spec.get("candidate") or "loops"
    roi = spec.get("roi", "full")
    say("panel {} (candidate {})".format(spec["panel_id"], candidate))

    paths = [str(p) for p in snakemake.input]
    stats = {}

    def _absorb_scores(path):
        """Pile-up scores; for Xa_vs_Xi keep the loops-called-on-Xa set."""
        table = C.read_table(path)
        cooler_col = "cooler" if "cooler" in table.columns else "file"
        for _, row in table.iterrows():
            name = str(row["file"])
            if name.startswith("Xi_"):
                continue
            key, allele = _pair_key(row[cooler_col])
            if allele is None:
                continue
            stats.setdefault(key, {})["Loop strength {}".format(allele)] = \
                float(row["mean"])

    def _absorb(path, column, target):
        table = C.read_table(path)
        cols = list(table.columns)
        key_col = "file" if "file" in cols else cols[0]
        val_col = column if column in cols else cols[1]
        for _, row in table.iterrows():
            key, allele = _pair_key(row[key_col])
            if allele is None:
                continue
            stats.setdefault(key, {})["{} {}".format(target, allele)] = \
                float(row[val_col])

    for path in paths:
        base = os.path.basename(path)
        if base == "scores.tsv":
            # Pile-up scores, NOT loop_stats: 01_03's writes of loop_stats*.tsv
            # are commented out, and the files Fig 3c was made from came from
            # 01_04 and hold pile-up scores under that name.
            _absorb_scores(path)
        elif "saddle_strength_selected" in base:
            _absorb(path, "value", "Saddle strength")

    comp_path = next((p for p in paths
                      if "compartmentalization" in os.path.basename(p)), None)
    if comp_path and os.path.exists(comp_path):
        comp = C.read_table(comp_path)
        if "dataset" in comp.columns and "compartmentalization" in comp.columns:
            for key, group in comp.groupby(
                    comp["dataset"].map(lambda d: _pair_key(d)[0])):
                stats.setdefault(key, {})["Compartmentalization"] = \
                    float(group["compartmentalization"].mean())

    for path in paths:
        if not path.endswith(".csv"):
            continue
        table = pd.read_csv(path, index_col=0)
        if "mean_all" not in table.index:
            continue
        for clone, value in table.loc["mean_all"].items():
            value = pd.to_numeric(value, errors="coerce")
            if not np.isfinite(value):
                continue
            for key in stats:
                if "_{}_".format(str(clone).split(".")[0]) in "_{}_".format(key):
                    stats[key].setdefault("Allelic ratio", float(value))

    table = pd.DataFrame(stats).T
    for base in ("Loop strength", "Saddle strength"):
        if "{} Xi".format(base) in table and "{} Xa".format(base) in table:
            table["{} Xi - Xa".format(base)] = (
                table["{} Xi".format(base)] - table["{} Xa".format(base)])
    say("  {} clone-experiments, columns: {}".format(
        len(table), ", ".join(map(str, table.columns))))

    xcol, ycol = AXES.get(candidate, AXES["loops"])
    if xcol not in table.columns or ycol not in table.columns:
        pl.empty_panel(out_path, "candidate {} needs columns {!r} and {!r};\n"
                                 "available: {}".format(
                                     candidate, xcol, ycol,
                                     ", ".join(map(str, table.columns))))
        C.write_n_items(snakemake, 0)
        raise SystemExit(0)

    data = table[[c for c in {xcol, ycol, "Allelic ratio"} if c in table]].dropna(
        subset=[xcol, ycol])
    data["Label"] = data.index
    say("  plotting {} points".format(len(data)))

    fig, ax = plt.subplots(figsize=(8, 6))
    colours = data.get("Allelic ratio")
    if colours is not None and colours.notna().any():
        norm = plt.Normalize(float(colours.min()), float(colours.max()))
        sc = ax.scatter(data[xcol], data[ycol], c=colours, cmap="Greens",
                        s=100, edgecolor="black", norm=norm)
        sm = plt.cm.ScalarMappable(cmap="Greens", norm=norm)
        sm.set_array([])
        fig.colorbar(sm, ax=ax, label="Mean allelic ratio")
    else:
        ax.scatter(data[xcol], data[ycol], s=100, edgecolor="black",
                   color="#9ecae1")

    for _, row in data.iterrows():
        ax.text(row[xcol], row[ycol], str(row["Label"]), fontsize=7,
                ha="center", va="bottom")

    ax.set_xlabel(xcol)
    ax.set_ylabel(ycol)
    ax.grid(False)
    if candidate in ("loops", "saddle"):
        ax.set_aspect("equal", adjustable="box")
    if candidate == "xi_xa":
        ax.axhline(0, color="gray", ls="--", lw=0.8)
        ax.axvline(0, color="gray", ls="--", lw=0.8)
    ax.set_title("Fig 3c candidate: {}   (ROI {})".format(candidate, roi),
                 fontsize=10)
    fig.tight_layout()

    pl.save(fig, out_path, dpi=spec.get("dpi", 300))
    C.write_n_items(snakemake, len(data))
    say("wrote {}".format(out_path))
