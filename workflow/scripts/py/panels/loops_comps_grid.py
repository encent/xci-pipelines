"""Fig 3a/3b, Fig 5b/5c/5g/5h, EFig 8c/8d, EFig 10c/10d.

The workhorse. Ten of the twenty-five target panels are the same figure with
different inputs: a row of loop pile-ups and/or a row of saddle heatmaps, one
column per clone-allele, with each column's loop or compartment strength
printed above it and each column's mean RNA-seq allelic ratio shown as a green
circle.

    kind=loops  pile-up row + the loop-strength numbers
    kind=comps  saddle row  + the saddle-strength numbers
    kind=both   both rows          (Fig 3a / 3b)

3a is Kdm5c (`Jarid` in the code) and 3b is Mecp2. `figure_legends.txt` has
these the other way round; it uses an older lettering scheme and `figure.pdf`
(2026-07-25) is authoritative.

THE GREEN CIRCLES
-----------------
The published legend calls them a METALoci D-score. They are not: they are the
mean RNA-seq allelic ratio from `01_07_get_allelic_ratio_gtf.py`. The colour
scale is fixture F-8 -- a truncated `Greens` under
`Normalize(0.0424769728421052, 0.40909353528)`. Those bounds are the WT range,
hard-coded in the original; recomputing them from the data would rescale a
degron-only run and quietly break comparability between figures.

PILE-UP RENDERING
-----------------
Verbatim from `01_04`: `cmap='coolwarm', scale='log', sym=True, vmax=2,
vmin=None, height=4, plot_ticks=True, center=3`. `scale='log'` + `sym=True`
means a symmetric log colour scale about 1.0, i.e. `LogNorm(1/2, 2)`. The score
is the mean of the central 3x3 (`center=3`), which `loop_pileup` already
computed and stored beside the matrix.
"""

import os
import re
import sys

sys.path.insert(0, os.path.join(snakemake.params.spec["libdir"],
                                "scripts", "py", "panels"))
sys.path.insert(0, snakemake.params.spec["libdir"])

import numpy as np                                             # noqa: E402
import _common as C                                            # noqa: E402
from lib import plotting as pl                                 # noqa: E402

spec = C.bootstrap(snakemake)
out_path = C.outputs(snakemake)


def _clone_of_cooler(name):
    parts = str(name).split("_")
    return parts[1] if len(parts) > 1 else str(name)


def _allele_of(name):
    return "Xa" if str(name).endswith("Xa") else (
        "Xi" if str(name).endswith("Xi") else "")


def _order_key(name, clones):
    clone = _clone_of_cooler(name)
    rank = clones.index(clone) if clone in clones else len(clones)
    return (rank, _allele_of(name), name)


def _load_allelic_ratio(paths, locus):
    """Mean allelic ratio per clone, from the AR_stats CSV (index x clone)."""
    import pandas as pd
    for path in paths:
        if not path.endswith(".csv") or locus not in os.path.basename(path):
            continue
        table = pd.read_csv(path, index_col=0)
        for key in ("mean_all", "mean_non_zero", "mean_>0.1"):
            if key in table.index:
                return {str(k): float(v)
                        for k, v in table.loc[key].items()
                        if np.isfinite(pd.to_numeric(v, errors="coerce"))}
    return {}


def _ar_for(clone, ar_map):
    """AR tables are keyed on the RNA-seq column names (`CL30.7`, `B1`)."""
    if clone in ar_map:
        return ar_map[clone]
    for key, value in ar_map.items():
        if str(key).split(".")[0] == clone:
            return value
    return None


with C.logging(snakemake) as say:
    pl.style()
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    locus = spec.get("locus", "Mecp2")
    roi = spec.get("roi", "full")
    comparison = spec.get("comparison", "Xa_vs_Xi")
    kind = spec.get("kind", "both")
    clones = list(spec.get("clones") or [])
    display = C.display_locus(locus, spec.get("display_names"))
    say("panel {}: {} / {} / {} / kind={}".format(
        spec["panel_id"], display, roi, comparison, kind))

    paths = [str(p) for p in snakemake.input]
    pileups = sorted(p for p in paths
                     if p.endswith(".npz") and os.sep + "pileups" + os.sep in p)
    saddles = sorted(p for p in paths
                     if p.endswith(".npz") and os.sep + "saddle" + os.sep in p)
    stats_tsvs = [p for p in paths if p.endswith(".tsv")]
    ar_map = _load_allelic_ratio(paths, locus)

    # --- the numbers printed above each column ---------------------------
    loop_strength, saddle_strength = {}, {}
    for path in stats_tsvs:
        table = C.read_table(path)
        base = os.path.basename(path)
        if base == "scores.tsv" or base.startswith("loop_stats"):
            # `scores.tsv` is the pile-up score table -- the quantity the
            # original's on-disk `loop_stats*.tsv` actually held. Key on both
            # the sample id and the cooler name so the `Xa_`/`Xi_` called-on
            # prefix does not have to be guessed at the lookup site.
            for _, row in table.iterrows():
                value = float(row["mean"])
                loop_strength[str(row["file"])] = value
                if "cooler" in table.columns:
                    loop_strength.setdefault(str(row["cooler"]), value)
        elif "saddle_strength_selected" in base:
            for _, row in table.iterrows():
                key = str(row[table.columns[0]])
                saddle_strength[key] = float(row[table.columns[1]])

    rows = []
    if kind in ("loops", "both") and pileups:
        rows.append("loops")
    if kind in ("comps", "both") and saddles:
        rows.append("comps")
    if not rows:
        pl.empty_panel(out_path, "no pile-up or saddle matrices for {} {} {}"
                       .format(display, roi, comparison))
        C.write_n_items(snakemake, 0)
        raise SystemExit(0)

    def _sample_of_pileup(path):
        return re.sub(r"^pileup_", "", os.path.basename(path)[: -len(".npz")])

    def _name_of_saddle(path):
        return re.sub(r"^Saddle_E\d+_", "",
                      os.path.basename(path)[: -len(".npz")])

    columns = {}
    if "loops" in rows:
        for path in pileups:
            columns.setdefault(_sample_of_pileup(path), {})["pileup"] = path
    if "comps" in rows:
        for path in saddles:
            columns.setdefault(_name_of_saddle(path), {})["saddle"] = path

    names = sorted(columns, key=lambda n: _order_key(n, clones))
    if clones:
        keep = [n for n in names if _clone_of_cooler(n) in clones]
        names = keep or names
    say("  {} columns: {}".format(len(names), ", ".join(names)))

    n_col = max(1, len(names))
    fig, axes = plt.subplots(len(rows), n_col,
                             figsize=(2.6 * n_col, 3.4 * len(rows)),
                             squeeze=False)
    cmap = pl.allelic_ratio_cmap()
    norm = pl.allelic_ratio_norm()

    for ri, row_kind in enumerate(rows):
        for ci, name in enumerate(names):
            ax = axes[ri][ci]
            path = columns[name].get(
                "pileup" if row_kind == "loops" else "saddle")
            if path is None:
                ax.axis("off")
                continue
            npz = C.read_npz(path)

            if row_kind == "loops":
                data = np.asarray(npz["data"], dtype=float)
                # scale='log', sym=True, vmax=2 -> LogNorm(1/2, 2)
                im = ax.imshow(data, cmap="coolwarm",
                               norm=LogNorm(vmin=1.0 / 2, vmax=2),
                               interpolation="none")
                score = loop_strength.get(name)
                if score is None:
                    centre = data.shape[0] // 2
                    score = float(np.nanmean(
                        data[centre - 1:centre + 2, centre - 1:centre + 2]))
                label = "loop {:.2f}".format(score)
            else:
                data = np.asarray(
                    npz["saddle"] if "saddle" in npz else
                    npz["interaction_sum"] / npz["interaction_count"],
                    dtype=float)
                im = ax.imshow(data, cmap="coolwarm",
                               norm=LogNorm(vmin=0.5, vmax=2.0),
                               interpolation="none")
                score = saddle_strength.get(name)
                if score is None and "score" in npz:
                    score = float(npz["score"])
                label = "saddle {:.2f}".format(score) if score is not None \
                    else "saddle n/a"

            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title("{}\n{}".format(name.replace("_", " "), label),
                         fontsize=6)
            say("    {:34s} {:10s} {}".format(name, row_kind, label))

            if ri == 0:
                ar = _ar_for(_clone_of_cooler(name), ar_map)
                if ar is not None:
                    ax.scatter([data.shape[1] * 0.88], [data.shape[0] * 0.12],
                               s=90, c=[cmap(norm(ar))], edgecolor="black",
                               linewidth=0.5, zorder=5, clip_on=False)

        axes[ri][0].set_ylabel(
            "loop pile-up" if row_kind == "loops" else "saddle", fontsize=8)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes.ravel().tolist(), fraction=0.02, pad=0.02)
    cbar.set_label("Mean RNA-seq allelic ratio\n(the figure legend calls this a "
                   "METALoci D-score)", fontsize=6)

    fig.suptitle("{} - {} ({}, {})".format(
        display, "loops + compartments" if kind == "both" else kind,
        comparison, roi), fontsize=11)

    pl.save(fig, out_path, dpi=spec.get("dpi", 300))
    C.write_n_items(snakemake, len(names), extra=dict(
        list(loop_strength.items())[:20] + list(saddle_strength.items())[:20]))
    say("wrote {}".format(out_path))
