"""Four Hi-C panels, dispatched on the producing rule.

    panel_pileup       one coolpuppy pile-up
    panel_saddle       saddle heatmap + the (AA+BB)/(AB+BA) strength profile
    panel_eigenvector  E1/E2/E3 on one axis
    panel_loop_overlay loops before and after the manual refinement
    panel_cooler_qc    the balanced matrix at the locus view window

None of these is a published panel on its own; together they are how a
published one gets debugged. In particular:

* **panel_eigenvector.** `compartments*.json` records WHICH eigenvector was
  selected but never its ORIENTATION. `eigs_cis` canonicalises the sign to
  GC-rich-positive and records any flip in `orientation.tsv`, which the report
  surfaces. A silent flip mirrors every saddle plot and would quietly corrupt
  fifteen of the twenty-five panels, so being able to see all three
  eigenvectors side by side is not a luxury.
* **panel_loop_overlay.** The refined loop sets are a frozen fixture --
  human judgement that cannot be recomputed. Raw counts are 55 / 32 / 41 / 26
  (Jarid_Xa, Jarid_Xi, Mecp2_Xa, Mecp2_Xi) and refined 14 / 13 / 10 / 12; this
  panel is the visual record of what was removed.

Pile-up rendering is verbatim from `01_04`: `cmap='coolwarm', scale='log',
sym=True, vmax=2, height=4, plot_ticks=True, center=3`. `scale='log'` with
`sym=True` and `vmax=2` is a log colour scale symmetric about 1, i.e.
`LogNorm(1/2, 2)`; `center=3` is the 3x3 the score is taken over.
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


def _pileup(say, paths):
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    npz_path = next(p for p in paths if p.endswith(".npz"))
    npz = C.read_npz(npz_path)
    data = np.asarray(npz["data"], dtype=float)
    n_loops = int(npz["n_loops"]) if "n_loops" in npz else 0
    flank = int(npz["flank"]) if "flank" in npz else 100000

    centre = data.shape[0] // 2
    score = float(np.nanmean(data[centre - 1:centre + 2, centre - 1:centre + 2]))

    fig, ax = plt.subplots(figsize=(4, 4))
    im = ax.imshow(data, cmap="coolwarm", norm=LogNorm(vmin=0.5, vmax=2.0),
                   interpolation="none")
    ticks = [0, data.shape[0] // 2, data.shape[0] - 1]
    ax.set_xticks(ticks)
    ax.set_xticklabels(["-{}k".format(flank // 1000), "0",
                        "+{}k".format(flank // 1000)], fontsize=7)
    ax.set_yticks(ticks)
    ax.set_yticklabels(["-{}k".format(flank // 1000), "0",
                        "+{}k".format(flank // 1000)], fontsize=7)
    ax.set_title("{}\nn = {} loops   score = {:.3f}".format(
        spec.get("sample", ""), n_loops, score), fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.046, label="obs / exp")
    fig.tight_layout()
    say("  {} loops, central 3x3 mean {:.4f}".format(n_loops, score))
    return fig, n_loops


def _saddle(say, paths):
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    npz_path = next(p for p in paths if p.endswith(".npz"))
    npz = C.read_npz(npz_path)
    if "saddle" in npz:
        S = np.asarray(npz["saddle"], dtype=float)
    else:
        with np.errstate(invalid="ignore", divide="ignore"):
            S = (np.asarray(npz["interaction_sum"], dtype=float)
                 / np.asarray(npz["interaction_count"], dtype=float))
    profile = (np.asarray(npz["saddle_strength_profile"], dtype=float)
               if "saddle_strength_profile" in npz else None)
    score = float(npz["score"]) if "score" in npz else float("nan")

    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    im = axes[0].imshow(S, cmap="coolwarm", norm=LogNorm(vmin=0.5, vmax=2.0),
                        interpolation="none", origin="lower")
    axes[0].set_xlabel("eigenvector rank")
    axes[0].set_ylabel("eigenvector rank")
    axes[0].set_title("{}\n{}   saddle score {:.2f}".format(
        spec.get("name", ""), spec.get("eig", ""), score), fontsize=8)
    fig.colorbar(im, ax=axes[0], fraction=0.046,
                 label="average obs / exp contact frequency")

    if profile is not None:
        x = np.arange(len(profile))
        axes[1].step(x, profile, where="pre")
        axes[1].axhline(0, c="grey", ls="--", lw=1)
        axes[1].set_xlim(0, len(x) // 2)
        axes[1].set_xlabel("extent")
        axes[1].set_ylabel("(AA + BB) / (AB + BA)")
        axes[1].set_title("saddle strength profile", fontsize=8)
        pl.despine(axes[1])
    else:
        axes[1].axis("off")
    fig.tight_layout()
    say("  saddle {}x{}, score {:.4f}".format(S.shape[0], S.shape[1], score))
    return fig, int(S.shape[0])


def _eigenvector(say, paths):
    import matplotlib.pyplot as plt

    path = next(p for p in paths if p.endswith(".tsv"))
    table = C.read_table(path)
    eigs = [c for c in table.columns if c.startswith("E") and c[1:].isdigit()]
    if not eigs:
        return None, "no E1/E2/E3 columns in " + os.path.basename(path)

    x = (C.need(table, "start", path) + C.need(table, "end", path)) / 2.0
    fig, ax = plt.subplots(figsize=(11, 4))
    for eig in eigs:
        ax.plot(x, table[eig].astype(float), lw=0.9, label=eig)
    ax.axhline(0, color="k", lw=0.7, ls="--", alpha=0.5)
    ax.set_xlabel("chrX position (bp)")
    ax.set_ylabel("eigenvector value")
    ax.set_title(
        "{}  {}\nsign is canonicalised GC-rich-positive; any flip is recorded "
        "in orientation.tsv".format(spec.get("name", ""), spec.get("roi", "")),
        fontsize=9)
    ax.legend(fontsize=7, ncol=len(eigs))
    pl.despine(ax)
    fig.tight_layout()
    say("  {} bins, {} eigenvectors".format(len(table), len(eigs)))
    return fig, len(table)


def _loop_overlay(say, paths):
    import matplotlib.pyplot as plt

    raw_path = next((p for p in paths
                     if os.sep + "raw" + os.sep in p and p.endswith(".tsv")), None)
    ref_path = next((p for p in paths
                     if os.sep + "refined" + os.sep in p and p.endswith(".tsv")),
                    None)
    if not raw_path or not ref_path:
        return None, "need both the raw and the refined loop table"

    raw = C.read_table(raw_path)
    refined = C.read_table(ref_path)
    say("  raw {} loops -> refined {} loops".format(len(raw), len(refined)))

    fig, ax = plt.subplots(figsize=(6, 6))
    for table, colour, label, size in ((raw, "#bbbbbb", "raw", 45),
                                       (refined, "#d62728", "refined", 20)):
        if {"start1", "start2"} <= set(table.columns):
            ax.scatter(table["start1"], table["start2"], s=size, c=colour,
                       label="{} (n={})".format(label, len(table)),
                       edgecolor="black", linewidth=0.3, alpha=0.85)
    ax.set_xlabel("anchor 1 (bp)")
    ax.set_ylabel("anchor 2 (bp)")
    ax.set_title("{}\nchromosight calls before and after manual refinement"
                 .format(spec.get("name", "")), fontsize=9)
    ax.legend(fontsize=8)
    pl.despine(ax)
    fig.tight_layout()
    return fig, len(refined)


def _cooler_qc(say, paths):
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    try:
        import cooler
    except ImportError:
        return None, "cooler is not importable in this environment"

    mcool = next(p for p in paths if p.endswith(".mcool"))
    name = spec.get("cooler", os.path.basename(mcool))
    locus = next((k for k in (snakemake.params.windows or {}) if k in name),
                 None)
    resolutions = cooler.fileops.list_coolers(mcool)
    clr = cooler.Cooler(mcool + "::" + resolutions[0])
    if locus:
        start, end = (int(v) for v in snakemake.params.windows[locus])
        region = "chrX:{}-{}".format(start, end)
    else:
        region = "chrX"
    matrix = clr.matrix(balance=True).fetch(region)

    fig, ax = plt.subplots(figsize=(6, 6))
    im = ax.matshow(matrix, cmap="YlOrRd",
                    norm=LogNorm(vmax=np.nanmax(matrix) or 1.0))
    ax.set_title("{}\n{}  ({})".format(name, region, resolutions[0]),
                 fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.046)
    fig.tight_layout()
    say("  {} at {}: {}x{}".format(name, region, *matrix.shape))
    return fig, int(matrix.shape[0])


with C.logging(snakemake) as say:
    pl.style()
    producer = spec.get("producer", "")
    paths = [str(p) for p in snakemake.input]
    say("panel {} ({})".format(spec["panel_id"], producer))

    handlers = {
        "panel_pileup": _pileup,
        "panel_saddle": _saddle,
        "panel_eigenvector": _eigenvector,
        "panel_loop_overlay": _loop_overlay,
        "panel_cooler_qc": _cooler_qc,
    }
    handler = handlers.get(producer)
    if handler is None:
        fig, n = None, "hic_panels.py has no branch for " + producer
    else:
        fig, n = handler(say, paths)

    if fig is None:
        pl.empty_panel(out_path, str(n))
        C.write_n_items(snakemake, 0)
        say("empty panel: {}".format(n))
    else:
        pl.save(fig, out_path, dpi=spec.get("dpi", 300))
        C.write_n_items(snakemake, n)
        say("wrote {}".format(out_path))
