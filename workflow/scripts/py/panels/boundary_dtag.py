"""Fig 4f + Fig 4g + EFig 6b — one figure per degron clone.

    [0,0]  CTCF / no-CTCF count heatmap at escaping boundaries   \
    [0,1]  pie: gene category of those boundaries (all escaping)  > EFig 6b
    [0,2]  pie: CTCF status                                      /
    rows   NodTAG (solid, blue) vs +dTAG (dashed, amber) metaplots
             H3K27me3 row  = Fig 4f
             H3K27ac  row  = Fig 4g
    col 1  all escaping boundaries
    col 2  CTCF+ escaping boundaries (NodTAG peaks)

Published n: E6A7 110 escaping boundaries, 80 CTCF+ / 30 CTCF-; F3 36, 23 / 13.
Fig 4f/4g quote the CTCF+ subsets: E6A7 n = 80, F3 n = 23.

The boundaries are ALWAYS the NodTAG escaping-valley boundaries, in both
conditions. That is not a shortcut: the question the panel asks is what happens
to the signal at a fixed set of positions after degradation. Re-calling
boundaries per condition would compare two different position sets and answer a
different question.

B1621 is the Rad21-degron clone, so its "CTCF" row is Rad21 and it is not one
of the published panels; it is emitted as an extra.
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
ESCAPING = "escaping-gene-valley"
SUBSETS = [("escaping", "All escaping boundaries"),
           ("escaping_ctcf", "CTCF <= 50kbp boundaries (NodTAG)")]


with C.logging(snakemake) as say:
    pl.style()
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    clone = spec.get("clone", "?")
    say("panel {} (clone {})".format(spec["panel_id"], clone))

    ctcf_paths = [p for p in map(str, snakemake.input)
                  if os.sep + "boundary_ctcf" + os.sep in p]
    npz_paths = {}
    for path in map(str, snakemake.input):
        if path.endswith(".npz"):
            stem = os.path.basename(path)[: -len(".npz")]
            if stem.startswith(str(clone) + "_"):
                npz_paths[stem[len(str(clone)) + 1:]] = path

    if not ctcf_paths:
        pl.empty_panel(out_path, "no boundary_ctcf table for " + str(clone))
        C.write_n_items(snakemake, 0)
        raise SystemExit(0)

    ctcf = C.read_table(ctcf_paths[0])
    escaping = ctcf[C.need(ctcf, "valley_class", ctcf_paths[0]) == ESCAPING]
    n_bnd = int(len(escaping))
    n_ctcf = int(C.need(escaping, "has_ctcf", ctcf_paths[0]).astype(bool).sum())
    n_noctcf = n_bnd - n_ctcf
    say("  {} escaping boundaries: {} CTCF+, {} CTCF-".format(
        n_bnd, n_ctcf, n_noctcf))

    marks = [m for m in pl.MARKS if m in npz_paths]
    marks += sorted(m for m in npz_paths if m not in pl.MARKS)
    say("  metaplot rows: {}".format(", ".join(marks) or "(none)"))

    n_rows = max(1, len(marks))
    fig = plt.figure(figsize=(13, 4 + 3.5 * n_rows))
    fig.suptitle("{} - escaping boundaries +/-dTAG   (Fig 4f/4g, EFig 6b)"
                 .format(clone), fontsize=13, y=1.01, fontweight="bold")
    gs = fig.add_gridspec(n_rows + 1, 3, height_ratios=[3] + [2] * n_rows,
                          hspace=0.7, wspace=0.45)

    ax_hm = fig.add_subplot(gs[0, 0])
    counts = np.array([[n_ctcf], [n_noctcf]], dtype=float)
    ax_hm.imshow(counts, cmap="YlOrBr", aspect="auto")
    for i, value in enumerate((n_ctcf, n_noctcf)):
        ax_hm.text(0, i, str(int(value)), ha="center", va="center",
                   fontsize=14, fontweight="bold")
    ax_hm.set_xticks([])
    ax_hm.set_yticks([0, 1])
    ax_hm.set_yticklabels(["CTCF <= 50kbp", "No CTCF <= 50kbp"], fontsize=9)
    ax_hm.set_title("CTCF at escaping boundaries\n(NodTAG peaks)\n"
                    "n = {} boundaries".format(n_bnd), fontsize=9)

    ax_p1 = fig.add_subplot(gs[0, 1])
    if n_bnd:
        ax_p1.pie([n_bnd], colors=[pl.CAT_COLORS[ESCAPING]], startangle=90,
                  wedgeprops={"linewidth": 0.5, "edgecolor": "white"})
    ax_p1.set_title("Boundary gene category\nn = {}".format(n_bnd), fontsize=9)
    ax_p1.legend(handles=[Patch(facecolor=pl.CAT_COLORS[ESCAPING],
                                label="Escaping ({})".format(n_bnd))],
                 loc="upper center", bbox_to_anchor=(0.5, -0.03), ncol=1,
                 fontsize=8, frameon=False)

    ax_p2 = fig.add_subplot(gs[0, 2])
    if n_bnd:
        wedges, _, texts = ax_p2.pie(
            [n_ctcf, n_noctcf],
            colors=[pl.CTCF_PIE_COLOR, pl.NO_CTCF_PIE_COLOR],
            autopct=pl.make_autopct(n_bnd), startangle=90, counterclock=False,
            pctdistance=0.72,
            wedgeprops={"linewidth": 0.5, "edgecolor": "white"})
        for t in texts:
            t.set_fontsize(8)
        ax_p2.legend(wedges, ["CTCF <= 50kbp", "No CTCF <= 50kbp"],
                     loc="upper center", bbox_to_anchor=(0.5, -0.03), ncol=1,
                     fontsize=8, frameon=False)
    ax_p2.set_title("CTCF status (NodTAG peaks)\nn = {}".format(n_bnd),
                    fontsize=9)

    n_used = {}
    rna_axes = []
    for ri, mark in enumerate(marks):
        npz = C.read_npz(npz_paths[mark])
        ylim = pl.MARK_YLIM.get(mark, (None, None))
        is_bottom = ri == len(marks) - 1

        for col, (subset, title) in enumerate(SUBSETS, start=1):
            ax = fig.add_subplot(gs[ri + 1, col])
            key_n = "profile_{}_nodtag".format(subset)
            key_d = "profile_{}_dtag".format(subset)
            n_key = "n_{}".format(subset)
            n = int(npz[n_key]) if n_key in npz else 0
            n_used[subset] = n

            drew = False
            if key_n in npz:
                prof = np.asarray(npz[key_n], dtype=float)
                ax.plot(prof, color=pl.NODTAG_LINE_COLOR, lw=2.0, zorder=2,
                        label="NodTAG (n={})".format(n))
                pl.boundary_xaxis(ax, nbins=len(prof))
                drew = True
            if key_d in npz:
                prof = np.asarray(npz[key_d], dtype=float)
                ax.plot(prof, color=pl.DTAG_LINE_COLOR, lw=1.6, zorder=1,
                        ls="--", label="+dTAG (n={})".format(n))
                if not drew:
                    pl.boundary_xaxis(ax, nbins=len(prof))
                drew = True
            if not drew:
                ax.text(0.5, 0.5, "no profile", ha="center", va="center",
                        transform=ax.transAxes, fontsize=8)

            if ylim != (None, None):
                ax.set_ylim(*ylim)
            ax.tick_params(axis="y", labelsize=6)
            ax.set_ylabel(mark, fontsize=8)
            if ri == 0:
                ax.set_title(title, fontsize=9)
            if is_bottom:
                ax.set_xlabel("Position relative to boundary", fontsize=7)
            pl.despine(ax)
            ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.25), ncol=2,
                      fontsize=6, frameon=False)
            if mark == "RNA-Seq":
                rna_axes.append(ax)

    if rna_axes:
        top = max(ax.get_ylim()[1] for ax in rna_axes)
        for ax in rna_axes:
            ax.set_ylim(top=top)

    pl.save(fig, out_path, dpi=spec.get("dpi", 300))
    C.write_n_items(snakemake, n_bnd, extra={
        "n_escaping": n_bnd, "n_escaping_ctcf": n_ctcf,
        "n_escaping_no_ctcf": n_noctcf,
        "n_profile_escaping": n_used.get("escaping", 0),
        "n_profile_escaping_ctcf": n_used.get("escaping_ctcf", 0),
    })
    say("wrote {}".format(out_path))
