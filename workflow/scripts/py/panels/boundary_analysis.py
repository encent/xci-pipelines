"""EFig 2e + 2f + 2g — one figure, three published panels, one per clone.

    left   [0,0]  2x2 contingency Escapee x CTCF + two-sided Fisher p   = EFig 2g
    middle [*,1]  H3K27me3 metaplot, escaping vs other boundaries       = EFig 2e
    right  [*,2]  H3K27ac  metaplot, same split                         = EFig 2f

They are three panels of ONE figure produced by ONE code path
(`xci_valleys_check_LAST.ipynb`). Splitting them into three rules would let
them drift; they are modelled as one.

THE TRAP (R-2 / P-1)
--------------------
The profiles come from ``P.boundary_profile(track, signal)`` — the
NO-background-subtraction path. ``P.stackup(track, signal, variant)`` reads the
same boundary anchors and the same 20 bp bigWigs but subtracts a
100-iteration random circular-shift background (`default_rng(seed=k+42)`).
Feeding this panel from the stackups yields smooth, plausible, WRONG curves:
EFig 2e/2f would not reproduce and nothing would look broken.

``P.boundary_profile`` carries the R-2 note in ``paths.py`` for the same
reason. If you are editing this file and reach for ``P.stackup``, stop.

Layout is verbatim from the original: figsize (13, 18), 5x3 gridspec with
height_ratios [3,2,2,2,2], hspace 0.7, wspace 0.45; rows are H3K27me3, CTCF,
H3K27ac, RNA-Seq; column 1 is escaping vs other, column 2 is CTCF+ vs CTCF-
among escaping boundaries. RNA-Seq has no fixed y-limit and is synchronised
across clones by the caller reading `MARK_YLIM`.
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


def _profile_for(npz, mask):
    """Mean profile over the rows of the stored matrix selected by `mask`.

    The .npz stores the full per-boundary matrix as well as the pooled
    profile, so the panel can re-split by gene class without recomputing
    anything from bigWigs -- which is the whole point of the data/picture
    split.
    """
    if "matrix" not in npz:
        return None, 0
    matrix = np.asarray(npz["matrix"], dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != len(mask):
        # The profile is still usable even if the row order cannot be matched.
        return np.asarray(npz["profile"], dtype=float), int(matrix.shape[0])
    sel = matrix[np.asarray(mask, dtype=bool)]
    if sel.size == 0:
        return None, 0
    with np.errstate(invalid="ignore"):
        return np.nanmean(sel, axis=0), int(sel.shape[0])


with C.logging(snakemake) as say:
    pl.style()
    import matplotlib.pyplot as plt
    from scipy.stats import fisher_exact

    clone = spec.get("clone", "?")
    say("panel {} (clone {})".format(spec["panel_id"], clone))

    ctcf_path = C.find_input(snakemake.input, ".tsv")
    ctcf_paths = [p for p in map(str, snakemake.input)
                  if os.sep + "boundary_ctcf" + os.sep in p]
    if not ctcf_paths:
        pl.empty_panel(out_path, "no boundary_ctcf table for " + str(clone))
        C.write_n_items(snakemake, 0)
        raise SystemExit(0)
    ctcf = C.read_table(ctcf_paths[0])
    track = os.path.basename(ctcf_paths[0])[: -len(".tsv")]

    is_escaping = (C.need(ctcf, "valley_class", ctcf_paths[0]) == ESCAPING).values
    has_ctcf = C.need(ctcf, "has_ctcf", ctcf_paths[0]).astype(bool).values

    # --- EFig 2g: the 2x2 and its Fisher p -------------------------------
    a = int((is_escaping & has_ctcf).sum())        # escaping, CTCF+
    b = int((~is_escaping & has_ctcf).sum())       # other,    CTCF+
    c = int((is_escaping & ~has_ctcf).sum())       # escaping, CTCF-
    d = int((~is_escaping & ~has_ctcf).sum())      # other,    CTCF-
    odds, p_fisher = fisher_exact([[a, b], [c, d]])
    n_total = len(ctcf)
    n_esc = int(is_escaping.sum())
    say("  contingency  escaping/CTCF+ {}  other/CTCF+ {}  "
        "escaping/CTCF- {}  other/CTCF- {}".format(a, b, c, d))
    say("  Fisher two-sided p = {:.3e}  (OR {:.3f})  n = {} boundaries, "
        "{} escaping".format(p_fisher, odds, n_total, n_esc))

    # --- the metaplots ----------------------------------------------------
    npz_paths = {}
    for path in map(str, snakemake.input):
        if not path.endswith(".npz"):
            continue
        stem = os.path.basename(path)[: -len(".npz")]
        if stem.startswith(track + "_"):
            npz_paths[stem[len(track) + 1:]] = path

    marks = [m for m in pl.MARKS if m in npz_paths]
    extra = [m for m in npz_paths if m not in pl.MARKS]
    marks += sorted(extra)
    say("  metaplot rows: {}".format(", ".join(marks) or "(none)"))

    n_rows = max(1, len(marks))
    fig = plt.figure(figsize=(13, 4 + 3.5 * n_rows))
    fig.suptitle("{}   (EFig 2e/2f/2g)".format(C.track_label(track)),
                 fontsize=13, y=1.01, fontweight="bold")
    gs = fig.add_gridspec(n_rows + 1, 3, height_ratios=[3] + [2] * n_rows,
                          hspace=0.7, wspace=0.45)

    # [0,0] contingency heatmap
    ax_hm = fig.add_subplot(gs[0, 0])
    matrix = np.array([[a, b], [c, d]], dtype=float)
    im = ax_hm.imshow(matrix, cmap="YlOrBr")
    for i in range(2):
        for j in range(2):
            ax_hm.text(j, i, "{:d}".format(int(matrix[i, j])), ha="center",
                       va="center", fontsize=13, fontweight="bold")
    ax_hm.set_xticks([0, 1])
    ax_hm.set_xticklabels(["Escaping", "Other"], fontsize=9)
    ax_hm.set_yticks([0, 1])
    ax_hm.set_yticklabels(["CTCF <= 50kbp", "No CTCF <= 50kbp"], fontsize=8)
    ax_hm.set_xlabel("Gene class", fontsize=9)
    ax_hm.set_title("CTCF x gene class\nFisher p = {:.2e}   (n = {} bnd)"
                    .format(p_fisher, n_total), fontsize=9)

    # [0,1] all-boundary gene-class pie
    ax_p1 = fig.add_subplot(gs[0, 1])
    cls = C.need(ctcf, "valley_class", ctcf_paths[0])
    counts = [int((cls == k).sum()) for k in pl.CAT_PIE_ORDER]
    wedges, _, texts = ax_p1.pie(
        counts, colors=[pl.CAT_COLORS[k] for k in pl.CAT_PIE_ORDER],
        explode=[0.15, 0, 0, 0], autopct=pl.make_autopct(n_total),
        startangle=90, counterclock=False, pctdistance=0.72,
        wedgeprops={"linewidth": 0.5, "edgecolor": "white"})
    for t in texts:
        t.set_fontsize(8)
    ax_p1.set_title("All boundaries\nn = {}".format(n_total), fontsize=9)
    ax_p1.legend(wedges, [pl.CAT_LABEL[k] for k in pl.CAT_PIE_ORDER],
                 loc="upper center", bbox_to_anchor=(0.5, -0.03), ncol=2,
                 fontsize=8, frameon=False)

    # [0,2] escaping-boundary CTCF pie  (also the Fig 2h content)
    ax_p2 = fig.add_subplot(gs[0, 2])
    if n_esc:
        wedges2, _, texts2 = ax_p2.pie(
            [a, c], colors=[pl.CTCF_PIE_COLOR, pl.NO_CTCF_PIE_COLOR],
            autopct=pl.make_autopct(n_esc), startangle=90, counterclock=False,
            pctdistance=0.72,
            wedgeprops={"linewidth": 0.5, "edgecolor": "white"})
        for t in texts2:
            t.set_fontsize(8)
        ax_p2.legend(wedges2, ["CTCF <= 50kbp", "No CTCF <= 50kbp"],
                     loc="upper center", bbox_to_anchor=(0.5, -0.03), ncol=1,
                     fontsize=8, frameon=False)
    else:
        ax_p2.text(0.5, 0.5, "No escaping\nboundaries", ha="center",
                   va="center", fontsize=10, transform=ax_p2.transAxes)
    ax_p2.set_title("Escaping boundaries\nn = {}".format(n_esc), fontsize=9)

    rna_axes = []
    for ri, mark in enumerate(marks):
        npz = C.read_npz(npz_paths[mark])
        colour = pl.MARK_COLOR.get(mark, "#333333")
        ylim = pl.MARK_YLIM.get(mark, (None, None))
        is_bottom = ri == len(marks) - 1

        for col, (primary, secondary, plabel, slabel) in enumerate((
            (is_escaping, ~is_escaping, "Escaping", "Other"),
            (is_escaping & has_ctcf, is_escaping & ~has_ctcf,
             "CTCF <= 50kbp", "No CTCF"),
        ), start=1):
            ax = fig.add_subplot(gs[ri + 1, col])
            prof_s, n_s = _profile_for(npz, secondary)
            prof_p, n_p = _profile_for(npz, primary)
            nbins = len(prof_p) if prof_p is not None else (
                len(prof_s) if prof_s is not None else 100)
            if prof_s is not None:
                ax.plot(prof_s, color=pl.OTHER_COLOR, lw=1.4, zorder=1,
                        label="{} (n={})".format(slabel, n_s))
            if prof_p is not None:
                ax.plot(prof_p, color=colour, lw=2.0, zorder=2,
                        label="{} (n={})".format(plabel, n_p))
            pl.boundary_xaxis(ax, nbins=nbins)
            if ylim != (None, None):
                ax.set_ylim(*ylim)
            ax.tick_params(axis="y", labelsize=6)
            ax.set_ylabel(mark, fontsize=8)
            if ri == 0:
                ax.set_title("Escaping vs other" if col == 1
                             else "CTCF <= 50kbp vs not (escaping)",
                             fontsize=9)
            if is_bottom:
                ax.set_xlabel("Position relative to boundary", fontsize=7)
            pl.despine(ax)
            ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.25), ncol=2,
                      fontsize=6, frameon=False)
            if mark == "RNA-Seq":
                rna_axes.append(ax)

    # RNA-Seq has no published y-limit; share one within the figure so the
    # two columns are comparable. Across clones, MARK_YLIM would be needed.
    if rna_axes:
        top = max(ax.get_ylim()[1] for ax in rna_axes)
        for ax in rna_axes:
            ax.set_ylim(top=top)

    pl.save(fig, out_path, dpi=spec.get("dpi", 300))
    C.write_n_items(snakemake, n_total, extra={
        "n_escaping": n_esc, "escaping_ctcf": a, "escaping_no_ctcf": c,
        "other_ctcf": b, "other_no_ctcf": d, "fisher_p": p_fisher,
    })
    say("wrote {}".format(out_path))
