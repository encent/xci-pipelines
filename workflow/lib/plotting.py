"""House style for every figure the pipeline emits.

One import, one look. Every ``panel_*`` script starts with::

    from lib import plotting as pl
    pl.style()

and finishes with ``pl.save(fig, snakemake.output[0])``.

THE ONE RULE THAT IS NOT COSMETIC
---------------------------------
``svg.fonttype = "none"`` — **SVG text stays text**, not outlines.

Decision D-07 keeps the final figure assembly manual: the pipeline stops at
panel level and a human composes Figure 3 in Illustrator, exactly as the
original did. That only works if the SVGs are editable, which means the text
must arrive as ``<text>`` elements with a ``font-family``, not as a soup of
``<path>``. matplotlib's default (``fonttype = "path"``) converts every glyph
to a path and makes the panel unusable for typesetting.

The same reasoning gives ``pdf.fonttype = 42`` (embed TrueType, keep the text
selectable) rather than the default Type 3.

The cost is that the renderer must have the font. We ask for a stack that
degrades to DejaVu Sans, which ships with matplotlib, so a missing Arial
changes the metrics but never the editability.

REPRODUCIBILITY
---------------
``svg.hashsalt`` is pinned and ``metadata={"Date": None}`` is passed on save,
so two runs of the same panel on the same data produce byte-identical SVGs.
Without this, matplotlib stamps a creation date and randomises element ids and
every figure looks changed to git and to the comparison harness.

COLOURS
-------
The Okabe-Ito palette below is transcribed from the original notebooks, not
chosen. `CAT_COLORS` in particular is load-bearing: the gene-class colours in
Fig 2g, Fig 4d, Fig 2h and EFig 2e-g are the same four hexes everywhere and a
reader compares them across figures.

The allelic-ratio colour map is fixture F-8
(``resources/fixtures/figure_colours.yaml``): a truncated ``Greens`` under
``Normalize(0.0424769728421052, 0.40909353528)``. Note the legend of Fig 3a/b,
5b/c/g/h and EFig 8c/d, 10c/d calls these circles a **METALoci D-score**. They
are not. They are the RNA-seq allelic ratio computed by
``01_07_get_allelic_ratio_gtf.py``. The published wording is wrong; the numbers
are right; we reproduce the numbers.
"""

from __future__ import annotations

import os
import warnings

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                       # noqa: E402
from matplotlib.colors import (           # noqa: E402
    LinearSegmentedColormap,
    Normalize,
)

__all__ = [
    "style", "save", "despine", "figure", "CAT_COLORS", "CAT_LABEL",
    "CAT_PIE_ORDER", "MARK_COLOR", "MARK_YLIM", "MARKS", "OTHER_COLOR",
    "CTCF_PIE_COLOR", "NODTAG_LINE_COLOR", "DTAG_LINE_COLOR",
    "allelic_ratio_cmap", "allelic_ratio_norm", "truncate_colormap",
    "boundary_xaxis", "make_autopct", "annotate_pairwise", "venn2_counts",
    "paired_violin", "empty_panel", "load_colours",
]

# ---------------------------------------------------------------------------
# palette — transcribed from xci_valleys_check_LAST.ipynb, do not re-choose
# ---------------------------------------------------------------------------
#: Okabe-Ito, colour-blind safe. Same four hexes in every gene-class figure.
CAT_COLORS = {
    "escaping-gene-valley": "#E69F00",      # orange
    "silent-gene-valley": "#56B4E9",        # sky blue
    "no-expressed-gene-valley": "#0072B2",  # blue
    "no-gene": "#CCCCCC",                   # light grey
}
CAT_LABEL = {
    "escaping-gene-valley": "Escaping",
    "silent-gene-valley": "Silent",
    "no-expressed-gene-valley": "No expressed",
    "no-gene": "No gene",
}
#: Legend / stack order. Bars stack no-gene -> no-expressed -> silent -> escaping.
CAT_STACK_ORDER = [
    "no-gene", "no-expressed-gene-valley", "silent-gene-valley",
    "escaping-gene-valley",
]
CAT_PIE_ORDER = [
    "escaping-gene-valley", "silent-gene-valley", "no-expressed-gene-valley",
    "no-gene",
]

#: Rows of the boundary metaplot figures, in order.
MARKS = ["H3K27me3", "CTCF", "H3K27ac", "RNA-Seq"]
MARK_COLOR = {
    "H3K27me3": "#0072B2", "CTCF": "#D55E00",
    "H3K27ac": "#009E73", "RNA-Seq": "#CC79A7",
}
#: Fixed y-limits so the five clone figures are comparable at a glance.
#: RNA-Seq is (None, None) and synchronised across clones after drawing.
MARK_YLIM = {
    "H3K27me3": (0, 1.2), "CTCF": (0, 0.16),
    "H3K27ac": (0, 0.16), "RNA-Seq": (None, None),
}
OTHER_COLOR = "#AAAAAA"
CTCF_PIE_COLOR = "#E69F00"
NO_CTCF_PIE_COLOR = "#CCCCCC"
NODTAG_LINE_COLOR = "#2166AC"   # dark blue, solid
DTAG_LINE_COLOR = "#E69F00"     # Okabe-Ito amber, dashed

#: Loop / saddle violins.
VIOLIN_FACE = "0.9"
HIGHLIGHT_COLOR = "red"
PLAIN_COLOR = "black"

_FONT_STACK = ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"]

_STYLED = False


def style(font_size=8, force=False):
    """Install the house style. Idempotent; safe to call from every script."""
    global _STYLED
    if _STYLED and not force:
        return
    matplotlib.rcParams.update({
        # --- THE point of this module -----------------------------------
        "svg.fonttype": "none",     # keep SVG text as text (D-07)
        "pdf.fonttype": 42,         # TrueType in PDF, text stays selectable
        "ps.fonttype": 42,
        # --- deterministic output ---------------------------------------
        "svg.hashsalt": "xci-pipeline",
        # --- type --------------------------------------------------------
        "font.family": "sans-serif",
        "font.sans-serif": _FONT_STACK,
        "font.size": font_size,
        "axes.titlesize": font_size + 1,
        "axes.labelsize": font_size,
        "xtick.labelsize": font_size - 1,
        "ytick.labelsize": font_size - 1,
        "legend.fontsize": font_size - 1,
        "figure.titlesize": font_size + 3,
        # --- ink ---------------------------------------------------------
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": False,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "lines.linewidth": 1.2,
        "patch.linewidth": 0.5,
        "legend.frameon": False,
        # --- output ------------------------------------------------------
        "figure.dpi": 100,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.transparent": False,
        "image.cmap": "viridis",
    })
    # seaborn, when present, must not clobber the rcParams above.
    try:
        import seaborn as sns
        sns.set_theme(context="paper", style="ticks", rc=matplotlib.rcParams)
    except Exception:                                          # pragma: no cover
        pass
    _STYLED = True


def figure(*args, **kwargs):
    """``plt.subplots`` with the house style guaranteed to be installed."""
    style()
    return plt.subplots(*args, **kwargs)


def despine(ax, top=True, right=True, left=False, bottom=False):
    for side, drop in (("top", top), ("right", right),
                       ("left", left), ("bottom", bottom)):
        if drop:
            ax.spines[side].set_visible(False)
    return ax


def save(fig, path, dpi=300, close=True):
    """Write exactly one file, in the format implied by its extension.

    Each output format is its own Snakemake job (``{ext}`` is a wildcard), so a
    panel script writes one file and never guesses at siblings.
    """
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    ext = os.path.splitext(path)[1].lower().lstrip(".")
    kwargs = dict(dpi=dpi, bbox_inches="tight")
    if ext == "svg":
        # Date=None strips the <dc:date> stamp so two runs are byte-identical.
        kwargs["metadata"] = {"Date": None}
    elif ext == "pdf":
        kwargs["metadata"] = {"CreationDate": None}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fig.savefig(path, **kwargs)
    if close:
        plt.close(fig)
    return path


def empty_panel(path, message, size=(6, 4)):
    """Write a legible placeholder instead of a zero-byte file.

    A panel with no data is a fact about the run, not a crash. The manifest
    records it as ``status: empty`` and the report shows this card, which is far
    easier to act on than an empty file or a traceback at hour 30.
    """
    style()
    fig, ax = plt.subplots(figsize=size)
    ax.axis("off")
    ax.text(0.5, 0.5, message, ha="center", va="center", fontsize=10,
            wrap=True, transform=ax.transAxes)
    return save(fig, path)


# ---------------------------------------------------------------------------
# colour helpers — fixture F-8
# ---------------------------------------------------------------------------
def load_colours(path=None):
    """Read ``resources/fixtures/figure_colours.yaml`` (F-8)."""
    import yaml
    if path is None:
        path = os.path.join(
            os.path.dirname(os.path.dirname(os.path.dirname(
                os.path.abspath(__file__)))),
            "resources", "fixtures", "figure_colours.yaml",
        )
    with open(path) as fh:
        return yaml.safe_load(fh)


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=256):
    """Slice a colormap. ``truncate: [0.0, 1.0]`` in F-8 is the identity."""
    if isinstance(cmap, str):
        cmap = plt.get_cmap(cmap)
    if (minval, maxval) == (0.0, 1.0):
        return cmap
    return LinearSegmentedColormap.from_list(
        "trunc({},{:.2f},{:.2f})".format(cmap.name, minval, maxval),
        cmap(np.linspace(minval, maxval, n)),
    )


def allelic_ratio_cmap(colours=None):
    """The truncated ``Greens`` of fixture F-8."""
    colours = colours or load_colours()
    ar = colours["allelic_ratio"]
    lo, hi = ar.get("truncate", [0.0, 1.0])
    return truncate_colormap(ar.get("cmap", "Greens"), lo, hi)


def allelic_ratio_norm(colours=None):
    """``Normalize(0.0424769728421052, 0.40909353528)`` — fixture F-8.

    Hard-coded in the original and impossible to recover from the data (it is
    the WT clone range, so a degron-only run would rescale and the green
    circles would stop being comparable between figures).

    The legend calls this a METALoci D-score. It is the RNA-seq allelic ratio.
    """
    colours = colours or load_colours()
    lo, hi = colours["allelic_ratio"]["norm"]
    return Normalize(lo, hi)


# ---------------------------------------------------------------------------
# boundary metaplot helpers
# ---------------------------------------------------------------------------
def boundary_xaxis(ax, nbins=100, flank=100000, rotation=30, fontsize=6):
    """The five-tick ``-100k .. bnd .. +100k`` axis used by every metaplot."""
    half = nbins // 2
    ticks = [0, half // 2, half, half + half // 2, nbins - 1]
    labels = [
        "-{}k".format(flank // 1000), "-{}k".format(flank // 2000), "bnd",
        "+{}k".format(flank // 2000), "+{}k".format(flank // 1000),
    ]
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, rotation=rotation, fontsize=fontsize)
    ax.axvline(half, color="black", ls="--", lw=0.8, alpha=0.6)
    return ax


def make_autopct(total_n):
    """Pie labels showing COUNTS, not percentages — the original's convention."""
    def autopct(pct):
        n = int(round(pct / 100.0 * total_n))
        return str(n) if n > 0 else ""
    return autopct


# ---------------------------------------------------------------------------
# statistics on plots
# ---------------------------------------------------------------------------
def annotate_pairwise(ax, data, x, y, pairs, test="Wilcoxon", pvalues=None,
                      **kwargs):
    """Draw significance brackets, with the p-values computed by the caller.

    The statistic is ALWAYS computed by the panel script (with scipy) and
    written into ``paper_numbers.tsv``; this function only draws it. That split
    is deliberate:

      * the numbers in the table and the numbers on the figure cannot disagree;
      * ``statannotations`` is an optional dependency — if it is missing the
        bracket is drawn by hand and the figure still carries the right p.

    It also stops us from re-introducing the original's ``annotate_violin_mwu``
    bug, which announced "Mann-Whitney" in the title and called
    ``scipy.stats.ttest_rel`` (deviation D-5). Honour
    ``config.legacy.stackup_violin_test`` when choosing the test, and pass the
    resulting p here.
    """
    if pvalues is None:
        raise ValueError("compute the p-values in the panel script and pass them")
    try:
        from statannotations.Annotator import Annotator
        annot = Annotator(ax, pairs, data=data, x=x, y=y, **kwargs)
        annot.configure(test=None, text_format="star", loc="inside", verbose=0)
        annot.set_pvalues(list(pvalues))
        annot.annotate()
        return ax
    except Exception:
        pass

    # --- manual fallback, identical numbers ------------------------------
    lo, hi = ax.get_ylim()
    span = hi - lo
    categories = list(data[x].unique()) if hasattr(data, "columns") else []
    for i, ((a, b), p) in enumerate(zip(pairs, pvalues)):
        try:
            x1, x2 = categories.index(a), categories.index(b)
        except ValueError:
            x1, x2 = 0, 1
        yy = hi - span * 0.05 * (i + 1)
        ax.plot([x1, x1, x2, x2], [yy, yy + span * 0.01, yy + span * 0.01, yy],
                lw=0.8, c="black")
        ax.text((x1 + x2) / 2.0, yy + span * 0.015, _stars(p),
                ha="center", va="bottom", fontsize=7)
    return ax


def _stars(p):
    if p is None or not np.isfinite(p):
        return "n/a"
    for thr, s in ((1e-4, "****"), (1e-3, "***"), (1e-2, "**"), (5e-2, "*")):
        if p < thr:
            return s
    return "ns"


def paired_violin(ax, table, group1, group2, value_label, highlight=None,
                  ylim=(0, 4.0), alpha=0.2):
    """The paired violin of ``02_00_analysis_*``: violin + joined per-clone dots.

    ``table`` is wide: one row per clone-allele, columns ``group1`` / ``group2``.
    Dots of the highlighted group are red, the other black, joined by a thin
    grey line — that line is the whole point, since the test is paired.
    """
    style()
    try:
        import seaborn as sns
        melted = table.melt(value_vars=[group1, group2],
                            var_name="group", value_name=value_label)
        sns.violinplot(data=melted, x="group", y=value_label, inner="quart",
                       color=VIOLIN_FACE, ax=ax)
    except Exception:
        ax.violinplot([table[group1].dropna().values,
                       table[group2].dropna().values], positions=[0, 1])
        ax.set_xticks([0, 1])
        ax.set_xticklabels([group1, group2])
        ax.set_ylabel(value_label)

    for idx in table.index:
        y_vals = [table.loc[idx, group1], table.loc[idx, group2]]
        colors = [
            HIGHLIGHT_COLOR if group1 == highlight else PLAIN_COLOR,
            HIGHLIGHT_COLOR if group2 == highlight else PLAIN_COLOR,
        ]
        ax.scatter([0, 1], y_vals, c=colors, zorder=5, alpha=alpha)
        ax.plot([0, 1], y_vals, c="gray", linewidth=0.7, alpha=0.5)

    if ylim is not None:
        ax.set_ylim(*ylim)
    despine(ax)
    return ax


def venn2_counts(ax, only_a, shared, only_b, label_a, label_b,
                 colors=("#E69F00", "#0072B2")):
    """Two-set Venn from counts, with a dependency-free fallback.

    ``matplotlib_venn`` is in ``environment.yml`` and is used when present. The
    fallback draws two circles and the three numbers, which is enough for the
    structural test and for a reader to check EFig 6c's 1/44/11 and 1/12/6.
    """
    try:
        from matplotlib_venn import venn2
        v = venn2(subsets=(int(only_a), int(only_b), int(shared)),
                  set_labels=(label_a, label_b), ax=ax)
        for patch, colour in zip(v.patches, colors):
            if patch is not None:
                patch.set_color(colour)
                patch.set_alpha(0.55)
        return ax
    except Exception:
        pass

    from matplotlib.patches import Circle
    ax.add_patch(Circle((0.38, 0.5), 0.30, color=colors[0], alpha=0.45, lw=0))
    ax.add_patch(Circle((0.62, 0.5), 0.30, color=colors[1], alpha=0.45, lw=0))
    ax.text(0.22, 0.5, str(int(only_a)), ha="center", va="center", fontsize=11)
    ax.text(0.50, 0.5, str(int(shared)), ha="center", va="center", fontsize=11)
    ax.text(0.78, 0.5, str(int(only_b)), ha="center", va="center", fontsize=11)
    ax.text(0.30, 0.86, label_a, ha="center", va="center", fontsize=9)
    ax.text(0.70, 0.86, label_b, ha="center", va="center", fontsize=9)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.axis("off")
    return ax
