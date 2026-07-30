# =====================================================================
#  90_visualise.smk — THE SINGLE VISUALISATION STEP
#
#  In the original, figures came out of the R normalization scripts, the
#  valley caller, six Stackups notebooks, three notebooks in a *different*
#  project, eight `01_0*` scripts, `metaloci figure`, four `02_0*` notebook
#  families -- and one of four near-identical `02_03` variants, only one of
#  which made Fig 6. There was no manifest and no single place to look.
#
#  Here there is exactly one: every image the pipeline emits is written by
#  a rule in THIS file, and nowhere else.
#  `tests/lint_figures_only_in_viz.py` fails the build if that stops being
#  true. Every other module emits data -- .tsv, .npz, .bed, .bw -- and this
#  module turns data into pictures.
#
#  ------------------------------------------------------------------
#  ONE LOGICAL STEP, ~40 PARALLEL JOBS
#  ------------------------------------------------------------------
#  `rule visualise` is a bodyless aggregator: no shell, no script, just an
#  input list and a touch file. To the user and to the DAG it is one
#  terminal step that runs last. To the scheduler it is a wide, cheap
#  fan-in over ~40 independent panel rules, so a tweak to one panel
#  redraws one panel.
#
#  ------------------------------------------------------------------
#  HOW ~40 RULES SHARE ONE OUTPUT PATTERN WITHOUT AMBIGUITY
#  ------------------------------------------------------------------
#  Every panel rule writes results/figures/{group}/{name}.{ext} and
#  constrains {name} to an alternation of exactly the registry ids it
#  owns. Two rules can therefore never match the same path, and a typo in
#  panels.yaml fails at DAG build instead of at hour 30.
#
#  ------------------------------------------------------------------
#  PANELS WHOSE UPSTREAM MODULE IS NOT BUILT YET
#  ------------------------------------------------------------------
#  This module was written against the frozen path contract in
#  workflow/lib/paths.py while several feature modules were still
#  placeholders. `_has_producer()` asks the already-registered rule graph
#  whether anything can make a given file; a panel none of whose inputs
#  can be produced (and which is not on disk) is marked `blocked` and left
#  out of the request, with the reason recorded in figure_manifest.tsv.
#  As each module lands its panels light up with no edit here.
#
#  ------------------------------------------------------------------
#  THREE TRAPS, ALL LIVE IN THIS FILE
#  ------------------------------------------------------------------
#  P-2  Fig 2g / Fig 4d plot the GENE-FILTERED valleys (375/301/258/300/321),
#       not the raw HMM calls (377/304/261/302/324). Both look plausible.
#       panel_valley_gene_content reads P.gene_content, which is built from
#       `filter_valleys_genes`, and paper_numbers.tsv asserts the totals.
#  R-2  EFig 2e/2f come from P.boundary_profile -- NO background
#       subtraction. P.stackup subtracts a 100-shift random circular
#       background. Same inputs, different arithmetic, different curves.
#       panel_stackup_composite exists for the non-paper extras only.
#  F-6  The Fig 6 grid order is hard-coded in the original notebook and is
#       not derivable. It comes from resources/fixtures/fig6_layout.yaml.
#
#  Replaces: 12_normalization_*.r (MA plots), TEST_03_00/01/02 notebooks,
#            Stackups_wt*/Stackups_degrons_FINAL*, xci_valleys_check_LAST,
#            xci_valley_overlap, xci_density_plots, 01_01/01_02/01_08,
#            00_08_viz_coolers, `metaloci figure`, 02_00/02_01/02_02/02_03.
# =====================================================================

import os as _os
import re as _re
import sys as _sys

from lib import panels as _panels

FIG = config.get("figures") or {}
FORMATS = FIG.get("formats", ["svg", "pdf"])
FIG_DPI = int(FIG.get("dpi", 300))

GROUP_RE = r"figure_\d+|extended_\d+|qc|extras"
EXT_RE = r"svg|pdf|png"

LIBDIR = _os.path.join(workflow.basedir, "workflow")
REGISTRY_PATH = config["paths"].get("panels", "config/panels.yaml")
if not _os.path.isabs(REGISTRY_PATH):
    REGISTRY_PATH = _os.path.join(workflow.basedir, REGISTRY_PATH)


# ---------------------------------------------------------------------
# Is there a rule that can make this file?
#
# 90_visualise.smk is included LAST, so every other module's rules are
# already registered and can be interrogated. This is what lets the
# visualisation stage be built in parallel with the feature stages
# instead of after them.
# ---------------------------------------------------------------------
def _registered_rules():
    for attr in ("_rules", "rules"):
        obj = getattr(workflow, attr, None)
        if isinstance(obj, dict):
            return list(obj.values())
        if obj is not None:
            try:
                return list(obj)
            except TypeError:
                pass
    return None


def _build_output_matchers():
    rules = _registered_rules()
    if rules is None:
        logger.warning(
            "90_visualise: cannot inspect the rule graph on this Snakemake "
            "version; every panel will be requested and a missing upstream "
            "will be a normal missing-input error."
        )
        return None
    from snakemake.io import regex as _sm_regex

    matchers = []
    for rule in rules:
        for out in getattr(rule, "output", []) or []:
            try:
                matchers.append(_re.compile("^" + _sm_regex(str(out))))
            except Exception:
                continue
    return matchers


_OUTPUT_MATCHERS = _build_output_matchers()


def _has_producer(path):
    """True if a rule can build `path`, or it is already on disk (a fixture)."""
    if _OUTPUT_MATCHERS is None:
        return True
    if _os.path.exists(path):
        return True
    return any(m.match(str(path)) for m in _OUTPUT_MATCHERS)


# ---------------------------------------------------------------------
# The sample / cooler universe, sliced the way panels need it
# ---------------------------------------------------------------------
VIZ_ME3 = SS.tracks_of(mark="H3K27me3")
VIZ_ME3_XI = SS.tracks_of(mark="H3K27me3", allele="Xi")


def _cond(track):
    return naming.parse_track(track)["condition"]


def _clone_of(track):
    return naming.parse_track(track)["clone"]


def _wt_xi_tracks():
    order = ["E6", "C5", "B1", "JTG", "CL30"]           # published Fig 2g order
    hit = [t for t in VIZ_ME3_XI if _cond(t) == "WT"]
    return sorted(hit, key=lambda t: (order.index(_clone_of(t))
                                      if _clone_of(t) in order else 99, t))


def _degron_xi_tracks(which):
    """`which` in {dTAG, NodTAG}: the degron Xi H3K27me3 tracks."""
    suffix = "-dTAG" if which == "dTAG" else "-NodTAG"
    return sorted(t for t in VIZ_ME3_XI if _cond(t).endswith(suffix))


def _xi_track(clone, which):
    for t in VIZ_ME3_XI:
        if _clone_of(t) != clone:
            continue
        cond = _cond(t)
        if which == "WT" and cond == "WT":
            return t
        if which == "dTAG" and cond.endswith("-dTAG"):
            return t
        if which == "NodTAG" and cond.endswith("-NodTAG"):
            return t
    raise WorkflowError(
        "no {} Xi H3K27me3 track for clone {!r}".format(which, clone))


def _dtag_pair_clones():
    """Clones with BOTH a dTAG and a NodTAG Xi H3K27me3 track."""
    out = []
    for clone in SS.clones:
        conds = {_cond(t) for t in VIZ_ME3_XI if _clone_of(t) == clone}
        if any(c.endswith("-dTAG") for c in conds) and \
           any(c.endswith("-NodTAG") for c in conds):
            out.append(clone)
    return sorted(out)


def _shared(name, fallback):
    """Reuse a global that an earlier module already computed, if it is there.

    31_hic_features.smk is included before this file, so its job tables
    (PILEUP_JOBS, COMP_NAMES, EIGS, ...) are already in this namespace.
    Deriving them a second time here would be two sources of truth for one
    enumeration, and the copy would rot the first time RD-3 changed a pairing
    rule -- silently, because a wrong-but-plausible sample id still names a
    file. The fallback keeps this module loadable on its own.
    """
    value = globals().get(name)
    return value if value else fallback()


def _comp_cooler_names():
    """Coolers that have compartments called on them (RD-3's COMP_NAMES).

    Fallback: the merge_comps memberships in the sheet, which are comma-joined.
    """
    def _fallback():
        out = set()
        for name in SS.merged_cooler_names("merged_comps"):
            out.update(n for n in str(name).split(",") if n)
        return sorted(out)

    return sorted(_shared("COMP_NAMES", _fallback))


def _metaloci_datasets(run):
    try:
        return SS.metaloci_datasets(run)
    except Exception:
        return []


def _roi_window(locus, roi="full"):
    cfg = config["loci"][locus]["regions"][roi]
    return int(cfg[0]), int(cfg[1])


ENRICH_ANCHOR = {"CTCF": "motif", "Rad21": "motif"}
ENRICH_PAPER_REF = {"CTCF": ["Fig 1e"], "H3K27ac": ["Fig 2a"],
                    "H3K27me3": ["Fig 2d"]}
ENRICH_GROUP = {"CTCF": "figure_01", "H3K27ac": "figure_02",
                "H3K27me3": "figure_02"}


# ---------------------------------------------------------------------
# Selectors — the registry says "@wt_clones", this says what that is
# ---------------------------------------------------------------------
def _enrichment_rows():
    rows = []
    for mark in ("CTCF", "H3K27ac", "H3K27me3"):
        if not SS.tracks_of(mark=mark):
            continue
        rows.append({
            "mark": mark,
            "anchor": ENRICH_ANCHOR.get(mark, "tss"),
            "_paper_ref": ENRICH_PAPER_REF.get(mark, []),
            "_group": ENRICH_GROUP.get(mark, "qc"),
        })
    return rows


def _fig6_composites():
    """Fig 6: one file per locus, ids taken from the F-6 layout fixture."""
    slug = {"Jarid": "kdm5c", "Mecp2": "mecp2"}
    return [
        {"locus": locus, "locus_slug": slug[locus], "run": "wt"}
        for locus in LOCI if locus in slug
    ]


def _eig_rows():
    return [{"roi": roi, "name": name}
            for roi in ROI_TYPES for name in _comp_cooler_names()]


def _eigs():
    return list(_shared("EIGS", lambda: ["E1", "E2", "E3"]))


def _saddle_rows():
    """One saddle figure per (roi, cooler, eigenvector).

    All three eigenvectors, deliberately: compartments*.json records WHICH
    eigenvector was selected but never its orientation, so being able to put
    E1/E2/E3 side by side is how a selection or sign problem gets diagnosed.
    """
    return [{"roi": roi, "name": name, "eig": eig}
            for roi in ROI_TYPES
            for name in _comp_cooler_names()
            for eig in _eigs()]


def _pileup_jobs():
    """(comparison, exp, sample) triples, from 31_hic_features.smk.

    For `Xa_vs_Xi` the sample id is NOT the cooler name: it carries the allele
    the loops were CALLED on (`Xa_Jarid_B1_WT_G1_Xa`), reproducing the
    original's `pileup_{Xa|Xi}_{cooler}.svg`. Re-deriving that convention here
    is exactly the sort of duplicated rule that drifts, so we read RD-3's job
    table instead of guessing.
    """
    jobs = globals().get("PILEUP_JOBS")
    return sorted(jobs) if jobs else []


def _pileup_rows():
    """Extras only, and only for `full`.

    All three ROIs would be three times more files than anyone opens, and
    nothing published depends on this fan-out: the paper pile-ups are drawn by
    panel_loops_comps_grid from the same .npz.
    """
    return [{"comparison": c, "roi": "full", "exp": e, "sample": s}
            for (c, e, s) in _pileup_jobs()]


def _stackup_rows():
    """Stackup panels, one row per (track, signal, variant) that actually exists.

    The third signal is per-clone, NOT globally CTCF: B1621 is the Rad21-degron
    clone and has no CTCF CUT&RUN at all, so a Cartesian product over a fixed
    ("H3K27me3", "H3K27ac", "CTCF") asks `boundary_stackup` for
    CTCF_B1621_Rad21-NodTAG_Xi, which is not in samples.tsv and aborts the DAG.

    `_boundary_signals` applies the original's selector from
    Stackups_degrons_FINAL.ipynb:

        if (signal=="Rad21" and "B1621" not in ...) or
           (signal=="CTCF"  and "B1621" in ...): continue

    RNA-Seq is excluded here deliberately: these extras only need the three
    chromatin rows, and none of them reached the paper (the archaeology notes 5B).
    """
    rows = []
    for track in VIZ_ME3_XI:
        for signal in _boundary_signals(track):
            if signal == "RNA-Seq":
                continue
            for variant in ("all", "motif_yes", "motif_no"):
                rows.append({"track": track, "signal": signal, "variant": variant})
    return rows


def _density_rows():
    rows = []
    for clone in _dtag_pair_clones():
        for allele in ("Xi", "Xa", "Gall"):
            if any(naming.parse_track(t)["allele"] == allele
                   and _clone_of(t) == clone for t in VIZ_ME3):
                rows.append({"clone": clone, "allele": allele})
    return rows


def _metaloci_rows():
    rows = []
    for run in config["metaloci"]["runs"]:
        for ds in _metaloci_datasets(run):
            for sig in SS.metaloci_signals_for(ds):
                rows.append({"run": run, "dataset": ds, "signal": sig})
    return rows


def _metaloci_kk_rows():
    return [{"run": run, "dataset": ds}
            for run in config["metaloci"]["runs"]
            for ds in _metaloci_datasets(run)]


def _metaloci_signal_rows():
    return [{"locus": locus, "signal": sig}
            for locus in LOCI
            for sig in config["metaloci"]["signals"]]


SELECTORS = {
    "wt_clones": lambda: [_clone_of(t) for t in _wt_xi_tracks()],
    "degron_clones": _dtag_pair_clones,
    "dtag_pair_clones": _dtag_pair_clones,
    "degron_nodtag_clones": lambda: [_clone_of(t)
                                     for t in _degron_xi_tracks("NodTAG")],
    "me3_tracks": lambda: list(VIZ_ME3),
    "me3_xi_tracks": lambda: list(VIZ_ME3_XI),
    "loci": lambda: list(LOCI),
    "rois": lambda: list(ROI_TYPES),
    "normgroups": lambda: list(SS.normgroups),
    "coolers": lambda: list(SS.cooler_names),
    "merged_loop_coolers": lambda: list(SS.merged_cooler_names("merged_loops")),
    "merged_comp_coolers": _comp_cooler_names,
    "enrichment_rows": _enrichment_rows,
    "fig6_composites": _fig6_composites,
    "eig_rows": _eig_rows,
    "saddle_rows": _saddle_rows,
    "pileup_rows": _pileup_rows,
    "stackup_rows": _stackup_rows,
    "density_rows": _density_rows,
    "metaloci_rows": _metaloci_rows,
    "metaloci_kk_rows": _metaloci_kk_rows,
    "metaloci_signal_rows": _metaloci_signal_rows,
}


# ---------------------------------------------------------------------
# Input functions — one per producer.
#
# These are the ONLY place a panel learns where its data lives, and every
# path comes from workflow/lib/paths.py. They are used twice: to build the
# rule's input list, and to decide whether the panel is buildable at all.
# ---------------------------------------------------------------------
def _in_valley_gene_content(panel):
    tracks = (_wt_xi_tracks() if panel.wildcards.get("cohort") == "WT"
              else _degron_xi_tracks("NodTAG") + _degron_xi_tracks("dTAG"))
    return [P.gene_content(t) for t in tracks]


def _in_boundary_ctcf_pies(panel):
    cohort = panel.wildcards.get("cohort")
    if cohort == "WT":
        tracks = _wt_xi_tracks()
    elif cohort == "degron":
        tracks = _degron_xi_tracks("NodTAG")
    else:                                   # both_sides: everything NodTAG-ish
        tracks = _wt_xi_tracks() + _degron_xi_tracks("NodTAG")
    return ([P.boundary_ctcf(t) for t in tracks]
            + [P.gene_content(t) for t in tracks])


def _in_boundary_analysis(panel):
    clone = panel.wildcards["clone"]
    which = "WT" if panel.series == "boundary_analysis_wt" else "NodTAG"
    track = _xi_track(clone, which)
    signals = _boundary_signals(track)
    # R-2: boundary_profile, NEVER stackup.
    return ([P.boundary_profile(track, s) for s in signals]
            + [P.boundary_ctcf(track), P.gene_content(track)])


def _boundary_signals(track):
    """Metaplot rows. B1621 is the Rad21 clone: its CTCF row is Rad21."""
    marks = ["H3K27me3", "H3K27ac", "RNA-Seq"]
    marks.append("Rad21" if _clone_of(track) == "B1621" else "CTCF")
    have = set(SS.marks)
    return [m for m in marks if m in have]


def _in_boundary_dtag(panel):
    clone = panel.wildcards["clone"]
    nodtag = _xi_track(clone, "NodTAG")
    signals = _boundary_signals(nodtag)
    return ([P.boundary_dtag(clone, s) for s in signals]
            + [P.boundary_ctcf(nodtag), P.gene_content(nodtag)])


def _in_loops_comps_grid(panel):
    w = panel.wildcards
    locus, roi, comparison, kind = (w["locus"], w["roi"], w["comparison"],
                                    w["kind"])
    out = []
    if kind in ("loops", "both"):
        # The number printed above each pile-up is the pile-up score, which is
        # what the on-disk `loop_stats*.tsv` of the original actually held.
        out += [P.pileup_score(comparison, roi)]
        for (cmp_, exp, sample) in _pileup_jobs():
            if cmp_ == comparison and exp.startswith(locus):
                out.append(P.pileup(comparison, roi, exp, sample))
    if kind in ("comps", "both"):
        scope = "refined_merged" if comparison == "Xa_vs_Xi" else "refined"
        out.append(P.saddle_strength_selected(scope, roi))
        for name in _comp_cooler_names():
            if name.startswith(locus):
                out.append(P.saddle(roi, name, "E1"))
    out.append(P.allelic_ratio_stats(locus, roi))
    return out


def _in_scatter_delta(panel):
    roi = panel.wildcards.get("roi", "full")
    out = [P.pileup_score("Xa_vs_Xi", roi),      # NOT loop_stats -- see above
           P.saddle_strength_selected("refined_merged", roi),
           P.ml_compartmentalization_final()]
    out += [P.allelic_ratio_stats(locus, roi) for locus in LOCI]
    return out


#: WHY THE VIOLINS DO NOT READ `loop_stats*.tsv`.
#
# Plan section 2.2 assigns `loop_stats{,_Xa,_Xi}.tsv` to `loop_strength`, the
# port of 01_03. RD-3 read 01_03 and found its writes of that filename are
# COMMENTED OUT (lines 210, 319-320). The `loop_stats*.tsv` files that exist on
# disk -- the ones Fig 3c and EFig 3a were actually made from -- were written by
# 01_04 and hold PILE-UP SCORES, a different quantity under the same name.
#
# So the paper panels read `P.pileup_score(comparison, roi)`. `P.loop_stats` is
# still produced and still compared against ground truth; it just is not what
# these figures plot. This is the kind of divergence that reproduces a
# plausible number, which is why it is spelled out here rather than fixed
# silently.
def _in_violin_loops(panel):
    w = panel.wildcards
    roi = w.get("roi", "full")
    comparison = w.get("comparison", "Xa_vs_Xi")
    return [P.pileup_score(comparison, roi)]


def _in_violin_saddle(panel):
    w = panel.wildcards
    return [P.saddle_strength_selected(w.get("scope", "refined_merged"),
                                       w.get("roi", "full"))]


def _in_density_ma(panel):
    w = panel.wildcards
    win = w.get("win", "100kb")
    masks = ("antivalley", "allcoverage")
    out = []
    for clone in _dtag_pair_clones():
        for mask in masks:
            out.append(P.density(clone, w.get("allele", "Xi"), win, mask))
    return out


def _in_valley_density(panel):
    w = panel.wildcards
    return [P.density(w["clone"], w["allele"], win, mask)
            for win in ("100kb",)
            for mask in ("antivalley", "allcoverage", "valley")]


def _in_valley_venn(panel):
    scope = panel.wildcards.get("scope", "escaping")
    return [P.valley_overlap(c, scope) for c in _dtag_pair_clones()]


def _in_valley_locus(panel):
    return [P.valley_states("chrX", panel.wildcards["track"])]


def _in_valley_sizes(panel):
    return [P.valley_sizes()]


def _in_valley_xa_xi_venn(panel):
    return [P.valley_xa_xi(panel.wildcards["track"])]


def _in_ma_plot(panel):
    return [P.scalefactors(panel.wildcards["normgroup"], "counts"),
            P.scalefactors(panel.wildcards["normgroup"], "factors")]


def _in_frip(panel):
    return [P.qc("frip.tsv")]


def _in_peak_venn(panel):
    return [P.qc("peak_overlap.tsv")]


def _in_correlation(panel):
    return [P.qc("correlation", "{}.npz".format(panel.wildcards["normgroup"]))]


def _in_fragment_sizes(panel):
    return [P.qc("fragment_sizes", "{}.tsv".format(s)) for s in SS.all_samples]


def _in_signal_enrichment(panel):
    w = panel.wildcards
    return [P.qc("enrichment", "{}_{}.npz".format(w["mark"], w["anchor"]))]


def _in_cooler_qc(panel):
    return [P.mcool("individual", panel.wildcards["cooler"])]


def _in_loop_overlay(panel):
    name = panel.wildcards["name"]
    return [P.loops_raw(name), P.loops_refined(name),
            P.mcool("merged_loops", name)]


def _in_eigenvector(panel):
    w = panel.wildcards
    return [P.eigs(w["roi"], w["name"])]


def _in_saddle(panel):
    w = panel.wildcards
    return [P.saddle(w["roi"], w["name"], w["eig"]),
            P.saddle_values(w["roi"], w["name"])]


def _in_pileup(panel):
    w = panel.wildcards
    return [P.pileup(w["comparison"], w["roi"], w["exp"], w["sample"]),
            P.pileup_score(w["comparison"], w["roi"])]


def _in_stackup_composite(panel):
    w = panel.wildcards
    return [P.stackup(w["track"], w["signal"], w["variant"])]


def _in_metaloci_gaudi(panel):
    w = panel.wildcards
    return [P.ml_signals(w["run"], w["dataset"]),
            P.ml_moran(w["run"], w["dataset"])]


def _in_metaloci_kk(panel):
    w = panel.wildcards
    return [P.ml_signals(w["run"], w["dataset"])]


def _in_metaloci_violin(panel):
    return [P.ml_compartmentalization(run, roi)
            for run in config["metaloci"]["runs"] for roi in ("full",)]


def _in_metaloci_scatter(panel):
    return [P.ml_compartmentalization_final()]


def _in_metaloci_composite(panel):
    """Fig 6: the montage of `_gtp` Gaudi plots and the violins.

    Its inputs are OTHER PANELS' outputs. That is the one place the
    visualisation stage depends on itself, and it is faithful: the original
    `02_03` notebook opened the already-written PDFs with PyMuPDF.
    """
    layout = _fig6_layout()
    locus = panel.wildcards["locus"]
    slug = panel.wildcards["locus_slug"]
    out = [P.fixture("fig6_layout.yaml"), P.ml_compartmentalization_final()]
    out += [P.ml_compartmentalization("consensus", "full")]
    for signal in layout["rows"]:
        for column in layout["columns"].get(slug, []):
            if column == "violin":
                out.append(P.figures(
                    "extras", "metaloci_violin_{}_{}".format(locus, signal),
                    "pdf"))
            else:
                run, dataset = _fig6_dataset(locus, column)
                if dataset is None:
                    continue
                out.append(P.figures(
                    "extras",
                    "metaloci_gaudi_{}_{}_{}".format(run, dataset, signal),
                    "pdf"))
    return out


def _fig6_layout():
    import yaml
    with open(P.fixture("fig6_layout.yaml")) as fh:
        return yaml.safe_load(fh)


#: The Fig 6 "Xa_consensus" column is the NodTAG-or-WT consensus pool, not any
#: Xa consensus. `metaloci_datasets("consensus")` also contains the CTCF-dTAG
#: pools, and `Mecp2_CTCF-dTAG_Xa` sorts first -- picking it would montage the
#: degron Xa into a WT figure and look entirely normal.
FIG6_CONSENSUS_GROUP = "NodTAG-or-WT"


def _fig6_dataset(locus, column):
    """Map a fixture column label onto a METALoci (run, dataset).

    The fixture speaks in figure labels (`Xa_consensus`, `E6`); METALoci speaks
    in dataset ids (`Mecp2_NodTAG-or-WT_Xa`, `Mecp2_E6_WT_G1_Xi`). This is the
    only place the two vocabularies meet, so a mismatch here silently montages
    the wrong panel into Fig 6.
    """
    if column == "Xa_consensus":
        for ds in _metaloci_datasets("consensus"):
            if (ds.startswith(locus + "_") and ds.endswith("_Xa")
                    and FIG6_CONSENSUS_GROUP in ds):
                return "consensus", ds
        return "consensus", None
    for ds in _metaloci_datasets("wt"):
        if ds.startswith("{}_{}_".format(locus, column)) and ds.endswith("_Xi"):
            return "wt", ds
    return "wt", None


def _in_coolbox_browser(panel):
    w = panel.wildcards
    locus = w["locus"]
    comparison = w.get("comparison", "dTAG_vs_NodTAG")
    out = [P.fixture("coolbox_display.yaml")]
    for name in SS.cooler_names:
        if not name.startswith(locus):
            continue
        d = naming.parse_cooler(name)
        is_degron = "TAG" in d["condition"]
        if comparison == "dTAG_vs_NodTAG" and not is_degron:
            continue
        out.append(P.mcool("individual", name))
    for t in VIZ_ME3_XI:
        out.append(P.staged_bigwig(t))
        out.append(P.staged_boundary(t))
    # Refined loops as BEDPE -- coolbox's HiCPeaksCoverage overlay.
    for name in SS.merged_cooler_names("merged_loops"):
        if name.startswith(locus):
            out.append(P.loops_refined(name, ext="bedpe"))
    # The GFFNC non-coding annotation track. 01_08 adds it at FOUR call sites
    # (lines 354, 367, 459, 472), always at TrackHeight(12) -- a fixed 12,
    # distinct from the ALLELIC GTF track's 12-for-Mecp2 / 7-for-Jarid, which
    # lives in the F-7 fixture as `gtf_track_height`. Omitting it renders a
    # plausible Fig 4e that is missing a published track, and no DAG error
    # would ever point at it.
    out.append(P.gtf_noncoding())
    # `allelic_gtf` is keyed on {clone}_{condition} (`B1_WT`, `E6A7_CTCF-dTAG`),
    # not on the bare clone -- 31_hic_features.smk calls that set AR_CLONES and
    # derives it from the verbatim column maps in hic.allelic_ratio. Asking for
    # `GRCm38.102_B1.gtf` names a file no rule declares.
    out.extend(P.allelic_gtf(ar_clone) for ar_clone in _ar_clones())
    return out


def _ar_clones():
    """Clone keys the allelic-ratio GTFs are built for (RD-3's AR_CLONES)."""
    def _fallback():
        ar = (config.get("hic") or {}).get("allelic_ratio") or {}
        return sorted(set(ar.get("wt_clone_map", {}).values())
                      | set(ar.get("degron_clone_map", {}).values()))

    return sorted(_shared("AR_CLONES", _fallback))


INPUTS = {
    "panel_valley_gene_content": _in_valley_gene_content,
    "panel_boundary_ctcf_pies": _in_boundary_ctcf_pies,
    "panel_boundary_analysis": _in_boundary_analysis,
    "panel_boundary_dtag": _in_boundary_dtag,
    "panel_loops_comps_grid": _in_loops_comps_grid,
    "panel_scatter_delta": _in_scatter_delta,
    "panel_violin_loops": _in_violin_loops,
    "panel_violin_saddle": _in_violin_saddle,
    "panel_density_ma": _in_density_ma,
    "panel_valley_density": _in_valley_density,
    "panel_valley_venn": _in_valley_venn,
    "panel_valley_locus": _in_valley_locus,
    "panel_valley_sizes": _in_valley_sizes,
    "panel_valley_xa_xi_venn": _in_valley_xa_xi_venn,
    "panel_ma_plot": _in_ma_plot,
    "panel_frip": _in_frip,
    "panel_peak_venn": _in_peak_venn,
    "panel_correlation": _in_correlation,
    "panel_fragment_sizes": _in_fragment_sizes,
    "panel_signal_enrichment": _in_signal_enrichment,
    "panel_cooler_qc": _in_cooler_qc,
    "panel_loop_overlay": _in_loop_overlay,
    "panel_eigenvector": _in_eigenvector,
    "panel_saddle": _in_saddle,
    "panel_pileup": _in_pileup,
    "panel_stackup_composite": _in_stackup_composite,
    "panel_metaloci_gaudi": _in_metaloci_gaudi,
    "panel_metaloci_kk": _in_metaloci_kk,
    "panel_metaloci_violin": _in_metaloci_violin,
    "panel_metaloci_heatmap": _in_metaloci_violin,
    "panel_metaloci_scatter": _in_metaloci_scatter,
    "panel_metaloci_composite": _in_metaloci_composite,
    "panel_coolbox_browser": _in_coolbox_browser,
}


def _panel_input_paths(panel):
    fn = INPUTS.get(panel.producer)
    if fn is None:
        raise WorkflowError(
            "panels.yaml routes {} to producer {!r}, which has no input "
            "function in 90_visualise.smk".format(panel.id, panel.producer))
    return [str(p) for p in fn(panel)]


# ---------------------------------------------------------------------
# The registry
# ---------------------------------------------------------------------
_FIGURE_ROOT = _os.path.join(P.results, "figures") + _os.sep


def _panel_output_stem(path):
    """`results/figures/{group}/{stem}.{ext}` -> `stem`, else None.

    Fig 6 montages PDFs that other panels in this module draw, so the
    visualisation stage has exactly one internal dependency. It cannot be
    resolved by asking the rule graph: this runs while 90_visualise.smk is
    still being parsed, so its own rules are not registered yet and every
    figure path would look unproducible. `panels.load` resolves these against
    the registry instead.
    """
    path = str(path)
    if not path.startswith(_FIGURE_ROOT):
        return None
    return _os.path.splitext(_os.path.basename(path))[0]


REG = _panels.load(
    REGISTRY_PATH,
    config,
    selectors=SELECTORS,
    is_available=_has_producer,
    input_paths=_panel_input_paths,
    panel_output_stem=_panel_output_stem,
)

PANEL_FILES = REG.files(P.figures)


def panel_inputs(w):
    return _panel_input_paths(REG.by_stem(w.name))


def panel_spec(w):
    """Everything the panel script needs that is not a file path."""
    panel = REG.by_stem(w.name)
    candidate = w.name.split("__", 1)[1] if "__" in w.name else None
    spec = dict(panel.wildcards or {})
    spec.update({
        "panel_id": panel.id,
        "stem": w.name,
        "series": panel.series,
        "paper_ref": list(panel.paper_ref),
        "group": panel.group,
        "producer": panel.producer,
        "confidence": panel.confidence,
        "candidate": candidate,
        "candidates": list(panel.candidates),
        "notes": panel.notes,
        "libdir": LIBDIR,
        "dpi": FIG_DPI,
        "colours": P.fixture("figure_colours.yaml"),
        "display_names": config["loci"].get("display_names", {}),
        "legacy": dict(config.get("legacy") or {}),
    })
    return spec


def _panel_out():
    return P.figures("{group}", "{name}", ext="{ext}")


def _name_re(producer):
    return REG.stem_alternation(producer)


# =====================================================================
#  THE SINGLE LOGICAL STEP
#
#  Bodyless by design: no shell:, no script:, no run:. Everything it
#  needs already exists when it fires, so it only has to record that the
#  stage completed.
# =====================================================================
rule visualise:
    """Every enabled panel, plus the manifest, the numbers and the report."""
    input:
        panels=PANEL_FILES,
        manifest=P.manifest(),
        numbers=P.paper_numbers(),
        report=P.report(),
    output:
        touch(P.visualise_done()),


# =====================================================================
#  PAPER PANELS
# =====================================================================
rule panel_valley_gene_content:
    """Fig 2g (WT) and Fig 4d (degron) — stacked gene-class bars per clone.

    TRAP P-2. Reads P.gene_content, which `valley_gene_content` builds from
    P.valleys_filtered -- the GENE-FILTERED set. The published bars are
    E6 375, C5 301, B1 258, JTG 300, CL30 321; the raw HMM calls give
    377/304/261/302/324, which is the single most likely silent divergence
    in the pipeline because the figure looks entirely correct either way.
    """
    input:
        panel_inputs,
    output:
        _panel_out(),
    params:
        spec=panel_spec,
    wildcard_constraints:
        name=_name_re("panel_valley_gene_content"),
        group=GROUP_RE,
        ext=EXT_RE,
    log:
        P.log("panel_valley_gene_content", "{group}_{name}_{ext}"),
    resources:
        mem_mb=4000,
    script:
        "../scripts/py/panels/valley_gene_content.py"


rule panel_boundary_ctcf_pies:
    """Fig 2h — CTCF status of escapee valley boundaries, per clone.

    Published n per clone: 110, 96, 74, 76, 76.
    """
    input:
        panel_inputs,
    output:
        _panel_out(),
    params:
        spec=panel_spec,
    wildcard_constraints:
        name=_name_re("panel_boundary_ctcf_pies"),
        group=GROUP_RE,
        ext=EXT_RE,
    log:
        P.log("panel_boundary_ctcf_pies", "{group}_{name}_{ext}"),
    resources:
        mem_mb=4000,
    script:
        "../scripts/py/panels/boundary_ctcf_pies.py"


rule panel_boundary_analysis:
    """EFig 2e + 2f + 2g — ONE figure, three published panels, per clone.

    TRAP R-2 / P-1. Fed from P.boundary_profile: NO background subtraction.
    P.stackup subtracts a 100-iteration random circular-shift background and
    would produce plausible curves that do not reproduce. The two rules read
    the same anchors and the same bigWigs and differ only in arithmetic.
    """
    input:
        panel_inputs,
    output:
        _panel_out(),
    params:
        spec=panel_spec,
    wildcard_constraints:
        name=_name_re("panel_boundary_analysis"),
        group=GROUP_RE,
        ext=EXT_RE,
    log:
        P.log("panel_boundary_analysis", "{group}_{name}_{ext}"),
    resources:
        mem_mb=8000,
    script:
        "../scripts/py/panels/boundary_analysis.py"


rule panel_boundary_dtag:
    """Fig 4f + Fig 4g + EFig 6b — ONE figure, three published panels.

    Rows are marks (H3K27me3 = 4f, H3K27ac = 4g); the top row of pies is
    EFig 6b. Boundaries are the NodTAG escaping-valley boundaries in both
    conditions -- a fixed position set, because the question is what happens
    to the signal there, not where the boundaries move to.
    """
    input:
        panel_inputs,
    output:
        _panel_out(),
    params:
        spec=panel_spec,
    wildcard_constraints:
        name=_name_re("panel_boundary_dtag"),
        group=GROUP_RE,
        ext=EXT_RE,
    log:
        P.log("panel_boundary_dtag", "{group}_{name}_{ext}"),
    resources:
        mem_mb=8000,
    script:
        "../scripts/py/panels/boundary_dtag.py"


rule panel_loops_comps_grid:
    """Fig 3a/3b, 5b/5c/5g/5h, EFig 8c/8d, 10c/10d — pile-ups + saddles.

    3a is Kdm5c (code name `Jarid`) and 3b is Mecp2. figure_legends.txt has
    them the other way round; figure.pdf (2026-07-25) is authoritative.
    Green circles are the RNA-seq allelic ratio under fixture F-8's
    Normalize(0.0424769728421052, 0.40909353528), despite the published
    legend calling them a METALoci D-score.
    """
    input:
        panel_inputs,
    output:
        _panel_out(),
    params:
        spec=panel_spec,
    wildcard_constraints:
        name=_name_re("panel_loops_comps_grid"),
        group=GROUP_RE,
        ext=EXT_RE,
    log:
        P.log("panel_loops_comps_grid", "{group}_{name}_{ext}"),
    resources:
        mem_mb=8000,
    script:
        "../scripts/py/panels/loops_comps_grid.py"


rule panel_scatter_delta:
    """Fig 3c — delta-compartment vs delta-loop strength. MEDIUM confidence.

    The original wrote about ten candidate scatters and none is marked as the
    published one; all four short-listed candidates are rendered.
    """
    input:
        panel_inputs,
    output:
        _panel_out(),
    params:
        spec=panel_spec,
    wildcard_constraints:
        name=_name_re("panel_scatter_delta"),
        group=GROUP_RE,
        ext=EXT_RE,
    log:
        P.log("panel_scatter_delta", "{group}_{name}_{ext}"),
    resources:
        mem_mb=4000,
    script:
        "../scripts/py/panels/scatter_delta.py"


rule panel_violin_loops:
    """EFig 3a — paired Xa/Xi loop-strength violin, N = 15 clone-alleles."""
    input:
        panel_inputs,
    output:
        _panel_out(),
    params:
        spec=panel_spec,
    wildcard_constraints:
        name=_name_re("panel_violin_loops"),
        group=GROUP_RE,
        ext=EXT_RE,
    log:
        P.log("panel_violin_loops", "{group}_{name}_{ext}"),
    resources:
        mem_mb=4000,
    script:
        "../scripts/py/panels/violin_loops.py"


rule panel_violin_saddle:
    """EFig 3b — paired Xa/Xi compartment (saddle) strength violin."""
    input:
        panel_inputs,
    output:
        _panel_out(),
    params:
        spec=panel_spec,
    wildcard_constraints:
        name=_name_re("panel_violin_saddle"),
        group=GROUP_RE,
        ext=EXT_RE,
    log:
        P.log("panel_violin_saddle", "{group}_{name}_{ext}"),
    resources:
        mem_mb=4000,
    script:
        "../scripts/py/panels/violin_saddle.py"


rule panel_density_ma:
    """EFig 6a — MA / density scatters of H3K27me3, Ctrl vs dTAG. MEDIUM."""
    input:
        panel_inputs,
    output:
        _panel_out(),
    params:
        spec=panel_spec,
    wildcard_constraints:
        name=_name_re("panel_density_ma"),
        group=GROUP_RE,
        ext=EXT_RE,
    log:
        P.log("panel_density_ma", "{group}_{name}_{ext}"),
    resources:
        mem_mb=8000,
    script:
        "../scripts/py/panels/density_ma.py"


rule panel_valley_venn:
    """EFig 6c — escapee-valley Venns, dTAG vs Ctrl. E6A7 1/44/11, F3 1/12/6."""
    input:
        panel_inputs,
    output:
        _panel_out(),
    params:
        spec=panel_spec,
    wildcard_constraints:
        name=_name_re("panel_valley_venn"),
        group=GROUP_RE,
        ext=EXT_RE,
    log:
        P.log("panel_valley_venn", "{group}_{name}_{ext}"),
    resources:
        mem_mb=2000,
    script:
        "../scripts/py/panels/valley_venn.py"


rule panel_metaloci_composite:
    """Fig 6 — the Gaudi montage.

    TRAP F-6. The grid order comes from resources/fixtures/fig6_layout.yaml
    and from nowhere else: it is hard-coded in
    "02_03_make_composite_metaloci_figures copy 2.ipynb", is not alphabetical
    and is not derivable. Wrong order = Fig 6 is silently wrong.
    """
    input:
        panel_inputs,
    output:
        _panel_out(),
    params:
        spec=panel_spec,
    wildcard_constraints:
        name=_name_re("panel_metaloci_composite"),
        group=GROUP_RE,
        ext=EXT_RE,
    log:
        P.log("panel_metaloci_composite", "{group}_{name}_{ext}"),
    resources:
        mem_mb=8000,
    script:
        "../scripts/py/panels/metaloci_composite.py"


rule panel_coolbox_browser:
    """Fig 4e (and the EFig 6d / Fig 1c-d style views).

    coolbox and the main environment cannot coexist -- see plan section 5 --
    so this is one of the two panels that run in their own conda env.
    Display ranges are fixture F-7; wrong values change the contrast of a
    published panel and nothing complains.
    """
    input:
        panel_inputs,
    output:
        _panel_out(),
    params:
        spec=panel_spec,
        display=P.fixture("coolbox_display.yaml"),
    wildcard_constraints:
        name=_name_re("panel_coolbox_browser"),
        group=GROUP_RE,
        ext=EXT_RE,
    log:
        P.log("panel_coolbox_browser", "{group}_{name}_{ext}"),
    threads: 4
    resources:
        mem_mb=16000,
    conda:
        "../envs/coolbox.yaml"
    script:
        "../scripts/py/panels/coolbox_browser.py"


# =====================================================================
#  EXTRAS (the archaeology notes F.2 / plan section 6.5)
# =====================================================================
rule panel_valley_locus:
    """Jarid / Mecp2 valley tracks: line plot, valley bins shaded red.

    Deviation D-2: 36 figures where the ground truth has 34, because we do
    not swallow the empty F3_CTCF-NodTAG_Gall valley set.
    """
    input:
        panel_inputs,
    output:
        _panel_out(),
    params:
        spec=panel_spec,
        windows=lambda w: {
            k: config["loci"][k]["view_window"] for k in LOCI
        },
    wildcard_constraints:
        name=_name_re("panel_valley_locus"),
        group=GROUP_RE,
        ext=EXT_RE,
    log:
        P.log("panel_valley_locus", "{group}_{name}_{ext}"),
    resources:
        mem_mb=4000,
    script:
        "../scripts/py/panels/valley_locus.py"


rule panel_valley_sizes:
    """Valley size distributions, one histogram per clone."""
    input:
        panel_inputs,
    output:
        _panel_out(),
    params:
        spec=panel_spec,
    wildcard_constraints:
        name=_name_re("panel_valley_sizes"),
        group=GROUP_RE,
        ext=EXT_RE,
    log:
        P.log("panel_valley_sizes", "{group}_{name}_{ext}"),
    resources:
        mem_mb=2000,
    script:
        "../scripts/py/panels/valley_sizes.py"


rule panel_valley_xa_xi_venn:
    """Xa-vs-Xi valley overlap. Deviation D-4: L/R, not left/right."""
    input:
        panel_inputs,
    output:
        _panel_out(),
    params:
        spec=panel_spec,
    wildcard_constraints:
        name=_name_re("panel_valley_xa_xi_venn"),
        group=GROUP_RE,
        ext=EXT_RE,
    log:
        P.log("panel_valley_xa_xi_venn", "{group}_{name}_{ext}"),
    resources:
        mem_mb=2000,
    script:
        "../scripts/py/panels/valley_xa_xi_venn.py"


rule panel_valley_density:
    """Density-scatter counterparts of EFig 6a, all alleles and masks."""
    input:
        panel_inputs,
    output:
        _panel_out(),
    params:
        spec=panel_spec,
    wildcard_constraints:
        name=_name_re("panel_valley_density"),
        group=GROUP_RE,
        ext=EXT_RE,
    log:
        P.log("panel_valley_density", "{group}_{name}_{ext}"),
    resources:
        mem_mb=8000,
    script:
        "../scripts/py/panels/density_ma.py"


rule panel_stackup_composite:
    """The `Stackups_*` lineage — WITH random-shift background subtraction.

    None of these reached the paper. The rule exists because that pipeline
    branch also writes the boundaries/{all,motif_yes,motif_no} sets, which
    ARE ground truth. Never route a paper panel here (correction R-2).
    Violin significance honours config.legacy.stackup_violin_test; we do not
    port the original's annotate_violin_mwu bug (D-5).
    """
    input:
        panel_inputs,
    output:
        _panel_out(),
    params:
        spec=panel_spec,
    wildcard_constraints:
        name=_name_re("panel_stackup_composite"),
        group=GROUP_RE,
        ext=EXT_RE,
    log:
        P.log("panel_stackup_composite", "{group}_{name}_{ext}"),
    resources:
        mem_mb=8000,
    script:
        "../scripts/py/panels/stackup_composite.py"


rule panel_ma_plot:
    """csaw normalisation QC. The R scripts emit counts; Python draws them."""
    input:
        panel_inputs,
    output:
        _panel_out(),
    params:
        spec=panel_spec,
    wildcard_constraints:
        name=_name_re("panel_ma_plot"),
        group=GROUP_RE,
        ext=EXT_RE,
    log:
        P.log("panel_ma_plot", "{group}_{name}_{ext}"),
    resources:
        mem_mb=8000,
    script:
        "../scripts/py/panels/qc_panels.py"


rule panel_frip:
    """Fraction of reads in peaks."""
    input:
        panel_inputs,
    output:
        _panel_out(),
    params:
        spec=panel_spec,
    wildcard_constraints:
        name=_name_re("panel_frip"),
        group=GROUP_RE,
        ext=EXT_RE,
    log:
        P.log("panel_frip", "{group}_{name}_{ext}"),
    resources:
        mem_mb=2000,
    script:
        "../scripts/py/panels/qc_panels.py"


rule panel_peak_venn:
    """Replicate peak overlaps."""
    input:
        panel_inputs,
    output:
        _panel_out(),
    params:
        spec=panel_spec,
    wildcard_constraints:
        name=_name_re("panel_peak_venn"),
        group=GROUP_RE,
        ext=EXT_RE,
    log:
        P.log("panel_peak_venn", "{group}_{name}_{ext}"),
    resources:
        mem_mb=2000,
    script:
        "../scripts/py/panels/qc_panels.py"


rule panel_correlation:
    """PCA plus Spearman and Pearson heatmaps, per normalisation group."""
    input:
        panel_inputs,
    output:
        _panel_out(),
    params:
        spec=panel_spec,
    wildcard_constraints:
        name=_name_re("panel_correlation"),
        group=GROUP_RE,
        ext=EXT_RE,
    log:
        P.log("panel_correlation", "{group}_{name}_{ext}"),
    resources:
        mem_mb=8000,
    script:
        "../scripts/py/panels/qc_panels.py"


rule panel_fragment_sizes:
    """Fragment-size distributions across every sample."""
    input:
        panel_inputs,
    output:
        _panel_out(),
    params:
        spec=panel_spec,
    wildcard_constraints:
        name=_name_re("panel_fragment_sizes"),
        group=GROUP_RE,
        ext=EXT_RE,
    log:
        P.log("panel_fragment_sizes", "{group}_{name}_{ext}"),
    resources:
        mem_mb=4000,
    script:
        "../scripts/py/panels/qc_panels.py"


rule panel_signal_enrichment:
    """Fig 1e / 2a / 2d — enrichment over gene bodies or motifs +/- 200 kb.

    Published panels, so NOT gated on qc.enabled.
    """
    input:
        panel_inputs,
    output:
        _panel_out(),
    params:
        spec=panel_spec,
    wildcard_constraints:
        name=_name_re("panel_signal_enrichment"),
        group=GROUP_RE,
        ext=EXT_RE,
    log:
        P.log("panel_signal_enrichment", "{group}_{name}_{ext}"),
    resources:
        mem_mb=8000,
    script:
        "../scripts/py/panels/qc_panels.py"


rule panel_cooler_qc:
    """Balanced contact matrix at the view window, per cooler."""
    input:
        panel_inputs,
    output:
        _panel_out(),
    params:
        spec=panel_spec,
        windows=lambda w: {k: config["loci"][k]["view_window"] for k in LOCI},
    wildcard_constraints:
        name=_name_re("panel_cooler_qc"),
        group=GROUP_RE,
        ext=EXT_RE,
    log:
        P.log("panel_cooler_qc", "{group}_{name}_{ext}"),
    resources:
        mem_mb=16000,
    script:
        "../scripts/py/panels/hic_panels.py"


rule panel_loop_overlay:
    """Loops before and after the manual refinement fixture.

    Raw 55 / 32 / 41 / 26, refined 14 / 13 / 10 / 12 -- the visual record of
    the human judgement that `resources/fixtures/loops_refined/` freezes.
    """
    input:
        panel_inputs,
    output:
        _panel_out(),
    params:
        spec=panel_spec,
        windows=lambda w: {k: config["loci"][k]["view_window"] for k in LOCI},
    wildcard_constraints:
        name=_name_re("panel_loop_overlay"),
        group=GROUP_RE,
        ext=EXT_RE,
    log:
        P.log("panel_loop_overlay", "{group}_{name}_{ext}"),
    resources:
        mem_mb=16000,
    script:
        "../scripts/py/panels/hic_panels.py"


rule panel_eigenvector:
    """All three eigenvectors on one axis.

    Read this together with work/features/compartments/eigenvector_orientation.tsv,
    which the report surfaces: a silent sign flip mirrors every saddle plot
    and would corrupt 15 of the 25 panels.
    """
    input:
        panel_inputs,
    output:
        _panel_out(),
    params:
        spec=panel_spec,
    wildcard_constraints:
        name=_name_re("panel_eigenvector"),
        group=GROUP_RE,
        ext=EXT_RE,
    log:
        P.log("panel_eigenvector", "{group}_{name}_{ext}"),
    resources:
        mem_mb=4000,
    script:
        "../scripts/py/panels/hic_panels.py"


rule panel_saddle:
    """Saddle heatmap plus the (AA+BB)/(AB+BA) strength profile, E1/E2/E3."""
    input:
        panel_inputs,
    output:
        _panel_out(),
    params:
        spec=panel_spec,
        saddle_cfg=config["hic"]["compartments"]["saddle"],
    wildcard_constraints:
        name=_name_re("panel_saddle"),
        group=GROUP_RE,
        ext=EXT_RE,
    log:
        P.log("panel_saddle", "{group}_{name}_{ext}"),
    resources:
        mem_mb=4000,
    script:
        "../scripts/py/panels/hic_panels.py"


rule panel_pileup:
    """coolpuppy pile-ups: coolwarm, log, sym, vmax 2, height 4, center 3."""
    input:
        panel_inputs,
    output:
        _panel_out(),
    params:
        spec=panel_spec,
        pileup_cfg=config["hic"]["pileup"],
    wildcard_constraints:
        name=_name_re("panel_pileup"),
        group=GROUP_RE,
        ext=EXT_RE,
    log:
        P.log("panel_pileup", "{group}_{name}_{ext}"),
    resources:
        mem_mb=8000,
    script:
        "../scripts/py/panels/hic_panels.py"


rule panel_metaloci_gaudi:
    """`metaloci figure` — the only external tool that emits images.

    Its quantitative output comes from `metaloci lm` upstream (moran_info.txt);
    the picture-maker lives in the picture layer. `_gtp` is the Gaudi TYPE
    plot (LMI quadrant colouring) -- that is the Fig 6 panel type.
    """
    input:
        panel_inputs,
    output:
        _panel_out(),
    params:
        spec=panel_spec,
        mldir=lambda w: P.ml_dir(*_gaudi_run_dataset(w.name)),
        window=lambda w: _roi_window(_gaudi_locus(w.name)),
        figure_cfg=config["metaloci"]["figure"],
    wildcard_constraints:
        name=_name_re("panel_metaloci_gaudi"),
        group=GROUP_RE,
        ext=EXT_RE,
    log:
        P.log("panel_metaloci_gaudi", "{group}_{name}_{ext}"),
    threads: 8
    resources:
        mem_mb=16000,
    conda:
        "../envs/metaloci.yaml"
    script:
        "../scripts/py/panels/metaloci_figure.py"


rule panel_metaloci_kk:
    """Kamada-Kawai and mixed-matrix plots. Same tool, same env."""
    input:
        panel_inputs,
    output:
        _panel_out(),
    params:
        spec=panel_spec,
        mldir=lambda w: P.ml_dir(*_gaudi_run_dataset(w.name)),
        window=lambda w: _roi_window(_gaudi_locus(w.name)),
        figure_cfg=config["metaloci"]["figure"],
    wildcard_constraints:
        name=_name_re("panel_metaloci_kk"),
        group=GROUP_RE,
        ext=EXT_RE,
    log:
        P.log("panel_metaloci_kk", "{group}_{name}_{ext}"),
    threads: 8
    resources:
        mem_mb=16000,
    conda:
        "../envs/metaloci.yaml"
    script:
        "../scripts/py/panels/metaloci_figure.py"


rule panel_metaloci_violin:
    """Xa-vs-Xi compartmentalization violins. Column 0 of the Fig 6 grid."""
    input:
        panel_inputs,
    output:
        _panel_out(),
    params:
        spec=panel_spec,
    wildcard_constraints:
        name=_name_re("panel_metaloci_violin"),
        group=GROUP_RE,
        ext=EXT_RE,
    log:
        P.log("panel_metaloci_violin", "{group}_{name}_{ext}"),
    resources:
        mem_mb=4000,
    script:
        "../scripts/py/panels/metaloci_panels.py"


rule panel_metaloci_heatmap:
    """The heatmap companion of panel_metaloci_violin, all six signals."""
    input:
        panel_inputs,
    output:
        _panel_out(),
    params:
        spec=panel_spec,
    wildcard_constraints:
        name=_name_re("panel_metaloci_heatmap"),
        group=GROUP_RE,
        ext=EXT_RE,
    log:
        P.log("panel_metaloci_heatmap", "{group}_{name}_{ext}"),
    resources:
        mem_mb=4000,
    script:
        "../scripts/py/panels/metaloci_panels.py"


rule panel_metaloci_scatter:
    """METALoci compartment-like strength against loop / saddle strength."""
    input:
        panel_inputs,
    output:
        _panel_out(),
    params:
        spec=panel_spec,
    wildcard_constraints:
        name=_name_re("panel_metaloci_scatter"),
        group=GROUP_RE,
        ext=EXT_RE,
    log:
        P.log("panel_metaloci_scatter", "{group}_{name}_{ext}"),
    resources:
        mem_mb=4000,
    script:
        "../scripts/py/panels/metaloci_panels.py"


def _gaudi_run_dataset(name):
    spec = REG.by_stem(name).wildcards
    return spec["run"], spec["dataset"]


def _gaudi_locus(name):
    dataset = REG.by_stem(name).wildcards["dataset"]
    for locus in LOCI:
        if dataset.startswith(locus):
            return locus
    return LOCI[0]


# =====================================================================
#  THE THREE SUMMARY ARTEFACTS
# =====================================================================
def _paper_numbers_inputs(_w=None):
    """Every table any published number is computed from, if it can be built."""
    wanted = []
    for t in _wt_xi_tracks() + _degron_xi_tracks("NodTAG") + _degron_xi_tracks("dTAG"):
        wanted += [P.gene_content(t), P.boundary_ctcf(t),
                   P.valleys("chrX", t), P.valleys_filtered("chrX", t)]
    for clone in _dtag_pair_clones():
        wanted += [P.valley_overlap(clone, "escaping"),
                   P.valley_overlap(clone, "all")]
    for roi in ("full",):
        wanted += [P.pileup_score("Xa_vs_Xi", roi),
                   P.pileup_score("dTAG_vs_NodTAG", roi),
                   P.loop_stats("Xa_vs_Xi", roi, "_Xa"),
                   P.loop_stats("Xa_vs_Xi", roi, "_Xi"),
                   P.loop_stats("dTAG_vs_NodTAG", roi),
                   P.saddle_strength_selected("refined_merged", roi),
                   P.saddle_strength_selected("refined", roi)]
    for name in SS.merged_cooler_names("merged_loops"):
        wanted += [P.loops_raw(name), P.loops_refined(name)]
    wanted.append(P.ml_compartmentalization_final())
    return sorted({p for p in wanted if _has_producer(p)})


PAPER_NUMBER_INPUTS = _paper_numbers_inputs()


rule figure_manifest:
    """panel -> paper panel, one row per emitted file. Deliverable D-07.

    Every registry entry appears, built or not: `status` is ok / empty /
    missing / skipped / blocked / disabled with the reason in `notes`, so
    "why is Fig 5c not in my results" is answerable from one file.
    """
    input:
        panels=PANEL_FILES,
    output:
        tsv=P.manifest(),
    params:
        # The EXPANDED, GATED registry, not the YAML path: the manifest has to
        # describe what the DAG actually did, and a fresh load() inside the job
        # can neither expand (the selectors live here) nor re-derive why a
        # panel was blocked (the rule graph is gone by then).
        registry=REG.to_records(),
        settings=dict(REG.settings or {}),
        spec=lambda w: {"libdir": LIBDIR, "results": P.results},
    log:
        P.log("figure_manifest"),
    resources:
        mem_mb=2000,
    script:
        "../scripts/py/panels/figure_manifest.py"


rule paper_numbers:
    """Every quantity quoted in the paper, with its value and its source.

    Includes the assertions that catch the traps: the gene-filtered Fig 2g
    totals (375/301/258/300/321, NOT 377/304/261/302/324), the raw Xi valley
    counts, the Fig 2h and EFig 6b pie n-values, the EFig 6c Venn
    cardinalities, the raw and refined loop counts, and the Wilcoxon
    statistics behind every violin.
    """
    input:
        tables=PAPER_NUMBER_INPUTS,
    output:
        tsv=P.paper_numbers(),
    params:
        registry=REGISTRY_PATH,
        spec=lambda w: {
            "libdir": LIBDIR,
            "strict": bool((REG.settings or {}).get("strict_paper_numbers", True)),
            "wt_tracks": _wt_xi_tracks(),
            "degron_nodtag_tracks": _degron_xi_tracks("NodTAG"),
            "degron_dtag_tracks": _degron_xi_tracks("dTAG"),
            "dtag_clones": _dtag_pair_clones(),
            "gene_content": {t: P.gene_content(t) for t in
                             _wt_xi_tracks() + _degron_xi_tracks("NodTAG")
                             + _degron_xi_tracks("dTAG")},
            "boundary_ctcf": {t: P.boundary_ctcf(t) for t in
                              _wt_xi_tracks() + _degron_xi_tracks("NodTAG")},
            "valleys": {t: P.valleys("chrX", t) for t in VIZ_ME3_XI},
            "valleys_filtered": {t: P.valleys_filtered("chrX", t)
                                 for t in VIZ_ME3_XI},
            "valley_overlap": {c: P.valley_overlap(c, "escaping")
                               for c in _dtag_pair_clones()},
            # Both quantities: the pile-up scores are what the paper plots,
            # `loop_stats` is the 01_03 metric kept for the comparison harness.
            "pileup_scores": {
                "Xa_vs_Xi": P.pileup_score("Xa_vs_Xi", "full"),
                "dTAG_vs_NodTAG": P.pileup_score("dTAG_vs_NodTAG", "full"),
            },
            "loop_stats": {
                "Xa": P.loop_stats("Xa_vs_Xi", "full", "_Xa"),
                "Xi": P.loop_stats("Xa_vs_Xi", "full", "_Xi"),
                "dTAG": P.loop_stats("dTAG_vs_NodTAG", "full"),
            },
            "loops_raw": {n: P.loops_raw(n)
                          for n in SS.merged_cooler_names("merged_loops")},
            "loops_refined": {n: P.loops_refined(n)
                              for n in SS.merged_cooler_names("merged_loops")},
            "saddle_strength": {
                "Xa_vs_Xi": P.saddle_strength_selected("refined_merged", "full"),
                "dTAG_vs_NodTAG": P.saddle_strength_selected("refined", "full"),
            },
            "compartmentalization": P.ml_compartmentalization_final(),
        },
    log:
        P.log("paper_numbers"),
    resources:
        mem_mb=4000,
    script:
        "../scripts/py/panels/paper_numbers.py"


rule report:
    """One self-contained HTML page with every panel inlined.

    No external assets: SVG and PNG go in as data URIs, so the file can be
    copied anywhere and still render. Sections follow the paper, and the
    eigenvector orientation table is surfaced because a silent sign flip
    there mirrors every saddle plot.
    """
    input:
        panels=PANEL_FILES,
        manifest=P.manifest(),
        numbers=P.paper_numbers(),
    output:
        html=P.report(),
    params:
        registry=REG.to_records(),
        settings=dict(REG.settings or {}),
        spec=lambda w: {
            "libdir": LIBDIR,
            "results": P.results,
            "pipeline_version": config.get("pipeline_version", "?"),
            "max_inline_bytes": int(
                (REG.settings or {}).get("report_max_inline_bytes", 6000000)),
            # RD-3 writes one orientation table per (roi, cooler) beside the
            # eigenvectors and rolls them into this one. Surfacing it is a
            # requirement, not a nicety: compartments*.json fixes WHICH
            # eigenvector but never its sign, and a silent flip mirrors every
            # saddle plot in fifteen of the twenty-five panels.
            "eigenvector_orientation": P.work(
                "features", "compartments", "eigenvector_orientation.tsv"),
            "provenance": P.provenance(),
        },
    log:
        P.log("report"),
    resources:
        mem_mb=4000,
    script:
        "../scripts/py/panels/report.py"


# ---------------------------------------------------------------------
# A one-line summary at DAG build, so a blocked panel is never a surprise.
# ---------------------------------------------------------------------
_counts = {}
for _p in REG:
    _counts[_p.status] = _counts.get(_p.status, 0) + 1
logger.info(
    "visualise: {} panels ({}) -> {} files".format(
        len(REG),
        ", ".join("{} {}".format(v, k) for k, v in sorted(_counts.items())),
        len(PANEL_FILES),
    )
)
