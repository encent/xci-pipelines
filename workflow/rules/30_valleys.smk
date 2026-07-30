# =====================================================================
#  30_valleys.smk — LAYER 2a
#
#  H3K27me3 valley calling and everything built on it. This is a leaf on
#  the signal branch, not a peer of it: it consumes exactly one artefact
#  class (merged5k H3K27me3) and contributes exactly one back (valley
#  BEDs, staged for coolbox).
#
#  It is also the primary exact-equality gate for testing -- 10 of the 25
#  target panels come from here, and the valley BEDs are compared byte for
#  byte. A 2-state HMM's output is discrete, so a small change in the input
#  bigWig can flip an entire valley rather than perturb it; the valley test
#  is therefore gated behind a passing bigWig test.
#
#  ------------------------------------------------------------------
#  THE TWO BOUNDARY-METAPLOT CODE PATHS ARE NOT INTERCHANGEABLE (R-2)
#  ------------------------------------------------------------------
#  `boundary_stackup` and `boundary_profile` take the SAME inputs -- Xi
#  valley boundaries and merged20 bigWigs -- and compute DIFFERENT
#  arithmetic. Routing a panel to the wrong one produces plausible-looking
#  curves that silently do not reproduce.
#
#      boundary_stackup      Stackups_*.ipynb lineage.
#                            APPLIES a 100-iteration circular random-shift
#                            background: numpy.random.default_rng(
#                            seed = clone_index + 42), averaged and
#                            subtracted. Feeds NO paper panel. It stays in
#                            the pipeline because it writes
#                            boundaries/{all,motif_yes,motif_no}, which ARE
#                            ground truth, and because it produces the
#                            non-paper composite extras.
#
#      boundary_profile      xci_valleys_check_LAST.ipynb lineage.
#      boundary_dtag_profile NO background subtraction, at all.
#                            These feed EFig 2e/2f/2g, Fig 4f/4g, EFig 6b.
#
#  ------------------------------------------------------------------
#  RAW vs GENE-FILTERED VALLEYS
#  ------------------------------------------------------------------
#  `call_valleys` emits the raw HMM set; `filter_valleys_genes` removes
#  valleys touching Mid1 / Tmem29 / Firre. Both are kept, and which one a
#  rule reads is load-bearing (correction P-2):
#
#      RAW       valley_boundaries, boundary_motif_split -- those files are
#                ground truth from a lineage that predates the filter.
#      FILTERED  valley_gene_content and everything gene-related below it:
#                boundary_ctcf_status, valley_overlap_dtag,
#                valley_xa_xi_overlap, valley_coverage_density.
#
#  The `chrX:0-3,285,000` call is the mm10 assembly gap, visible to the
#  HMM only because it sits outside the first blacklist interval (which
#  starts at 3,286,700) and `fillna(0.0)` makes it look like the low
#  state. It is a genuine artefact -- but it overlaps none of the three
#  excluded genes, so it survives the filter and IS counted in the
#  published Fig 2g bars (ruling O-7). Reproduce it;
#  `valleys.drop_leading_gap` exists to revisit that on future data.
#
#  ------------------------------------------------------------------
#  ACCEPTANCE NUMBERS
#  ------------------------------------------------------------------
#  Raw Xi valleys:       B1 261, C5 304, CL30 324, E6 377, JTG 302,
#                        B1621-NodTAG 284, B1621-dTAG 281,
#                        E6A7-NodTAG 321, E6A7-dTAG 281,
#                        F3-NodTAG 351, F3-dTAG 205.
#  Gene-FILTERED Xi valleys -- what Fig 2g actually plots:
#                        E6 375, C5 301, B1 258, JTG 300, CL30 321.
#  Escaping boundaries (Fig 2h pies):  110, 96, 74, 76, 76.
#  Rows in boundaries/all/ (raw set):
#                        B1 522, C5 608, CL30 648, E6 754, JTG 604,
#                        B1621-NodTAG 568, E6A7-NodTAG 642, F3-NodTAG 702.
#
#  Expected deviations: D-08 (18 valley BEDs where ground truth has 17, 36
#  locus figures where it has 34) and D-4 (side vocabulary unified on L/R).
#  Both are documented in the module notes.
#
#  Replaces: TEST_03_00_H3K27me3_boundary_calling.ipynb,
#            TEST_03_01_* / Stackups_* notebooks,
#            xci_valleys_check_LAST.ipynb, xci_valley_overlap.ipynb,
#            xci_density_plots.ipynb, TEST_03_02_valley_sizes.ipynb
# =====================================================================

VAL = config["valleys"]
BOUND = VAL["boundaries"]
STACK = VAL["stackups"]
BANAL = VAL["boundary_analysis"]

ME3_TRACKS = SS.tracks_of(mark="H3K27me3")
ME3_XI_TRACKS = SS.tracks_of(mark="H3K27me3", allele="Xi")

# Wildcards this module introduces. `track` and `signal` are not in the global
# set (Snakefile / naming.WC) because other modules use `track` with a wider
# vocabulary; constraining them per rule keeps an ambiguous match a DAG-build
# error rather than an hour-30 surprise.
TRACK_WC = r"[A-Za-z0-9\-]+_[A-Za-z0-9\.]+_[A-Za-z0-9\-]+_(Gall|Xa|Xi)"
SIGNAL_WC = "|".join(naming.ALL_MARKS)


# ---------------------------------------------------------------------
# CLONE ALIAS MAP -- the fixtures do not speak our grammar
#
# The external deliveries under resources/fixtures/ were named by the
# collaborators and use clone aliases that do NOT match the pipeline's
# {mark}_{clone}_{condition}_{allele} grammar. They are enumerated here, once,
# and never re-derived at a call site. Note in particular:
#
#   * the CTCF peak files call CL30 "CL30.7" (the subclone id), so a naive
#     f"NPC_{clone}" builds a filename that does not exist;
#   * the escapee/silent tables call the 4-day dTAG condition "4DdTAG" while
#     the CTCF peak files call the same condition "dTAG";
#   * WT clones carry no condition token in either delivery.
#
# Verified against the fixture filenames actually on disk. A missing key is a
# loud error, not a guessed filename.
# ---------------------------------------------------------------------
CTCF_PEAK_ALIAS = {
    ("B1", "WT"): "NPC_B1",
    ("C5", "WT"): "NPC_C5",
    ("CL30", "WT"): "NPC_CL30.7",
    ("E6", "WT"): "NPC_E6",
    ("JTG", "WT"): "NPC_JTG",
    ("B1621", "Rad21-NodTAG"): "NPC_B1621_NodTAG",
    ("B1621", "Rad21-dTAG"): "NPC_B1621_dTAG",
    ("E6A7", "CTCF-NodTAG"): "NPC_E6A7_NodTAG",
    ("E6A7", "CTCF-dTAG"): "NPC_E6A7_dTAG",
    ("F3", "CTCF-NodTAG"): "NPC_F3_NodTAG",
    ("F3", "CTCF-dTAG"): "NPC_F3_dTAG",
}

GENE_TABLE_ALIAS = {
    ("B1", "WT"): "B1",
    ("C5", "WT"): "C5",
    ("CL30", "WT"): "CL30",
    ("E6", "WT"): "E6",
    ("JTG", "WT"): "JTG",
    ("B1621", "Rad21-NodTAG"): "B1621_NodTAG",
    ("B1621", "Rad21-dTAG"): "B1621_4DdTAG",
    ("E6A7", "CTCF-NodTAG"): "E6A7_NodTAG",
    ("E6A7", "CTCF-dTAG"): "E6A7_4DdTAG",
    ("F3", "CTCF-NodTAG"): "F3_NodTAG",
    ("F3", "CTCF-dTAG"): "F3_4DdTAG",
}

# Recorded because the originals used them, NOT used as a filename anywhere:
# these are column keys inside the allelic-ratio CSVs, which this module does
# not read. Kept here so nobody re-invents "F3_0h" as a BED stem.
ALLELIC_RATIO_COLUMN_ALIAS = {
    ("B1621", "Rad21-NodTAG"): "Rad21B1621_nodTAG",
    ("E6A7", "CTCF-NodTAG"): "E6A7_NodTAG",
    ("F3", "CTCF-NodTAG"): "F3_0h",
    ("F3", "CTCF-dTAG"): "F3_4d",
    ("CL30", "WT"): "CL30.7",
}


def _alias(table, track, what):
    d = naming.parse_track(track)
    key = (d["clone"], d["condition"])
    if key not in table:
        raise WorkflowError(
            f"No {what} fixture alias for clone={d['clone']} "
            f"condition={d['condition']} (track {track}).\n"
            f"Known keys: {sorted(table)}\n"
            "Add the alias to the map at the top of "
            "workflow/rules/30_valleys.smk -- never guess the filename at the "
            "call site."
        )
    return table[key]


def _ctcf_peak_file(track):
    """Absolute path to the per-clone CTCF consensus-peak fixture for `track`."""
    stem = _alias(CTCF_PEAK_ALIAS, track, "CTCF peak")
    return P.fixture(
        "CTCFpeak_per_clone", f"{stem}_consensusPeaks_withMeanRatioandDirection.bed"
    )


def _escapee_file(track):
    stem = _alias(GENE_TABLE_ALIAS, track, "escapee")
    return P.fixture("escapees_per_clone", f"{stem}_escapees_overlap.bed")


def _silent_file(track):
    stem = _alias(GENE_TABLE_ALIAS, track, "silent")
    return P.fixture("silent_per_clone", f"{stem}_escapees_no_overlap.bed")


def _signal_track(track, signal):
    """Companion track of another mark: swap the leading token, then verify.

    Rad21 exists for B1621 only and CTCF for the other degron clones, so an
    unchecked `swap_mark` produces a path no rule can build and Snakemake
    reports it as a missing input file rather than as a bad request.
    """
    companion = naming.swap_mark(track, signal)
    if companion not in SS.tracks_of(mark=signal):
        raise WorkflowError(
            f"{signal} was requested at the boundaries of {track}, but "
            f"{companion} is not in config/samples.tsv.\n"
            "Ask for a signal this clone actually has."
        )
    return companion


def _clone_track(clone, which):
    """The Xi H3K27me3 track of `clone` in the dTAG or NodTAG condition."""
    for t in ME3_XI_TRACKS:
        d = naming.parse_track(t)
        if d["clone"] != clone:
            continue
        cond = d["condition"]
        if which == "dTAG" and cond.endswith("-dTAG"):
            return t
        if which == "NodTAG" and cond.endswith("-NodTAG"):
            return t
    raise WorkflowError(f"no {which} Xi H3K27me3 track for clone {clone!r}")


def _track_with_allele(track, allele):
    d = naming.parse_track(track)
    return naming.track_id(d["mark"], d["clone"], d["condition"], allele)


def _clone_track_allele(clone, which, allele):
    """The `allele` H3K27me3 track of `clone` in the dTAG or NodTAG condition."""
    return _track_with_allele(_clone_track(clone, which), allele)


def _clone_index(track):
    """The original's `k` -- the clone's position in its generation's file list.

    `Stackups_*` seeded the random-shift background with
    `np.random.default_rng(seed=k+42)`, where k enumerated a sorted glob of Xa/Xi
    bigWig PAIRS. The two generations were globbed separately, so WT clones are
    numbered 0..4 and degron clones 0..2 -- k is NOT unique across the pipeline
    and must not be replaced by a global index, or every stackup background
    changes. Reconstructed here as the position of this track's clone in the
    sorted clone list of its own generation.
    """
    d = naming.parse_track(track)
    generation = [
        naming.parse_track(t)["clone"]
        for t in ME3_XI_TRACKS
        if (naming.parse_track(t)["condition"] == "WT") == (d["condition"] == "WT")
    ]
    ordered = sorted(set(generation))
    return ordered.index(d["clone"])


def _degron_clones():
    """Clones having both a dTAG and a NodTAG Xi H3K27me3 track."""
    out = []
    for clone in SS.clones:
        conds = {
            naming.parse_track(t)["condition"]
            for t in ME3_XI_TRACKS
            if naming.parse_track(t)["clone"] == clone
        }
        has_dtag = any(c.endswith("-dTAG") for c in conds)
        has_nodtag = any(c.endswith("-NodTAG") for c in conds)
        if has_dtag and has_nodtag:
            out.append(clone)
    return sorted(out)


DEGRON_CLONES = _degron_clones()
# An empty alternation would match everything; `$^` matches nothing instead.
DEGRON_WC = "|".join(DEGRON_CLONES) if DEGRON_CLONES else "$^"



# ---------------------------------------------------------------------
# Calling
# ---------------------------------------------------------------------
rule call_valleys:
    """2-state GaussianHMM over 5 kb H3K27me3 signal. THE exact-equality target.

    Fitted on chromosome-minus-blacklist intervals of >= 500 kb, concatenated
    with lengths=; Viterbi restarted inside each interval; bins outside every
    training interval keep state -1 and can never be valleys. No smoothing, no
    signal transform, no threshold, no minimum valley length.

    Writes BED3 -- chrom/start/end, no name, no score, no strand. This is the
    file the testing team byte-compares. `valley_boundaries` writes BED4; the
    two formats differ deliberately, so do not "fix" either to match the other.

    D-08: the original wrapped `bioframe.merge` in a bare try/except that
    printed "No valleys found!" and continued, which is why
    H3K27me3_F3_CTCF-NodTAG_Gall_valleys.* is absent from the ground truth. We
    do NOT swallow it -- an empty valley set produces a valid empty BED, an
    empty bigWig and a loud warning. A correct run therefore produces 18 degron
    valley BEDs where ground truth has 17, and 36 locus figures where ground
    truth has 34. That is expected, not a failure.

    Emits the per-bin state table alongside the BED/bigWig. That table is not
    an original output -- it exists so that when a valley boundary moves, the
    testing team can localise it to a bin instead of diffing two BED files and
    guessing.
    """
    input:
        bw=P.merged5k("H3K27me3", "{track}"),
        chrom_sizes=P.chrom_sizes(),
        blacklist=P.blacklist(),
    output:
        bed=P.valleys("{chrom}", "{track}"),
        bw=P.valleys("{chrom}", "{track}", ext="bw"),
        states=P.valley_states("{chrom}", "{track}"),
    params:
        resolution=VAL["resolution"],
        hmm=VAL["hmm"],
        valley_state=VAL["valley_state"],
        merge_gap=VAL["merge_gap"],
        drop_leading_gap=VAL["drop_leading_gap"],
        blas_threads=VAL["blas_threads"],
        region_size_threshold=VAL["region_size_threshold"],
    wildcard_constraints:
        track=TRACK_WC,
    log:
        P.log("call_valleys", "{chrom}_{track}"),
    benchmark:
        P.benchmark("call_valleys", "{chrom}_{track}")
    threads: 1
    resources:
        mem_mb=8000,
    script:
        "../scripts/py/call_valleys.py"


# ---------------------------------------------------------------------
# Fig 4d valley source (D-18)
#
# Fig 2g reproduces FROM COMPUTATION -- 375/301/258/300/321, exact -- so it uses
# no fixture. Fig 4d does not: our four CTCF-degron Xi tracks come out at
# 293/291/343/204 against the published 321/281/351/205, a spread of more than
# an order of magnitude (F3-dTAG misses by ONE valley; E6A7-NodTAG by 28).
#
# ONLY those four tracks may be served from the fixture. The other seven Xi
# tracks -- B1, C5, CL30, E6, JTG and BOTH B1621 Rad21 conditions -- reproduce
# count-exactly, and a fixture there would mask a real regression.
# ---------------------------------------------------------------------
FIG4D_FIXTURE_TRACKS = (
    "H3K27me3_E6A7_CTCF-NodTAG_Xi",
    "H3K27me3_E6A7_CTCF-dTAG_Xi",
    "H3K27me3_F3_CTCF-NodTAG_Xi",
    "H3K27me3_F3_CTCF-dTAG_Xi",
)


def _valley_source(chrom, track):
    """The valley BED the gene-content filter should read for this track.

    Falls back to computation for everything except the four Fig 4d tracks, and
    only when `valleys.fig4d_source: fixture`. The pipeline's own call is
    produced either way, so the comparison is always available.
    """
    if (
        chrom == "chrX"
        and track in FIG4D_FIXTURE_TRACKS
        and VAL.get("fig4d_source", "fixture") == "fixture"
    ):
        return P.valleys_fig4d_fixture(track)
    return P.valleys(chrom, track)


rule filter_valleys_genes:
    """The 3-gene filter behind Fig 2g / 4d.

    The paper's bar totals are the FILTERED counts (E6 375, C5 301, B1 258,
    JTG 300, CL30 321), not the raw HMM calls (377 / 304 / 261 / 302 / 324).
    Both sets are emitted; the panels draw from this one.
    """
    input:
        bed=lambda w: _valley_source(w.chrom, w.track),
        # always built, so the fixture and the computation can be compared even
        # when the fixture is the one being filtered
        computed=P.valleys("{chrom}", "{track}"),
    output:
        bed=P.valleys_filtered("{chrom}", "{track}"),
    params:
        exclude=VAL["exclude_genes"],
        source=lambda w: (
            "fixture" if _valley_source(w.chrom, w.track) != P.valleys(w.chrom, w.track)
            else "computed"
        ),
    wildcard_constraints:
        track=TRACK_WC,
    log:
        P.log("filter_valleys_genes", "{chrom}_{track}"),
    benchmark:
        P.benchmark("filter_valleys_genes", "{chrom}_{track}")
    threads: 1
    resources:
        mem_mb=2000,
    script:
        "../scripts/py/filter_valleys_genes.py"


# ---------------------------------------------------------------------
# Boundaries
# ---------------------------------------------------------------------
rule valley_boundaries:
    """Valley edges as points, labelled L/R.

    The original used two different flip conventions across notebooks; the
    unified one is by side (`legacy.stackup_flip_rule: side`). `parity`
    reproduces the plain variant.
    """
    input:
        bed=P.valleys("chrX", "{track}"),
    output:
        bed=P.boundaries("all", "{track}"),
    params:
        side_vocabulary=BOUND["side_vocabulary"],
        flip_rule=config["legacy"]["stackup_flip_rule"],
    wildcard_constraints:
        track=TRACK_WC,
    log:
        P.log("valley_boundaries", "{track}"),
    benchmark:
        P.benchmark("valley_boundaries", "{track}")
    threads: 1
    resources:
        mem_mb=2000,
    script:
        "../scripts/py/valley_boundaries.py"


rule boundary_motif_split:
    """Split boundaries by presence of a CTCF motif within `motif_distance`."""
    input:
        bed=P.boundaries("all", "{track}"),
        motifs=P.ctcf_motifs(chrx_only=True),
    output:
        yes=P.boundaries("motif_yes", "{track}"),
        no=P.boundaries("motif_no", "{track}"),
    params:
        distance=BOUND["motif_distance"],
    wildcard_constraints:
        track=TRACK_WC,
    log:
        P.log("boundary_motif_split", "{track}"),
    benchmark:
        P.benchmark("boundary_motif_split", "{track}")
    threads: 1
    resources:
        mem_mb=4000,
    script:
        "../scripts/py/boundary_motif_split.py"


rule boundary_ctcf_status:
    """CTCF-peak status per boundary -- the input to Fig 2h and EFig 6b.

    Uses the per-clone CTCF consensus-peak delivery (a fixture), not MACS2
    output. The asymmetric flanks are the original's: a peak counts if it lies
    within 50 kb on the valley-interior side or 10 kb on the exterior side.

    Boundaries are re-derived from `gene_content` rather than read from
    `boundaries/all/`, because the paper panels need each boundary's gene class
    AND the 3-gene filter applied. `boundaries/all/` is the pre-filter ground
    truth from the `Stackups_*` lineage and carries no class; using it here
    would put the boundaries of the three excluded valleys back into the Fig 2h
    denominators. Both sets exist on purpose -- see the module notes.
    """
    input:
        gene_content=P.gene_content("{track}"),
        peaks=lambda w: _ctcf_peak_file(w.track),
    output:
        tsv=P.boundary_ctcf("{track}"),
    params:
        inside=BOUND["ctcf_inside_flank"],
        outside=BOUND["ctcf_outside_flank"],
    wildcard_constraints:
        track=TRACK_WC,
    log:
        P.log("boundary_ctcf_status", "{track}"),
    benchmark:
        P.benchmark("boundary_ctcf_status", "{track}")
    threads: 1
    resources:
        mem_mb=4000,
    script:
        "../scripts/py/boundary_ctcf_status.py"


# ---------------------------------------------------------------------
# Gene content and overlaps
# ---------------------------------------------------------------------
rule valley_gene_content:
    """Gene-class composition per valley set — the data behind Fig 2g / 4d."""
    input:
        bed=P.valleys_filtered("chrX", "{track}"),
        gtf=P.gtf(),
        escapees=lambda w: _escapee_file(w.track),
        silent=lambda w: _silent_file(w.track),
    output:
        tsv=P.gene_content("{track}"),
    params:
        classes=VAL["gene_classes"],
    wildcard_constraints:
        track=TRACK_WC,
    log:
        P.log("valley_gene_content", "{track}"),
    benchmark:
        P.benchmark("valley_gene_content", "{track}")
    threads: 2
    resources:
        mem_mb=12000,
    script:
        "../scripts/py/valley_gene_content.py"


rule valley_overlap_dtag:
    """Valley overlap between dTAG and NodTAG for one clone — EFig 6c Venns.

    Each condition is classified with ITS OWN escapee delivery
    (`E6A7_NodTAG` vs `E6A7_4DdTAG`), because escape status is what the degron
    changes -- classifying both sides with the NodTAG list would silently answer
    a different question. `xci_valley_overlap.ipynb` cell 6 does the same.
    """
    input:
        nodtag=lambda w: P.valleys_filtered("chrX", _clone_track(w.clone, "NodTAG")),
        dtag=lambda w: P.valleys_filtered("chrX", _clone_track(w.clone, "dTAG")),
        escapees=lambda w: _escapee_file(_clone_track(w.clone, "NodTAG")),
        escapees_dtag=lambda w: _escapee_file(_clone_track(w.clone, "dTAG")),
    output:
        tsv=P.valley_overlap("{clone}", "{scope}"),
    params:
        scope=lambda w: w.scope,
    wildcard_constraints:
        clone=DEGRON_WC,
        scope="all|escaping",
    log:
        P.log("valley_overlap_dtag", "{clone}_{scope}"),
    benchmark:
        P.benchmark("valley_overlap_dtag", "{clone}_{scope}")
    threads: 1
    resources:
        mem_mb=2000,
    script:
        "../scripts/py/valley_overlap.py"


rule valley_xa_xi_overlap:
    """Xa-vs-Xi valley overlap (secondary analysis, not in the paper).

    Reads the **RAW** valley set, not the gene-filtered one.

    This is the same pre-filter lineage as `valley_boundaries` and
    `boundary_motif_split`: `TEST_03_01_..._final_degrons_overlap.ipynb`
    predates the 3-gene filter, which was introduced later in
    `xci_valleys_check_LAST.ipynb`. Correction P-2 says the pre-filter
    generation reads raw and only the gene-content lineage reads filtered; that
    argument was applied to the other two rules and missed here.

    Measured on B1621-NodTAG, after normalising the original's `left`/`right`
    to `L`/`R` per deviation D-4:

        RAW       231 / 831 / 337   <- exact match to ground truth
        FILTERED  230 / 822 / 334   <- 13 boundaries short

    and every one of the 13 missing boundaries falls inside Firre or Mid1 --
    precisely the genes the filter removes. That is not a tolerance question,
    it is the wrong input.
    """
    input:
        xa=lambda w: P.valleys("chrX", _track_with_allele(w.track, "Xa")),
        xi=lambda w: P.valleys("chrX", _track_with_allele(w.track, "Xi")),
    output:
        tsv=P.valley_xa_xi("{track}"),
    wildcard_constraints:
        track=TRACK_WC,
    log:
        P.log("valley_xa_xi_overlap", "{track}"),
    benchmark:
        P.benchmark("valley_xa_xi_overlap", "{track}")
    threads: 1
    resources:
        mem_mb=2000,
    script:
        "../scripts/py/valley_xa_xi_overlap.py"


rule valley_sizes:
    """Valley size distribution across every called set."""
    input:
        beds=[P.valleys("chrX", t) for t in ME3_TRACKS],
    output:
        tsv=P.valley_sizes(),
    log:
        P.log("valley_sizes"),
    benchmark:
        P.benchmark("valley_sizes")
    threads: 1
    resources:
        mem_mb=2000,
    script:
        "../scripts/py/valley_sizes.py"


# ---------------------------------------------------------------------
# Signal around boundaries
#
# TWO DISTINCT RULES. See correction R-2 -- do not merge them.
# ---------------------------------------------------------------------
rule boundary_stackup:
    """Signal matrix around boundaries, WITH random-shift background subtraction.

    100 circular shifts, seeded `background_seed_base + k`. This is the
    `Stackups_*` lineage. It does NOT feed EFig 2e/2f.
    """
    input:
        boundaries=lambda w: P.boundaries(w.variant, w.track),
        valleys=P.valleys("chrX", "{track}"),
        bw=lambda w: P.merged20(w.signal, _signal_track(w.track, w.signal)),
        chrom_sizes=P.chrom_sizes(),
    output:
        npz=P.stackup("{track}", "{signal}", "{variant}"),
    params:
        flank=STACK["flank"],
        nbins=STACK["nbins"],
        shifts=STACK["background_shifts"],
        seed_base=STACK["background_seed_base"],
        clone_index=lambda w: _clone_index(w.track),
        aggregate=STACK["aggregate"],
        subtract_background=True,
    wildcard_constraints:
        track=TRACK_WC,
        signal=SIGNAL_WC,
    log:
        P.log("boundary_stackup", "{track}_{signal}_{variant}"),
    benchmark:
        P.benchmark("boundary_stackup", "{track}_{signal}_{variant}")
    threads: 8
    resources:
        mem_mb=16000,
    script:
        "../scripts/py/boundary_signal_matrix.py"


rule boundary_profile:
    """Signal matrix around boundaries, WITHOUT background subtraction.

    This is what EFig 2e/2f actually plot (`xci_valleys_check_LAST.ipynb`).
    Identical inputs to `boundary_stackup`, different arithmetic -- routing a
    panel to the wrong one of these produces a plausible, wrong figure.
    """
    input:
        boundaries=P.boundaries("all", "{track}"),
        valleys=P.valleys("chrX", "{track}"),
        bw=lambda w: P.merged20(w.signal, _signal_track(w.track, w.signal)),
        chrom_sizes=P.chrom_sizes(),
    output:
        npz=P.boundary_profile("{track}", "{signal}"),
    params:
        flank=BANAL["flank"],
        nbins=BANAL["nbins"],
        shifts=0,
        seed_base=STACK["background_seed_base"],
        clone_index=lambda w: _clone_index(w.track),
        aggregate=STACK["aggregate"],
        subtract_background=False,
    wildcard_constraints:
        track=TRACK_WC,
        signal=SIGNAL_WC,
    log:
        P.log("boundary_profile", "{track}_{signal}"),
    benchmark:
        P.benchmark("boundary_profile", "{track}_{signal}")
    threads: 8
    resources:
        mem_mb=16000,
    script:
        "../scripts/py/boundary_signal_matrix.py"


rule boundary_dtag_profile:
    """Signal at NodTAG boundaries, compared between dTAG and NodTAG.

    Behind Fig 4f/4g and the EFig 6b pies. Boundaries always come from the
    NodTAG condition -- the comparison asks what happens to the signal at a
    fixed set of positions, so re-calling boundaries per condition would
    change the question.
    """
    input:
        # `ctcf` already carries the NodTAG boundary coordinates, their gene
        # class and their CTCF status, so it IS the boundary set here -- there
        # is deliberately no second, unfiltered boundary input that could
        # disagree with it.
        ctcf=lambda w: P.boundary_ctcf(_clone_track(w.clone, "NodTAG")),
        chrom_sizes=P.chrom_sizes(),
        nodtag_bw=lambda w: P.merged20(
            w.signal, _signal_track(_clone_track(w.clone, "NodTAG"), w.signal)
        ),
        dtag_bw=lambda w: P.merged20(
            w.signal, _signal_track(_clone_track(w.clone, "dTAG"), w.signal)
        ),
    output:
        npz=P.boundary_dtag("{clone}", "{signal}"),
    params:
        flank=BANAL["flank"],
        nbins=BANAL["nbins"],
    wildcard_constraints:
        clone=DEGRON_WC,
        signal=SIGNAL_WC,
    log:
        P.log("boundary_dtag_profile", "{clone}_{signal}"),
    benchmark:
        P.benchmark("boundary_dtag_profile", "{clone}_{signal}")
    threads: 8
    resources:
        mem_mb=16000,
    script:
        "../scripts/py/boundary_dtag_profile.py"


rule valley_coverage_density:
    """Binned H3K27me3 coverage inside / outside valleys — EFig 6a MA plots.

    A NodTAG-vs-dTAG comparison over one fixed bin partition, so it takes BOTH
    conditions of a degron clone and is defined only for degron clones. The bins
    come from the NodTAG valley call, gene-FILTERED, which is the set
    `xci_density_plots.ipynb` reads (`VALLEY_DIR = .../valleys_filtered`).

    The autosome control is part of the figure, not a QC afterthought: every
    density/MA panel has three columns, chr7 / chrX / both. So the valley input
    is a LIST -- the filtered chrX set plus the raw set of each
    `chromosomes.valley_control` chromosome (the 3-gene filter is a no-op off
    chrX, and no filtered control set is produced). Rows carry their chrom, and
    the statistics are computed per chromosome, so the panel slices whichever
    column it wants.
    """
    input:
        nodtag_bw=lambda w: P.merged20(
            "H3K27me3", _clone_track_allele(w.clone, "NodTAG", w.allele)
        ),
        dtag_bw=lambda w: P.merged20(
            "H3K27me3", _clone_track_allele(w.clone, "dTAG", w.allele)
        ),
        valleys=lambda w: [
            P.valleys_filtered("chrX", _clone_track_allele(w.clone, "NodTAG", w.allele))
        ]
        + [
            P.valleys(c, _clone_track_allele(w.clone, "NodTAG", w.allele))
            for c in config["chromosomes"].get("valley_control", [])
        ],
        chrom_sizes=P.chrom_sizes(),
    output:
        tsv=P.density("{clone}", "{allele}", "{win}", "{mask}"),
    params:
        window=lambda w: {"10kb": 10000, "100kb": 100000}[w.win],
        mask=lambda w: w.mask,
    wildcard_constraints:
        clone=DEGRON_WC,
        win="10kb|100kb",
        mask="valley|antivalley|allcoverage",
    log:
        P.log("valley_coverage_density", "{clone}_{allele}_{win}_{mask}"),
    benchmark:
        P.benchmark("valley_coverage_density", "{clone}_{allele}_{win}_{mask}")
    threads: 4
    resources:
        mem_mb=12000,
    script:
        "../scripts/py/valley_coverage_density.py"


# ---------------------------------------------------------------------
# WHAT THIS MODULE PRODUCES — the `valleys` convenience target
#
# Assembled here rather than in the Snakefile so that the eligibility rules
# (which clone has which mark, which track has a fixture alias) live next to
# the rules that depend on them. `rule valleys` in the Snakefile is a thin
# alias over VALLEY_TARGETS.
# ---------------------------------------------------------------------
GENE_CONTENT_TRACKS = [
    t
    for t in ME3_XI_TRACKS
    if (naming.parse_track(t)["clone"], naming.parse_track(t)["condition"])
    in GENE_TABLE_ALIAS
]


def _boundary_signals(track):
    """The marks this clone actually has, in the original's plotting order.

    Rad21 exists for B1621 only and CTCF for every other clone, so the signal
    list is per track, not global.
    """
    out = []
    for signal in ("H3K27me3", "H3K27ac", "CTCF", "Rad21", "RNA-Seq"):
        if naming.swap_mark(track, signal) in SS.tracks_of(mark=signal):
            out.append(signal)
    return out


def _valley_targets():
    out = []
    # calling, both chromosomes
    for chrom in VALLEY_CHROMS:
        for track in ME3_TRACKS:
            out.append(P.valleys(chrom, track))
            out.append(P.valley_states(chrom, track))
    # the gene filter, chrX only (that is the only chromosome the paper uses)
    out += [P.valleys_filtered("chrX", t) for t in ME3_TRACKS]
    # boundaries and their motif split, Xi only -- matching the ground truth
    for track in ME3_XI_TRACKS:
        out += [P.boundaries(v, track) for v in ("all", "motif_yes", "motif_no")]
    # gene content and CTCF status, Xi tracks that have the fixture deliveries
    out += [P.gene_content(t) for t in GENE_CONTENT_TRACKS]
    out += [P.boundary_ctcf(t) for t in GENE_CONTENT_TRACKS]
    # Xa-vs-Xi overlap (secondary), and the size table
    out += [P.valley_xa_xi(t) for t in ME3_XI_TRACKS]
    out.append(P.valley_sizes())
    # the two boundary-metaplot paths, kept visibly separate (R-2)
    for track in ME3_XI_TRACKS:
        for signal in _boundary_signals(track):
            out.append(P.boundary_profile(track, signal))          # no background
            for variant in ("all", "motif_yes", "motif_no"):
                out.append(P.stackup(track, signal, variant))      # background
    # degron-only products
    for clone in DEGRON_CLONES:
        out += [P.valley_overlap(clone, s) for s in ("all", "escaping")]
        nodtag = _clone_track(clone, "NodTAG")
        for signal in _boundary_signals(nodtag):
            if naming.swap_mark(_clone_track(clone, "dTAG"), signal) in SS.tracks_of(
                mark=signal
            ):
                out.append(P.boundary_dtag(clone, signal))
        for allele in ("Xa", "Xi"):
            for win in ("10kb", "100kb"):
                for mask in ("valley", "antivalley", "allcoverage"):
                    out.append(P.density(clone, allele, win, mask))
    return sorted(set(out))


VALLEY_TARGETS = _valley_targets()
