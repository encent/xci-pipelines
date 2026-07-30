import re

# =====================================================================
#  31_hic_features.smk — LAYER 2b
#
#  mcools -> loops, loop strengths, pile-ups, eigenvectors, saddles,
#  allelic-ratio GTFs. Fifteen of the twenty-five target panels are drawn
#  from data produced here.
#
#  NO FIGURES. The original 01_02 / 01_03 / 01_04 / 01_05 / 01_06 wrote
#  SVG and PNG directly -- saddle plots, pile-up images, eigenvector
#  panels, loop overlays. Every one of those is now a `panel_*` rule in
#  90_visualise.smk drawing from the `.npz` / `.tsv` this module emits.
#  The original's plot arguments are preserved verbatim, as comments, next
#  to the rule that produces the matrix.
#
#  Four things here are not what they look like:
#
#  * **Refined loops are frozen input; raw loops are reproducible.**
#    chromosight's raw calls (55 / 32 / 41 / 26) come back from code. The
#    published sets (14 / 13 / 10 / 12) are those minus hand-deleted false
#    positives, and that judgement lives only in the files. `paper` mode
#    uses the fixtures and fails loudly without them; `auto` mode falls
#    through to raw and says so in writing.
#
#  * **Two different "loop strength" numbers exist and they are not
#    interchangeable.** `loop_strength` is 01_03: per loop, the mean of the
#    central 3x3 observed/expected pixels, compared with a paired Wilcoxon.
#    `loop_pileup` is 01_04: one coolpuppy pile-up per sample, scored on
#    the central 3x3 of the AGGREGATE. The `loop_stats*.tsv` that exist in
#    the ground truth -- the ones the paper's violins and scatters read --
#    are the PILE-UP scores; 01_03's writes of that same filename are
#    commented out. Here they are separate files in separate directories so
#    they can never be confused again. Section 4 of
#    the module notes says which panel reads which.
#
#  * **Eigenvector sign is a real reproduction risk (plan 9.5, R-3).**
#    `cooltools.eigs_cis` returns an arbitrary global sign, and the
#    `compartments*.json` fixture records only WHICH eigenvector, never its
#    ORIENTATION. We canonicalise to GC-rich-positive and record every flip
#    in `eigenvector_orientation.tsv`, because an unrecorded flip mirrors
#    the saddle plot and looks perfectly plausible.
#
#  * **D-9: `loop_anchor_stackup` reads `merged20`**, the csaw-normalised
#    tracks, where the original read `data/bigwigs/merged_bw/` (deepTools
#    means of the originally delivered bigWigs, never csaw-normalised).
#    Scope, so nobody panics: this touches ONLY the stackups. `loop_pileup`
#    reads coolers, so every paper pile-up -- Fig 3a/b, 5b/g, EFig 8c,
#    EFig 10c -- is unaffected, and the 01_04 stackups never reached the
#    paper at all.
#
#  Replaces: 01_00_loops_chromosight.py, 01_02_loops_chromosight_viz.py
#            (data half), 01_03_loops_chromosight_stats.py,
#            01_04_loops_pileups_and_stackups.py,
#            01_05_compartments_cooltools.py, 01_06_compartments_refine.py,
#            01_07_get_allelic_ratio_gtf.py
# =====================================================================

import os as _os

CS = HIC["chromosight"]
LS = HIC["loop_strength"]
PU = HIC["pileup"]
COMP = HIC["compartments"]
STK = HIC["stackup"]
AR = HIC["allelic_ratio"]

COMPARISONS = ["Xa_vs_Xi", "dTAG_vs_NodTAG"]
CONSENSUS = "consensus_loops"
PILEUP_SETS = COMPARISONS + [CONSENSUS]
EIGS = [f"E{i + 1}" for i in range(COMP["n_eigs"])]

# Compartments are called on two universes, mirroring the original's two runs:
# every individual cooler, and every merged_for_compartments pool.
COMP_SCOPES = {"refined": SS.cooler_names, "refined_merged": MERGED_COMPS}
COMP_NAMES = sorted(set(SS.cooler_names) | set(MERGED_COMPS))
COMP_FIXTURE = {"refined": "compartments_{roi}.json",
                "refined_merged": "compartments_merged_{roi}.json"}

AR_CLONES = sorted(set(AR["wt_clone_map"].values()) | set(AR["degron_clone_map"].values()))
# The loop sets that were ever HAND-CURATED, and therefore the only ones with a
# refined fixture: {locus}_{allele}, i.e. Jarid_Xa, Jarid_Xi, Mecp2_Xa, Mecp2_Xi.
#
# `MERGED_LOOPS` is wider -- it also holds the pooled `*_dTAG_*` merges, and
# `call_loops` DOES run on all 8, because the original did. But the original
# never curated the dTAG sets: $AE/results/loops_chromosight_merged_refined/ has
# exactly these four subdirectories, the ground-truth inventory section 9 lists exactly these
# four, and nothing downstream reads a dTAG refined set -- `loop_strength` and
# the original's `01_04` both consume only {locus}_{Xa|Xi}.
#
# EVERYTHING BELOW THIS LINE THAT CONCERNS REFINED LOOPS MUST USE `LOOP_SETS`.
# A bare `MERGED_LOOPS` below here is a bug on sight: `stage_refined_loops`,
# `loops_bedpe` and two convenience targets all reached past this narrowing to
# the wider set, which made `hic` and `hic_features` unsatisfiable in the
# shipped `fixtures.mode: paper` default -- four jobs demanding fixtures that
# do not and should not exist.
LOOP_SETS = [n for n in MERGED_LOOPS if n.count("_") == 1]


# ---------------------------------------------------------------------
# Experiment enumeration
#
# The original derived its experiment lists by string surgery on mcool
# filenames. Here they come off the sample sheet, but the ids are
# byte-identical to the original's so the output tree still matches.
# ---------------------------------------------------------------------
def _rep_suffix(row):
    return f"_rep{row['replicate']}" if row["replicate"] else ""


def _dtag_experiments():
    """`{locus}_{clone}_{tag}_{allele}[_rep{n}]` -> (dTAG cooler, NodTAG cooler).

    Reproduces `'_'.join(x.split('_')[:2] + x.split('_')[3:])`: the condition
    token is dropped, so the two members of a pair collapse onto one id.
    Pairing is on (locus, clone, snpsplit_tag, allele, replicate), which is
    what that string surgery was approximating.
    """
    out = {}
    for _, d in SS.coolers.iterrows():
        if not d["condition"].endswith("-dTAG"):
            continue
        partner = SS.coolers[
            (SS.coolers["locus"] == d["locus"])
            & (SS.coolers["clone"] == d["clone"])
            & (SS.coolers["snpsplit_tag"] == d["snpsplit_tag"])
            & (SS.coolers["allele"] == d["allele"])
            & (SS.coolers["replicate"] == d["replicate"])
            & (SS.coolers["condition"].str.endswith("-NodTAG"))
        ]
        if partner.empty:
            continue
        exp = (f"{d['locus']}_{d['clone']}_{d['snpsplit_tag']}_{d['allele']}"
               f"{_rep_suffix(d)}")
        out[exp] = (d["cooler"], partner.iloc[0]["cooler"])
    return dict(sorted(out.items()))


def _xaxi_experiments():
    """`{locus}_{clone}_{condition}[_rep{n}]` -> (Xa cooler, Xi cooler)."""
    out = {}
    for _, d in SS.coolers.iterrows():
        if d["allele"] != "Xa":
            continue
        partner = SS.coolers[
            (SS.coolers["locus"] == d["locus"])
            & (SS.coolers["clone"] == d["clone"])
            & (SS.coolers["condition"] == d["condition"])
            & (SS.coolers["replicate"] == d["replicate"])
            & (SS.coolers["allele"] == "Xi")
        ]
        if partner.empty:
            continue
        out[f"{d['locus']}_{d['clone']}_{d['condition']}{_rep_suffix(d)}"] = (
            d["cooler"], partner.iloc[0]["cooler"])
    return dict(sorted(out.items()))


DTAG_EXP = _dtag_experiments()
XAXI_EXP = _xaxi_experiments()


def _pileup_jobs():
    """(comparison, exp, sample) -> (mcool scope, cooler, loop set).

    The loop set is always a merged WT/NodTAG cooler -- the original never
    piles up per-clone loop calls. For `Xa_vs_Xi` the sample id carries the
    allele the loops were CALLED on, reproducing `pileup_{Xa|Xi}_{cooler}`.
    """
    jobs = {}
    for exp, (dtag, nodtag) in DTAG_EXP.items():
        locus, allele = exp.split("_")[0], exp.split("_")[3]
        for c in (dtag, nodtag):
            jobs[("dTAG_vs_NodTAG", exp, c)] = ("individual", c, f"{locus}_{allele}")
    for exp, (xa, xi) in XAXI_EXP.items():
        locus = exp.split("_")[0]
        for called_on in ("Xa", "Xi"):
            for c in (xa, xi):
                jobs[("Xa_vs_Xi", exp, f"{called_on}_{c}")] = (
                    "individual", c, f"{locus}_{called_on}")
    for locus in LOCI:
        for allele in ("Xa", "Xi"):
            jobs[(CONSENSUS, locus, f"consensus_{locus}_{allele}")] = (
                "merged_loops", f"{locus}_{allele}", f"{locus}_{allele}")
            jobs[(CONSENSUS, locus, f"consensus_{locus}_CTCF-dTAG_{allele}")] = (
                "merged_comps", f"{locus}_CTCF-dTAG_{allele}", f"{locus}_{allele}")
    return {k: v for k, v in sorted(jobs.items()) if v[2] in MERGED_LOOPS}


PILEUP_JOBS = _pileup_jobs()


def _stackup_jobs():
    """(comparison, signal, pair) -> (loop set, track1, track2).

    A job exists only when BOTH tracks are in the sample sheet, reproducing the
    original's `if os.path.exists(bw1) and os.path.exists(bw2)`.
    """
    have = {m: set(SS.tracks_of(mark=m)) for m in naming.MARKS}
    have["AcMe3"] = set(SS.acme3_tracks())
    jobs = {}

    def add(comparison, signal, pair, loopset, t1, t2):
        if loopset in MERGED_LOOPS and t1 in have.get(signal, ()) and t2 in have.get(signal, ()):
            jobs[(comparison, signal, pair)] = (loopset, t1, t2)

    for signal in STK["signals"]:
        for exp, (dtag, nodtag) in DTAG_EXP.items():
            a1, a2 = naming.dataset_adj(dtag), naming.dataset_adj(nodtag)
            locus, allele = exp.split("_")[0], exp.split("_")[3]
            add("dTAG_vs_NodTAG", signal, f"{a1}_{a2}", f"{locus}_{allele}",
                f"{signal}_{a1}", f"{signal}_{a2}")
        for exp, (xa, xi) in XAXI_EXP.items():
            a1, a2 = naming.dataset_adj(xa), naming.dataset_adj(xi)
            locus = exp.split("_")[0]
            for called_on in ("Xa", "Xi"):
                add("Xa_vs_Xi", signal, f"{called_on}_{a1}_{a2}",
                    f"{locus}_{called_on}", f"{signal}_{a1}", f"{signal}_{a2}")
    return dict(sorted(jobs.items()))


STACKUP_JOBS = _stackup_jobs()


def _signal_bigwig(signal, track):
    return P.acme3_bw20(track) if signal == "AcMe3" else P.merged20(signal, track)


def _comp_mcool(name):
    """Individual coolers and merged_comps pools share one flat name space."""
    return P.mcool("individual" if name in SS.cooler_names else "merged_comps", name)


def _orientation_tsv(roi, name):
    return _os.path.join(_os.path.dirname(P.eigs(roi, name)), "orientation.tsv")


# ---------------------------------------------------------------------
# Loop calling
# ---------------------------------------------------------------------
rule call_loops:
    """chromosight loop detection on a merged cooler at 5 kb (01_00).

    Reproducible, and an exact-equality test gate: the raw call counts are
    Jarid_Xa 55, Jarid_Xi 32, Mecp2_Xa 41, Mecp2_Xi 26.
    `--min-separation 10000` is `2 * resolution` in the original.
    """
    input:
        mcool=P.mcool("merged_loops", "{name}"),
    output:
        tsv=P.loops_raw("{name}"),
        json=P.loops_raw("{name}", ext="json"),
    params:
        prefix=lambda w: P.loops_raw(w.name)[: -len(".tsv")],
        resolution=HIC["resolution"],
        pattern=CS["pattern"],
        norm=CS["norm"],
        min_dist=CS["min_dist"],
        max_dist=CS["max_dist"],
        pearson=CS["pearson"],
        perc_undetected=CS["perc_undetected"],
        perc_zero=CS["perc_zero"],
        n_mads=CS["n_mads"],
        min_separation=CS["min_separation"],
    wildcard_constraints:
        name=_alt(MERGED_LOOPS),
    log:
        P.log("call_loops", "{name}"),
    benchmark:
        P.benchmark("call_loops", "{name}")
    threads: 40
    resources:
        mem_mb=32000,
    shell:
        r"""
        mkdir -p "$(dirname {output.tsv})" "$(dirname {log})"
        chromosight detect \
            --threads {threads} \
            --norm {params.norm} \
            --pattern {params.pattern} \
            --min-dist {params.min_dist} \
            --max-dist {params.max_dist} \
            --pearson {params.pearson} \
            --perc-undetected {params.perc_undetected} \
            --perc-zero {params.perc_zero} \
            --n-mads {params.n_mads} \
            --min-separation {params.min_separation} \
            --no-plotting \
            {input.mcool}::resolutions/{params.resolution} \
            {params.prefix} > {log} 2>&1
        """


rule stage_refined_loops:
    """The manual false-positive deletion, as a versioned input (D-05).

    `results/loops_chromosight_merged_refined/` was made by copying the raw
    calls and deleting rows by eye in a PDF viewer. Not derivable from code, so
    `paper` mode ships it: 14 / 13 / 10 / 12 rows for Jarid_Xa / Jarid_Xi /
    Mecp2_Xa / Mecp2_Xi. A merged cooler with no fixture is a loud error
    carrying a real curation recipe, never a silent fallthrough.

    `auto` mode uses the raw calls unrefined and writes the refinement note
    telling the user which overlays to inspect and what to edit.
    """
    input:
        raw=P.loops_raw("{name}"),
    output:
        tsv=P.loops_refined("{name}"),
        note=P.loops_refined("{name}", ext="refinement.md"),
    params:
        mode=FIX["mode"],
        fixture=lambda w: P.fixture("loops_refined", f"{w.name}.tsv"),
    wildcard_constraints:
        name=_alt(LOOP_SETS),
    log:
        P.log("stage_refined_loops", "{name}"),
    benchmark:
        P.benchmark("stage_refined_loops", "{name}")
    script:
        "../scripts/py/stage_refined_loops.py"


rule loops_bedpe:
    """Refined loops as BEDPE, bgzipped and pairix-indexed (01_02, data half).

    Six columns, no header -- the format coolbox's `HiCPeaksCoverage` reads.
    The `.px2` index is what makes it random-access at plot time; the original
    produced it outside any script.
    """
    input:
        tsv=P.loops_refined("{name}"),
    output:
        bedpe=P.loops_refined("{name}", ext="bedpe"),
        bgz=P.loops_refined("{name}", ext="bedpe.bgz"),
        px2=P.loops_refined("{name}", ext="bedpe.bgz.px2"),
    wildcard_constraints:
        name=_alt(LOOP_SETS),
    log:
        P.log("loops_bedpe", "{name}"),
    benchmark:
        P.benchmark("loops_bedpe", "{name}")
    resources:
        mem_mb=2000,
    shell:
        r"""
        mkdir -p "$(dirname {output.bedpe})" "$(dirname {log})"
        {{
            tail -n +2 {input.tsv} | cut -f1-6 > {output.bedpe}
            sort -k1,1 -k2,2n -k4,4 -k5,5n {output.bedpe} | bgzip -c > {output.bgz}
            pairix -f -p bedpe {output.bgz}
        }} > {log} 2>&1
        """


# ---------------------------------------------------------------------
# Loop strength  (01_03)
# ---------------------------------------------------------------------
def _loop_strength_inputs(wildcards):
    pairs = DTAG_EXP if wildcards.comparison == "dTAG_vs_NodTAG" else XAXI_EXP
    coolers = sorted({c for pair in pairs.values() for c in pair})
    return {
        "mcools": [P.mcool("individual", c) for c in coolers],
        "loops": [P.loops_refined(n) for n in LOOP_SETS],
        "roi_table": P.roi_table(),
    }


rule loop_strength:
    """Per-loop observed/expected contact frequency and its paired test (01_03).

    For each loop: `expected_cis(clr, nproc=N, ignore_diags=0, intra_only=True,
    view_df=<balanced span>, chunksize=1e6)`, divide `balanced` by
    `balanced.avg` at that distance, then `np.nanmean` of the 3x3 pixel block
    centred on the loop. Groups are compared with a Wilcoxon signed-rank test
    paired by loop index -- two-sided and both one-sided, as the original
    printed into its violin titles.

    `{suffix}` is "" for dTAG_vs_NodTAG and `_Xa` / `_Xi` for Xa_vs_Xi, naming
    which allele's loop set was used. The per-loop long table, which the
    original only ever held in memory, is emitted alongside the summary.

    NOT the same number as `loop_pileup`'s score -- see the module docstring.
    """
    input:
        unpack(_loop_strength_inputs),
    output:
        stats=P.loop_stats("{comparison}", "{roi}", "{suffix}"),
        per_loop=P.loop_stats("{comparison}", "{roi}", "{suffix}").replace(
            ".tsv", "_per_loop.tsv"),
        wilcoxon=P.loop_stats("{comparison}", "{roi}", "{suffix}").replace(
            ".tsv", "_wilcoxon.tsv"),
    params:
        experiments=lambda w: DTAG_EXP if w.comparison == "dTAG_vs_NodTAG" else XAXI_EXP,
        center=LS["center"],
        expected=LS["expected"],
        resolution=HIC["resolution"],
    wildcard_constraints:
        comparison="|".join(COMPARISONS),
        suffix="|_Xa|_Xi",
    log:
        P.log("loop_strength", "{comparison}_{roi}{suffix}"),
    benchmark:
        P.benchmark("loop_strength", "{comparison}_{roi}{suffix}")
    threads: 40
    resources:
        mem_mb=32000,
    script:
        "../scripts/py/loop_strength.py"


# ---------------------------------------------------------------------
# Pile-ups  (01_04, pile-up half)
# ---------------------------------------------------------------------
def _pileup_spec(wildcards):
    key = (wildcards.comparison, wildcards.exp, wildcards.sample)
    if key not in PILEUP_JOBS:
        raise WorkflowError(
            f"no pile-up job {key}; valid samples for "
            f"{wildcards.comparison}/{wildcards.exp} are "
            + ", ".join(s for c, e, s in PILEUP_JOBS
                        if (c, e) == (wildcards.comparison, wildcards.exp))
        )
    return PILEUP_JOBS[key]


rule loop_pileup:
    """coolpuppy pile-up of one cooler over one refined loop set (01_04).

    Parameters, all verified against the original call:
        features_format='bedpe', view_df=<full balanced span>, flank=100_000,
        min_diag=0, nproc=N, seed=42, store_stripes=False,
        expected_df=expected_cis(ignore_diags=0, intra_only=True), nshifts=0,
        minshift=1e5, maxshift=1e6, clr_weight_name='weight', local=False,
        trans=False, by_strand=False, by_distance=False, mindist='auto',
        maxdist=None, flip_negative_strand=False
    Score = `np.nanmean` of the central 3x3; median and std are stored too,
    because the original recorded all three.

    `view_df` is the span of BALANCED bins -- first to last bin with a non-null
    weight -- not the ROI. The ROI only filters which loops enter.

    NO PLOT. `panel_pileup` / `panel_loops_comps_grid` draw it. The original's
    arguments, preserved verbatim for that stage:
        plotpup.plot(pup, score=True, cmap='coolwarm', scale='log', sym=True,
                     vmax=2, vmin=None, height=4, plot_ticks=True, center=3)
    """
    input:
        mcool=lambda w: P.mcool(_pileup_spec(w)[0], _pileup_spec(w)[1]),
        loops=lambda w: P.loops_refined(_pileup_spec(w)[2]),
        roi_table=P.roi_table(),
    output:
        npz=P.pileup("{comparison}", "{roi}", "{exp}", "{sample}"),
        score=P.pileup("{comparison}", "{roi}", "{exp}", "{sample}").replace(
            ".npz", ".score.tsv"),
    params:
        flank=PU["flank"],
        min_diag=PU["min_diag"],
        nshifts=PU["nshifts"],
        minshift=PU["minshift"],
        maxshift=PU["maxshift"],
        seed=PU["seed"],
        score=PU["score"],
        expected=LS["expected"],
        resolution=HIC["resolution"],
        cooler=lambda w: _pileup_spec(w)[1],
        locus=lambda w: _pileup_spec(w)[2].split("_")[0],
    wildcard_constraints:
        comparison="|".join(PILEUP_SETS),
        exp=r"[^/]+",
        sample=r"[^/]+",
    log:
        P.log("loop_pileup", "{comparison}_{roi}_{exp}_{sample}"),
    benchmark:
        P.benchmark("loop_pileup", "{comparison}_{roi}_{exp}_{sample}")
    threads: 40
    resources:
        mem_mb=40000,
    script:
        "../scripts/py/loop_pileup.py"


rule pileup_scores:
    """Collect every pile-up score of one comparison + ROI into one table.

    ADDED beyond the plan's twelve rules: `P.pileup_score()` is a
    per-(comparison, roi) path while `loop_pileup` is per-sample, so something
    has to join them. Columns are the ground truth's -- `file, median, mean,
    std` -- which makes THIS file, not `loops/stats/loop_stats*.tsv`, the
    drop-in replacement for the `loop_stats*.tsv` that 01_04 actually wrote and
    that Fig 3c and EFig 3a read.
    """
    input:
        lambda w: [
            P.pileup(w.comparison, w.roi, exp, sample).replace(".npz", ".score.tsv")
            for (comparison, exp, sample) in PILEUP_JOBS
            if comparison == w.comparison
        ],
    output:
        P.pileup_score("{comparison}", "{roi}"),
    wildcard_constraints:
        comparison="|".join(PILEUP_SETS),
    log:
        P.log("pileup_scores", "{comparison}_{roi}"),
    benchmark:
        P.benchmark("pileup_scores", "{comparison}_{roi}")
    resources:
        mem_mb=2000,
    shell:
        r"""
        mkdir -p "$(dirname {output})" "$(dirname {log})"
        {{
            head -n 1 {input[0]}
            for f in {input}; do tail -n +2 "$f"; done
        }} > {output} 2> {log}
        """


# ---------------------------------------------------------------------
# Loop-anchor stackups  (01_04, stackup half)  — DECLARED DEVIATION D-9
# ---------------------------------------------------------------------
def _stackup_spec(wildcards):
    key = (wildcards.comparison, wildcards.signal, wildcards.pair)
    if key not in STACKUP_JOBS:
        raise WorkflowError(f"no stackup job {key}")
    return STACKUP_JOBS[key]


rule loop_anchor_stackup:
    """Signal around both feet of every loop (`bbi.stackup`, 01_04).

    `bbi.stackup(bw, chrom, mid - 100_000, mid + 100_000, bins=100)` over the
    concatenation of (chrom1, start1, end1) and (chrom2, start2, end2),
    aggregated across anchors with `np.nanmean`.

    **D-9.** The original read `data/bigwigs/merged_bw/` -- deepTools means of
    the bigWigs as originally delivered, never csaw/TMM-normalised. This rule
    reads `merged20`, the csaw-normalised tracks, on the archaeologist's
    recommendation (the archaeology notes section J.5).

    Scope, so nobody panics: this is the ONLY rule that changes. `loop_pileup`
    reads coolers, so every paper pile-up (Fig 3a/b, 5b/g, EFig 8c, EFig 10c)
    is unaffected, and the 01_04 stackups never reached the paper. Structural
    test only; no ground-truth diff is expected to pass.
    """
    input:
        loops=lambda w: P.loops_refined(_stackup_spec(w)[0]),
        bw1=lambda w: _signal_bigwig(w.signal, _stackup_spec(w)[1]),
        bw2=lambda w: _signal_bigwig(w.signal, _stackup_spec(w)[2]),
        roi_table=P.roi_table(),
    output:
        npz=P.loop_anchor_stackup("{comparison}", "{roi}", "{signal}", "{pair}"),
    params:
        flank=STK["flank"],
        nbins=STK["nbins"],
        tracks=lambda w: list(_stackup_spec(w)[1:]),
        locus=lambda w: _stackup_spec(w)[0].split("_")[0],
    wildcard_constraints:
        comparison="|".join(COMPARISONS),
        signal=naming.WC["mark"],
        pair=r"[^/]+",
    log:
        P.log("loop_anchor_stackup", "{comparison}_{roi}_{signal}_{pair}"),
    benchmark:
        P.benchmark("loop_anchor_stackup", "{comparison}_{roi}_{signal}_{pair}")
    threads: 8
    resources:
        mem_mb=16000,
    script:
        "../scripts/py/loop_anchor_stackup.py"


# ---------------------------------------------------------------------
# Compartments  (01_05 / 01_06)
# ---------------------------------------------------------------------
# ---------------------------------------------------------------------
# GC phasing track
#
# `bioframe.frac_gc` sees only the bin table and the FASTA -- never the contact
# matrix, never the ROI. So every cooler sharing a bin table has an IDENTICAL
# GC track. For this dataset that is all 204 of them (chrX only, 5 kb): one
# distinct file, md5 4ed6660643f07d805e74f4308c8623a3, matching ground truth.
#
# The original recomputed it 204 times, reloading a 2.8 GB FASTA each time.
# Same defect class as `index_bam`: recomputing what already exists, invisible
# to a correctness test because the answer is right -- just paid for 204 times.
#
# Computed once per bin table, then reused with an ASSERTION. The identity is a
# property of this dataset, not of the rule, so a blind copy would be silently
# wrong on coolers at another resolution or spanning several chromosomes.
# `gc_track` therefore checks each consumer's bin table against the canonical
# one and fails loudly on a mismatch.
# ---------------------------------------------------------------------
def _gc_bin_key():
    """Identifies a bin table: what `frac_gc` actually depends on.

    Deliberately NOT a cooler name. Naming the shared file after one arbitrary
    cooler would hide that the sharing is a property of the BINNING.
    """
    chroms = "-".join(config["chromosomes"]["downstream"])
    return f"{chroms}_{HIC['resolution']}"


GC_BIN_KEY = _gc_bin_key()
# Representative cooler: a source of bins only, never part of the output name.
GC_REPRESENTATIVE = sorted(COMP_NAMES)[0]


rule gc_track_canonical:
    """GC fraction per bin, computed ONCE for this bin table (01_05)."""
    input:
        mcool=_comp_mcool(GC_REPRESENTATIVE),
        fasta=P.fasta(),
        fai=P.fasta() + ".fai",
    output:
        tsv=P.gc_track_canonical("{binkey}"),
        bins=P.gc_track_canonical("{binkey}") + ".bins",
    params:
        resolution=HIC["resolution"],
    wildcard_constraints:
        binkey=re.escape(GC_BIN_KEY),
    log:
        P.log("gc_track_canonical", "{binkey}"),
    benchmark:
        P.benchmark("gc_track_canonical", "{binkey}")
    threads: 4
    resources:
        mem_mb=16000,
    script:
        "../scripts/py/gc_track_canonical.py"


rule gc_track:
    """Per-(cooler, ROI) GC track: assert the bin table, then reuse (01_05).

    The `{roi}` in the path is the original's directory layout; the CONTENT is
    roi-independent by construction, and the rule says so rather than leaving a
    reader to wonder.
    """
    input:
        mcool=lambda w: _comp_mcool(w.name),
        gc=P.gc_track_canonical(GC_BIN_KEY),
        bins=P.gc_track_canonical(GC_BIN_KEY) + ".bins",
    output:
        tsv=P.gc_track("{roi}", "{name}"),
    params:
        resolution=HIC["resolution"],
    wildcard_constraints:
        name=_alt(COMP_NAMES),
    log:
        P.log("gc_track", "{roi}_{name}"),
    benchmark:
        P.benchmark("gc_track", "{roi}_{name}")
    threads: 1
    resources:
        mem_mb=4000,
    script:
        "../scripts/py/gc_track.py"


rule eigs_cis:
    """First three cis eigenvectors over one ROI (01_05).

    `cooltools.eigs_cis(clr, gc_cov, view_df=view, n_eigs=3)` where `view` is a
    SINGLE-ROW DataFrame `chrX:<roi.start>-<roi.end>` named `chrX`. There is no
    arms file: the "view" is the capture window itself.

    Plan 9.5 lives here. cooltools returns an arbitrary global sign per
    eigenvector; with `hic.compartments.canonicalise_sign: true` each E is
    oriented so Spearman(E, GC) > 0 and the flip is written to
    `orientation.tsv` beside the track. The correlations are recorded even when
    nothing flips, because they are also what `auto` mode selects on.
    """
    input:
        mcool=lambda w: _comp_mcool(w.name),
        gc=P.gc_track("{roi}", "{name}"),
        roi_table=P.roi_table(),
        chrom_sizes=P.chrom_sizes(),
    output:
        tsv=P.eigs("{roi}", "{name}"),
        bw=[P.eigs("{roi}", "{name}", which=e) for e in EIGS],
        orientation=_orientation_tsv("{roi}", "{name}"),
    params:
        n_eigs=COMP["n_eigs"],
        resolution=HIC["resolution"],
        canonicalise=COMP["canonicalise_sign"],
        eigs=EIGS,
        locus=lambda w: w.name.split("_")[0],
    wildcard_constraints:
        name=_alt(COMP_NAMES),
    log:
        P.log("eigs_cis", "{roi}_{name}"),
    benchmark:
        P.benchmark("eigs_cis", "{roi}_{name}")
    threads: 8
    resources:
        mem_mb=24000,
    script:
        "../scripts/py/eigs_cis.py"


rule saddle:
    """Saddle matrix and strength for all three eigenvectors (01_05).

    `cooltools.saddle(clr, cvd, eig_track, 'cis', n_bins=38,
    qrange=(0.025, 0.975), view_df)` on `cooltools.expected_cis(clr, view_df)`.

    TWO strength numbers are stored, and they are different quantities:
      * `Saddle_values_{name}.tsv` -- the original's own score,
            (mean(S[:8,:8]) + mean(S[32:,32:])) / (mean(S[:8,32:]) + mean(S[32:,:8]))
        This is what `refine_compartments` selects from, and what EFig 3b plots.
      * `saddle_strength_{E}.tsv` -- `cooltools.api.saddle.saddle_strength(...)[8]`,
        which the original computed and then commented out of its TSV.

    All three eigenvectors in one job, because `Saddle_values_{name}.tsv` is a
    single row with an E1/E2/E3 column each.

    NO PLOT. `panel_saddle` draws it from the `.npz`.
    """
    input:
        mcool=lambda w: _comp_mcool(w.name),
        eigs=P.eigs("{roi}", "{name}"),
        roi_table=P.roi_table(),
    output:
        npz=[P.saddle("{roi}", "{name}", e) for e in EIGS],
        values=P.saddle_values("{roi}", "{name}"),
        strength=[
            _os.path.join(_os.path.dirname(P.saddle_values("{roi}", "{name}")),
                          f"saddle_strength_{e}.tsv")
            for e in EIGS
        ],
    params:
        n_bins=COMP["saddle"]["n_bins"],
        qrange=COMP["saddle"]["qrange"],
        extent=COMP["saddle"]["strength_extent"],
        eigs=EIGS,
        resolution=HIC["resolution"],
        locus=lambda w: w.name.split("_")[0],
    wildcard_constraints:
        name=_alt(COMP_NAMES),
    log:
        P.log("saddle", "{roi}_{name}"),
    benchmark:
        P.benchmark("saddle", "{roi}_{name}")
    threads: 8
    resources:
        mem_mb=24000,
    script:
        "../scripts/py/saddle.py"


rule refine_compartments:
    """Select the curated eigenvector for one cooler and stage its outputs (01_06).

    `compartments{,_merged}_{roi}.json` -- 6 committed files, e.g.
    `"Mecp2_NodTAG-or-WT_Xi": "E3"` -- names the eigenvector a human chose after
    looking at the maps. `01_06` then copied that eigenvector's bigWig and
    saddle out into a `refined` tree.

    **The JSON records WHICH eigenvector, never its ORIENTATION** (plan 9.5).
    The sign canonicalisation applied by `eigs_cis` is carried through into
    `selection.tsv` so the choice and the flip are visible in one place.

    Under `eigenvector_selection: auto` the choice is max |Spearman(E, GC)| and
    the three correlations are recorded next to it. The selected bigWig always
    gets the SAME stable filename, `Comp_selected_{name}.bw`, in both modes --
    the ground truth's `Comp_E2_{name}.bw` embeds a choice that only exists in
    `paper` mode, and a DAG cannot depend on a filename it will not know until
    the job runs.

    `{cscope}` is `refined` (per-clone coolers) or `refined_merged` (the
    merged_for_compartments pools); the wildcard is not called `scope` because
    that name is globally constrained to the mcool scopes.
    """
    input:
        eigs=P.eigs("{roi}", "{name}"),
        orientation=_orientation_tsv("{roi}", "{name}"),
        values=P.saddle_values("{roi}", "{name}"),
        gc=P.gc_track("{roi}", "{name}"),
        fixture=lambda w: P.fixture(
            "compartments", COMP_FIXTURE[w.cscope].format(roi=w.roi)),
    output:
        tsv=P.comps_refined("{cscope}", "{roi}", "{name}", "Comp_{name}.tsv"),
        bw=P.comps_refined("{cscope}", "{roi}", "{name}", "Comp_selected_{name}.bw"),
        selection=P.comps_refined("{cscope}", "{roi}", "{name}", "selection.tsv"),
    params:
        selection=COMP["eigenvector_selection"],
        eigs=EIGS,
        eig_bw=lambda w: {e: P.eigs(w.roi, w.name, which=e) for e in EIGS},
    wildcard_constraints:
        cscope="refined|refined_merged",
        name=_alt(COMP_NAMES),
    log:
        P.log("refine_compartments", "{cscope}_{roi}_{name}"),
    benchmark:
        P.benchmark("refine_compartments", "{cscope}_{roi}_{name}")
    resources:
        mem_mb=4000,
    script:
        "../scripts/py/refine_compartments.py"


rule saddle_strength_selected:
    """`saddle_strength_selected_{roi}.tsv` -- `file`, `value` (01_06 tail).

    ADDED beyond the plan's twelve rules, for the same reason as
    `pileup_scores`: 01_06 wrote one per-ROI table while looping over coolers,
    and per-cooler outputs plus a per-ROI table cannot be one Snakemake job.
    Behind EFig 3b and Fig 3c. Emits the eigenvector-choice table beside it, so
    a reader of the strengths can see which eigenvector each came from.
    """
    input:
        selection=lambda w: [
            P.comps_refined(w.cscope, w.roi, n, "selection.tsv")
            for n in COMP_SCOPES[w.cscope]
        ],
    output:
        selected=P.saddle_strength_selected("{cscope}", "{roi}"),
        choice=_os.path.join(
            _os.path.dirname(P.saddle_strength_selected("{cscope}", "{roi}")),
            "eigenvector_selection_{roi}.tsv"),
    wildcard_constraints:
        cscope="refined|refined_merged",
    log:
        P.log("saddle_strength_selected", "{cscope}_{roi}"),
    benchmark:
        P.benchmark("saddle_strength_selected", "{cscope}_{roi}")
    resources:
        mem_mb=2000,
    script:
        "../scripts/py/saddle_strength_selected.py"


rule eigenvector_orientation:
    """One table of every eigenvector sign flip the pipeline applied (plan 9.5).

    ADDED beyond the plan's twelve rules: the flip record is one global file
    while `eigs_cis` is per (roi, name). It matters far more than its size
    suggests. `compartments*.json` fixes only WHICH eigenvector, so if our
    orientation differs from the ground truth's, every saddle plot is mirrored
    and fifteen of the twenty-five panels are quietly wrong (risk R-3). **The
    visualisation stage must surface this table, not merely produce it.**

    Lives under `work/` because nothing outside 90_visualise.smk may write to
    `results/`; a `panel_*` rule copies it to `results/tables/`.
    """
    input:
        [_orientation_tsv(roi, name) for roi in ROI_TYPES for name in COMP_NAMES],
    output:
        P.work("features", "compartments", "eigenvector_orientation.tsv"),
    log:
        P.log("eigenvector_orientation"),
    benchmark:
        P.benchmark("eigenvector_orientation")
    resources:
        mem_mb=2000,
    shell:
        r"""
        mkdir -p "$(dirname {output})" "$(dirname {log})"
        {{
            head -n 1 {input[0]}
            for f in {input}; do tail -n +2 "$f"; done
        }} > {output} 2> {log}
        """


# ---------------------------------------------------------------------
# Allelic ratio  (01_07)
# ---------------------------------------------------------------------
rule allelic_gtf:
    """Per-clone GTF with `allelic_ratio "<value>";` injected (01_07).

    Thirteen clones, 5 WT + 8 degron, keyed by the verbatim column maps in
    `hic.allelic_ratio`. This is the "D-score" source for the green colour
    scale in Fig 3a/b, 5b/c/g/h, EFig 8c/d and 10c/d.

    The two joins are NOT symmetric, and the asymmetry is the original's:
    degron rows join on `gene_id`, WT rows on `gene_name` OR an exact
    (start, end) match. Reproduced as-is.

    `GRCm38.102_NC.gtf`, the other output of 01_07, is NOT produced here --
    `make_noncoding_gtf` in 00_reference.smk owns that path. Section 8 of
    the module notes records that the two derive it
    differently and that the discrepancy is open.
    """
    input:
        gtf=P.gtf(),
        wt=config["paths"]["inputs"]["rnaseq_allelic_ratio"]["wt"],
        degron=config["paths"]["inputs"]["rnaseq_allelic_ratio"]["degron"],
    output:
        gtf=P.allelic_gtf("{ar_clone}"),
        bgz=P.allelic_gtf("{ar_clone}", ext="gtf.bgz"),
        tbi=P.allelic_gtf("{ar_clone}", ext="gtf.bgz.tbi"),
    params:
        wt_map=AR["wt_clone_map"],
        degron_map=AR["degron_clone_map"],
        chrom="X",
    wildcard_constraints:
        ar_clone=_alt(AR_CLONES),
    log:
        P.log("allelic_gtf", "{ar_clone}"),
    benchmark:
        P.benchmark("allelic_gtf", "{ar_clone}")
    threads: 2
    resources:
        mem_mb=16000,
    script:
        "../scripts/py/allelic_gtf.py"


rule allelic_ratio_stats:
    """Allelic-ratio summary per locus and ROI (01_07, AR_stats half).

    Nine rows -- sum / mean / median x all / non_zero / >0.1 -- by clone, over
    the genes fully inside the ROI. Two files per job: WT and degron.

    Written with `sep='\\t'` despite the `.csv` extension, because that is what
    01_07 did and what 02_02 reads back.

    The degron table carries no coordinates of its own; 01_07 back-fills them
    by joining `gene_id` against the GTF, which is why the GTF is an input.
    """
    input:
        gtf=P.gtf(),
        wt=config["paths"]["inputs"]["rnaseq_allelic_ratio"]["wt"],
        degron=config["paths"]["inputs"]["rnaseq_allelic_ratio"]["degron"],
        roi_table=P.roi_table(),
    output:
        wt=P.allelic_ratio_stats("{locus}", "{roi}"),
        degron=P.allelic_ratio_stats("{locus}", "{roi}", degron=True),
    params:
        wt_map=AR["wt_clone_map"],
        degron_map=AR["degron_clone_map"],
        thresholds=AR["stats_thresholds"],
        chrom="X",
    log:
        P.log("allelic_ratio_stats", "{locus}_{roi}"),
    benchmark:
        P.benchmark("allelic_ratio_stats", "{locus}_{roi}")
    resources:
        mem_mb=4000,
    script:
        "../scripts/py/allelic_ratio_stats.py"


# ---------------------------------------------------------------------
# Convenience target
# ---------------------------------------------------------------------
rule hic_features:
    """Everything 31_hic_features.smk produces.

    The Snakefile's `hic` alias stops at refined loops plus METALoci; this is
    the full feature surface the visualisation stage draws from.
    """
    input:
        [P.loops_refined(n, ext=e) for n in LOOP_SETS
         for e in ("tsv", "bedpe", "bedpe.bgz", "bedpe.bgz.px2")],
        [P.loop_stats(c, roi, s)
         for roi in ROI_TYPES
         for c, sfx in (("dTAG_vs_NodTAG", [""]), ("Xa_vs_Xi", ["_Xa", "_Xi"]))
         for s in sfx],
        [P.pileup(c, roi, e, s) for (c, e, s) in PILEUP_JOBS for roi in ROI_TYPES],
        [P.pileup_score(c, roi) for c in PILEUP_SETS for roi in ROI_TYPES],
        [P.loop_anchor_stackup(c, roi, sig, pair)
         for (c, sig, pair) in STACKUP_JOBS for roi in ROI_TYPES],
        [P.saddle_values(roi, n) for roi in ROI_TYPES for n in COMP_NAMES],
        [P.saddle_strength_selected(cscope, roi)
         for cscope in COMP_SCOPES for roi in ROI_TYPES],
        P.work("features", "compartments", "eigenvector_orientation.tsv"),
        [P.allelic_gtf(c, ext=e) for c in AR_CLONES
         for e in ("gtf", "gtf.bgz", "gtf.bgz.tbi")],
        [P.allelic_ratio_stats(locus, roi, degron=d)
         for locus in LOCI for roi in ROI_TYPES for d in (False, True)],
