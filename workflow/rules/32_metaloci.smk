# =====================================================================
#  32_metaloci.smk — LAYER 2c.  THIS MODULE PRODUCES FIG 6.
#
#  mcool + 5 kb bedGraphs -> METALoci prep / layout / lm -> moran_info.txt
#  -> compartment-like strength tables.
#
#  Replaces: 01_09_metaloci_run_WT.ipynb, 01_09_metaloci_run_dtag.ipynb
#            (both minus their `metaloci figure` cell), the undocumented
#            `subset_mlo_from_file2_to_file1` patch, and
#            01_10_metaloci_compartments.ipynb.
#
#  NO FIGURES HERE. `metaloci figure` -- the Gaudi/KK/composite PDFs and
#  the `_gtp.pdf` panels Fig 6 is assembled from -- belongs to
#  90_visualise.smk. This module emits `.mlo`, `moran_info.txt`,
#  `.signals` and `compartmentalization_*.tsv` only.
#
#  ==================================================================
#  THE FOUR THINGS THAT WILL SILENTLY CORRUPT FIG 6
#  ==================================================================
#
#  1. THE `-pl 15` ARGPARSE TRAP.  `metaloci layout ... -pl 15` does NOT
#     mean "persistence length 15". argparse splits it into `-p`
#     (`--plot`, store_true) plus `-l 15` (`--pl`). The original's KK and
#     mixed-matrix PDFs exist only because `--plot` was set by accident.
#     A re-developer writing `--pl 15` loses them with no error.
#     ==> we always emit `-p -l <value>` explicitly. 15 for `full`,
#         10 for `escape` / `non_escape`.
#
#  2. THE KK TRANSPLANT (`metaloci_transplant_kk`). Between `layout` and
#     `lm` the FULL-region Kamada-Kawai embedding is sliced into the
#     `escape` and `non_escape` objects. `shift = abs(full.start -
#     sub.start) // resolution` -- **the abs() is mandatory**; without it
#     Mecp2 escape gets -66 and Jarid escape -280, and a negative index
#     slices from the end of the array. The full derivation, the inverted
#     function name, and the four replaced fields are in
#     `workflow/scripts/py/metaloci_transplant_kk.py`. Read it before
#     touching anything here.
#     CONSEQUENCE: the sub-region `-l 10` layout is DISCARDED. We compute
#     it anyway because the target object's start/end/poi/matrix and the
#     `_old.mlo` sidecar come from it.
#
#  3. SIGNAL NAMES ARE FILE BASENAMES. `metaloci prep` takes a signal's
#     name from the basename of its bedGraph. A mismatch does not error --
#     the signal is silently dropped and the Gaudi grid comes back one row
#     short. `metaloci_signals` stages symlinks under exactly
#     `naming.metaloci_signal_name(mark, cooler)` and `metaloci_lm` hard-
#     fails if any name in `.signals` is missing from the computed
#     `lmi_info`.
#
#  4. `moran_info.txt` IS NOT WRITTEN BY `metaloci lm`. In METALoci 1.3.2
#     it is appended by `metaloci figure`. Plan section 2.2 attributes it to
#     `lm`; the plan is wrong about the producer. Since `figure` is RD-5's
#     and `moran_info.txt` is data, `metaloci_lm` recomputes it from the
#     `.mlo` pickles using figure.py's arithmetic verbatim, including two
#     stray whitespace runs that are genuinely in the ground-truth bytes.
#     See `workflow/scripts/py/metaloci_lm.py`.
#
#  ==================================================================
#  DETERMINISM (plan section 9.3) — WHAT TESTERS MAY AND MAY NOT DIFF
#  ==================================================================
#  `metaloci lm` runs a 9999-permutation Local Moran test through
#  `esda.moran.Moran_Local(y, w, permutations=n, n_jobs=1)`. METALoci
#  1.3.2 exposes **no seed**: not on the CLI (`metaloci lm --help` offers
#  -p/-v/-a/-i/-m/-t/-f/-b/-q/-po and nothing else) and not in the call.
#  esda 2.7.0's `Moran_Local.__init__` DOES accept `seed=None` -- so a
#  seed is one keyword away upstream, and if METALoci ever passes it the
#  tolerance below can be tightened to exact. Today it cannot.
#
#    deterministic  sq1..sq4, q1..q4, r_value, and everything derived from
#                   the Kamada-Kawai layout (L-BFGS from a fixed circular
#                   start, reproducible for a pinned networkx)
#    stochastic     LMI_pvalue, and therefore the significance calls that
#                   move sq1..sq4 at the boundary
#    tolerance      |p_new - p_gt| <= 3*sqrt(p(1-p)/9999); the alpha=0.05
#                   significance call must agree for >= 99% of bins
#    NEVER          byte-diff a `.mlo`. They are pickles (the ground-truth inventory 8).
#
#  One further precision trap, documented at length in
#  `metaloci_compartmentalization.py`: METALoci stores Sig/Lag/ZSig/ZLag/
#  LMI_score/LMI_pvalue as `np.half`. `pearsonr` on float16 returns a
#  **float16 r** under numpy 1.26, which is why the ground-truth
#  `pearsonr` column carries 3-5 significant digits and `pearsonp` carries
#  17. Do not cast to float64 to "improve" it.
#
#  ==================================================================
#  WILDCARDS
#  ==================================================================
#  Plan section 2.2 lists `roi` as a wildcard of `metaloci_layout_sub` and
#  `metaloci_transplant_kk`. It cannot be: METALoci names its objects
#  `chrX_{start}_{end}_0.mlo`, with no roi token anywhere in the path, and
#  the transplant's outputs sit beside that object. So those three rules
#  carry `{run}/{dataset}/{start}/{end}` and the roi is *derived* from the
#  coordinates by `_roi_of()`, which also rejects coordinates that are not
#  a region of the dataset's locus. `{roi}` survives as a real wildcard
#  only where we name the file ourselves --
#  `metaloci_compartmentalization`.
#
#  Every rule's `log:` and `benchmark:` carry the SAME wildcard set as its
#  `output:`. Four rules in 10_signal.smk had to be fixed for exactly this;
#  it is not optional.
#
#  `{start}` and `{end}` are constrained to the coordinate sets that
#  actually exist in `config.loci`, and `{dataset}` to `[^/]+` so it cannot
#  swallow a path separator. `full` and the sub-regions share Mecp2's start
#  and the full end, so `ruleorder: metaloci_layout_full >
#  metaloci_layout_sub` resolves the one genuine overlap.
# =====================================================================

import os

ML = config["metaloci"]
ML_LAYOUT = ML["layout"]
ML_LM = ML["lm"]
ML_RESO = int(ML["resolution"])
ML_RUNS = list(ML["runs"])
ML_CHROM = "chrX"                       # METALoci is chrX-only here (all coolers are)
ML_THREADS = int(config["resources"]["max_threads_per_job"])
ML_MIN_POLYGON_AREA = 1e-2              # 01_10 stats_polygons: drop area <= 1e-2

# Region order is the original's `regions = ['full', 'escape', 'non_escape']`.
# `moran_info.txt` is append-ordered, so this is not cosmetic.
ML_REGION_ORDER = ("full", "escape", "non_escape")
ML_SUB_ROIS = ("escape", "non_escape")


def _locus_of(dataset):
    """`Mecp2_E6_WT_G1_Xi` / `Jarid_NodTAG-or-WT_Xa` -> the locus token."""
    for locus in LOCI:
        if dataset.startswith(f"{locus}_"):
            return locus
    raise WorkflowError(
        f"METALoci dataset {dataset!r} does not start with a known locus "
        f"({', '.join(LOCI)}). Cooler names are {{locus}}_... by construction."
    )


def _region(dataset, roi):
    start, end = config["loci"][_locus_of(dataset)]["regions"][roi]
    return int(start), int(end)


def _mlo(run, dataset, roi):
    start, end = _region(dataset, roi)
    return P.ml_mlo(run, dataset, start, end)


def _transplant_shift(dataset, roi):
    """abs(full.start - sub.start) // resolution.

    THE abs() IS MANDATORY -- see the module header and the design review R-1.
    Computed here from config so the assertion in the script is data-driven,
    and asserted again below so a config edit that breaks it fails at DAG
    build rather than three hours into a run.
    """
    full_start, _ = _region(dataset, "full")
    sub_start, _ = _region(dataset, roi)
    return abs(full_start - sub_start) // ML_RESO


def _mcool_of(run, dataset):
    """wt / degron datasets are individual coolers; consensus are merged_comps.

    Verbatim from the notebooks: the WT and DEGRON runs read
    `../data/mcoolers/{dataset}.mcool`, the CONSENSUS run reads
    `../data/mcoolers/merged_for_compartments/{dataset}.mcool`.
    """
    scope = "merged_comps" if run == "consensus" else "individual"
    return P.mcool(scope, dataset)


def _signal_source(run, mark, dataset):
    """Where one signal's 5 kb bedGraph comes from, per run.

    section_reports/hic_metaloci.md section 7. The consensus run takes every
    signal from the cross-clone consensus tree; wt/degron take AcMe3 from its
    own tree and the rest from the per-track bedGraphs. All three are the
    csaw/TMM-normalised tracks, never the raw deepTools merges.
    """
    adj = naming.dataset_adj(dataset)          # e.g. `E6_WT_Xi`, `NodTAG-or-WT_Xa`
    name = f"{mark}_{adj}"                     # == naming.metaloci_signal_name()
    if run == "consensus":
        return P.consensus5k(name)
    if mark == "AcMe3":
        return P.acme3_bg5k(name)
    return P.bedgraph5k(mark, name)


def _signals_of(run, dataset):
    """[(mark, metaloci_signal_name, source_path), ...] for one dataset."""
    out = []
    for mark in SS.metaloci_signals_for(dataset):
        out.append(
            (mark, naming.metaloci_signal_name(mark, dataset), _signal_source(run, mark, dataset))
        )
    return out


def _staged_signal_dir(run, dataset):
    """Symlink farm inside the dataset's METALoci working directory.

    Derived from `P.ml_dir`, not a new path root. METALoci itself only ever
    touches `signal/`, `tmp/`, `{chrom}/`, `*.signals`, `moran_info.txt`,
    `bad_regions.txt` and `version_log.txt` in a working directory, so a
    sibling `signal_input/` is inert.
    """
    return os.path.join(P.ml_dir(run, dataset), "signal_input")


def _prep_signal_tsv(run, dataset):
    """`metaloci prep`'s output: `{work_dir}signal/{chrom}/{chrom}_signal.tsv`."""
    return os.path.join(P.ml_dir(run, dataset), "signal", ML_CHROM, f"{ML_CHROM}_signal.tsv")


def _transplant_flag(run, dataset, roi):
    """Provenance stamp for one transplant; the DAG edge lm depends on.

    The patched `.mlo` cannot be a Snakemake output: the transplant rewrites
    the very file `metaloci_layout_sub` produced, in place, which is what the
    original did and what `metaloci lm` requires (it opens the object by its
    coordinate-derived name). So the declared outputs are the `_old.mlo`
    backup -- a real artefact of the original -- and this stamp, which records
    shift, size and the two file paths.
    """
    return _mlo(run, dataset, roi)[: -len(".mlo")] + ".transplanted"


def _old_mlo(run, dataset, roi):
    return _mlo(run, dataset, roi)[: -len(".mlo")] + "_old.mlo"


def _mlo_stem_pattern(suffix):
    """`.../metaloci/{run}/{dataset}/chrX/objects/chrX_{start}_{end}_0<suffix>`."""
    base = P.ml_mlo("{run}", "{dataset}", "{start}", "{end}")
    return base[: -len(".mlo")] + suffix


def _roi_of(dataset, start, end):
    """Resolve a sub-region roi from its coordinates, and validate it.

    The transplant carries `{start}`/`{end}` rather than `{roi}` because
    METALoci names its objects by coordinates only. This is where the two
    representations are reconciled -- and where a request for a region that is
    not an `escape`/`non_escape` of this dataset's locus fails at DAG build.
    """
    locus = _locus_of(dataset)
    for roi, (s, e) in config["loci"][locus]["regions"].items():
        if int(s) == int(start) and int(e) == int(end):
            if roi not in ML_SUB_ROIS:
                raise WorkflowError(
                    f"transplant requested for {locus}/{roi} "
                    f"({start}-{end}). Only {', '.join(ML_SUB_ROIS)} are "
                    "transplanted; the `full` layout is the SOURCE."
                )
            return roi
    raise WorkflowError(
        f"{start}-{end} is not a region of locus {locus!r} "
        f"(dataset {dataset!r}). Known: {config['loci'][locus]['regions']}"
    )


ML_DATASETS = {run: SS.metaloci_datasets(run) for run in ML_RUNS}

# --- the transplant shifts, asserted at DAG build ------------------------
# Expected, from config.loci: Mecp2 escape +66, Jarid escape +280, both
# non_escape 0. A negative or unexpected value here means the `abs()` was lost
# or the region coordinates moved; either way nothing downstream is meaningful.
for _locus in LOCI:
    for _roi in ML_SUB_ROIS:
        _fs, _ = config["loci"][_locus]["regions"]["full"]
        _ss, _ = config["loci"][_locus]["regions"][_roi]
        _shift = abs(int(_fs) - int(_ss)) // ML_RESO
        if _shift < 0:
            raise WorkflowError(
                f"METALoci transplant shift for {_locus}/{_roi} is {_shift}. "
                "A negative shift slices from the END of the layout array and "
                "silently corrupts every escape panel. See the design review R-1."
            )

# Constrain the coordinate wildcards to the coordinates that actually exist, so
# a typo cannot invent a region. `full` and `sub` overlap on Mecp2's start and
# on the full end, which is what the ruleorder below is for.
ML_FULL_STARTS = sorted({str(config["loci"][l]["regions"]["full"][0]) for l in LOCI})
ML_FULL_ENDS = sorted({str(config["loci"][l]["regions"]["full"][1]) for l in LOCI})
ML_SUB_STARTS = sorted(
    {str(config["loci"][l]["regions"][r][0]) for l in LOCI for r in ML_SUB_ROIS}
)
ML_SUB_ENDS = sorted(
    {str(config["loci"][l]["regions"][r][1]) for l in LOCI for r in ML_SUB_ROIS}
)


# ---------------------------------------------------------------------
# Signals
# ---------------------------------------------------------------------
rule metaloci_signals:
    """Stage this dataset's signal bedGraphs under METALoci-legal basenames
    and write `{dataset}.signals`.

    See trap 3 in the module header: the basename IS the signal name, and a
    mismatch is silent. Per-dataset subsetting (`Rad21` only for B1621, `CTCF`
    never for B1621, C5C10 CUT&RUN dropped except RNA-Seq) comes from
    `SS.metaloci_signals_for()`, which encodes `01_10` cell 4.
    """
    input:
        beds=lambda w: [s[2] for s in _signals_of(w.run, w.dataset)],
    output:
        signals=P.ml_signals("{run}", "{dataset}"),
        staged=directory(_staged_signal_dir("{run}", "{dataset}")),
    params:
        signal_names=lambda w: [s[1] for s in _signals_of(w.run, w.dataset)],
        dataset=lambda w: w.dataset,
    wildcard_constraints:
        dataset=r"[^/]+",
    log:
        P.log("metaloci_signals", "{run}_{dataset}"),
    resources:
        mem_mb=2000,
    conda:
        "../envs/metaloci.yaml"
    script:
        "../scripts/py/metaloci_signals.py"


# ---------------------------------------------------------------------
# prep
# ---------------------------------------------------------------------
rule metaloci_prep:
    """metaloci prep -w ... -c ... -r 5000 -s chrom.sizes -d sig1.bed ... sigN.bed

    `-t/--summarize_type` is left unset, so bins are summarised by the default
    **median** -- verified against `metaloci prep --help` in the original
    `metaloci2` env. The `-d` list is given in config signal order, which fixes
    the column order of `chrX_signal.tsv`.
    """
    input:
        signals=P.ml_signals("{run}", "{dataset}"),
        staged=_staged_signal_dir("{run}", "{dataset}"),
        beds=lambda w: [s[2] for s in _signals_of(w.run, w.dataset)],
        mcool=lambda w: _mcool_of(w.run, w.dataset),
        chrom_sizes=P.chrom_sizes(),
    output:
        signal_tsv=_prep_signal_tsv("{run}", "{dataset}"),
    params:
        work_dir=lambda w: P.ml_dir(w.run, w.dataset),
        resolution=ML_RESO,
        staged_beds=lambda w: " ".join(
            os.path.join(_staged_signal_dir(w.run, w.dataset), f"{s[1]}.bed")
            for s in _signals_of(w.run, w.dataset)
        ),
    wildcard_constraints:
        dataset=r"[^/]+",
    log:
        P.log("metaloci_prep", "{run}_{dataset}"),
    benchmark:
        P.benchmark("metaloci_prep", "{run}_{dataset}")
    threads: 4
    resources:
        mem_mb=16000,
    conda:
        "../envs/metaloci.yaml"
    shell:
        r"""
        mkdir -p "{params.work_dir}" "$(dirname {log})"
        metaloci prep \
            -w {params.work_dir} \
            -c {input.mcool} \
            -r {params.resolution} \
            -s {input.chrom_sizes} \
            -d {params.staged_beds} > {log} 2>&1
        test -s {output.signal_tsv} \
            || {{ echo "metaloci prep produced no signal table -- check that every" \
                       "-d basename matches a name in {input.signals}" >&2; exit 1; }}
        """


# ---------------------------------------------------------------------
# layout
#
# Two rules, one output pattern. `full` wins the overlap by ruleorder; the
# wildcard constraints keep every other combination unambiguous.
# ---------------------------------------------------------------------
ruleorder: metaloci_layout_full > metaloci_layout_sub


rule metaloci_layout_full:
    """Kamada-Kawai layout of the FULL region, persistence length 15.

        metaloci layout -w .. -c .. -r 5000 -g chrX:{s}-{e}_0 \
                        -a -o 1.5 -p -l 15 -m -t N -f

    `-a` makes `-o` an absolute minimum interaction strength rather than a
    fraction of top interactions. `-p -l 15` is written out because `-pl 15`
    is the argparse trap (module header, trap 1).

    This is the layout the sub-regions inherit, so it must exist before either
    transplant. That ordering is structural: `metaloci_transplant_kk` takes
    this file as an input.
    """
    input:
        signal_tsv=_prep_signal_tsv("{run}", "{dataset}"),
        mcool=lambda w: _mcool_of(w.run, w.dataset),
    output:
        mlo=P.ml_mlo("{run}", "{dataset}", "{start}", "{end}"),
    params:
        work_dir=lambda w: P.ml_dir(w.run, w.dataset),
        resolution=ML_RESO,
        chrom=ML_CHROM,
        cutoff=ML_LAYOUT["cutoff"],
        pl=ML_LAYOUT["persistence_length"]["full"],
        absolute="-a" if ML_LAYOUT["absolute"] else "",
        plot="-p" if ML_LAYOUT["plot"] else "",
        mp="-m" if ML_LAYOUT["multiprocess"] else "",
        force="-f" if ML_LAYOUT["force"] else "",
    wildcard_constraints:
        dataset=r"[^/]+",
        start="|".join(ML_FULL_STARTS),
        end="|".join(ML_FULL_ENDS),
    log:
        P.log("metaloci_layout_full", "{run}_{dataset}_{start}_{end}"),
    benchmark:
        P.benchmark("metaloci_layout_full", "{run}_{dataset}_{start}_{end}")
    threads: ML_THREADS
    resources:
        mem_mb=32000,
    conda:
        "../envs/metaloci.yaml"
    shell:
        r"""
        mkdir -p "{params.work_dir}" "$(dirname {log})"
        metaloci layout \
            -w {params.work_dir} \
            -c {input.mcool} \
            -r {params.resolution} \
            -g {params.chrom}:{wildcards.start}-{wildcards.end}_0 \
            {params.absolute} -o {params.cutoff} \
            {params.plot} -l {params.pl} \
            {params.mp} -t {threads} {params.force} > {log} 2>&1
        test -s {output.mlo} || {{ echo "no .mlo produced -- see {log}" >&2; exit 1; }}
        """


rule metaloci_layout_sub:
    """Kamada-Kawai layout of an `escape` / `non_escape` region, pl 10.

    THIS LAYOUT IS THEN THROWN AWAY. `metaloci_transplant_kk` overwrites its
    four `kk_*` fields with a slice of the full-region layout, so `-l 10`
    never reaches a figure or a statistic (module header, trap 2). It is still
    computed, because the object's own `start`/`end`/`poi`/`matrix`/
    `subset_matrix` come from here and because the `_old.mlo` sidecar the
    original left behind is this file.
    """
    input:
        signal_tsv=_prep_signal_tsv("{run}", "{dataset}"),
        mcool=lambda w: _mcool_of(w.run, w.dataset),
    output:
        mlo=P.ml_mlo("{run}", "{dataset}", "{start}", "{end}"),
    params:
        work_dir=lambda w: P.ml_dir(w.run, w.dataset),
        resolution=ML_RESO,
        chrom=ML_CHROM,
        cutoff=ML_LAYOUT["cutoff"],
        pl=ML_LAYOUT["persistence_length"]["escape"],
        absolute="-a" if ML_LAYOUT["absolute"] else "",
        plot="-p" if ML_LAYOUT["plot"] else "",
        mp="-m" if ML_LAYOUT["multiprocess"] else "",
        force="-f" if ML_LAYOUT["force"] else "",
    wildcard_constraints:
        dataset=r"[^/]+",
        start="|".join(ML_SUB_STARTS),
        end="|".join(ML_SUB_ENDS),
    log:
        P.log("metaloci_layout_sub", "{run}_{dataset}_{start}_{end}"),
    benchmark:
        P.benchmark("metaloci_layout_sub", "{run}_{dataset}_{start}_{end}")
    threads: ML_THREADS
    resources:
        mem_mb=32000,
    conda:
        "../envs/metaloci.yaml"
    shell:
        r"""
        mkdir -p "{params.work_dir}" "$(dirname {log})"
        metaloci layout \
            -w {params.work_dir} \
            -c {input.mcool} \
            -r {params.resolution} \
            -g {params.chrom}:{wildcards.start}-{wildcards.end}_0 \
            {params.absolute} -o {params.cutoff} \
            {params.plot} -l {params.pl} \
            {params.mp} -t {threads} {params.force} > {log} 2>&1
        test -s {output.mlo} || {{ echo "no .mlo produced -- see {log}" >&2; exit 1; }}
        """


# ---------------------------------------------------------------------
# The transplant
# ---------------------------------------------------------------------
rule metaloci_transplant_kk:
    """Slice the FULL-region Kamada-Kawai layout into a sub-region `.mlo`.

    `subset_mlo_from_file2_to_file1` from `01_09_*` cells 13-14, AFTER layout
    and BEFORE lm. Runs twice per dataset: `full -> non_escape`, then
    `full -> escape`.

        shift = abs(full.start - sub.start) // resolution   <- abs() MANDATORY
        size  = len(sub.kk_nodes)

    Replaces exactly four fields, sliced `[shift : shift+size]`:
    `kk_nodes`, `kk_restraints_matrix`, `kk_coords`, `kk_distances`. The
    target keeps its own start / end / resolution / signal. The original
    object is renamed `*_old.mlo`.

    Expected shifts, asserted by the script against `expected_shift`:
    Mecp2 non_escape 0, Mecp2 escape +66, Jarid non_escape 0,
    Jarid escape +280.

    The patched `.mlo` is written in place and is therefore NOT a declared
    output -- see `_transplant_flag()`. The declared outputs are the
    `_old.mlo` backup and the provenance stamp.
    """
    input:
        full=lambda w: _mlo(w.run, w.dataset, "full"),
        sub=P.ml_mlo("{run}", "{dataset}", "{start}", "{end}"),
    output:
        old=_mlo_stem_pattern("_old.mlo"),
        flag=_mlo_stem_pattern(".transplanted"),
    params:
        expected_shift=lambda w: _transplant_shift(
            w.dataset, _roi_of(w.dataset, w.start, w.end)
        ),
        resolution=ML_RESO,
        roi=lambda w: _roi_of(w.dataset, w.start, w.end),
        dataset=lambda w: w.dataset,
    wildcard_constraints:
        dataset=r"[^/]+",
        start="|".join(ML_SUB_STARTS),
        end="|".join(ML_SUB_ENDS),
    log:
        P.log("metaloci_transplant_kk", "{run}_{dataset}_{start}_{end}"),
    resources:
        mem_mb=8000,
    conda:
        "../envs/metaloci.yaml"
    script:
        "../scripts/py/metaloci_transplant_kk.py"


# ---------------------------------------------------------------------
# lm
# ---------------------------------------------------------------------
rule metaloci_lm:
    """Local Moran's I over all three regions of one dataset, then
    `moran_info.txt`.

        metaloci lm -w .. -s ..signals -g chrX:{s}-{e}_0 -m -t N -f

    Defaults, unchanged: 9999 permutations, p 0.05, quadrants [1,3].
    `-f` clears only `lmi_info` / `lmi_geometry`, never the transplanted
    `kk_*` fields.

    One job per dataset rather than per region, because `moran_info.txt` is a
    single per-dataset file and because METALoci does a read-modify-write of
    `version_log.txt` on every invocation.

    `moran_info.txt` is derived here, not by `metaloci figure` -- module
    header, trap 4.
    """
    input:
        signals=P.ml_signals("{run}", "{dataset}"),
        signal_tsv=_prep_signal_tsv("{run}", "{dataset}"),
        mlos=lambda w: [_mlo(w.run, w.dataset, r) for r in ML_REGION_ORDER],
        transplants=lambda w: [
            _transplant_flag(w.run, w.dataset, r) for r in ML_SUB_ROIS
        ],
    output:
        moran=P.ml_moran("{run}", "{dataset}"),
    params:
        work_dir=lambda w: P.ml_dir(w.run, w.dataset),
        dataset=lambda w: w.dataset,
        chrom=ML_CHROM,
        regions=lambda w: [
            (r, *_region(w.dataset, r), _mlo(w.run, w.dataset, r)) for r in ML_REGION_ORDER
        ],
        permutations=ML_LM["permutations"],
        pvalue=ML_LM["pvalue"],
        quadrants=ML_LM["quadrants"],
        multiprocess=ML_LAYOUT["multiprocess"],
        force=ML_LAYOUT["force"],
    wildcard_constraints:
        dataset=r"[^/]+",
    log:
        P.log("metaloci_lm", "{run}_{dataset}"),
    benchmark:
        P.benchmark("metaloci_lm", "{run}_{dataset}")
    threads: ML_THREADS
    resources:
        mem_mb=24000,
    conda:
        "../envs/metaloci.yaml"
    script:
        "../scripts/py/metaloci_lm.py"


# ---------------------------------------------------------------------
# Compartment-like strength
# ---------------------------------------------------------------------
rule metaloci_compartmentalization:
    """`01_10_metaloci_compartments.ipynb` -- one table per (run, region).

        compartmentalization (%) = (sq1 + sq3) / total_bins * 100,  p < 0.05

    plus `pearsonr(ZSig, ZLag)` and per-quadrant polygon statistics. 3 regions
    x 3 runs = 9 tables, where the original hand-edited two variables 9 times.

    Depends on `moran_info.txt`, not on the `.mlo` alone: the `.mlo` files
    exist after `layout`, but only `metaloci_lm` puts `lmi_info` in them.
    """
    input:
        moran=lambda w: [P.ml_moran(w.run, d) for d in ML_DATASETS[w.run]],
        mlos=lambda w: [_mlo(w.run, d, w.roi) for d in ML_DATASETS[w.run]],
    output:
        tsv=P.ml_compartmentalization("{run}", "{roi}"),
    params:
        run=lambda w: w.run,
        roi=lambda w: w.roi,
        pvalue=ML_LM["pvalue"],
        min_polygon_area=ML_MIN_POLYGON_AREA,
        jobs=lambda w: [
            (
                d,
                _mlo(w.run, d, w.roi),
                [(mark, name) for mark, name, _src in _signals_of(w.run, d)],
            )
            for d in ML_DATASETS[w.run]
        ],
    log:
        P.log("metaloci_compartmentalization", "{run}_{roi}"),
    benchmark:
        P.benchmark("metaloci_compartmentalization", "{run}_{roi}")
    threads: 4
    resources:
        mem_mb=16000,
    conda:
        "../envs/metaloci.yaml"
    script:
        "../scripts/py/metaloci_compartmentalization.py"


rule metaloci_merge_tables:
    """`compartmentalization_full_final.tsv`.

    NOTE the deviation, spelled out in the script: the ground-truth final table
    holds wt + degron only (188 rows, 44 datasets) even though the consensus
    per-run table exists. Plan section 2.2 says "the three runs", so all runs in
    `config.metaloci.runs` are merged here. Compare on the wt + degron subset.
    """
    input:
        tables=[P.ml_compartmentalization(r, "full") for r in ML_RUNS],
    output:
        tsv=P.ml_compartmentalization_final(),
    params:
        runs=ML_RUNS,
    log:
        P.log("metaloci_merge_tables"),
    resources:
        mem_mb=4000,
    conda:
        "../envs/metaloci.yaml"
    script:
        "../scripts/py/metaloci_merge_tables.py"


# ---------------------------------------------------------------------
# Convenience target
# ---------------------------------------------------------------------
rule metaloci:
    """Everything this module produces. Fig 6's data layer."""
    input:
        P.ml_compartmentalization_final(),
        [
            P.ml_compartmentalization(run, roi)
            for run in ML_RUNS
            for roi in ML_REGION_ORDER
        ],
        [P.ml_moran(run, d) for run in ML_RUNS for d in ML_DATASETS[run]],
