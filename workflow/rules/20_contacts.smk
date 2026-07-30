# =====================================================================
#  20_contacts.smk — LAYER 1b
#
#  Capture Hi-C coolers -> corrected -> merged -> balanced/zoomified mcool.
#  This is the entire Hi-C input side; 31_hic_features.smk and
#  32_metaloci.smk both start from the `.mcool` files produced here.
#
#  Four things about these coolers are not obvious and are load-bearing:
#
#  * They are **allele-specific Capture Hi-C, chrX only** (`nchroms = 1`),
#    with real coverage in exactly two captured windows. That is why ICE
#    is run with a blacklist that masks everything OUTSIDE the capture
#    region: balancing across the untargeted 99% of chrX would divide the
#    signal by noise.
#
#  * A handful of pixels are ligation artefacts that were repaired by hand
#    at 5 kb. Which pixels, and how each is replaced, is human judgement
#    that exists nowhere but inside `00_01`; it ships as
#    `resources/fixtures/pixel_fixes.tsv` (D-10).
#
#  * The merged coolers are a hand-enumerated pooling, not a groupby. Two
#    scopes exist and they are NOT nested: `merged_loops` pools WT+NodTAG
#    into `{locus}_{allele}` (9 Mecp2, 6 Jarid -> the "N = 15" of EFig 3),
#    while `merged_comps` also splits the degrons by target. A dTAG cooler
#    belongs to TWO compartment groups at once. Ships as
#    `resources/fixtures/cooler_membership.tsv` (D-10).
#
#  * F3 is the only clone with two Capture Hi-C replicates. They are summed
#    with `cooler merge` into a parent cooler that then behaves exactly like
#    a single-replicate one -- 48 raw + 4 merges = the 52 individual mcools
#    of the ground truth, and the 52 keys of the compartment fixtures.
#
#  NO FIGURES. `00_08_viz_coolers.py`'s per-cooler QC PDFs are drawn by
#  `panel_cooler_qc` in 90_visualise.smk, from the mcools produced here.
#
#  Replaces: 00_00_copy_and_rename_input_files.py,
#            00_01_correct_pixels_in_coolers.py,
#            00_03_merge_F3_coolers.sh + `00_03_merge_F3_coolers copy.sh`,
#            00_04_balance_and_zoomify.sh,
#            00_09_mcool_to_hic.py (off by default)
# =====================================================================

import re as _re

HIC = config["hic"]
BAL = HIC["balance"]
FIX = config["fixtures"]

# Cooler name universe, split by how each one comes into being.
RAW_COOLERS = [n for n in SS.cooler_names if SS.cooler_path(n)]
REPMERGE_COOLERS = [n for n in SS.cooler_names if not SS.cooler_path(n)]
MERGED_LOOPS = SS.merged_cooler_names("merged_loops")
MERGED_COMPS = SS.merged_cooler_names("merged_comps")
ALL_COOLER_NAMES = sorted(set(SS.cooler_names) | set(MERGED_LOOPS) | set(MERGED_COMPS))

MCOOL_SCOPES = {
    "individual": SS.cooler_names,
    "merged_loops": MERGED_LOOPS,
    "merged_comps": MERGED_COMPS,
}


def _alt(names):
    """An alternation regex over `names`, or one that can never match."""
    return "|".join(_re.escape(n) for n in names) if names else r"(?!x)x"


def _locus_of(name):
    """Leading token of any cooler name -- merged or not -- is the locus."""
    locus = name.split("_", 1)[0]
    if locus not in LOCI:
        raise WorkflowError(
            f"cooler {name!r} does not start with a locus from config.loci "
            f"({', '.join(LOCI)}); the ICE blacklist is chosen by that prefix."
        )
    return locus


def _mcool_source(wildcards):
    """The `.cool` that `balance_zoomify` reads, per scope."""
    if wildcards.scope == "individual":
        return P.cool_corrected(wildcards.name)
    return P.cool_merged(wildcards.scope, wildcards.name)


# ---------------------------------------------------------------------
# Raw input staging
# ---------------------------------------------------------------------
rule link_cooler:
    """Symlink the user's cooler in under its frozen name.

    The original `00_00` renamed physical copies from a hand-maintained TSV.
    Here `config/coolers.tsv` carries the mapping and nothing is copied, so
    the source directory stays read-only by construction.
    """
    input:
        lambda w: SS.cooler_path(w.name),
    output:
        P.cool_raw("{name}"),
    wildcard_constraints:
        name=_alt(RAW_COOLERS),
    log:
        P.log("link_cooler", "{name}"),
    shell:
        r"""
        mkdir -p "$(dirname {output})" "$(dirname {log})"
        ln -sfn "$(readlink -f {input})" {output} 2> {log}
        """


# ---------------------------------------------------------------------
# Artefact-pixel repair
# ---------------------------------------------------------------------
rule correct_cooler_pixels:
    """Replace the hand-identified artefact pixels (00_01).

    Reads `pixel_fixes.tsv`; each row names one 5 kb x 5 kb pixel and the rule
    for replacing it (mean of its existing neighbours, or a fixed 5-neighbour
    directional mean). Under `fixtures.mode: auto` the step degrades to a copy
    and says so in the log -- the pixel list is specific to these two capture
    regions and means nothing on other data.

    NEVER run the original in place: `00_01` writes into the ground-truth
    `coolers_corrected/` directory.
    """
    input:
        cool=P.cool_raw("{name}"),
        fixes=P.fixture("pixel_fixes.tsv"),
    output:
        P.cool_corrected("{name}"),
    params:
        mode=FIX["mode"],
        resolution=HIC["resolution"],
    wildcard_constraints:
        name=_alt(RAW_COOLERS),
    log:
        P.log("correct_cooler_pixels", "{name}"),
    benchmark:
        P.benchmark("correct_cooler_pixels", "{name}")
    threads: 2
    resources:
        mem_mb=8000,
    script:
        "../scripts/py/correct_cooler_pixels.py"


# ---------------------------------------------------------------------
# Merging
# ---------------------------------------------------------------------
rule merge_f3_replicates:
    """Sum the Capture Hi-C replicates of one clone-allele (`cooler merge`).

    Named for F3 because F3 is the only clone that has replicates, but the
    membership is derived from `config/coolers.tsv`, not hard-coded: any cooler
    row with an empty `path` and `_rep{{n}}` siblings is a replicate merge.
    Writes into the SAME directory as `correct_cooler_pixels`, matching the
    original -- the merge happens after correction, so the corrected replicates
    are its inputs and the product needs no further pixel repair.
    """
    input:
        lambda w: [P.cool_corrected(c) for c in SS.replicate_members(w.name)],
    output:
        P.cool_corrected("{name}"),
    wildcard_constraints:
        name=_alt(REPMERGE_COOLERS),
    log:
        P.log("merge_f3_replicates", "{name}"),
    benchmark:
        P.benchmark("merge_f3_replicates", "{name}")
    threads: 2
    resources:
        mem_mb=8000,
    shell:
        r"""
        mkdir -p "$(dirname {output})" "$(dirname {log})"
        cooler merge {output} {input} > {log} 2>&1
        """


rule merge_coolers_loops:
    """Pool WT+NodTAG (or all dTAG) clone-alleles for loop calling.

    `cooler merge` sums pixels and does NOT re-balance; the pooled cooler is
    balanced afterwards by `balance_zoomify`, exactly as in the original.
    """
    input:
        lambda w: [P.cool_corrected(c) for c in SS.merged_cooler_members("merged_loops", w.name)],
    output:
        P.cool_merged("merged_loops", "{name}"),
    wildcard_constraints:
        name=_alt(MERGED_LOOPS),
    log:
        P.log("merge_coolers_loops", "{name}"),
    benchmark:
        P.benchmark("merge_coolers_loops", "{name}")
    threads: 2
    resources:
        mem_mb=16000,
    shell:
        r"""
        mkdir -p "$(dirname {output})" "$(dirname {log})"
        cooler merge {output} {input} > {log} 2>&1
        """


rule merge_coolers_comps:
    """Pool clone-alleles by condition group for compartment calling.

    Sixteen groups, and they overlap: `Mecp2_F3_CTCF-dTAG_G2_Xa` is a member of
    both `Mecp2_CTCF-dTAG_Xa` and the pooled `Mecp2_dTAG_Xa`. That is why the
    `merge_comps` column of `config/coolers.tsv` is comma-separated.
    """
    input:
        lambda w: [P.cool_corrected(c) for c in SS.merged_cooler_members("merged_comps", w.name)],
    output:
        P.cool_merged("merged_comps", "{name}"),
    wildcard_constraints:
        name=_alt(MERGED_COMPS),
    log:
        P.log("merge_coolers_comps", "{name}"),
    benchmark:
        P.benchmark("merge_coolers_comps", "{name}")
    threads: 2
    resources:
        mem_mb=16000,
    shell:
        r"""
        mkdir -p "$(dirname {output})" "$(dirname {log})"
        cooler merge {output} {input} > {log} 2>&1
        """


# ---------------------------------------------------------------------
# Balancing + multi-resolution
# ---------------------------------------------------------------------
rule balance_zoomify:
    """ICE-balance and zoomify to the 11 published resolutions (00_04).

    Every parameter here was read back out of the `bins/weight` HDF5
    attributes of the ground-truth mcool, so these are what actually ran and
    not what a script claims: ICE, cis_only, mad_max 3, min_nnz 10, min_count
    0, ignore_diags 2, tol 1e-5, divisive_weights False. (The last four are
    `cooler balance` defaults and are therefore not passed explicitly, exactly
    as in the original command line.)

    `--cis-only` and the locus blacklist together are the whole trick: the
    blacklist masks everything outside the capture window, so ICE converges on
    the ~1 Mb that has data. The blacklist is chosen by the cooler's locus
    prefix, exactly as the original's `if [[ $filename =~ ^Jarid ]]`.

    The original had to be edited and re-run three times, once per input
    directory. Here `{scope}` does that.
    """
    input:
        cool=_mcool_source,
        blacklist=lambda w: P.ice_blacklist(_locus_of(w.name)),
    output:
        P.mcool("{scope}", "{name}"),
    params:
        resolutions=",".join(str(r) for r in HIC["zoomify_resolutions"]),
        mad_max=BAL["mad_max"],
        force="--force" if BAL["force"] else "",
        cis_only="--cis-only" if BAL["cis_only"] else "",
    wildcard_constraints:
        name=_alt(ALL_COOLER_NAMES),
    log:
        P.log("balance_zoomify", "{scope}_{name}"),
    benchmark:
        P.benchmark("balance_zoomify", "{scope}_{name}")
    threads: 24
    resources:
        mem_mb=24000,
    shell:
        r"""
        mkdir -p "$(dirname {output})" "$(dirname {log})"
        cooler zoomify \
            -n {threads} \
            -r {params.resolutions} \
            -o {output} \
            --balance \
            --balance-args "{params.force} -p {threads} --blacklist {input.blacklist} --mad-max {params.mad_max} {params.cis_only}" \
            {input.cool} > {log} 2>&1
        """


# ---------------------------------------------------------------------
# Optional .hic export
# ---------------------------------------------------------------------
rule mcool_to_hic:
    """`.hic` export for Juicebox (00_09). OFF by default.

    Nothing in the paper depends on it -- it exists so a reader can open the
    contact maps in Juicebox. HiCLift pins its own numpy/cooler versions, so it
    lives in its own conda env rather than being forced into the solved core
    env; that env is only ever built if this rule is actually requested.

    Turn it on with  hic.mcool_to_hic.enabled: true  and run  ./run.sh contacts_hic .
    """
    input:
        P.mcool("{scope}", "{name}"),
    output:
        P.hic("{scope}", "{name}"),
    params:
        chrom_sizes=P.chrom_sizes(),
        resolution=HIC["resolution"],
        memory=HIC["mcool_to_hic"]["memory"],
        assembly=config["genome"]["build"],
        prefix=lambda w: P.hic(w.scope, w.name)[: -len(".hic")],
    wildcard_constraints:
        name=_alt(ALL_COOLER_NAMES),
    log:
        P.log("mcool_to_hic", "{scope}_{name}"),
    benchmark:
        P.benchmark("mcool_to_hic", "{scope}_{name}")
    conda:
        "../envs/hiclift.yaml"
    threads: 8
    resources:
        mem_mb=40000,
    shell:
        r"""
        mkdir -p "$(dirname {output})" "$(dirname {log})"
        HiCLift \
            --input {input}::resolutions/{params.resolution} \
            --input-format cooler \
            --out-pre {params.prefix} \
            --output-format hic \
            --out-chromsizes {params.chrom_sizes} \
            --in-assembly {params.assembly} \
            --out-assembly {params.assembly} \
            --memory {params.memory} > {log} 2>&1
        """


# ---------------------------------------------------------------------
# Convenience target (the `contacts` alias itself lives in the Snakefile)
# ---------------------------------------------------------------------
rule contacts_hic:
    """Every `.hic` export. Requires  hic.mcool_to_hic.enabled: true ."""
    input:
        [P.hic(scope, name) for scope, names in MCOOL_SCOPES.items() for name in names]
        if HIC["mcool_to_hic"]["enabled"]
        else [],
