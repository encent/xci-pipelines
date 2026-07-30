# =====================================================================
#  10_signal.smk — LAYER 1a
#
#  BAMs -> normalised bigWigs / bedGraphs. This is the dominant cost in
#  the whole pipeline: 17 csaw jobs gate ~250 bamCoverage jobs over
#  ~763 GB of BAM.
#
#  The normalisation is two-stage, exactly as the original:
#    1. csaw counts 10 kb bins on AUTOSOMES ONLY and edgeR-TMM derives a
#       per-sample NormFactor. chrX is excluded by design -- allelic and
#       clone-genotype differences corrupt TMM (D-04).
#    2. the Gall scale factor is applied to that sample's Xa and Xi
#       tracks too, because the three allele splits share one library.
#
#  Step 2 is where the original had a real bug: it broadcast Gall factors
#  to Xa/Xi *positionally*, over an alphabetically sorted file list, and
#  the WT H3K27ac `selected` run leaves an odd count -- so 5 of 27
#  bigWigs carry another sample's factor. We join on
#  (mark, clone, condition, replicate) instead. See D-03: the original
#  behaviour is still reachable, bit-for-bit, via
#  `legacy.h3k27ac_positional_scalefactors: true`.
#
#  Replaces: 12_normalization_*.r, 13_merge_normalized_bigwigs*.py,
#            15_bedgraph_all_FINAL.sh, 16_acme3_dtag_NEW.py,
#            17_consensus_bw_and_bedgraph.py
# =====================================================================

NORM = config["normalization"]
BC = NORM["bam_coverage"]
MERGE = NORM["merge"]
ACME3 = NORM["acme3"]
CONSENSUS = NORM["consensus"]


def _mark_of(track_or_sample):
    return track_or_sample.split("_", 1)[0]


def _bam_coverage_read_mode(mark):
    """--extendReads / --centerReads / neither, per mark.

    RNA-Seq is never extended (spliced reads would be bridged across introns);
    point-source factors are centered; broad marks are extended.
    """
    if mark in BC.get("no_extend_marks", []):
        return ""
    if mark in BC.get("center_reads_marks", []):
        return "--centerReads --extendReads"
    return "--extendReads"


# ---------------------------------------------------------------------
# Raw input staging
# ---------------------------------------------------------------------
rule link_bam:
    """Symlink the user's BAM into the pipeline namespace under the frozen
    sample id. Never copies -- the BAMs are ~763 GB and read-only."""
    input:
        lambda w: SS.bam_of(w.sample),
    output:
        P.bam("{sample}"),
    log:
        P.log("link_bam", "{sample}"),
    shell:
        r"""
        mkdir -p "$(dirname {output})" "$(dirname {log})"
        ln -sfn "$(readlink -f {input})" {output} 2> {log}
        """


rule index_bam:
    """Make an index available next to the symlink, REUSING the existing one.

    Prefer an existing sibling index; only build one if genuinely absent.

    This used to run `samtools index` unconditionally, which on a production
    run meant a full read of ~763 GB of BAM to regenerate indexes that were
    already there. All 273 in-use BAMs in the reference dataset ship with a
    `.bai`, so the entire cost was waste, and it landed in front of the first
    `bam_coverage` job -- the single largest cost of the valley chain. It is a
    performance defect only, so no correctness test would ever have caught it.

    It also made the module internally inconsistent: `link_bam` deliberately
    symlinks rather than copies, for exactly the same reason.

    Both naming conventions occur in the wild and both are accepted:

        /path/sample.bam.bai      (samtools default)
        /path/sample.bai          (Picard / older GATK)

    An existing index is COPIED, not symlinked, and that is deliberate.

    A symlink here is a write-through hazard: anything that later opens the
    link for writing -- `samtools index`, or Snakemake touching its own output
    to update an mtime -- writes through into the source tree. When that tree
    is mode r--r--r--, as the reference data is, the write fails and Snakemake
    deletes the output it thinks it just made. Phase 4 hit exactly this.

    **This is not a niche case; it is our target audience.** Someone who
    downloads the published data and mounts it read-only, or receives it on a
    read-only volume, hits it on their first run. A `.bai` is ~2.5 MB against a
    ~2 GB BAM, so copying all 273 costs well under a gigabyte -- a rounding
    error against the ~763 GB of reads it avoids.

    Nothing is ever written next to the source BAM. A genuinely missing index
    is built INTO our namespace by pointing `samtools index` at our own symlink
    and giving it an explicit output path.
    """
    input:
        bam=P.bam("{sample}"),
    output:
        bai=P.bam("{sample}") + ".bai",
    params:
        # The real file, not the symlink: its siblings are where an index lives.
        src=lambda w: SS.bam_of(w.sample),
    log:
        P.log("index_bam", "{sample}"),
    threads: 4
    resources:
        mem_mb=4000,
    shell:
        r"""
        mkdir -p "$(dirname {output.bai})" "$(dirname {log})"
        {{
          SRC="$(readlink -f "{params.src}")"
          FOUND=""
          for CAND in "$SRC.bai" "${{SRC%.bam}}.bai"; do
            if [ -s "$CAND" ]; then FOUND="$CAND"; break; fi
          done
          if [ -n "$FOUND" ]; then
            echo "reusing existing index (copied, not linked): $FOUND"
            cp -f "$FOUND" {output.bai}
            chmod u+w {output.bai}
          else
            echo "no sibling .bai for $SRC -- building one (this reads the whole BAM)"
            # Read through OUR symlink and write to an explicit path inside our
            # namespace, so the source directory is never a write target.
            samtools index -@ {threads} {input.bam} {output.bai}
          fi
        }} > {log} 2>&1
        """


# ---------------------------------------------------------------------
# Stage 1 — csaw / edgeR-TMM scale factors  (AUTOSOMES ONLY)
# ---------------------------------------------------------------------
rule csaw_counts:
    """10 kb bin counts over the csaw chromosome set, for the Gall (unsplit)
    BAMs of one normgroup. Emits the count matrix as data, not a plot -- the
    MA plots are drawn later, in the visualisation stage."""
    input:
        bams=lambda w: [P.bam(s) for s in SS.normgroup_samples(w.normgroup, allele="Gall")],
        bais=lambda w: [P.bam(s) + ".bai" for s in SS.normgroup_samples(w.normgroup, allele="Gall")],
        blacklist=P.blacklist(),
        chroms=P.chrom_subset("csaw"),
    output:
        counts=P.scalefactors("{normgroup}", "counts"),
        totals=P.scalefactors("{normgroup}", "totals"),
    params:
        bin_width=NORM["csaw"]["bin_width"],
        pe=NORM["csaw"]["pe"],
        max_frag=NORM["csaw"]["max_frag"],
        dedup=str(NORM["csaw"]["dedup"]).upper(),
        minq="NA" if NORM["csaw"]["minq"] is None else NORM["csaw"]["minq"],
        r_lib=NORM.get("r_library_override") or "",
        samples=lambda w: ",".join(SS.normgroup_samples(w.normgroup, allele="Gall")),
    log:
        P.log("csaw_counts", "{normgroup}"),
    benchmark:
        P.benchmark("csaw_counts", "{normgroup}")
    threads: 8
    resources:
        mem_mb=48000,
    script:
        "../scripts/R/csaw_counts.R"


rule csaw_norm_factors:
    """edgeR TMM on the bin counts -> LibSize, NormFactor, and the single
    derived quantity everything downstream uses:

        scaleFactor = 1e6 / (LibSize * NormFactor)
    """
    input:
        counts=P.scalefactors("{normgroup}", "counts"),
        totals=P.scalefactors("{normgroup}", "totals"),
    output:
        P.scalefactors("{normgroup}", "factors"),
    params:
        norm_method=NORM["csaw"]["norm_method"],
        formula=NORM["scale_factor_formula"],
        r_lib=NORM.get("r_library_override") or "",
    log:
        P.log("csaw_norm_factors", "{normgroup}"),
    resources:
        mem_mb=16000,
    script:
        "../scripts/R/csaw_norm_factors.R"


rule import_scale_factors:
    """D-09 alternative to computing the factors.

    csaw's restrict set is autosomal, so a chrX-only BAM subset cannot
    reproduce the ground-truth factors *at all* -- not approximately, but by
    construction. Rather than make everything downstream untestable, this rule
    injects factors recovered from the original bigWigs (the minimum positive
    value on an autosome is exactly the applied scale factor).

    The recovery is bit-exact for REPLICATE-LEVEL tracks only. Merged, 5 kb and
    AcMe3 tracks are always recomputed from the corrected replicates -- their
    values are averages and log-ratios, from which no single factor can be
    recovered (measured: 28% / 0.8% / 0.015% of bins are integer multiples).
    """
    input:
        fixture=lambda w: P.fixture("scalefactors", f"{w.normgroup}.tsv"),
    output:
        P.scalefactors("{normgroup}", "factors"),
    log:
        P.log("import_scale_factors", "{normgroup}"),
    shell:
        r"""
        mkdir -p "$(dirname {output})" "$(dirname {log})"
        cp -f {input.fixture} {output} 2> {log}
        echo "imported ground-truth scale factors from {input.fixture}" >> {log}
        """


# Both rules can produce `{normgroup}.factors.tsv`; config decides which wins.
if NORM["scale_factor_source"] == "fixture":

    ruleorder: import_scale_factors > csaw_norm_factors

else:

    ruleorder: csaw_norm_factors > import_scale_factors


rule propagate_scale_factors:
    """Broadcast each Gall scale factor onto that sample's Xa and Xi splits.

    D-03 lives here. Default: an explicit join on
    (mark, clone, condition, replicate). Legacy: the original's positional
    broadcast over an alphabetically sorted file list, reproduced exactly so
    the testing team can prove equivalence against the ground truth.
    """
    input:
        factors=P.scalefactors("{normgroup}", "factors"),
        sheet=config["paths"]["sheets"]["samples"],
    output:
        P.scalefactors("{normgroup}", "allsamples"),
    params:
        legacy=config["legacy"]["h3k27ac_positional_scalefactors"],
        # BUG-B: the legacy broadcast applies ONLY to these normgroups.
        legacy_normgroups=config["legacy"].get("legacy_normgroups", ["WT_H3K27ac"]),
        normgroup=lambda w: w.normgroup,
    log:
        P.log("propagate_scale_factors", "{normgroup}"),
    script:
        "../scripts/py/propagate_scale_factors.py"


# ---------------------------------------------------------------------
# Stage 2 — coverage tracks
# ---------------------------------------------------------------------
def _sf_table(wildcards):
    return P.scalefactors(SS.normgroup_of(wildcards.sample), "allsamples")


rule bam_coverage:
    """Replicate-level bigWig at 10 bp. NOT 20 bp -- the archaeologist verified
    the grid against the bigWig interval boundaries, and the briefs' "@20bp"
    refers to the merged tracks."""
    input:
        bam=P.bam("{sample}"),
        bai=P.bam("{sample}") + ".bai",
        factors=_sf_table,
    output:
        P.rep10("{mark}", "{sample}"),
    params:
        bin_size=BC["bin_size"],
        read_mode=lambda w: _bam_coverage_read_mode(w.mark),
        normalize=BC["normalize_using"],
    log:
        P.log("bam_coverage", "{mark}_{sample}"),
    benchmark:
        P.benchmark("bam_coverage", "{mark}_{sample}")
    threads: 8
    resources:
        mem_mb=12000,
    shell:
        r"""
        mkdir -p "$(dirname {output})" "$(dirname {log})"
        SF=$(awk -F'\t' -v s="{wildcards.sample}" \
                 'NR==1{{for(i=1;i<=NF;i++) if($i=="scaleFactor") c=i; next}}
                  $1==s {{print $c; found=1; exit}}
                  END{{if(!found) exit 1}}' {input.factors}) \
            || {{ echo "no scaleFactor row for {wildcards.sample} in {input.factors}" >&2; exit 1; }}
        echo "scaleFactor={wildcards.sample} -> $SF" > {log}
        bamCoverage \
            --bam "$(readlink -f {input.bam})" \
            --outFileName {output} \
            --binSize {params.bin_size} \
            --scaleFactor "$SF" \
            --normalizeUsing {params.normalize} \
            {params.read_mode} \
            --numberOfProcessors {threads} >> {log} 2>&1
        """


def _reps_of(wildcards):
    """Replicate bigWigs behind one merged track.

    `single_replicate_self_average: true` reproduces the original's handling of
    CL30 H3K27ac, which has only rep1: it was averaged with itself rather than
    passed through, which is a no-op numerically but keeps the code path -- and
    the output header -- identical.
    """
    reps = [P.rep10(_mark_of(s), s) for s in SS.reps_of_track(wildcards.track)]
    if len(reps) == 1 and MERGE.get("single_replicate_self_average", True):
        return reps * 2
    return reps


rule merge_reps_20:
    """Replicate mean at 20 bp -- the browser/stackup grid."""
    input:
        _reps_of,
    output:
        P.merged20("{mark}", "{track}"),
    params:
        bin_size=MERGE["bin_sizes"]["merged"],
        operation=MERGE["operation"],
    log:
        P.log("merge_reps_20", "{mark}_{track}"),
    threads: 8
    resources:
        mem_mb=12000,
    script:
        "../scripts/py/merge_bigwigs.py"


rule merge_reps_5k:
    """Replicate mean at 5 kb -- the valley-calling and METALoci grid."""
    input:
        _reps_of,
    output:
        P.merged5k("{mark}", "{track}"),
    params:
        bin_size=MERGE["bin_sizes"]["merged_5k"],
        operation=MERGE["operation"],
    log:
        P.log("merge_reps_5k", "{mark}_{track}"),
    threads: 8
    resources:
        mem_mb=8000,
    script:
        "../scripts/py/merge_bigwigs.py"


rule rebin_rep_5k:
    """A single replicate rebinned to 5 kb, without merging."""
    input:
        lambda w: [P.rep10(w.mark, w.sample)],
    output:
        P.rep5k("{mark}", "{sample}"),
    params:
        bin_size=MERGE["bin_sizes"]["rep_5k"],
        operation=MERGE["operation"],
    log:
        P.log("rebin_rep_5k", "{mark}_{sample}"),
    threads: 8
    resources:
        mem_mb=8000,
    script:
        "../scripts/py/merge_bigwigs.py"


rule bigwig_to_bedgraph:
    """5 kb bedGraph — METALoci's signal input format (15_bedgraph_all_FINAL.sh)."""
    input:
        P.merged5k("{mark}", "{track}"),
    output:
        P.bedgraph5k("{mark}", "{track}"),
    log:
        P.log("bigwig_to_bedgraph", "{mark}_{track}"),
    resources:
        mem_mb=4000,
    shell:
        r"""
        mkdir -p "$(dirname {output})" "$(dirname {log})"
        bigWigToBedGraph {input} {output} 2> {log}
        """


# ---------------------------------------------------------------------
# AcMe3 = log2(H3K27ac / H3K27me3)
# ---------------------------------------------------------------------
def _acme3_inputs(wildcards, resolution):
    """AcMe3_{clone}_{cond}_{allele} -> its H3K27ac and H3K27me3 partners."""
    suffix = wildcards.track.split("_", 1)[1]
    fn = P.merged20 if resolution == 20 else P.merged5k
    return {
        "ac": fn("H3K27ac", f"H3K27ac_{suffix}"),
        "me3": fn("H3K27me3", f"H3K27me3_{suffix}"),
    }


rule acme3_bw:
    """AcMe3 bigWig at 20 bp.

    The 20 bp / 5 kb asymmetry is deliberate and matches the ground truth: the
    bigWig is 20 bp (browser use) while the bedGraph is 5 kb (METALoci use).
    """
    input:
        unpack(lambda w: _acme3_inputs(w, 20)),
    output:
        P.acme3_bw20("{track}"),
    params:
        bin_size=ACME3["bw_bin_size"],
    log:
        P.log("acme3_bw", "{track}"),
    threads: 8
    resources:
        mem_mb=12000,
    shell:
        r"""
        mkdir -p "$(dirname {output})" "$(dirname {log})"
        bigwigCompare \
            --bigwig1 {input.ac} --bigwig2 {input.me3} \
            --operation log2 --binSize {params.bin_size} \
            --outFileName {output} --outFileFormat bigwig \
            --numberOfProcessors {threads} > {log} 2>&1
        """


rule acme3_bedgraph:
    """AcMe3 bedGraph at 5 kb — the METALoci AcMe3 signal."""
    input:
        unpack(lambda w: _acme3_inputs(w, 5000)),
    output:
        P.acme3_bg5k("{track}"),
    params:
        bin_size=ACME3["bedgraph_bin_size"],
    log:
        P.log("acme3_bedgraph", "{track}"),
    threads: 8
    resources:
        mem_mb=8000,
    shell:
        r"""
        mkdir -p "$(dirname {output})" "$(dirname {log})"
        bigwigCompare \
            --bigwig1 {input.ac} --bigwig2 {input.me3} \
            --operation log2 --binSize {params.bin_size} \
            --outFileName {output} --outFileFormat bedgraph \
            --numberOfProcessors {threads} > {log} 2>&1
        """


# ---------------------------------------------------------------------
# Cross-clone consensus
# ---------------------------------------------------------------------
def _consensus_members(wildcards):
    """Member tracks of one consensus, as bigWig paths.

    AcMe3 consensus is built from the 20 bp AcMe3 bigWigs, not from the 5 kb
    bedGraphs, which reproduces the original's mean-of-log2 (rather than
    log2-of-mean) ordering.
    """
    mark, group, allele = _parse_consensus(wildcards.name)
    members = SS.consensus_members(mark, group, allele)
    if mark == "AcMe3":
        return [P.acme3_bw20(t) for t in members]
    return [P.merged5k(mark, t) for t in members]


def _parse_consensus(name):
    """`H3K27me3_NodTAG-or-WT_Xi` -> (mark, group, allele)."""
    mark, rest = name.split("_", 1)
    group, allele = rest.rsplit("_", 1)
    return mark, group, allele


rule consensus_bedgraph:
    """Mean across the clones of a group, at 5 kb (17_consensus_bw_and_bedgraph.py).

    Uses `bigwigAverage`, which exists only from deepTools 3.5.5 — independent
    confirmation that 3.5.5 is the version that actually ran, despite the
    original env reporting 3.5.2.
    """
    input:
        _consensus_members,
    output:
        bedgraph=P.consensus5k("{name}"),
    params:
        bin_size=CONSENSUS["bin_size"],
    log:
        P.log("consensus_bedgraph", "{name}"),
    threads: 8
    resources:
        mem_mb=16000,
    shell:
        r"""
        mkdir -p "$(dirname {output.bedgraph})" "$(dirname {log})"
        bigwigAverage \
            --bigwigs {input} \
            --binSize {params.bin_size} \
            --outFileName {output.bedgraph} \
            --outFileFormat bedgraph \
            --numberOfProcessors {threads} > {log} 2>&1
        """
