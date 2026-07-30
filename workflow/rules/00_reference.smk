# =====================================================================
#  00_reference.smk — LAYER 0
#
#  Genome assets every other layer depends on. Nothing here reads a BAM
#  or a cooler, so this is the first thing to run and the cheapest thing
#  to re-run.
#
#  Each asset is either downloaded (config `genome.<x>: auto`) or taken
#  from a local path the user supplied. `_source()` decides which, so no
#  rule below has to care.
#
#  Replaces: manual downloads, 00_02_blacklist.py, 07_get_CTCF_motifs.ipynb,
#            09_get_TSS.sh, 01_07 (noncoding-GTF part), a hand-written ROI TSV.
# =====================================================================


def _source(key):
    """Local path for a genome asset, or None when it must be downloaded."""
    value = config["genome"].get(key, "auto")
    return None if value in (None, "auto", "") else str(value)


def _url(key):
    return config["genome"]["urls"][key]


# ---------------------------------------------------------------------
# Downloads. A user-supplied local file takes the same code path as a
# downloaded one, so switching between them never changes a downstream rule.
# ---------------------------------------------------------------------
rule get_chrom_sizes:
    """UCSC mm10 chrom.sizes (65 sequences)."""
    output:
        protected(P.chrom_sizes()),
    params:
        src=_source("chrom_sizes") or "",
        url=_url("chrom_sizes"),
    log:
        P.log("get_chrom_sizes"),
    shell:
        r"""
        mkdir -p "$(dirname {output})" "$(dirname {log})"
        if [ -n "{params.src}" ]; then
            cp -f "{params.src}" {output} 2> {log}
        else
            curl -fsSL "{params.url}" -o {output} 2> {log}
        fi
        """


rule get_blacklist:
    """ENCODE mm10 blacklist v2. Masked out of csaw counting and valley calling."""
    output:
        protected(P.blacklist()),
    params:
        src=_source("blacklist") or "",
        url=_url("blacklist"),
    log:
        P.log("get_blacklist"),
    shell:
        r"""
        mkdir -p "$(dirname {output})" "$(dirname {log})"
        if [ -n "{params.src}" ]; then
            cp -f "{params.src}" {output} 2> {log}
        else
            curl -fsSL "{params.url}" | gzip -dc > {output} 2> {log}
        fi
        """


rule get_fasta:
    """mm10 FASTA — used only for GC phasing of the compartment eigenvectors."""
    output:
        fa=protected(P.fasta()),
        fai=protected(P.fasta() + ".fai"),
    params:
        src=_source("fasta") or "",
        url=_url("fasta"),
    log:
        P.log("get_fasta"),
    resources:
        mem_mb=4000,
    shell:
        r"""
        mkdir -p "$(dirname {output.fa})" "$(dirname {log})"
        {{
            if [ -n "{params.src}" ]; then
                cp -f "{params.src}" {output.fa}
            else
                curl -fsSL "{params.url}" | gzip -dc > {output.fa}
            fi
            samtools faidx {output.fa}
        }} 2> {log}
        """


rule get_gtf:
    """Ensembl GRCm38 release 102 — the release the original code actually read
    (the archaeology notes OQ-7; methods.pdf cites release 100, which is not what ran)."""
    output:
        protected(P.gtf()),
    params:
        src=_source("gtf") or "",
        url=_url("gtf"),
    log:
        P.log("get_gtf"),
    resources:
        mem_mb=4000,
    shell:
        r"""
        mkdir -p "$(dirname {output})" "$(dirname {log})"
        if [ -n "{params.src}" ]; then
            cp -f "{params.src}" {output} 2> {log}
        else
            curl -fsSL "{params.url}" | gzip -dc > {output} 2> {log}
        fi
        """


rule get_gencode_gtf:
    """GENCODE vM25 — TSS extraction only. Kept separate from the Ensembl GTF
    because the original used a different annotation for TSS than for genes."""
    output:
        protected(P.genome("gencode.vM25.annotation.gtf")),
    params:
        src=_source("gencode_gtf") or "",
        url=_url("gencode_gtf"),
    log:
        P.log("get_gencode_gtf"),
    resources:
        mem_mb=4000,
    shell:
        r"""
        mkdir -p "$(dirname {output})" "$(dirname {log})"
        if [ -n "{params.src}" ]; then
            cp -f "{params.src}" {output} 2> {log}
        else
            curl -fsSL "{params.url}" | gzip -dc > {output} 2> {log}
        fi
        """


rule get_gff3_x:
    """Ensembl chrX GFF3 — the allelic-GTF builder needs GFF3-only attributes."""
    output:
        protected(P.genome("Mus_musculus.GRCm38.102.chromosome.X.gff3")),
    params:
        src=_source("gff3_x") or "",
        url=_url("gff3_x"),
    log:
        P.log("get_gff3_x"),
    shell:
        r"""
        mkdir -p "$(dirname {output})" "$(dirname {log})"
        if [ -n "{params.src}" ]; then
            cp -f "{params.src}" {output} 2> {log}
        else
            curl -fsSL "{params.url}" | gzip -dc > {output} 2> {log}
        fi
        """


# ---------------------------------------------------------------------
# Derived assets
# ---------------------------------------------------------------------
rule make_noncoding_gtf:
    """Non-coding feature track for coolbox's GFFNC track (Fig 4e, EFig 6d).

    Ports the first block of `01_07_get_allelic_ratio_gtf.py`, which is the
    original's sole producer of `GRCm38.102_NC.gtf`.

    An earlier implementation of this rule derived the file from the
    genome-wide GTF by excluding `gene_biotype "protein_coding"`. That was
    wrong on four counts and RD-3 caught it while porting `allelic_gtf`:

      source     chrX GFF3, not the genome-wide GTF
      filter     a GFF3 `feature` WHITELIST, not a GTF biotype exclusion
      scope      chrX only
      naming     Ensembl `X`, never renamed to `chrX`

    The last one is the quiet killer: renaming to `chrX` here would leave
    coolbox with an empty non-coding track and no error.

    The rule keeps living in 00_reference.smk rather than moving into
    31_hic_features.smk with the rest of `01_07`, because it is a pure
    function of the reference annotation with no sample dependency. Two rules
    may not declare the same output, so `allelic_gtf` deliberately does not
    produce it -- see the module notes section 8.
    """
    input:
        gff3=P.genome("Mus_musculus.GRCm38.102.chromosome.X.gff3"),
    output:
        gtf=P.gtf_noncoding(),
        bgz=P.gtf_noncoding() + ".bgz",
        tbi=P.gtf_noncoding() + ".bgz.tbi",
    log:
        P.log("make_noncoding_gtf"),
    benchmark:
        P.benchmark("make_noncoding_gtf")
    resources:
        mem_mb=8000,
    script:
        "../scripts/py/make_noncoding_gtf.py"


rule make_ctcf_motifs:
    """CTCF motif occurrences (JASPAR MA0139.1), genome-wide and chrX-only.

    Prefers the committed fixture. The original built this in an R/AnnotationHub
    session that is not reproducible offline, and the motif set is an input to
    the boundary analysis, so it is versioned rather than recomputed. Set
    `fixtures.ctcf_motifs: recompute` to scan the FASTA instead.
    """
    input:
        fa=P.fasta(),
        fai=P.fasta() + ".fai",
    output:
        genome_wide=P.ctcf_motifs(),
        chrx=P.ctcf_motifs(chrx_only=True),
    params:
        fixture=P.fixture("CTCF_mm10_X_only.bed.gz"),
        mode=config.get("fixtures", {}).get("ctcf_motifs", "fixture"),
        jaspar=config["genome"]["jaspar_motif"],
        chrx=config["loci"]["Mecp2"]["chrom"],
    log:
        P.log("make_ctcf_motifs"),
    threads: 4
    resources:
        mem_mb=8000,
    script:
        "../scripts/py/make_ctcf_motifs.py"


rule make_tss:
    """TSS BED from GENCODE vM25 (09_get_TSS.sh)."""
    input:
        gtf=P.genome("gencode.vM25.annotation.gtf"),
    output:
        P.tss(),
    log:
        P.log("make_tss"),
    resources:
        mem_mb=4000,
    shell:
        r"""
        mkdir -p "$(dirname {output})" "$(dirname {log})"
        awk 'BEGIN{{OFS="\t"}}
             $3=="transcript" {{
                 if ($7=="+") {{ s=$4-1; e=$4 }} else {{ s=$5-1; e=$5 }}
                 match($0, /gene_name "[^"]+"/)
                 name = (RSTART>0) ? substr($0, RSTART+11, RLENGTH-12) : "."
                 print $1, s, e, name, 0, $7
             }}' {input.gtf} \
          | sort -k1,1 -k2,2n -u > {output} 2> {log}
        """


rule write_regions_of_interest:
    """The `loci:` config block as a TSV, so scripts read one file rather than
    re-parsing YAML. One row per (locus, roi)."""
    output:
        P.roi_table(),
    log:
        P.log("write_regions_of_interest"),
    run:
        import os

        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        os.makedirs(os.path.dirname(log[0]), exist_ok=True)
        n = 0
        with open(output[0], "w") as fh:
            fh.write("locus\troi\tchrom\tstart\tend\tdisplay_name\n")
            for locus in LOCI:
                block = config["loci"][locus]
                shown = config["loci"]["display_names"].get(locus, locus)
                for roi in ROI_TYPES:
                    start, end = block["regions"][roi]
                    fh.write(f"{locus}\t{roi}\t{block['chrom']}\t{start}\t{end}\t{shown}\n")
                    n += 1
        with open(log[0], "w") as fh:
            fh.write(f"wrote {n} region rows\n")


rule make_ice_blacklist:
    """The two flanks of the locus chromosome, so ICE balances on the captured
    region only (00_02_blacklist.py).

    **Exactly two rows, on the locus chromosome only.** An earlier version of
    this rule also emitted one whole-chromosome row for every other chromosome
    in mm10.chrom.sizes, on the reasoning that "everything outside the capture
    window" should be masked. That is wrong and it is fatal: these Capture Hi-C
    coolers contain a single chromosome (`nchroms == 1`), and
    `cooler balance --blacklist` calls `bedslice()` per region, so the first
    non-chrX row raises `ValueError: Unknown sequence label: chr1` and every
    `balance_zoomify` job dies -- taking mcools, loops, compartments, pile-ups
    and METALoci with it.

    The original writes two rows. So do we.
    """
    input:
        chrom_sizes=P.chrom_sizes(),
    output:
        P.ice_blacklist("{locus}"),
    log:
        P.log("make_ice_blacklist", "{locus}"),
    run:
        import os

        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        os.makedirs(os.path.dirname(log[0]), exist_ok=True)

        block = config["loci"][wildcards.locus]
        keep_chrom = block["chrom"]
        # The original masks the flanks of the `full` region, not of
        # `view_window`. They are currently equal for both loci, but `full` is
        # what 00_02_blacklist.py reads, so use it.
        win_start, win_end = block["regions"]["full"]

        sizes = {}
        with open(input.chrom_sizes) as fh:
            for line in fh:
                if line.strip():
                    name, length = line.split()[:2]
                    sizes[name] = int(length)

        if keep_chrom not in sizes:
            raise WorkflowError(
                f"locus {wildcards.locus} is on {keep_chrom}, which is not in "
                f"{input.chrom_sizes}"
            )

        length = sizes[keep_chrom]
        rows = []
        if win_start > 0:
            rows.append((keep_chrom, 0, win_start))
        if win_end < length:
            rows.append((keep_chrom, win_end, length))

        with open(output[0], "w") as fh:
            for chrom, start, end in rows:
                fh.write(f"{chrom}\t{start}\t{end}\n")
        with open(log[0], "w") as fh:
            fh.write(
                f"{wildcards.locus}: kept {keep_chrom}:{win_start}-{win_end}, "
                f"blacklisted {len(rows)} interval(s) on {keep_chrom} only\n"
            )


rule write_chrom_subsets:
    """The two independent chromosome sets (D-04).

    `csaw` must never contain chrX — allelic and clone-genotype differences
    corrupt TMM. lib/samples.py asserts that at load time; this rule only
    materialises both sets for the shell-level tools.
    """
    input:
        chrom_sizes=P.chrom_sizes(),
    output:
        P.chrom_subset("{scope}"),
    wildcard_constraints:
        scope="csaw|downstream",
    log:
        P.log("write_chrom_subsets", "{scope}"),
    run:
        import os

        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        os.makedirs(os.path.dirname(log[0]), exist_ok=True)

        wanted = list(config["chromosomes"][wildcards.scope])
        if wildcards.scope == "downstream":
            wanted += list(config["chromosomes"].get("valley_control") or [])

        available = set()
        with open(input.chrom_sizes) as fh:
            for line in fh:
                if line.strip():
                    available.add(line.split()[0])

        missing = [c for c in wanted if c not in available]
        if missing:
            raise WorkflowError(
                f"chromosomes {missing} appear in config.chromosomes."
                f"{wildcards.scope} but not in {input.chrom_sizes}"
            )

        with open(output[0], "w") as fh:
            for chrom in wanted:
                fh.write(f"{chrom}\n")
        with open(log[0], "w") as fh:
            fh.write(f"{wildcards.scope}: {' '.join(wanted)}\n")


# The `reference` convenience target lives in the Snakefile alongside the other
# section aliases, so all of them are visible in one place.
