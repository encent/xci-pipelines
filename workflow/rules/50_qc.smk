# =====================================================================
#  50_qc.smk — optional quality control.
#
#  OFF by default (`qc.enabled: false`). Turning it on roughly doubles the
#  runtime of the signal branch, because FastQC and MACS2 both re-read
#  763 GB of BAM.
#
#  Replaces: 00_do_fastqc.sh, 01_calculate_fragment_sizes.sh,
#            02_merge_bams.sh, 03_call_peaks.sh, 04_plot_frip.py,
#            05_plot_correlation.sh, 06_plot_venn.py,
#            08/10/11_plot_*_enrichment.py
#
#  ------------------------------------------------------------------
#  TWO THINGS TO KNOW BEFORE TOUCHING THIS MODULE
#  ------------------------------------------------------------------
#  1. THE PEAKS HERE ARE QC-ONLY. They are NOT the CTCF peaks used by the
#     valley/boundary analysis. Those come from the external delivery
#     resources/fixtures/CTCFpeak_per_clone/*_consensusPeaks_*.bed.
#     Never wire MACS2 output into 30_valleys. (Plan correction N-5.)
#
#  2. `signal_enrichment` is NOT gated on qc.enabled. It feeds Fig 1e,
#     Fig 2a and Fig 2d, which are pipeline-produced paper panels even
#     though they sit outside the 25-panel target list. (Open item O-12.)
#
#  ------------------------------------------------------------------
#  NO FIGURES HERE
#  ------------------------------------------------------------------
#  Every rule emits tables or .npz. The FRiP dot plot, the replicate
#  Venns, the correlation heatmaps, the fragment-size distributions and
#  the enrichment heatmaps are all drawn by 90_visualise.smk.
#  tests/lint_figures_only_in_viz.py enforces this.
# =====================================================================

QC = config.get("qc", {})
QC_ON = bool(QC.get("enabled", False))
MACS2 = QC.get("macs2", {})
ENRICH = QC.get("enrichment", {"flank": 200000, "nbins": 100})

# MACS2 parameters, verified directly from
# Antonia_Edith_01/normalization/03_call_peaks.sh.
# the archaeology notes section D.1 is WRONG here (it carried methods.pdf's
# "--shift -100 --extsize 200"); the design review section 1 (B-1) is the authority.
#   -f BAMPE  -g mm  --keep-dup all  --nomodel  --extsize 180   (NO --shift)
#   CTCF, Rad21 : -q 0.01                             -> narrowPeak
#   H3K27ac     : -q 0.05                             -> narrowPeak
#   H3K27me3    : -q 0.05 --broad --broad-cutoff 0.1  -> broadPeak
# `-g mm` is MACS2's alias for exactly 1.87e9, which is why config stores the
# number but the shell emits `mm`: the command line then matches the original
# byte for byte.
#
# Deviation D-12: the surviving on-disk WT H3K27ac peaks are 30 .broadPeak +
# 30 .gappedPeak, which the saved script cannot produce -- a commented-out
# `elif` branch was live at run time. Two generations exist. We emit what the
# saved script says. QC-only; do not chase it.


def _macs2_qvalue(mark):
    return MACS2.get("qvalue", {}).get(mark, 0.05)


def _macs2_broad(mark):
    return mark in MACS2.get("broad_marks", ["H3K27me3"])


def _qc_tracks():
    """Tracks peak calling runs on: Gall only (allelic peaks are not meaningful)."""
    return [t for t in SS.tracks if naming.parse_track(t)["allele"] == "Gall"]


def _enrich_anchor(mark):
    return "motif" if mark in ("CTCF", "Rad21") else "tss"


# ---------------------------------------------------------------- FastQC --
rule fastqc:
    """FastQC on one BAM. Replaces 00_do_fastqc.sh (`fastqc -q -t 40`)."""
    input:
        bam=lambda w: P.bam(w.sample),
    output:
        zip=P.qc("fastqc", "{sample}_fastqc.zip"),
        html=P.qc("fastqc", "{sample}_fastqc.html"),
    log:
        P.log("fastqc", "{sample}"),
    benchmark:
        P.benchmark("fastqc", "{sample}")
    threads: 4
    resources:
        mem_mb=config["resources"]["mem_mb"]["small"],
    shell:
        r"""
        fastqc -q -t {threads} -o "$(dirname {output.zip})" {input.bam} > {log} 2>&1
        """


# -------------------------------------------------------- fragment sizes --
rule fragment_sizes:
    """Fragment-size distribution as a TABLE.

    The original (01_calculate_fragment_sizes.sh -> ATACseqQC::fragSizeDist)
    emitted a PDF directly. We emit the counts; panel_fragment_sizes draws it.
    """
    input:
        bam=lambda w: P.bam(w.sample),
        bai=lambda w: P.bam(w.sample) + ".bai",
    output:
        tsv=P.qc("fragment_sizes", "{sample}.tsv"),
    log:
        P.log("fragment_sizes", "{sample}"),
    benchmark:
        P.benchmark("fragment_sizes", "{sample}")
    threads: 2
    resources:
        mem_mb=config["resources"]["mem_mb"]["small"],
    script:
        "../scripts/py/qc_fragment_sizes.py"


# ------------------------------------------------------------ peak calling --
rule merge_bams_for_peaks:
    """samtools merge over a track's replicates. Peak calling only."""
    input:
        bams=lambda w: [P.bam(s) for s in SS.reps_of_track(w.track)],
    output:
        bam=temp(P.qc("merged_bams", "{track}.bam")),
        bai=temp(P.qc("merged_bams", "{track}.bam.bai")),
    log:
        P.log("merge_bams_for_peaks", "{track}"),
    benchmark:
        P.benchmark("merge_bams_for_peaks", "{track}")
    threads: 8
    resources:
        mem_mb=config["resources"]["mem_mb"]["medium"],
    shell:
        r"""
        set -euo pipefail
        {{
          n=$(echo {input.bams} | wc -w)
          if [ "$n" -eq 1 ]; then
            cp {input.bams} {output.bam}
          else
            samtools merge -@ {threads} -f {output.bam} {input.bams}
          fi
          samtools index -@ {threads} {output.bam}
        }} > {log} 2>&1
        """


rule call_peaks:
    """MACS2, no control. QC ONLY -- see the module header."""
    input:
        bam=P.qc("merged_bams", "{track}.bam"),
        bai=P.qc("merged_bams", "{track}.bam.bai"),
    output:
        peaks=P.qc("peaks", "{track}_peaks.bed"),
    params:
        qvalue=lambda w: _macs2_qvalue(naming.parse_track(w.track)["mark"]),
        broad=lambda w: _macs2_broad(naming.parse_track(w.track)["mark"]),
        broad_cutoff=MACS2.get("broad_cutoff", 0.1),
        extsize=MACS2.get("extsize", 180),
        outdir=P.qc("peaks"),
    log:
        P.log("call_peaks", "{track}"),
    benchmark:
        P.benchmark("call_peaks", "{track}")
    threads: 2
    resources:
        mem_mb=config["resources"]["mem_mb"]["medium"],
    shell:
        r"""
        set -euo pipefail
        {{
          BROAD=""
          if [ "{params.broad}" = "True" ]; then
            BROAD="--broad --broad-cutoff {params.broad_cutoff}"
          fi
          # -g mm == 1.87e9; written as `mm` to match the original command line.
          macs2 callpeak \
            -t {input.bam} \
            -f BAMPE \
            -g mm \
            -n {wildcards.track} \
            --outdir {params.outdir} \
            -q {params.qvalue} \
            --keep-dup all \
            --nomodel \
            --extsize {params.extsize} \
            $BROAD
          # normalise the extension so downstream rules need not branch
          if [ -f "{params.outdir}/{wildcards.track}_peaks.broadPeak" ]; then
            cp "{params.outdir}/{wildcards.track}_peaks.broadPeak" {output.peaks}
          else
            cp "{params.outdir}/{wildcards.track}_peaks.narrowPeak" {output.peaks}
          fi
        }} > {log} 2>&1
        """


rule frip:
    """Fraction of reads in peaks, as a table.

    FRiP = (bedtools intersect -u -a bam -b peaks | wc -l) / (samtools view -c -F 0x4)
    """
    input:
        bams=[P.qc("merged_bams", f"{t}.bam") for t in _qc_tracks()],
        peaks=[P.qc("peaks", f"{t}_peaks.bed") for t in _qc_tracks()],
    output:
        tsv=P.qc("frip.tsv"),
    log:
        P.log("frip"),
    benchmark:
        P.benchmark("frip")
    threads: 4
    resources:
        mem_mb=config["resources"]["mem_mb"]["medium"],
    script:
        "../scripts/py/qc_frip.py"


rule peak_overlap:
    """Replicate peak-overlap counts, for the Venn the viz stage draws."""
    input:
        peaks=[P.qc("peaks", f"{t}_peaks.bed") for t in _qc_tracks()],
    output:
        tsv=P.qc("peak_overlap.tsv"),
    log:
        P.log("peak_overlap"),
    threads: 2
    resources:
        mem_mb=config["resources"]["mem_mb"]["small"],
    script:
        "../scripts/py/qc_peak_overlap.py"


# ----------------------------------------------------------- correlation --
rule bigwig_correlation:
    """multiBigwigSummary over one normalisation group -> .npz.

    The original (05_plot_correlation.sh) went straight to plotPCA /
    plotCorrelation. We stop at the matrix; panel_correlation draws it.
    """
    input:
        bws=lambda w: [
            P.rep10(SS.normgroup_mark(w.normgroup), s)
            for s in SS.normgroup_samples(w.normgroup)
        ],
    output:
        npz=P.qc("correlation", "{normgroup}.npz"),
    log:
        P.log("bigwig_correlation", "{normgroup}"),
    benchmark:
        P.benchmark("bigwig_correlation", "{normgroup}")
    threads: 16
    resources:
        mem_mb=config["resources"]["mem_mb"]["large"],
    shell:
        r"""
        multiBigwigSummary bins \
            -b {input.bws} \
            -o {output.npz} \
            --numberOfProcessors {threads} > {log} 2>&1
        """


# ------------------------------------------------------------ enrichment --
# NOT gated on qc.enabled -- feeds Fig 1e, Fig 2a, Fig 2d (open item O-12).
rule signal_enrichment:
    """Stackup of one mark over its anchor set. Replaces 08/10/11_plot_*.py.

    bbi.stackup(bw, chrom, centre +/- flank, bins=nbins), then a mean profile.
    Anchors: CTCF motifs for CTCF/Rad21, TSS for the histone marks.
    """
    input:
        bws=lambda w: [P.merged20(w.mark, t) for t in SS.tracks_of(mark=w.mark, allele="Gall")],
        anchors=lambda w: P.ctcf_motifs() if w.anchor == "motif" else P.tss(),
        chrom_sizes=P.chrom_sizes(),
    output:
        npz=P.qc("enrichment", "{mark}_{anchor}.npz"),
    params:
        flank=ENRICH.get("flank", 200000),
        nbins=ENRICH.get("nbins", 100),
    log:
        P.log("signal_enrichment", "{mark}_{anchor}"),
    benchmark:
        P.benchmark("signal_enrichment", "{mark}_{anchor}")
    threads: 8
    resources:
        mem_mb=config["resources"]["mem_mb"]["large"],
    script:
        "../scripts/py/qc_signal_enrichment.py"


# --------------------------------------------------------------- targets --
def qc_targets():
    """Everything 90_visualise may ask of this module.

    Only the enrichment stackups unless qc.enabled -- those are always on
    because they feed Fig 1e / 2a / 2d.
    """
    out = [
        P.qc("enrichment", f"{mark}_{_enrich_anchor(mark)}.npz")
        for mark in SS.marks
        if mark != "RNA-Seq"
    ]
    if not QC_ON:
        return out
    out += [P.qc("fastqc", f"{s}_fastqc.zip") for s in SS.all_samples]
    out += [P.qc("fragment_sizes", f"{s}.tsv") for s in SS.all_samples]
    out += [P.qc("frip.tsv"), P.qc("peak_overlap.tsv")]
    out += [P.qc("correlation", f"{ng}.npz") for ng in SS.normgroups]
    return out
