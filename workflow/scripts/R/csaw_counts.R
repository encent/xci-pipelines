# =====================================================================
#  csaw_counts.R — 10 kb bin counts for one normgroup
#
#  Counts are taken on AUTOSOMES ONLY. chrX is excluded deliberately:
#  the clones differ in X genotype and the allelic split makes X coverage
#  non-comparable between samples, so including it corrupts the TMM
#  estimate that the whole pipeline hangs off (decision D-04).
#  lib/samples.py refuses to start if chrX appears in the csaw set.
#
#  This rule emits the count matrix as DATA. The MA plots that the
#  original drew here are produced later, in the visualisation stage —
#  no rule outside 90_visualise.smk writes an image.
#
#  Replaces the counting half of 12_normalization_*.r
# =====================================================================

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

r_lib <- snakemake@params[["r_lib"]]
if (nzchar(r_lib)) {
  .libPaths(c(r_lib, .libPaths()))
  message("R library override: ", r_lib)
}

suppressPackageStartupMessages({
  library(csaw)
  library(rtracklayer)
})

message("csaw      ", as.character(packageVersion("csaw")))
message("edgeR     ", as.character(packageVersion("edgeR")))

bam_files <- unlist(snakemake@input[["bams"]])
sample_ids <- strsplit(snakemake@params[["samples"]], ",")[[1]]
stopifnot(length(bam_files) == length(sample_ids))

restrict <- readLines(snakemake@input[["chroms"]])
restrict <- restrict[nzchar(restrict)]
if ("chrX" %in% restrict) {
  stop("chrX is in the csaw restrict set. See decision D-04 — this corrupts TMM.")
}
message("restrict: ", paste(restrict, collapse = " "))

blacklist <- rtracklayer::import(snakemake@input[["blacklist"]])
message("blacklist intervals: ", length(blacklist))

minq_raw <- snakemake@params[["minq"]]
minq <- if (identical(as.character(minq_raw), "NA")) NA_integer_ else as.integer(minq_raw)

param <- csaw::readParam(
  pe       = snakemake@params[["pe"]],
  max.frag = as.integer(snakemake@params[["max_frag"]]),
  minq     = minq,
  dedup    = as.logical(snakemake@params[["dedup"]]),
  restrict = restrict,
  discard  = blacklist
)

bin_width <- as.integer(snakemake@params[["bin_width"]])
message("counting ", length(bam_files), " BAMs in ", bin_width, " bp bins")

binned <- csaw::windowCounts(
  bam.files = bam_files,
  bin       = TRUE,
  width     = bin_width,
  param     = param
)

counts <- SummarizedExperiment::assay(binned)
colnames(counts) <- sample_ids
regions <- SummarizedExperiment::rowRanges(binned)

out <- data.frame(
  chrom = as.character(GenomicRanges::seqnames(regions)),
  start = GenomicRanges::start(regions) - 1L,   # BED-style, 0-based
  end   = GenomicRanges::end(regions),
  counts,
  check.names = FALSE
)

gz <- gzfile(snakemake@output[["counts"]], "wt")
write.table(out, gz, sep = "\t", quote = FALSE, row.names = FALSE)
close(gz)
message("wrote ", nrow(out), " bins x ", length(sample_ids), " samples")

totals <- data.frame(
  sample  = sample_ids,
  LibSize = binned$totals,
  bam     = bam_files,
  stringsAsFactors = FALSE
)
write.table(
  totals, snakemake@output[["totals"]],
  sep = "\t", quote = FALSE, row.names = FALSE
)
message("library sizes: ", paste(binned$totals, collapse = " "))

sink(type = "message")
sink(type = "output")
close(log_file)
