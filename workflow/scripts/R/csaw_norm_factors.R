# =====================================================================
#  csaw_norm_factors.R — edgeR TMM on the binned counts
#
#  Produces the one quantity the rest of the pipeline consumes:
#
#      scaleFactor = 1e6 / (LibSize * NormFactor)
#
#  That product is also what can be *recovered* from a finished bigWig
#  (the minimum positive value on an autosome equals it exactly), which is
#  what makes decision D-09's fixture path possible. LibSize and
#  NormFactor are emitted separately as well, but only the product is
#  recoverable — so the fixture format carries `scaleFactor` alone.
#
#  Replaces the normFactors half of 12_normalization_*.r
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
  library(edgeR)
  library(SummarizedExperiment)
})

message("csaw   ", as.character(packageVersion("csaw")))
message("edgeR  ", as.character(packageVersion("edgeR")))

counts_df <- read.delim(gzfile(snakemake@input[["counts"]]), check.names = FALSE)
totals_df <- read.delim(snakemake@input[["totals"]], stringsAsFactors = FALSE)

meta_cols <- c("chrom", "start", "end")
sample_ids <- setdiff(colnames(counts_df), meta_cols)
stopifnot(setequal(sample_ids, totals_df$sample))
totals_df <- totals_df[match(sample_ids, totals_df$sample), ]

counts <- as.matrix(counts_df[, sample_ids, drop = FALSE])
mode(counts) <- "numeric"
message("bins: ", nrow(counts), "  samples: ", ncol(counts))

# Rebuild the minimal SummarizedExperiment csaw::normFactors expects.
se <- SummarizedExperiment(assays = list(counts = counts))
se$totals <- totals_df$LibSize

norm_method <- snakemake@params[["norm_method"]]
if (!identical(norm_method, "TMM")) {
  stop("normalization.csaw.norm_method is '", norm_method,
       "', but the published tracks were produced with TMM.")
}

norm_factors <- csaw::normFactors(se, se.out = FALSE)
message("NormFactor: ", paste(round(norm_factors, 6), collapse = " "))

lib_size <- totals_df$LibSize
scale_factor <- 1e6 / (lib_size * norm_factors)

out <- data.frame(
  sample      = sample_ids,
  LibSize     = lib_size,
  NormFactor  = norm_factors,
  scaleFactor = scale_factor,
  stringsAsFactors = FALSE
)
write.table(
  out, snakemake@output[[1]],
  sep = "\t", quote = FALSE, row.names = FALSE
)

message("scaleFactor:")
for (i in seq_len(nrow(out))) {
  message("  ", out$sample[i], "  ", format(out$scaleFactor[i], digits = 10))
}

sink(type = "message")
sink(type = "output")
close(log_file)
