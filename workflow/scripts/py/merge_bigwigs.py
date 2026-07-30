"""Average replicate bigWigs onto a target bin grid.

Used by three rules that differ only in bin size and input count:

  merge_reps_20   replicate mean @ 20 bp    (browser / stackup grid)
  merge_reps_5k   replicate mean @ 5 kb     (valley calling, METALoci)
  rebin_rep_5k    one replicate  @ 5 kb     (no merging, just rebinning)

Two behaviours worth knowing about, both inherited from the original:

*Single-replicate self-average.* CL30 H3K27ac has only rep1. The original
passed it to the averaging tool twice rather than special-casing it. Numerically
that is a no-op, but it keeps the output header and binning identical to the
two-replicate tracks, so downstream byte-comparisons behave the same. The rule
does the duplication; this script just notices and records it.

*Mean, for every mark including RNA-Seq.* The original never switched to sum or
max, so neither do we; `normalization.merge.operation` exists to make that
visible rather than to invite changing it.
"""

import os
import subprocess
import sys

inputs = list(snakemake.input)
output = snakemake.output[0]
log_path = snakemake.log[0]
bin_size = int(snakemake.params.bin_size)
operation = str(snakemake.params.operation)
threads = int(snakemake.threads)

os.makedirs(os.path.dirname(output), exist_ok=True)
os.makedirs(os.path.dirname(log_path), exist_ok=True)

if operation != "mean":
    raise SystemExit(
        f"normalization.merge.operation is {operation!r}, but the original "
        "pipeline only ever averaged. Changing this will not reproduce the "
        "published tracks."
    )

distinct = sorted(set(inputs))
self_averaged = len(distinct) == 1 and len(inputs) > 1

with open(log_path, "w") as log:
    log.write(f"inputs ({len(inputs)}):\n")
    for path in inputs:
        log.write(f"    {path}\n")
    if self_averaged:
        log.write(
            "NOTE: a single replicate averaged with itself -- reproduces the "
            "original's handling of single-replicate tracks (numerically a "
            "no-op, kept for header/binning parity).\n"
        )
    log.write(f"binSize={bin_size}  operation={operation}  threads={threads}\n")
    log.flush()

    cmd = [
        "bigwigAverage",
        "--bigwigs", *inputs,
        "--binSize", str(bin_size),
        "--outFileName", output,
        "--outFileFormat", "bigwig",
        "--numberOfProcessors", str(threads),
    ]
    log.write("$ " + " ".join(cmd) + "\n")
    log.flush()

    result = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
    if result.returncode != 0:
        log.write(f"\nbigwigAverage failed with exit code {result.returncode}\n")
        sys.exit(result.returncode)
