"""Fragment-size distribution as a table (QC).

The original ran ATACseqQC::fragSizeDist and wrote a PDF straight out. We emit
the counts; panel_fragment_sizes draws them. Same numbers, one producer.
"""
import collections

import pandas as pd
import pysam

bam_path = snakemake.input.bam
out_tsv = snakemake.output.tsv
sample = snakemake.wildcards.sample
MAX_LEN = 1000

counts: collections.Counter = collections.Counter()
with pysam.AlignmentFile(bam_path, "rb") as bam:
    for read in bam.fetch(until_eof=True):
        if not read.is_proper_pair or read.is_unmapped or read.is_secondary:
            continue
        if not read.is_read1:          # count each fragment once
            continue
        size = abs(read.template_length)
        if 0 < size <= MAX_LEN:
            counts[size] += 1

pd.DataFrame(
    {"sample": sample,
     "fragment_size": list(range(1, MAX_LEN + 1)),
     "count": [counts.get(i, 0) for i in range(1, MAX_LEN + 1)]}
).to_csv(out_tsv, sep="\t", index=False)
