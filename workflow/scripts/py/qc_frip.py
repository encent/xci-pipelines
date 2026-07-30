"""FRiP (fraction of reads in peaks) as a table (QC).

FRiP = reads overlapping a peak / total mapped reads, exactly as the original:
  bedtools intersect -u -a <bam> -b <peaks> | wc -l
  samtools view -c -F 0x4 <bam>
"""
import os
import subprocess

import pandas as pd

rows = []
for bam, peaks in zip(snakemake.input.bams, snakemake.input.peaks):
    track = os.path.basename(bam)[: -len(".bam")]
    in_peaks = int(
        subprocess.run(
            f"bedtools intersect -u -a {bam} -b {peaks} | wc -l",
            shell=True, capture_output=True, text=True, check=True,
        ).stdout.strip()
    )
    total = int(
        subprocess.run(
            ["samtools", "view", "-c", "-F", "0x4", bam],
            capture_output=True, text=True, check=True,
        ).stdout.strip()
    )
    n_peaks = sum(1 for ln in open(peaks) if ln.strip() and not ln.startswith("#"))
    rows.append(
        {"track": track, "reads_in_peaks": in_peaks, "total_reads": total,
         "n_peaks": n_peaks, "frip": (in_peaks / total) if total else float("nan")}
    )

pd.DataFrame(rows).sort_values("track").to_csv(snakemake.output.tsv, sep="\t", index=False)
