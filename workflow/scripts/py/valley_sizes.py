"""Valley size distribution across every called set.

Replaces ``TEST_03_02_valley_sizes.ipynb``, whose output directory
(``merged_5000/valleys_distribution``) never made it to disk, and the size
histograms later re-implemented in ``xci_valleys_check_LAST.ipynb``. Not a paper
panel; it is the QC view that makes the ``chrX:0-3,285,000`` assembly-gap valley
visible as a 3.3 Mb outlier in every clone.

One row per valley, tagged with the track it came from, so the visualisation
stage can facet however it likes without re-reading 30 BED files::

    track mark clone condition allele chrom start end size size_class
"""

import os
import re
import sys

import pandas as pd

log_path = snakemake.log[0]
os.makedirs(os.path.dirname(log_path), exist_ok=True)
os.makedirs(os.path.dirname(snakemake.output.tsv), exist_ok=True)

SIZE_CLASSES = [
    (25_000, "under 25 Kbp"),
    (50_000, "25-50 Kbp"),
    (100_000, "50-100 Kbp"),
    (500_000, "100-500 Kbp"),
]


def size_class(width: int) -> str:
    for limit, label in SIZE_CLASSES:
        if width <= limit:
            return label
    return "over 500 Kbp"


with open(log_path, "w") as log:

    def say(msg: str) -> None:
        log.write(msg + "\n")
        log.flush()
        print(msg, file=sys.stderr)

    frames = []
    for path in snakemake.input.beds:
        track = re.sub(r"_valleys\.bed$", "", os.path.basename(path))
        fields = track.split("_")
        df = pd.read_csv(
            path, sep="\t", header=None, names=["chrom", "start", "end"]
        )
        df.insert(0, "track", track)
        df.insert(1, "mark", fields[0] if len(fields) > 0 else "")
        df.insert(2, "clone", fields[1] if len(fields) > 1 else "")
        df.insert(3, "condition", fields[2] if len(fields) > 2 else "")
        df.insert(4, "allele", fields[3] if len(fields) > 3 else "")
        df["size"] = df["end"] - df["start"]
        df["size_class"] = [size_class(int(w)) for w in df["size"]]
        say(
            f"{track}: {len(df)} valleys, median {int(df['size'].median()) if len(df) else 0} bp, "
            f"max {int(df['size'].max()) if len(df) else 0} bp"
        )
        frames.append(df)

    out = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    out.to_csv(snakemake.output.tsv, sep="\t", index=False)
    say(f"wrote {len(out)} rows to {snakemake.output.tsv}")
