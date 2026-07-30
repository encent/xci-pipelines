"""Split valley boundaries by proximity to a CTCF motif.

From ``Stackups_*_with_motif.ipynb`` / ``_without_motif.ipynb``. Verbatim rule::

    boundaries_with_motifs = bioframe.closest(boundaries, motifs)
    motif_yes  <=>  distance <= 1000
    motif_no   <=>  distance >  1000

``bioframe.closest`` reports distance 0 for an overlap, so a boundary sitting
inside a motif lands in ``motif_yes``. The two outputs partition the input:
``len(motif_yes) + len(motif_no) == len(all)``, asserted below -- if that ever
fails, ``closest`` dropped rows (it does that when a chromosome has no motif at
all) and the split is silently lossy.

Motifs come from ``resources/genome/CTCF_mm10_X_only.bed``, the FIMO JASPAR
MA0139.1 scan (44,859 chrX motifs). The original read it with ``schema='bed3'``
and emitted a pandas ParserWarning about the extra columns; we read three
columns explicitly instead, which is the same data without the warning.

These files feed ``boundary_stackup`` only. They are ground truth, but no paper
panel depends on them (the archaeology notes 5B: the motif split was an exploration that did
not reach the paper).
"""

import os
import sys

import bioframe
import pandas as pd

log_path = snakemake.log[0]
os.makedirs(os.path.dirname(log_path), exist_ok=True)
for out in (snakemake.output.yes, snakemake.output.no):
    os.makedirs(os.path.dirname(out), exist_ok=True)

threshold = int(snakemake.params.distance)
COLUMNS = ["chrom", "start", "end", "side"]

with open(log_path, "w") as log:

    def say(msg: str) -> None:
        log.write(msg + "\n")
        log.flush()
        print(msg, file=sys.stderr)

    boundaries = pd.read_csv(
        snakemake.input.bed, sep="\t", header=None, names=COLUMNS
    )

    motifs = pd.read_csv(
        snakemake.input.motifs,
        sep="\t",
        header=None,
        usecols=[0, 1, 2],
        names=["chrom", "start", "end"],
        comment="#",
    )
    say(f"{len(boundaries)} boundaries, {len(motifs)} CTCF motifs")

    if len(boundaries) == 0:
        near = far = boundaries
        say("no boundaries; writing two empty files")
    else:
        closest = bioframe.closest(boundaries, motifs)
        near = closest[closest["distance"] <= threshold]
        far = closest[closest["distance"] > threshold]

        if len(near) + len(far) != len(boundaries):
            raise SystemExit(
                f"motif split lost rows: {len(near)} + {len(far)} != "
                f"{len(boundaries)}. bioframe.closest drops boundaries on "
                "chromosomes with no motif -- check that the motif BED covers "
                "every chromosome present in the boundary file."
            )

    say(
        f"motif_yes (<= {threshold} bp): {len(near)}   "
        f"motif_no (> {threshold} bp): {len(far)}"
    )

    near[COLUMNS].to_csv(snakemake.output.yes, sep="\t", index=False, header=False)
    far[COLUMNS].to_csv(snakemake.output.no, sep="\t", index=False, header=False)
