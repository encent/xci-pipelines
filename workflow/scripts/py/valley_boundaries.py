"""Turn valley intervals into their two 1-bp end points.

From the ``Stackups_*_without_motif.ipynb`` lineage, which is the last writer of
the ground-truth ``boundaries/all/`` files.

Format -- BED4, and deliberately NOT the BED3 that ``call_valleys`` writes:

    chrom   start   start+1   side

``side`` is ``L`` for a valley start and ``R`` for a valley end. Row order is
L,R,L,R,... following the valley order, so ``len(boundaries) == 2 * len(valleys)``
exactly. The original's plain stackup variant relied on that ordering and
flipped rows by ``i % 2 == 1``; the motif variants flipped by
``side == "R"``. Those agree only as long as the ordering holds, which is why
``legacy.stackup_flip_rule`` defaults to ``side`` -- it is the one that stays
correct if the set is ever subset (as ``boundary_motif_split`` does).

Vocabulary note (deviation D-4): the original used ``L``/``R`` here but
``left``/``right`` in the non-paper Xa/Xi overlap outputs. The pipeline unifies
on ``L``/``R`` everywhere. Normalise before diffing those older files.

Ground-truth row counts, ``boundaries/all/`` (the ground-truth inventory 4.4):
    WT      B1 522, C5 608, CL30 648, E6 754, JTG 604
    degron  B1621-NodTAG 568, E6A7-NodTAG 642, F3-NodTAG 702

Reads the RAW valley set: these boundary files are ground truth and were
produced before the 3-gene filter existed. ``boundary_ctcf_status`` re-derives
its own boundaries from the FILTERED set for the paper panels.
"""

import os
import sys

import bioframe
import pandas as pd

log_path = snakemake.log[0]
os.makedirs(os.path.dirname(log_path), exist_ok=True)
os.makedirs(os.path.dirname(snakemake.output.bed), exist_ok=True)

sides = list(snakemake.params.side_vocabulary)
if sides != ["L", "R"]:
    raise SystemExit(
        f"valleys.boundaries.side_vocabulary is {sides!r}; the pipeline is "
        "written against ['L', 'R'] (deviation D-4). Change the panels and the "
        "test comparators together, or leave it alone."
    )
left_tag, right_tag = sides

with open(log_path, "w") as log:

    def say(msg: str) -> None:
        log.write(msg + "\n")
        log.flush()
        print(msg, file=sys.stderr)

    valleys = bioframe.read_table(snakemake.input.bed, schema="bed3")

    rows = []
    for _, v in valleys.iterrows():
        rows.append((v["chrom"], int(v["start"]), left_tag))
        rows.append((v["chrom"], int(v["end"]), right_tag))

    boundaries = pd.DataFrame(rows, columns=["chrom", "start", "side"])
    boundaries["end"] = boundaries["start"] + 1

    say(
        f"{len(valleys)} valleys -> {len(boundaries)} boundaries "
        f"({int((boundaries['side'] == left_tag).sum())} {left_tag}, "
        f"{int((boundaries['side'] == right_tag).sum())} {right_tag}), "
        f"flip rule = {snakemake.params.flip_rule}"
    )

    boundaries[["chrom", "start", "end", "side"]].to_csv(
        snakemake.output.bed, sep="\t", index=False, header=False
    )
