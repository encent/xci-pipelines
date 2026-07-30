"""Match Xa valley boundaries to Xi valley boundaries within one clone.

Secondary analysis. It reproduces
``TEST_03_01_..._final_degrons_overlap.ipynb``, which is the only surviving
producer of ``valleys_{overlap,xa_unique,xi_unique}_annotation/``. **No paper
panel depends on it.**

The matcher is a greedy two-pointer over sorted positions, per chromosome and
per side, with ``|xa - xi| <= threshold`` (15 kb). Greedy means the result is
order-dependent: a boundary is consumed by the first partner within threshold,
not by the nearest one. That is what the original did and it is preserved.

A matched pair collapses to a single 1-bp interval at the MEAN of the two
coordinates, so these rows are NOT 5 kb-quantised even though the valleys are
(e.g. ``chrX 6777500 6777501``). Unmatched boundaries keep their own
coordinate.

Vocabulary: the original wrote ``left``/``right`` here while
``boundaries/all/`` used ``L``/``R``. The pipeline unifies on ``L``/``R``
(deviation D-4) -- normalise before diffing the old files.

Emits one row per boundary::

    chrom start end side category
    category in {overlap, xa_unique, xi_unique}
"""

import os
import sys

import bioframe
import numpy as np
import pandas as pd

log_path = snakemake.log[0]
os.makedirs(os.path.dirname(log_path), exist_ok=True)
os.makedirs(os.path.dirname(snakemake.output.tsv), exist_ok=True)

THRESHOLD = 15_000
SIDES = (("L", "start"), ("R", "end"))


def match_greedy(xa_pos, xi_pos, threshold):
    """Greedy two-pointer matching of sorted positions -> list of (i, j)."""
    pairs = []
    i = j = 0
    while i < len(xa_pos) and j < len(xi_pos):
        a, b = xa_pos[i], xi_pos[j]
        if abs(a - b) <= threshold:
            pairs.append((i, j))
            i += 1
            j += 1
        elif b < a - threshold:
            j += 1
        else:
            i += 1
    return pairs


with open(log_path, "w") as log:

    def say(msg: str) -> None:
        log.write(msg + "\n")
        log.flush()
        print(msg, file=sys.stderr)

    xa = bioframe.read_table(snakemake.input.xa, schema="bed3")
    xi = bioframe.read_table(snakemake.input.xi, schema="bed3")
    say(f"{len(xa)} Xa valleys, {len(xi)} Xi valleys, threshold {THRESHOLD} bp")

    records = []
    for chrom in sorted(set(xa["chrom"]) | set(xi["chrom"])):
        xa_c = xa[xa["chrom"] == chrom].reset_index(drop=True)
        xi_c = xi[xi["chrom"] == chrom].reset_index(drop=True)

        for side, column in SIDES:
            xa_pos = xa_c[column].astype(int).to_numpy()
            xi_pos = xi_c[column].astype(int).to_numpy()
            if xa_pos.size == 0 and xi_pos.size == 0:
                continue

            xa_order = np.argsort(xa_pos)
            xi_order = np.argsort(xi_pos)
            pairs = match_greedy(
                xa_pos[xa_order].tolist(), xi_pos[xi_order].tolist(), THRESHOLD
            )

            matched_xa, matched_xi = set(), set()
            for ia, ib in pairs:
                a_idx, b_idx = xa_order[ia], xi_order[ib]
                matched_xa.add(int(a_idx))
                matched_xi.add(int(b_idx))
                mean_pos = int((int(xa_pos[a_idx]) + int(xi_pos[b_idx])) / 2)
                mean_pos = max(0, mean_pos)
                records.append((chrom, mean_pos, mean_pos + 1, side, "overlap"))

            for k, pos in enumerate(xa_pos):
                if k not in matched_xa:
                    records.append((chrom, int(pos), int(pos) + 1, side, "xa_unique"))
            for k, pos in enumerate(xi_pos):
                if k not in matched_xi:
                    records.append((chrom, int(pos), int(pos) + 1, side, "xi_unique"))

    out = pd.DataFrame(
        records, columns=["chrom", "start", "end", "side", "category"]
    ).sort_values(["chrom", "start", "side"], kind="stable")

    for category, n in out["category"].value_counts().items():
        say(f"  {category:10s} {n}")

    out.to_csv(snakemake.output.tsv, sep="\t", index=False)
    say(f"wrote {snakemake.output.tsv}")
