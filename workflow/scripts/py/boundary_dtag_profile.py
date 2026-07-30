"""Signal at NodTAG escaping boundaries, compared NodTAG vs dTAG -- Fig 4f/4g, EFig 6b.

From ``xci_valleys_check_LAST.ipynb`` cell 19. NO background subtraction -- this
is the ``boundary_profile`` lineage, not the ``boundary_stackup`` one
(correction R-2).

The design question this rule answers is "what happens to the signal at a FIXED
set of positions when the degron is depleted", so:

* boundaries come from the NodTAG condition only, and are never re-called on the
  dTAG data;
* CTCF status is likewise the NodTAG peak set;
* the same two profiles are computed twice, once from the NodTAG bigWig and once
  from the dTAG bigWig of the same clone and signal.

Boundary subsets (both are plotted; EFig 6b is the pie over the second):

    escaping        every boundary of an escaping-gene valley
    escaping_ctcf   those with a NodTAG CTCF peak within the 50 kb inward /
                    10 kb outward window

Acceptance (EFig 6b): E6A7 n = 110 escaping boundaries, 80 CTCF+ / 30 CTCF-;
F3 n = 36, 23 / 13. Fig 4f/4g quote the CTCF+ subsets, E6A7 n = 80 and F3 n = 23.

Geometry is identical to ``boundary_signal_matrix.py``: ``bbi.stackup`` over
``pos +/- flank`` in ``nbins`` bins, rows with ``side == "R"`` reversed,
``np.nanmean`` down the rows, boundaries within ``flank`` of a chromosome end
dropped (the original did the same here).

Emits an ``.npz`` with, for each subset in {escaping, escaping_ctcf} and each
condition in {nodtag, dtag}::

    profile_{subset}_{condition}   (nbins,)
    matrix_{subset}_{condition}    (n_used, nbins)
    n_{subset}                     boundaries used
plus ``positions_{subset}``, ``sides_{subset}``, ``flank``, ``nbins``.
"""

import os
import sys

import bbi
import numpy as np
import pandas as pd

log_path = snakemake.log[0]
os.makedirs(os.path.dirname(log_path), exist_ok=True)
os.makedirs(os.path.dirname(snakemake.output.npz), exist_ok=True)

flank = int(snakemake.params.flank)
nbins = int(snakemake.params.nbins)
clone = snakemake.wildcards.clone
signal = snakemake.wildcards.signal
ESCAPING = "escaping-gene-valley"


def profile_at(bw_path, chroms, positions, sides):
    """Mean signal profile over `positions`, R-side rows reversed. No background."""
    if len(positions) == 0:
        return np.full(nbins, np.nan), np.empty((0, nbins), dtype=float)
    matrix = np.asarray(
        bbi.stackup(bw_path, chroms, positions - flank, positions + flank, bins=nbins),
        dtype=float,
    )
    right = sides == "R"
    matrix[right] = matrix[right][:, ::-1]
    return np.nanmean(matrix, axis=0), matrix


with open(log_path, "w") as log:

    def say(msg: str) -> None:
        log.write(msg + "\n")
        log.flush()
        print(msg, file=sys.stderr)

    ctcf = pd.read_csv(snakemake.input.ctcf, sep="\t")
    escaping = ctcf[ctcf["valley_class"] == ESCAPING].reset_index(drop=True)
    say(
        f"{clone} {signal}: {len(ctcf)} NodTAG boundaries, "
        f"{len(escaping)} at escaping valleys"
    )

    chrom_sizes = pd.read_csv(
        snakemake.input.chrom_sizes, sep="\t", header=None, names=["chrom", "size"]
    ).set_index("chrom")["size"]

    subsets = {
        "escaping": escaping,
        "escaping_ctcf": escaping[escaping["has_ctcf"]].reset_index(drop=True),
    }

    payload = {"flank": flank, "nbins": nbins}
    for name, subset in subsets.items():
        positions = subset["start"].to_numpy(dtype=np.int64)
        sides = subset["side"].to_numpy(dtype=object)
        chroms = subset["chrom"].to_numpy(dtype=object)

        sizes = np.array(
            [int(chrom_sizes.get(c, 0)) for c in chroms], dtype=np.int64
        )
        keep = (positions >= flank) & (positions <= sizes - flank)
        dropped = int((~keep).sum())
        positions, sides, chroms = positions[keep], sides[keep], chroms[keep]
        if dropped:
            say(f"  {name}: dropped {dropped} boundaries within {flank} bp of a chrom end")

        for condition, bw_path in (
            ("nodtag", snakemake.input.nodtag_bw),
            ("dtag", snakemake.input.dtag_bw),
        ):
            profile, matrix = profile_at(bw_path, chroms, positions, sides)
            payload[f"profile_{name}_{condition}"] = profile
            payload[f"matrix_{name}_{condition}"] = matrix

        payload[f"n_{name}"] = len(positions)
        payload[f"positions_{name}"] = positions
        payload[f"sides_{name}"] = np.asarray(sides, dtype=str)
        say(f"  {name}: n = {len(positions)}")

    np.savez_compressed(snakemake.output.npz, **payload)
    say(f"wrote {snakemake.output.npz}")
