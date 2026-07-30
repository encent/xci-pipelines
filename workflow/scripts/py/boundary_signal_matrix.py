"""Signal matrix around valley boundaries. Serves TWO rules with DIFFERENT maths.

    boundary_stackup   subtract_background = True    Stackups_*.ipynb lineage
    boundary_profile   subtract_background = False   xci_valleys_check_LAST.ipynb

Same inputs, same geometry, different arithmetic. Routing a panel to the wrong
one produces a plausible-looking curve that silently does not reproduce
(correction R-2). The rule that feeds EFig 2e/2f/2g, Fig 4f/4g and EFig 6b is
the one WITHOUT background subtraction. ``boundary_stackup`` feeds no paper
panel at all -- it exists because it is ground truth and because it produces the
non-paper composite extras.

Geometry (both paths)
---------------------
``bbi.stackup(bw, chrom, pos - flank, pos + flank, bins=nbins)``, one row per
boundary, then rows whose ``side == "R"`` are reversed so that "into the valley"
always points the same way. That is ``FLIPPED = True`` and
``legacy.stackup_flip_rule: side``. The original's plain stackup variant flipped
on ``i % 2 == 1`` instead, which is equivalent only while the boundary file is
strictly L,R,L,R -- it breaks on the motif subsets, so ``side`` is the default.

Background (stackup path only)
------------------------------
100 circular shifts of the whole boundary set,
``numpy.random.default_rng(seed = background_seed_base + clone_index)``,
``shift = rng.integers(0, chromsize)``, wrapped by subtracting ``chromsize``
where a coordinate runs past the end. A shift is REDRAWN if it would make any
valley straddle the origin -- roughly a quarter of draws, so the rejection loop
is load-bearing for the seed stream, not a formality. That is why this script
reads the valley BED as well as the boundary BED: the acceptance test is defined
on valleys, and the original applied it even when stacking a motif subset.
The 100 shifted matrices are averaged and subtracted.

Deviations, both documented in the module notes
----------------------------------------------------------
* Boundaries closer than ``flank`` to either chromosome end are dropped. The
  profile path always did this; the stackup path passed negative coordinates
  straight to ``bbi``. At most two boundaries per set are affected (the
  assembly-gap valley's edges) and the behaviour of ``bbi`` on a negative start
  is not contractual.
* ``clone_index`` is reconstructed as the track's position in the sorted Xi
  track list of its generation (WT or degron); the original took it from a glob
  enumeration order. It agrees for WT. See ``30_valleys.smk``.
* Only the boundary-centred matrix is emitted. The original also built a
  three-segment left-flank / valley-body / right-flank stackup at
  ``nbins_valley_body = 50`` for H3K27me3; that is a different anchor (valleys,
  not boundaries) and is not produced by this rule.

Emits an ``.npz``::

    matrix        (n_used, nbins)  flipped, raw signal
    background    (n_used, nbins)  mean of the 100 shifts (stackup path only)
    matrix_final  matrix, minus background on the stackup path
    profile       aggregate of matrix_final over rows (mean or median)
    profile_raw   aggregate of matrix over rows, never background-subtracted
    chroms, positions, sides, n_used, n_dropped, flank, nbins, seed,
    subtract_background, aggregate
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
n_shifts = int(snakemake.params.shifts)
seed = int(snakemake.params.seed_base) + int(snakemake.params.clone_index)
aggregate = str(snakemake.params.aggregate)
subtract_background = bool(snakemake.params.subtract_background)

BOUNDARY_COLUMNS = ["chrom", "start", "end", "side"]


def shifted_valleys(valleys: pd.DataFrame, chromsize: int, shift: int):
    """The original's acceptance test: None if any valley straddles the origin."""
    moved = valleys.copy()
    moved["start"] = moved["start"] + shift
    moved["end"] = moved["end"] + shift
    wrapped = (moved["start"] >= chromsize) | (moved["end"] >= chromsize)
    moved.loc[wrapped, "start"] -= chromsize
    moved.loc[wrapped, "end"] -= chromsize
    if ((moved["start"] < 0) | (moved["end"] < 0)).any():
        return None
    return moved.reset_index(drop=True)


def shifted_positions(positions: np.ndarray, chromsize: int, shift: int) -> np.ndarray:
    moved = positions + shift
    moved[moved >= chromsize] -= chromsize
    return moved


def flip_right(matrix: np.ndarray, sides: np.ndarray) -> np.ndarray:
    out = np.array(matrix, dtype=float, copy=True)
    right = sides == "R"
    out[right] = out[right][:, ::-1]
    return out


with open(log_path, "w") as log:

    def say(msg: str) -> None:
        log.write(msg + "\n")
        log.flush()
        print(msg, file=sys.stderr)

    boundaries = pd.read_csv(
        snakemake.input.boundaries, sep="\t", header=None, names=BOUNDARY_COLUMNS
    )
    chrom_sizes = pd.read_csv(
        snakemake.input.chrom_sizes, sep="\t", header=None, names=["chrom", "size"]
    ).set_index("chrom")["size"]

    say(
        f"{len(boundaries)} boundaries, flank={flank}, nbins={nbins}, "
        f"aggregate={aggregate}, subtract_background={subtract_background}, "
        f"seed={seed}"
    )

    chrom = boundaries["chrom"].iloc[0] if len(boundaries) else "chrX"
    chromsize = int(chrom_sizes[chrom]) if chrom in chrom_sizes.index else 0

    positions = boundaries["start"].to_numpy(dtype=np.int64)
    sides = boundaries["side"].to_numpy(dtype=object)
    chroms = boundaries["chrom"].to_numpy(dtype=object)

    keep = (positions >= flank) & (positions <= chromsize - flank)
    n_dropped = int((~keep).sum())
    positions, sides, chroms = positions[keep], sides[keep], chroms[keep]
    if n_dropped:
        say(f"dropped {n_dropped} boundaries closer than {flank} bp to a chromosome end")

    if len(positions) == 0:
        say("WARNING: no usable boundaries; writing an empty matrix")
        matrix = np.empty((0, nbins), dtype=float)
        background = np.empty((0, nbins), dtype=float)
        matrix_final = matrix
        profile = np.full(nbins, np.nan)
        profile_raw = np.full(nbins, np.nan)
    else:
        matrix = flip_right(
            bbi.stackup(
                snakemake.input.bw, chroms, positions - flank, positions + flank,
                bins=nbins,
            ),
            sides,
        )
        say(f"stacked {matrix.shape[0]} x {matrix.shape[1]}")

        background = np.empty((0, nbins), dtype=float)
        if subtract_background and n_shifts > 0:
            valleys = pd.read_csv(
                snakemake.input.valleys,
                sep="\t",
                header=None,
                names=["chrom", "start", "end"],
            )
            rng = np.random.default_rng(seed=seed)
            shifts, redraws = [], 0
            while len(shifts) < n_shifts:
                candidate = int(rng.integers(0, chromsize))
                if len(valleys) and shifted_valleys(valleys, chromsize, candidate) is None:
                    redraws += 1
                    continue
                shifts.append(candidate)
            say(f"{n_shifts} accepted shifts, {redraws} redrawn (valley crossed origin)")

            total = np.zeros_like(matrix, dtype=float)
            for shift in shifts:
                moved = shifted_positions(positions, chromsize, shift)
                total += flip_right(
                    bbi.stackup(
                        snakemake.input.bw, chroms, moved - flank, moved + flank,
                        bins=nbins,
                    ),
                    sides,
                )
            background = total / len(shifts)

        matrix_final = matrix - background if background.size else matrix
        reducer = np.nanmedian if aggregate == "median" else np.nanmean
        profile = reducer(matrix_final, axis=0)
        profile_raw = reducer(matrix, axis=0)

    np.savez_compressed(
        snakemake.output.npz,
        matrix=matrix,
        background=background,
        matrix_final=matrix_final,
        profile=profile,
        profile_raw=profile_raw,
        chroms=np.asarray(chroms, dtype=str),
        positions=np.asarray(positions, dtype=np.int64),
        sides=np.asarray(sides, dtype=str),
        n_used=len(positions),
        n_dropped=n_dropped,
        flank=flank,
        nbins=nbins,
        seed=seed,
        subtract_background=subtract_background,
        aggregate=aggregate,
    )
    say(f"wrote {snakemake.output.npz}")
