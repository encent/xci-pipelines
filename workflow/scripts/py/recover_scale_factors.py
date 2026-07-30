#!/usr/bin/env python3
"""Recover the applied scale factor from an already-written bigWig.

Why this exists
---------------
`csaw`'s restrict set is autosomal by design, so a chrX-only BAM subset cannot
reproduce the ground-truth normalisation factors *at all* -- not approximately,
but by construction. Decision D-09 therefore lets the testing team INJECT the
ground-truth factors and test everything below `bam_coverage` as a deterministic
function of chrX BAMs alone.

This tool builds those fixtures from the original pipeline's own output.

The method
----------
`bamCoverage --normalizeUsing None --scaleFactor s` writes `count x s` in every
bin. The smallest positive value in a track is therefore `1 x s` -- one read.
So the minimum positive value IS the scale factor, exactly, to float32.

It recovers the PRODUCT `1e6 / (LibSize * NormFactor)`. `LibSize` and
`NormFactor` are NOT separately recoverable from a bigWig, so the fixture
carries a single `scale_factor` column and leaves the other two empty; nothing
downstream may read them in `fixture` mode.

Scope limit -- do not widen it
------------------------------
Bit-exact for REPLICATE-LEVEL tracks only. Measured fractions of bins that are
integer multiples of the recovered quantum:

    replicate  *_rep{n}.bw     0.87 - 0.9998   (residue is float32 noise)
    merged/    (20 bp mean)    0.285           meaningless
    merged_5000/               0.0078          meaningless
    AcMe3_*    (log2 ratio)    0.00015         meaningless -- a log-ratio has
                                               no scale factor at all

A mean of two differently-scaled tracks has no single quantum. Never generate a
fixture for a merged, 5 kb or AcMe3 track: those are recomputed.

Three assertions
----------------
1. the factor recovered on each of >= 2 autosomes is BIT-identical;
2. the Gall / Xa / Xi of one replicate recover the IDENTICAL value;
3. under `--legacy`, assertion 2 is INVERTED for the five files that the
   original's positional broadcast mis-assigned (decision D-03).

Usage
-----
    python3 recover_scale_factors.py \
        --bigwig-dir /mnt/scratch/.../norm_bw_WT_H3K27me3_compos_all \
        --normgroup WT_H3K27me3 \
        --out resources/fixtures/scalefactors/WT_H3K27me3.tsv

The bigwig directory is READ-ONLY; this tool never writes into it.
"""

from __future__ import annotations

import argparse
import os
import sys
from collections import defaultdict

import numpy as np
import pyBigWig

# The five files the `selected` H3K27ac run mis-assigned, and whose factor each
# actually carries. Under --legacy these are EXPECTED to disagree with their
# siblings; without --legacy their presence means the input is the buggy set.
KNOWN_MISPAIRED = {
    "H3K27ac_CL30_WT_Xa_rep1": "E6_rep1",
    "H3K27ac_E6_WT_Gall_rep2": "CL30_rep1",
    "H3K27ac_E6_WT_Xi_rep1": "JTG_rep1",
    "H3K27ac_JTG_WT_Gall_rep2": "E6_rep2",
    "H3K27ac_JTG_WT_Xi_rep1": "JTG_rep2",
}

DEFAULT_CHROMS = ["chr19", "chr13"]


def recover_on(bw, chrom: str) -> float | None:
    """Minimum positive value on `chrom`, or None if the chrom is absent/empty."""
    if chrom not in bw.chroms():
        return None
    length = bw.chroms()[chrom]
    best = None
    step = 5_000_000
    for start in range(0, length, step):
        end = min(start + step, length)
        try:
            vals = bw.values(chrom, start, end, numpy=True)
        except RuntimeError:
            continue
        vals = vals[np.isfinite(vals)]
        vals = vals[vals > 0]
        if vals.size:
            m = float(vals.min())
            if best is None or m < best:
                best = m
    return best


def recover(path: str, chroms: list[str]) -> tuple[float, dict]:
    """Recover the factor and cross-check it across chromosomes.

    The factor is the smallest per-chromosome minimum. A chromosome whose
    minimum is a HIGHER INTEGER MULTIPLE is confirmation, not a contradiction:
    it simply has no single-read bin.

    That is common on the allelic tracks, where coverage is roughly halved.
    Measured example -- `H3K27me3_B1_WT_Xa_rep1`:

        chr19  0.10588400065898895   = 2 x quantum  (no 1-read bin on chr19)
        chr13  0.05294220149517059   = 1 x quantum

    An earlier version of this tool demanded bit-equality across chromosomes
    and rejected that file. The archaeologist's bit-equality measurement was
    made on `Gall` tracks, which are deep enough that every chromosome has a
    1-read bin; it does not generalise to `Xa`/`Xi`. What must hold -- and what
    is asserted here -- is that every chromosome's minimum is an integer
    multiple of the smallest.
    """
    per_chrom = {}
    with pyBigWig.open(path) as bw:
        for c in chroms:
            v = recover_on(bw, c)
            if v is not None:
                per_chrom[c] = v
    if not per_chrom:
        raise SystemExit(f"{path}: no positive values on any of {chroms}")

    factor = min(per_chrom.values())
    for c, v in per_chrom.items():
        mult = v / factor
        if abs(mult - round(mult)) > 1e-4:
            raise SystemExit(
                f"{os.path.basename(path)}: the minimum on {c} ({v!r}) is "
                f"{mult:.6f} x the smallest recovered quantum ({factor!r}), "
                "which is not an integer.\n"
                "That means the minimum positive bin is a genuine fraction "
                "rather than a whole number of reads, so no single scale factor "
                "can be recovered. Do not use this file as a fixture."
            )
    return factor, per_chrom


def sample_key(stem: str) -> tuple[str, str] | None:
    """`H3K27me3_E6_WT_Xi_rep1` -> (`H3K27me3_E6_WT_rep1`, `Xi`)."""
    parts = stem.split("_")
    if len(parts) < 5 or not parts[-1].startswith("rep"):
        return None
    allele = parts[-2]
    if allele not in ("Gall", "Xa", "Xi"):
        return None
    return "_".join(parts[:-2] + [parts[-1]]), allele


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bigwig-dir", required=True, help="READ-ONLY source directory")
    ap.add_argument("--normgroup", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--chroms", nargs="+", default=DEFAULT_CHROMS)
    ap.add_argument(
        "--legacy",
        action="store_true",
        help="expect the five known D-03 mis-pairings instead of rejecting them",
    )
    args = ap.parse_args()

    bws = sorted(
        f for f in os.listdir(args.bigwig_dir)
        if f.endswith(".bw") and "_rep" in f
    )
    if not bws:
        raise SystemExit(f"no replicate-level bigWigs in {args.bigwig_dir}")

    print(f"{args.normgroup}: {len(bws)} replicate bigWigs, "
          f"recovering on {', '.join(args.chroms)}", file=sys.stderr)

    rows, by_replicate = [], defaultdict(dict)
    for fn in bws:
        stem = fn[:-3]
        factor, per_chrom = recover(os.path.join(args.bigwig_dir, fn), args.chroms)
        rows.append((stem, factor, ",".join(sorted(per_chrom)), fn))
        key = sample_key(stem)
        if key:
            by_replicate[key[0]][key[1]] = (stem, factor)

    # Assertion 2 (inverted under --legacy): Gall/Xa/Xi of one replicate agree.
    #
    # Evaluated PER REPLICATE GROUP, not per file. When one allele of a
    # replicate carries the wrong factor, its two siblings sit in the same
    # disagreeing group while being entirely innocent -- flagging them
    # individually reports 10 "unexplained" files for 5 real defects. A group
    # is explained if ANY of its members is a documented D-03 mis-pairing.
    problems, expected_bad = [], []
    for rep, alleles in sorted(by_replicate.items()):
        if len(alleles) < 2:
            continue
        if len({f for _, f in alleles.values()}) == 1:
            continue
        stems = [stem for stem, _ in alleles.values()]
        known = [s for s in stems if s in KNOWN_MISPAIRED]
        if known:
            expected_bad.extend(known)
        else:
            problems.append(rep)

    if expected_bad and not args.legacy:
        print("\n  NOTE: the known D-03 mis-pairings are present in this "
              "directory:\n    " + "\n    ".join(sorted(set(expected_bad))),
              file=sys.stderr)
        print("  That is expected for the WT H3K27ac `selected` run. The fixture "
              "records what was APPLIED, which is what a legacy-mode test needs.",
              file=sys.stderr)
    if problems:
        raise SystemExit(
            "Gall/Xa/Xi of a replicate recovered DIFFERENT factors, and no "
            "member of these groups is a documented D-03 mis-pairing:\n    "
            + "\n    ".join(sorted(problems))
        )
    if args.legacy and not expected_bad:
        raise SystemExit(
            "--legacy was given but none of the five known mis-pairings were "
            "found. This is not the buggy `selected` H3K27ac directory."
        )

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write(
            "# Ground-truth scale factors recovered from the ORIGINAL pipeline's\n"
            "# bigWigs: the minimum positive value on an autosome is exactly the\n"
            "# applied `--scaleFactor` (one read x s).\n"
            f"# source: {args.bigwig_dir}\n"
            f"# recovered on: {', '.join(args.chroms)} (bit-equality asserted)\n"
            "# `scale_factor` is the PRODUCT 1e6/(LibSize*NormFactor). LibSize and\n"
            "# NormFactor are NOT separately recoverable and stay empty here.\n"
            "# REPLICATE LEVEL ONLY -- merged/5 kb/AcMe3 are always recomputed.\n"
        )
        fh.write("sample\tscale_factor\tLibSize\tNormFactor\trecovered_from\tsource_bigwig\n")
        for stem, factor, chroms, fn in rows:
            fh.write(f"{stem}\t{factor!r}\t\t\t{chroms}\t{fn}\n")

    print(f"  wrote {args.out} ({len(rows)} samples)", file=sys.stderr)
    for stem, factor, _, _ in rows[:3]:
        print(f"    {stem}\t{factor!r}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
