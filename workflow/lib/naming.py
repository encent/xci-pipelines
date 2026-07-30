"""Filename grammar for the XCI pipeline.

The original pipeline's filename grammar is load-bearing, not cosmetic:

  * METALoci joins a signal bedGraph to a dataset **by filename** -- a mismatch
    makes it silently drop the signal.
  * ``01_08_coolbox`` pairs a bigWig to its valley BED by stem.
  * ``Stackups_*`` swaps only the leading mark token to find a companion track.

So the grammar is frozen here, in one place, and every rule derives its paths
from these functions rather than formatting strings inline.

Grammar
-------
replicate sample   {mark}_{clone}_{condition}_{allele}_rep{n}
merged track       {mark}_{clone}_{condition}_{allele}
AcMe3 track        AcMe3_{clone}_{condition}_{allele}
consensus track    {mark}_{group}_{allele}
valley call        {merged track}_valleys.{bed,bw}
cooler             {locus}_{clone}_{condition}_{G1|G2}_{allele}[_rep{n}]
merged cooler      {locus}[_dTAG]_{allele}                (loops)
                   {locus}_{group}_{allele}               (compartments)
METALoci signal    {mark}_{dataset_adj}
"""

from __future__ import annotations

import re
from typing import Iterable

# --------------------------------------------------------------------------
# Controlled vocabularies. Anything outside these is a hard error at load time.
# --------------------------------------------------------------------------
MARKS = ("H3K27me3", "H3K27ac", "CTCF", "Rad21", "RNA-Seq")
DERIVED_MARKS = ("AcMe3",)
ALL_MARKS = MARKS + DERIVED_MARKS

WT_CLONES = ("E6", "JTG", "CL30", "B1", "C5")
DEGRON_CLONES = ("E6A7", "F3", "C5C10", "B1621")
CLONES = WT_CLONES + DEGRON_CLONES

CONDITIONS = ("WT", "CTCF-dTAG", "CTCF-NodTAG", "Rad21-dTAG", "Rad21-NodTAG")
ALLELES = ("Gall", "Xa", "Xi")
SNPSPLIT_TAGS = ("G1", "G2")          # G1 = C57BL/6J, G2 = Cast/EiJ
GROUPS = ("NodTAG-or-WT", "CTCF-dTAG", "Rad21-dTAG", "dTAG")
LOCI = ("Mecp2", "Jarid")             # Jarid == Kdm5c
ROI_TYPES = ("full", "escape", "non_escape")
SIDES = ("L", "R")                    # unified; the originals also used left/right

# Wildcard constraints, exported to Snakemake so an ambiguous match fails at
# DAG build rather than at hour 30.
WC = {
    "mark": "|".join(re.escape(m) for m in ALL_MARKS),
    "clone": "|".join(CLONES),
    "condition": "|".join(re.escape(c) for c in CONDITIONS),
    "allele": "|".join(ALLELES),
    "rep": r"[12]",
    "group": "|".join(re.escape(g) for g in GROUPS),
    "locus": "|".join(LOCI),
    "roi": "|".join(ROI_TYPES),
    "chrom": r"chr[0-9XY]+",
    "run": "wt|degron|consensus",
    "scope": "individual|merged_loops|merged_comps",
    "variant": "all|motif_yes|motif_no",
    "side": "L|R",
}


class NamingError(ValueError):
    """Raised when a name does not fit the frozen grammar."""


# --------------------------------------------------------------------------
# Constructors
# --------------------------------------------------------------------------
def sample_id(mark: str, clone: str, condition: str, allele: str, rep: int) -> str:
    """Replicate-level sample id, e.g. ``H3K27me3_E6_WT_Xi_rep1``."""
    _check(mark, ALL_MARKS, "mark")
    _check(clone, CLONES, "clone")
    _check(condition, CONDITIONS, "condition")
    _check(allele, ALLELES, "allele")
    if int(rep) not in (1, 2):
        raise NamingError(f"replicate must be 1 or 2, got {rep!r}")
    return f"{mark}_{clone}_{condition}_{allele}_rep{int(rep)}"


def track_id(mark: str, clone: str, condition: str, allele: str) -> str:
    """Replicate-merged track id, e.g. ``H3K27me3_E6_WT_Xi``."""
    _check(mark, ALL_MARKS, "mark")
    _check(clone, CLONES, "clone")
    _check(condition, CONDITIONS, "condition")
    _check(allele, ALLELES, "allele")
    return f"{mark}_{clone}_{condition}_{allele}"


def consensus_id(mark: str, group: str, allele: str) -> str:
    """Cross-clone consensus track, e.g. ``H3K27me3_NodTAG-or-WT_Xi``."""
    _check(mark, ALL_MARKS, "mark")
    _check(group, GROUPS, "group")
    _check(allele, ALLELES, "allele")
    return f"{mark}_{group}_{allele}"


def valley_id(track: str) -> str:
    """``H3K27me3_E6_WT_Xi`` -> ``H3K27me3_E6_WT_Xi_valleys``."""
    return f"{track}_valleys"


def swap_mark(track: str, new_mark: str) -> str:
    """Swap only the leading mark token, preserving clone/condition/allele.

    This reproduces the original's cross-mark lookup:
        f"{dir}/{signal}_{'_'.join(basename.split('_')[1:])}.bw"
    """
    _check(new_mark, ALL_MARKS, "mark")
    rest = track.split("_", 1)
    if len(rest) != 2:
        raise NamingError(f"cannot swap mark on {track!r}: no '_' separator")
    return f"{new_mark}_{rest[1]}"


def dataset_adj(cooler_name: str) -> str:
    """METALoci's join key between a cooler and its signal files.

    Verbatim from the original::

        dataset.replace('_rep1','').replace('_rep2','')
               .replace('_G1','').replace('_G2','')
               .replace('Jarid_','').replace('Mecp2_','')

    e.g. ``Jarid_E6_WT_G1_Xi`` -> ``E6_WT_Xi``.
    """
    out = cooler_name
    for token in ("_rep1", "_rep2", "_G1", "_G2"):
        out = out.replace(token, "")
    for token in ("Jarid_", "Mecp2_"):
        out = out.replace(token, "")
    return out


def metaloci_signal_name(mark: str, cooler_name: str) -> str:
    """``H3K27ac`` + ``Jarid_E6_WT_G1_Xi`` -> ``H3K27ac_E6_WT_Xi``.

    METALoci takes the signal name from the *file basename*, so this must be
    exactly the bedGraph stem staged for that dataset.
    """
    return f"{mark}_{dataset_adj(cooler_name)}"


# --------------------------------------------------------------------------
# Parsers
# --------------------------------------------------------------------------
_SAMPLE_RE = re.compile(
    rf"^(?P<mark>{WC['mark']})_(?P<clone>{WC['clone']})_"
    rf"(?P<condition>{WC['condition']})_(?P<allele>{WC['allele']})_rep(?P<rep>[12])$"
)
_TRACK_RE = re.compile(
    rf"^(?P<mark>{WC['mark']})_(?P<clone>{WC['clone']})_"
    rf"(?P<condition>{WC['condition']})_(?P<allele>{WC['allele']})$"
)
_COOLER_RE = re.compile(
    rf"^(?P<locus>{WC['locus']})_(?P<clone>{WC['clone']})_"
    rf"(?P<condition>{WC['condition']})_(?P<tag>G1|G2)_(?P<allele>Xa|Xi)"
    rf"(?:_rep(?P<rep>[12]))?$"
)


def parse_sample(sample: str) -> dict:
    m = _SAMPLE_RE.match(sample)
    if not m:
        raise NamingError(f"{sample!r} is not a valid replicate sample id")
    d = m.groupdict()
    d["rep"] = int(d["rep"])
    return d


def parse_track(track: str) -> dict:
    m = _TRACK_RE.match(track)
    if not m:
        raise NamingError(f"{track!r} is not a valid merged track id")
    return m.groupdict()


def parse_cooler(name: str) -> dict:
    m = _COOLER_RE.match(name)
    if not m:
        raise NamingError(f"{name!r} is not a valid cooler name")
    d = m.groupdict()
    d["rep"] = int(d["rep"]) if d["rep"] else None
    return d


def track_of(sample: str) -> str:
    """Drop the ``_rep{n}`` suffix: sample id -> merged track id."""
    d = parse_sample(sample)
    return track_id(d["mark"], d["clone"], d["condition"], d["allele"])


# --------------------------------------------------------------------------
def _check(value: str, allowed: Iterable[str], what: str) -> None:
    if value not in allowed:
        raise NamingError(
            f"unknown {what} {value!r}; allowed: {', '.join(sorted(allowed))}"
        )
