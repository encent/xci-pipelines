"""Replace hand-identified artefact pixels in a Capture Hi-C cooler.

Reimplements ``00_01_correct_pixels_in_coolers.py``. The original hard-coded
the pixel list inline and wrote in place over the ground-truth directory; here
the list is the fixture ``resources/fixtures/pixel_fixes.tsv`` and the output
is a new cooler in the pipeline's own work tree.

The replacement arithmetic is reproduced exactly, including its two distinct
behaviours:

``mode = mean``
    ``round(sum(v) / len(v))`` over whatever subset of the eight +/-5 kb
    neighbours actually exists in the pixel table. A sparse cooler can be
    missing several, and the original's ``try/except`` silently dropped them.

``mode in {left, right, top, bottom}``
    ``round((v_a + v_b + ...) / 5)`` over a FIXED five-neighbour subset. The
    original indexes those five slots by name, so a missing neighbour would
    raise; both directional dots used here are interior pixels with all eight
    neighbours present, and this script asserts that rather than quietly
    averaging fewer.

Neighbour numbering, verbatim from the original ``replace()``::

    v1 (-5k,-5k)   v2 (-5k, 0)   v3 (-5k,+5k)
    v4 ( 0 ,-5k)                 v5 ( 0 ,+5k)
    v6 (+5k,-5k)   v7 (+5k, 0)   v8 (+5k,+5k)
"""
import fnmatch
import os
import sys

import cooler
import pandas as pd

# Offsets of v1..v8 as (d_start1, d_start2), in bins of `resolution`.
NEIGHBOURS = {
    "v1": (-1, -1), "v2": (-1, 0), "v3": (-1, +1),
    "v4": (0, -1), "v5": (0, +1),
    "v6": (+1, -1), "v7": (+1, 0), "v8": (+1, +1),
}
DIRECTIONAL = {
    "left": ("v1", "v2", "v4", "v6", "v7"),
    "right": ("v2", "v3", "v5", "v7", "v8"),
    "top": ("v1", "v2", "v3", "v4", "v5"),
    "bottom": ("v4", "v5", "v6", "v7", "v8"),
}


def log(msg, fh):
    fh.write(msg + "\n")
    fh.flush()


def applicable(fixes: pd.DataFrame, name: str) -> pd.DataFrame:
    """Rows of the fixture whose pattern and snpsplit tag match this cooler."""
    keep = []
    for _, r in fixes.iterrows():
        if not fnmatch.fnmatch(name, r["cooler_pattern"]):
            continue
        tag = str(r["tag"])
        # The original tested `"G1" in basename(cf)`, a substring test.
        if tag != "any" and tag not in name:
            continue
        keep.append(r)
    return pd.DataFrame(keep, columns=fixes.columns)


def replacement(lookup: dict, s1: int, s2: int, res: int, mode: str):
    """The value that replaces pixel (s1, s2)."""
    values = {}
    for key, (d1, d2) in NEIGHBOURS.items():
        v = lookup.get((s1 + d1 * res, s2 + d2 * res))
        if v is not None:
            values[key] = v
    if mode == "mean":
        if not values:
            return None
        # Iteration order v1..v8 reproduces the original's append order; the
        # sum is order-independent for integers, so this is exact either way.
        return round(sum(values.values()) / len(values))
    wanted = DIRECTIONAL[mode]
    missing = [k for k in wanted if k not in values]
    if missing:
        raise ValueError(
            f"pixel ({s1}, {s2}) mode={mode} needs neighbours {list(wanted)} "
            f"but {missing} are absent from the pixel table. The original "
            f"00_01 would have raised here too; the fixture row is wrong for "
            f"this cooler, or the cooler is not the one it was written for."
        )
    return round(sum(values[k] for k in wanted) / 5)


def main():
    smk = snakemake  # noqa: F821 -- injected by Snakemake
    name = smk.wildcards.name
    res = int(smk.params.resolution)
    os.makedirs(os.path.dirname(smk.output[0]), exist_ok=True)
    os.makedirs(os.path.dirname(smk.log[0]), exist_ok=True)

    with open(smk.log[0], "w") as fh:
        clr = cooler.Cooler(smk.input.cool)
        pixels = clr.pixels()[:]
        bins = clr.bins()[:][["chrom", "start", "end"]]

        if smk.params.mode != "paper":
            log(f"fixtures.mode={smk.params.mode}: pixel repair SKIPPED for {name}.", fh)
            log("The artefact list is specific to the two published capture regions;", fh)
            log("applying it to other data would corrupt real contacts.", fh)
            cooler.create_cooler(smk.output[0], bins, pixels)
            return 0

        fixes = applicable(pd.read_csv(smk.input.fixes, sep="\t", comment="#"), name)
        log(f"{name}: {len(fixes)} pixel fix(es) from {smk.input.fixes}", fh)
        if fixes.empty:
            log("no fixture row matches this cooler -- copying unchanged.", fh)
            cooler.create_cooler(smk.output[0], bins, pixels)
            return 0

        joined = clr.pixels(join=True)[:]
        lookup = {
            (int(s1), int(s2)): c
            for s1, s2, c in zip(joined["start1"], joined["start2"], joined["count"])
        }
        index = {
            (int(s1), int(s2)): i
            for i, (s1, s2) in enumerate(zip(joined["start1"], joined["start2"]))
        }

        applied = 0
        for _, r in fixes.iterrows():
            s1, s2, mode = int(r["bin1_start"]), int(r["bin2_start"]), str(r["mode"])
            if (s1, s2) not in index:
                # The original used `.index[0]` and would have raised IndexError.
                log(f"  MISS  {s1}-{s1 + res} x {s2}-{s2 + res}: no such pixel; skipped", fh)
                continue
            new = replacement(lookup, s1, s2, res, mode)
            if new is None:
                log(f"  MISS  {s1} x {s2}: no neighbours at all; skipped", fh)
                continue
            i = index[(s1, s2)]
            old = pixels.loc[i, "count"]
            pixels.loc[i, "count"] = new
            applied += 1
            log(f"  FIX   {s1}-{s1 + res} x {s2}-{s2 + res}  mode={mode}  "
                f"{old} -> {new}", fh)

        log(f"{name}: applied {applied}/{len(fixes)} fixes", fh)
        cooler.create_cooler(smk.output[0], bins, pixels)
    return 0


if __name__ == "__main__":
    sys.exit(main())
