"""Transplant the `full`-region Kamada-Kawai layout into a sub-region `.mlo`.

This reproduces `subset_mlo_from_file2_to_file1` from
`01_09_metaloci_run_WT.ipynb` / `_dtag.ipynb` cells 13-14. It runs AFTER
`metaloci layout` and BEFORE `metaloci lm`. It is undocumented in the paper and
essential to reproduce: without it every `escape` / `non_escape` Gaudi plot and
every sub-region compartment-like strength is computed on a *different* 2D
embedding than the published one.

FOUR DETAILS. GET ANY OF THEM WRONG AND FIG 6 IS SILENTLY WRONG.
================================================================

1. `abs()` IS MANDATORY
-----------------------
    shift = abs(full.start - sub.start) // resolution
    size  = len(sub.kk_nodes)

the archaeology notes section J.5 originally dropped the `abs()` (corrected by amendment A-1 /
the design review R-1). Without it:

    locus   region       correct shift   un-abs()-ed
    Mecp2   non_escape        0               0      <- accidentally safe
    Mecp2   escape          +66             -66      <- WRONG
    Jarid   non_escape        0               0      <- accidentally safe
    Jarid   escape         +280            -280      <- WRONG

A negative index slices from the END of the array, so both `escape` regions
would be built from the wrong bins -- with no error, no shape mismatch, and a
perfectly plausible-looking plot. `non_escape` shares its start with `full`, so
it is shift 0 either way and was accidentally safe.

This script therefore (a) computes the shift with `abs()`, (b) asserts it
equals the value the rule derived from `config.loci`, and (c) refuses to run on
a negative shift at all. `expected_shift` is passed in as a params value so the
assertion is data-driven rather than a hard-coded 66/280.

2. THE ORIGINAL FUNCTION NAME IS INVERTED
-----------------------------------------
`subset_mlo_from_file2_to_file1(file1=<full>, file2=<sub>)` **reads file1 and
writes into file2**. Anyone following the name transplants backwards, which
overwrites the full-region layout with a 330 kb slice of itself. The names here
are `source_full` / `target_sub` precisely so that cannot happen.

3. IT RUNS TWICE PER DATASET
----------------------------
`full -> non_escape`, then `full -> escape` (cell 14, which is mislabelled
`# Calculate compartmentalization`). Both reads come from the same `full`
object. In this pipeline that is two independent jobs, one per `roi`.

4. ONLY FOUR FIELDS ARE REPLACED
--------------------------------
`kk_nodes`, `kk_restraints_matrix`, `kk_coords`, `kk_distances`, sliced
`[shift : shift+size]`. The target keeps its own `start`, `end`, `resolution`,
`region`, `poi`, `matrix`, `subset_matrix` and signal (`new_object =
mlobject_scape`). The original object is renamed `*_old.mlo`.

CONSEQUENCE, UNDOCUMENTED IN THE PAPER
--------------------------------------
The sub-region `metaloci layout -l 10` is COMPLETELY DISCARDED. The `-l 10`
layouts are computed, written, and then overwritten by a slice of the `-l 15`
full-region layout. We still run the sub-region layouts, because the target
object's own `start`/`end`/`poi`/`matrix` and the `_old.mlo` sidecar all come
from them, and because `metaloci lm` needs a valid object to patch.

`metaloci lm -f` is SAFE after this: it clears only `lmi_info` and
`lmi_geometry` (tools/ml.py), never the `kk_*` fields. `metaloci layout -f`
is NOT safe -- it removes and re-creates the object -- which is why the DAG
always orders layout before transplant.

FILE FORMAT NOTE
----------------
A `.mlo` is `pickle.dump(mlobject.__dict__)` (`metaloci/mlo.py: save()`), i.e.
a plain **dict**, not a class instance. That is why the original indexes it as
`mlobject['kk_nodes']`, and why this script never needs METALoci importable.
"""

import os
import pickle
import shutil
import sys

import numpy as np

full_path = snakemake.input.full
sub_path = snakemake.input.sub
flag_path = snakemake.output.flag
old_path = snakemake.output.old

expected_shift = int(snakemake.params.expected_shift)
resolution = int(snakemake.params.resolution)
roi = snakemake.params.roi
dataset = snakemake.params.dataset

os.makedirs(os.path.dirname(snakemake.log[0]), exist_ok=True)
log = open(snakemake.log[0], "w")


def say(msg):
    log.write(f"{msg}\n")


def fail(msg):
    say(f"ERROR {msg}")
    log.close()
    sys.exit(f"metaloci_transplant_kk[{dataset}/{roi}]: {msg}")


def load(path):
    with open(path, "rb") as fh:
        return pickle.load(fh)


KK_FIELDS = ("kk_nodes", "kk_restraints_matrix", "kk_coords", "kk_distances")


def slice_from_full(full, shift, size):
    """The four transplanted fields, exactly as the original slices them."""
    nodes = {k: v for k, v in full["kk_nodes"].items() if shift <= k < shift + size}
    nodes = {k - shift: v for k, v in nodes.items()}
    return {
        "kk_nodes": nodes,
        "kk_restraints_matrix": full["kk_restraints_matrix"][
            shift : shift + size, shift : shift + size
        ],
        "kk_coords": full["kk_coords"][shift : shift + size],
        "kk_distances": full["kk_distances"][shift : shift + size, shift : shift + size],
    }


full = load(full_path)

# --- pick the pristine target -------------------------------------------
# Normal path: the live `.mlo` is what `metaloci layout` just wrote. The one
# case where it is not is a manual `rm` of the flag file without re-running
# layout; then the live object is already a transplant product and re-slicing
# it would use the wrong `size` if `full` had dropped bins. Detect that and
# treat the existing `_old.mlo` as the pristine object instead.
target_path = sub_path
if os.path.exists(old_path):
    live = load(sub_path)
    probe_shift = abs(int(full["start"]) - int(live["start"])) // int(full["resolution"])
    probe_size = len(live["kk_nodes"])
    probe = slice_from_full(full, probe_shift, probe_size)
    already = np.array_equal(np.asarray(live["kk_coords"]), np.asarray(probe["kk_coords"]))
    if already:
        say(
            f"{sub_path} already carries the full-region layout and "
            f"{old_path} exists -- re-running from the pristine backup."
        )
        target_path = old_path

sub = load(target_path)

# --- the arithmetic ------------------------------------------------------
if int(full["resolution"]) != int(sub["resolution"]):
    fail(
        f"resolution mismatch: full={full['resolution']} sub={sub['resolution']}. "
        "The shift is expressed in bins; two resolutions make it meaningless."
    )
if int(full["resolution"]) != resolution:
    fail(f"object resolution {full['resolution']} != config {resolution}")

shift = abs(int(full["start"]) - int(sub["start"])) // int(full["resolution"])
size = len(sub["kk_nodes"])

say(f"dataset          {dataset}")
say(f"roi              {roi}")
say(f"full             {full_path}  start={full['start']} end={full['end']}")
say(f"sub              {target_path}  start={sub['start']} end={sub['end']}")
say(f"resolution       {full['resolution']}")
say(f"shift            {shift}   (expected {expected_shift})")
say(f"size             {size}   (len(full.kk_nodes)={len(full['kk_nodes'])})")

if shift < 0:
    fail(f"negative shift {shift} -- abs() was lost somewhere; refusing to slice")
if shift != expected_shift:
    fail(
        f"shift {shift} != expected {expected_shift} derived from config.loci. "
        "Either the region coordinates changed or the .mlo objects are stale."
    )
if size <= 0:
    fail(f"sub-region has {size} kk_nodes -- layout produced nothing to patch")
if shift + size > len(full["kk_coords"]):
    fail(
        f"slice [{shift}:{shift + size}] runs past the full layout "
        f"({len(full['kk_coords'])} coords). The sub-region is not contained "
        "in the full region."
    )

patch = slice_from_full(full, shift, size)

if len(patch["kk_nodes"]) != size:
    # Not fatal in the original, which silently accepted a short dict. Loud here:
    # it means the full layout dropped bins inside the sub-region window.
    say(
        f"WARNING transplanted kk_nodes has {len(patch['kk_nodes'])} entries, "
        f"expected {size} -- the full layout is missing bins in this window."
    )

# --- write ---------------------------------------------------------------
# `new_object = mlobject_scape`: the target keeps everything else it owns.
new_object = sub
for field in KK_FIELDS:
    new_object[field] = patch[field]

say(f"replaced fields  {', '.join(KK_FIELDS)}")
say(f"kk_restraints_matrix {np.shape(patch['kk_restraints_matrix'])}")
say(f"kk_distances         {np.shape(patch['kk_distances'])}")
say(f"kk_coords            {np.shape(patch['kk_coords'])}")

if target_path != old_path:
    if os.path.exists(old_path):
        os.remove(old_path)
    shutil.move(sub_path, old_path)
    say(f"backed up original layout to {old_path}")

with open(sub_path, "wb") as fh:
    pickle.dump(new_object, fh)
say(f"wrote {sub_path}")

# Verify the round trip before declaring success -- a truncated pickle here
# would only surface as a confusing `metaloci lm` crash much later.
check = load(sub_path)
for field in ("kk_coords", "kk_distances", "kk_restraints_matrix"):
    if np.shape(check[field]) != np.shape(patch[field]):
        fail(f"round-trip shape mismatch on {field}")
if len(check["kk_nodes"]) != len(patch["kk_nodes"]):
    fail("round-trip length mismatch on kk_nodes")
if int(check["start"]) != int(sub["start"]) or int(check["end"]) != int(sub["end"]):
    fail("the target lost its own start/end -- the transplant went the wrong way")

with open(flag_path, "w") as fh:
    fh.write(f"dataset\t{dataset}\n")
    fh.write(f"roi\t{roi}\n")
    fh.write(f"source_full\t{full_path}\n")
    fh.write(f"target_sub\t{sub_path}\n")
    fh.write(f"backup\t{old_path}\n")
    fh.write(f"shift\t{shift}\n")
    fh.write(f"size\t{size}\n")
    fh.write(f"resolution\t{resolution}\n")
    fh.write("fields\t" + ",".join(KK_FIELDS) + "\n")

say("done")
log.close()
