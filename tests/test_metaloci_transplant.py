#!/usr/bin/env python3
"""The highest-value test in the METALoci module: the KK transplant shifts.

Run directly (no pytest needed, no heavy data, ~1 s)::

    python3 tests/test_metaloci_transplant.py

WHAT IT GUARDS
--------------
`metaloci_transplant_kk` slices the FULL-region Kamada-Kawai layout into the
`escape` / `non_escape` `.mlo` objects at

    shift = abs(full.start - sub.start) // resolution

the archaeology notes section J.5 originally omitted the `abs()` (corrected by amendment A-1 /
the design review R-1). A negative shift does not raise: numpy happily slices from the
END of the array, so every `escape` Gaudi plot and every `escape`
compartment-like strength would be built from the wrong bins, silently. This
test:

1. recomputes the four shifts from `config/config.yaml` and asserts
   Mecp2 non_escape 0, Mecp2 escape +66, Jarid non_escape 0, Jarid escape +280;
2. asserts the un-`abs()`-ed formula really does go negative for both `escape`
   regions -- i.e. that the bug is reachable, not hypothetical;
3. demonstrates on synthetic arrays that a negative shift selects a *different,
   plausible-looking* slice, which is why it is silent;
4. runs the real `workflow/scripts/py/metaloci_transplant_kk.py` end to end on
   synthetic `.mlo` pickles and checks the four transplanted fields, the
   `_old.mlo` backup, the preserved `start`/`end`, and the provenance stamp;
5. checks that the script REFUSES a wrong or negative `expected_shift`.
"""

from __future__ import annotations

import os
import pickle
import sys
import tempfile
import types

import numpy as np
import yaml

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(REPO, "workflow", "scripts", "py", "metaloci_transplant_kk.py")
CONFIG = os.path.join(REPO, "config", "config.yaml")

EXPECTED = {
    ("Mecp2", "non_escape"): 0,
    ("Mecp2", "escape"): 66,
    ("Jarid", "non_escape"): 0,
    ("Jarid", "escape"): 280,
}

failures = []


def check(name, condition, detail=""):
    if condition:
        print(f"  PASS  {name}")
    else:
        print(f"  FAIL  {name}   {detail}")
        failures.append(name)


# ---------------------------------------------------------------------
# 1 + 2. the shifts, straight from config
# ---------------------------------------------------------------------
def test_shifts_from_config():
    print("[1] transplant shifts derived from config/config.yaml")
    cfg = yaml.safe_load(open(CONFIG))
    resolution = int(cfg["metaloci"]["resolution"])
    check("resolution is 5000", resolution == 5000, f"got {resolution}")

    for (locus, roi), expected in sorted(EXPECTED.items()):
        regions = cfg["loci"][locus]["regions"]
        full_start = int(regions["full"][0])
        sub_start = int(regions[roi][0])

        shift = abs(full_start - sub_start) // resolution
        check(f"{locus}/{roi} shift == {expected}", shift == expected, f"got {shift}")
        check(f"{locus}/{roi} shift is non-negative", shift >= 0, f"got {shift}")

        naive = (full_start - sub_start) // resolution
        if expected == 0:
            check(
                f"{locus}/{roi} is the accidentally-safe shift-0 case",
                naive == 0,
                f"got {naive}",
            )
        else:
            check(
                f"{locus}/{roi} WITHOUT abs() would be {-expected} (bug is reachable)",
                naive == -expected,
                f"got {naive}",
            )


# ---------------------------------------------------------------------
# 3. why a negative shift is silent
# ---------------------------------------------------------------------
def test_negative_shift_is_silent():
    """What the un-abs()-ed formula actually does, per locus.

    Both failure modes are worth naming, because they are different:

      Jarid escape   n_full 390, size 110, shift -280 -> coords[-280:-170]
                     == coords[110:220]. Right length, WRONG BINS, no error.
                     This is the fully silent corruption.
      Mecp2 escape   n_full 232, size 166, shift -66  -> coords[-66:100]
                     == coords[166:100] -> EMPTY. `kk_nodes` still comes back
                     with 100 plausible entries, so the object looks half
                     valid; `metaloci lm` then dies on an empty
                     kk_distances diagonal, far from the cause.
    """
    print("[3] what a negative shift really does, per locus")
    cfg = yaml.safe_load(open(CONFIG))
    reso = int(cfg["metaloci"]["resolution"])

    for locus, expect in (("Jarid", "wrong-bins"), ("Mecp2", "empty")):
        regions = cfg["loci"][locus]["regions"]
        n_full = (int(regions["full"][1]) - int(regions["full"][0])) // reso
        size = (int(regions["escape"][1]) - int(regions["escape"][0])) // reso
        shift = abs(int(regions["full"][0]) - int(regions["escape"][0])) // reso
        coords = np.arange(n_full, dtype=float).reshape(-1, 1)

        good = coords[shift : shift + size]
        bad = coords[-shift : -shift + size]

        check(
            f"{locus}: correct slice is bins {shift}..{shift + size - 1}",
            len(good) == size and int(good[0, 0]) == shift,
            f"len={len(good)} first={good[0, 0] if len(good) else None}",
        )
        if expect == "wrong-bins":
            check(
                f"{locus}: negative shift gives {size} bins from the WRONG offset",
                len(bad) == size and int(bad[0, 0]) != shift,
                f"len={len(bad)} first={bad[0, 0] if len(bad) else None}",
            )
        else:
            check(
                f"{locus}: negative shift gives an EMPTY slice",
                len(bad) == 0,
                f"len={len(bad)}",
            )
            # ... while kk_nodes still yields a plausible non-empty dict.
            nodes = {k - (-shift): k for k in range(n_full) if -shift <= k < -shift + size}
            check(
                f"{locus}: kk_nodes still looks plausible ({len(nodes)} entries)",
                0 < len(nodes) < size,
                f"got {len(nodes)}",
            )
        check(
            f"{locus}: the two slices genuinely differ",
            good.shape != bad.shape or not np.array_equal(good, bad),
            "they matched, which would make the bug untestable",
        )


# ---------------------------------------------------------------------
# 4 + 5. the real script, end to end
# ---------------------------------------------------------------------
def _make_mlo(start, end, n, resolution=5000, seed=0):
    rng = np.random.default_rng(seed)
    coords = rng.normal(size=(n, 2))
    dist = np.abs(np.subtract.outer(np.arange(n), np.arange(n))).astype(float)
    return {
        "region": f"chrX:{start}-{end}_0",
        "chrom": "chrX",
        "start": start,
        "end": end,
        "poi": 0,
        "resolution": resolution,
        "kk_nodes": {i: (float(i), float(i)) for i in range(n)},
        "kk_coords": coords,
        "kk_distances": dist,
        "kk_restraints_matrix": dist * 2.0,
        "lmi_info": {},
        "lmi_geometry": None,
    }


def _run_script(full_path, sub_path, old_path, flag_path, log_path, expected_shift):
    stub = types.SimpleNamespace(
        input=types.SimpleNamespace(full=full_path, sub=sub_path),
        output=types.SimpleNamespace(old=old_path, flag=flag_path),
        params=types.SimpleNamespace(
            expected_shift=expected_shift,
            resolution=5000,
            roi="escape",
            dataset="Mecp2_E6_WT_G1_Xi",
        ),
        log=[log_path],
    )
    src = open(SCRIPT).read()
    g = {"__name__": "__main__", "snakemake": stub}
    exec(compile(src, SCRIPT, "exec"), g)


def test_script_end_to_end():
    print("[4] workflow/scripts/py/metaloci_transplant_kk.py, end to end")
    # Mecp2 escape: full 73,315,000-74,475,000 (232 bins),
    #               escape 73,645,000-74,475,000 (166 bins), shift 66.
    full = _make_mlo(73_315_000, 74_475_000, 232, seed=1)
    sub = _make_mlo(73_645_000, 74_475_000, 166, seed=2)

    with tempfile.TemporaryDirectory() as tmp:
        full_path = os.path.join(tmp, "chrX_73315000_74475000_0.mlo")
        sub_path = os.path.join(tmp, "chrX_73645000_74475000_0.mlo")
        old_path = os.path.join(tmp, "chrX_73645000_74475000_0_old.mlo")
        flag_path = os.path.join(tmp, "chrX_73645000_74475000_0.transplanted")
        log_path = os.path.join(tmp, "logs", "transplant.log")
        for path, obj in ((full_path, full), (sub_path, sub)):
            with open(path, "wb") as fh:
                pickle.dump(obj, fh)

        pristine_sub_coords = sub["kk_coords"].copy()
        _run_script(full_path, sub_path, old_path, flag_path, log_path, 66)

        with open(sub_path, "rb") as fh:
            patched = pickle.load(fh)
        with open(old_path, "rb") as fh:
            backup = pickle.load(fh)

        check(
            "kk_coords == full.kk_coords[66:66+166]",
            np.array_equal(patched["kk_coords"], full["kk_coords"][66:232]),
        )
        check(
            "kk_distances == full[66:232, 66:232]",
            np.array_equal(patched["kk_distances"], full["kk_distances"][66:232, 66:232]),
        )
        check(
            "kk_restraints_matrix == full[66:232, 66:232]",
            np.array_equal(
                patched["kk_restraints_matrix"],
                full["kk_restraints_matrix"][66:232, 66:232],
            ),
        )
        check("kk_nodes re-keyed to 0..165", sorted(patched["kk_nodes"]) == list(range(166)))
        check(
            "kk_nodes[0] is full's node 66",
            patched["kk_nodes"][0] == full["kk_nodes"][66],
        )
        check("target kept its own start", patched["start"] == 73_645_000)
        check("target kept its own end", patched["end"] == 74_475_000)
        check("target kept its own resolution", patched["resolution"] == 5000)
        check(
            "backup holds the PRISTINE sub-region layout",
            np.array_equal(backup["kk_coords"], pristine_sub_coords),
        )
        check(
            "the transplant actually changed something",
            not np.array_equal(patched["kk_coords"], pristine_sub_coords),
        )

        stamp = dict(
            line.rstrip("\n").split("\t", 1) for line in open(flag_path) if "\t" in line
        )
        check("stamp records shift 66", stamp.get("shift") == "66", stamp.get("shift"))
        check("stamp records size 166", stamp.get("size") == "166", stamp.get("size"))
        check(
            "stamp names all four transplanted fields",
            stamp.get("fields")
            == "kk_nodes,kk_restraints_matrix,kk_coords,kk_distances",
            stamp.get("fields"),
        )

        # ---- idempotence: re-running must not corrupt anything -----------
        _run_script(full_path, sub_path, old_path, flag_path, log_path, 66)
        with open(sub_path, "rb") as fh:
            again = pickle.load(fh)
        check(
            "re-running is idempotent",
            np.array_equal(again["kk_coords"], patched["kk_coords"]),
        )
        with open(old_path, "rb") as fh:
            backup2 = pickle.load(fh)
        check(
            "re-running preserves the pristine backup",
            np.array_equal(backup2["kk_coords"], pristine_sub_coords),
        )


def test_script_rejects_wrong_shift():
    print("[5] the script refuses a wrong or negative expected_shift")
    full = _make_mlo(73_315_000, 74_475_000, 232, seed=1)
    sub = _make_mlo(73_645_000, 74_475_000, 166, seed=2)

    for bad_shift, label in ((-66, "negative (-66)"), (65, "off-by-one (65)")):
        with tempfile.TemporaryDirectory() as tmp:
            full_path = os.path.join(tmp, "full.mlo")
            sub_path = os.path.join(tmp, "sub.mlo")
            for path, obj in ((full_path, full), (sub_path, sub)):
                with open(path, "wb") as fh:
                    pickle.dump(obj, fh)
            raised = False
            try:
                _run_script(
                    full_path,
                    sub_path,
                    os.path.join(tmp, "sub_old.mlo"),
                    os.path.join(tmp, "sub.transplanted"),
                    os.path.join(tmp, "logs", "t.log"),
                    bad_shift,
                )
            except SystemExit:
                raised = True
            check(f"rejects {label}", raised, "the script accepted it")
            check(
                f"leaves the target untouched after rejecting {label}",
                np.array_equal(pickle.load(open(sub_path, "rb"))["kk_coords"], sub["kk_coords"]),
            )


if __name__ == "__main__":
    test_shifts_from_config()
    test_negative_shift_is_silent()
    test_script_end_to_end()
    test_script_rejects_wrong_shift()
    print()
    if failures:
        print(f"{len(failures)} FAILURE(S): {', '.join(failures)}")
        sys.exit(1)
    print("all transplant checks passed")
