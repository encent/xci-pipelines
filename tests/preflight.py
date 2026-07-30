#!/usr/bin/env python3
"""Pre-flight checks. Run by install.sh; safe to run any time.

Catches the problems a new user actually hits, before a 2-day run starts.
"""
from __future__ import annotations

import os
import shutil
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "workflow"))

import yaml  # noqa: E402

from lib import errors  # noqa: E402
from lib.paths import Paths  # noqa: E402

NEED_GB = 600
problems: list[str] = []


def check(label: str, ok: bool, detail: str = "") -> bool:
    print(f"  [{'ok' if ok else '!!'}] {label}" + (f" -- {detail}" if detail and not ok else ""))
    return ok


def main() -> int:
    cfg = yaml.safe_load(open("config/config.yaml"))
    P = Paths(cfg)

    # 1. conda
    if not check("conda on PATH", shutil.which("conda") is not None):
        problems.append(errors.conda_missing().render())

    # 2. cores
    n = os.cpu_count() or 1
    check(f"{n} CPU cores", True)
    if n < 4:
        print("       (fewer than 4 cores: expect a very long run)")

    # 3. disk
    os.makedirs(P.data, exist_ok=True)
    free_gb = shutil.disk_usage(P.data).free // (1024**3)
    if not check(f"{free_gb} GB free at {P.data}", free_gb >= NEED_GB):
        problems.append(errors.not_enough_disk(NEED_GB, free_gb, P.data).render())

    # 4. input directories
    for key, path in (cfg["paths"]["inputs"].get("bam_dirs") or {}).items():
        if not check(f"input dir {key}", os.path.isdir(path) and os.access(path, os.R_OK), path):
            problems.append(errors.missing_input_dir(key, path).render())

    cooldir = cfg["paths"]["inputs"].get("cooler_dir")
    if cooldir:
        check("cooler dir", os.path.isdir(cooldir), cooldir)

    # 5. sample sheet + BAM indexes
    sheet = cfg["paths"]["sheets"]["samples"]
    if os.path.exists(sheet):
        import csv

        missing_bam, missing_bai = [], []
        with open(sheet) as fh:
            for row in csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"):
                if row.get("use", "true").lower() not in ("true", "1", "yes"):
                    continue
                bam = row["bam"]
                if not os.path.exists(bam):
                    missing_bam.append(bam)
                elif not (os.path.exists(bam + ".bai") or os.path.exists(bam[:-4] + ".bai")):
                    missing_bai.append(bam)
        check(f"sample sheet ({sheet})", True)
        if not check("every BAM exists", not missing_bam, f"{len(missing_bam)} missing"):
            problems.append(
                f"\n  {len(missing_bam)} BAM file(s) in {sheet} do not exist, e.g.\n"
                f"    {missing_bam[0]}\n"
                "  Fix the paths, or set use=false on those rows.\n"
            )
        if not check("every BAM is indexed", not missing_bai, f"{len(missing_bai)} unindexed"):
            problems.append(errors.missing_bam_index(missing_bai[0]).render())
    else:
        check(f"sample sheet ({sheet})", False, "not generated yet")
        print("       run: python3 workflow/scripts/py/make_sheets.py")

    # 6. the chrX-in-csaw trap
    if not check("chrX not in the csaw chromosome set", "chrX" not in cfg["chromosomes"]["csaw"]):
        problems.append(errors.chrx_in_csaw_set(cfg["chromosomes"]["csaw"]).render())

    # 7. fixtures
    fx = "resources/fixtures/fixtures.sha256"
    check("fixtures present", os.path.exists(fx))

    if problems:
        print("\n" + "\n".join(problems))
        return 1
    print("\n  All pre-flight checks passed.\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
