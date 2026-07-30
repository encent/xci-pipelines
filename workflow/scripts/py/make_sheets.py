#!/usr/bin/env python3
"""Generate config/samples.tsv and config/coolers.tsv from the data on disk.

Run by install.sh. Re-runnable and idempotent; refuses to overwrite a sheet the
user has edited unless --force is given.

Sample sheet source of truth, in order:
  1. the BAM files themselves (their names already encode the full grammar)
  2. MetaDataSheet_CutaRun_allsamples_CTCFpaper.csv, for the per-clone
     G1/G2 -> Xa/Xi map and the parental genotype
"""
from __future__ import annotations

import argparse
import os
import re
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from lib import naming  # noqa: E402

# Normalisation groups: WT = one csaw run per mark over all clones;
# degron = one run PER CLONE-GROUP. 17 in total.
DEGRON_CLONE_GROUPS = {"F3", "E6A7", "B1621", "C5C10"}

# Per-clone Xi parental genome, from the paper (Methods, Table H.1).
XI_GENOTYPE = {
    "E6": "G1", "JTG": "G1", "CL30": "G1",     # C57BL/6J-Xi
    "B1": "G2", "C5": "G2",                     # Cast/EiJ-Xi
    "E6A7": "G1", "F3": "G1", "C5C10": "G2", "B1621": "G2",
}
GENOTYPE_NAME = {"G1": "C57BL-6J", "G2": "CAST-EiJ"}

BAM_RE = re.compile(
    r"^(?P<mark>H3K27me3|H3K27ac|CTCF|Rad21|RNA-Seq)_"
    r"(?P<clone>[A-Za-z0-9]+)_"
    r"(?P<condition>WT|CTCF-dTAG|CTCF-NodTAG|Rad21-dTAG|Rad21-NodTAG)_"
    r"(?P<allele>Gall|Xa|Xi)_rep(?P<rep>[12])\.bam$"
)
COOL_RE = re.compile(
    r"^(?P<locus>Jarid|Mecp2)_(?P<clone>[A-Za-z0-9]+)_"
    r"(?P<condition>WT|CTCF-dTAG|CTCF-NodTAG|Rad21-dTAG|Rad21-NodTAG)_"
    r"(?P<tag>G1|G2)_(?P<allele>Xa|Xi)(?:_rep(?P<rep>[12]))?\.cool$"
)


def normgroup(mark: str, condition: str, clone: str) -> str:
    if condition == "WT":
        return f"WT_{mark}"
    group = clone if clone in DEGRON_CLONE_GROUPS else clone
    return f"dTAG_{mark}_{group}"


def scan_bams(bam_dirs: dict) -> pd.DataFrame:
    rows = []
    for key, d in bam_dirs.items():
        if not os.path.isdir(d):
            print(f"  [skip] {key}: {d} not found", file=sys.stderr)
            continue
        n = 0
        for fn in sorted(os.listdir(d)):
            m = BAM_RE.match(fn)
            if not m:
                continue
            g = m.groupdict()
            tag = ""
            if g["allele"] in ("Xa", "Xi"):
                xi = XI_GENOTYPE.get(g["clone"])
                if xi:
                    tag = xi if g["allele"] == "Xi" else ("G2" if xi == "G1" else "G1")
            rows.append(dict(
                sample=fn[:-4],
                mark=g["mark"], clone=g["clone"], condition=g["condition"],
                allele=g["allele"], replicate=int(g["rep"]),
                genotype=GENOTYPE_NAME.get(tag, ""), snpsplit_tag=tag,
                normgroup=normgroup(g["mark"], g["condition"], g["clone"]),
                input_type="all",
                bam=os.path.join(d, fn), use="true",
            ))
            n += 1
        print(f"  {key}: {n} BAMs", file=sys.stderr)
    return pd.DataFrame(rows)


def apply_input_type_filters(df: pd.DataFrame, filters: dict) -> pd.DataFrame:
    """Reproduce the `selected` runs. Only WT_H3K27ac is live (drops CL30 rep2)."""
    for ng, spec in (filters or {}).items():
        mask = df["normgroup"] == ng
        if not mask.any():
            continue
        df.loc[mask, "input_type"] = "selected"
        for pat in spec.get("drop", []):
            drop = mask & df["sample"].str.contains(pat, regex=True)
            df.loc[drop, "use"] = "false"
            print(f"  {ng}: input_type=selected, dropped {int(drop.sum())} "
                  f"sample(s) matching {pat!r}", file=sys.stderr)
    return df


def scan_coolers(cooler_dir: str) -> pd.DataFrame:
    rows = []
    if not os.path.isdir(cooler_dir):
        print(f"  [skip] coolers: {cooler_dir} not found", file=sys.stderr)
        return pd.DataFrame(columns=[
            "cooler", "locus", "clone", "condition", "snpsplit_tag", "allele",
            "replicate", "path", "merge_loops", "merge_comps", "use"])
    for fn in sorted(os.listdir(cooler_dir)):
        m = COOL_RE.match(fn)
        if not m:
            continue
        g = m.groupdict()
        rows.append(dict(
            cooler=fn[:-5], locus=g["locus"], clone=g["clone"],
            condition=g["condition"], snpsplit_tag=g["tag"], allele=g["allele"],
            replicate=g["rep"] or "", path=os.path.join(cooler_dir, fn),
            merge_loops="", merge_comps="", use="true",
        ))
    df = pd.DataFrame(rows)
    df = add_replicate_merges(df)
    print(f"  coolers: {len(df)} rows", file=sys.stderr)
    return df


def add_replicate_merges(df: pd.DataFrame) -> pd.DataFrame:
    """Add the replicate-merged parent of every cooler that has replicates.

    `00_03_merge_F3_coolers.sh` sums `_rep1` + `_rep2` with `cooler merge`; F3 is
    the only clone with two Capture Hi-C replicates, so this adds exactly the
    four `Mecp2_F3_CTCF-{dTAG,NodTAG}_{G2_Xa,G1_Xi}` rows -- 48 raw + 4 = the 52
    individual mcools of the ground truth, and the 52 keys of
    `resources/fixtures/compartments/compartments_{roi}.json`.

    A merged row carries an EMPTY `path`: it is produced by
    `merge_f3_replicates`, not symlinked from the user's cooler directory.
    """
    if df.empty:
        return df
    have = set(df["cooler"])
    extra = []
    for stem, grp in df[df["replicate"] != ""].groupby(
        df["cooler"].str.replace(r"_rep[12]$", "", regex=True)
    ):
        if len(grp) < 2 or stem in have:
            continue
        r = grp.iloc[0].to_dict()
        r.update(cooler=stem, replicate="", path="")
        extra.append(r)
    if extra:
        names = ", ".join(sorted(r["cooler"] for r in extra))
        print(f"  replicate merges: {len(extra)} ({names})", file=sys.stderr)
        df = pd.concat([df, pd.DataFrame(extra)], ignore_index=True)
        df = df.sort_values("cooler").reset_index(drop=True)
    return df


def assign_merge_membership(df: pd.DataFrame, membership_fixture: str) -> pd.DataFrame:
    """Membership is a FIXTURE (D-10), not a rule -- it was enumerated by hand.

    A cooler can be in SEVERAL merged_comps groups: every dTAG cooler is in both
    its specific group (`*_CTCF-dTAG_*` / `*_Rad21-dTAG_*`) and the pooled
    `*_dTAG_*` group. The cell is therefore comma-separated, and
    `lib/samples.py` splits on commas.
    """
    if not os.path.exists(membership_fixture) or df.empty:
        return df
    mem = pd.read_csv(membership_fixture, sep="\t", dtype=str, comment="#")
    groups = {"merge_loops": set(), "merge_comps": set()}
    unknown = []
    for _, r in mem.iterrows():
        hit = df["cooler"] == r["cooler"]
        if not hit.any():
            unknown.append(r["cooler"])
            continue
        col = "merge_loops" if r["scope"] == "merged_loops" else "merge_comps"
        groups[col].add(r["merged_name"])
        cur = [v for v in str(df.loc[hit, col].iloc[0]).split(",") if v]
        if r["merged_name"] not in cur:
            df.loc[hit, col] = ",".join(cur + [r["merged_name"]])
    if unknown:
        print(f"  [warn] membership references {len(unknown)} unknown cooler(s): "
              f"{', '.join(sorted(set(unknown)))}", file=sys.stderr)
    print(f"  merged coolers: {len(groups['merge_loops'])} for loops, "
          f"{len(groups['merge_comps'])} for compartments", file=sys.stderr)
    return df


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/config.yaml")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    import yaml
    cfg = yaml.safe_load(open(args.config))
    inputs = cfg["paths"]["inputs"]

    print("Scanning input directories...", file=sys.stderr)
    samples = scan_bams(inputs.get("bam_dirs", {}))
    if samples.empty:
        print("\nNo BAM files matched the expected naming pattern.\n"
              "Expected:  {mark}_{clone}_{condition}_{allele}_rep{n}.bam\n"
              "Check paths.inputs.bam_dirs in config/config.yaml.", file=sys.stderr)
    samples = apply_input_type_filters(samples, cfg["normalization"].get("input_type_filters"))

    coolers = scan_coolers(inputs.get("cooler_dir", ""))
    coolers = assign_merge_membership(coolers, "resources/fixtures/cooler_membership.tsv")

    for path, df in ((cfg["paths"]["sheets"]["samples"], samples),
                     (cfg["paths"]["sheets"]["coolers"], coolers)):
        if os.path.exists(path) and not args.force:
            print(f"\n{path} already exists -- not overwriting (use --force).", file=sys.stderr)
            continue
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        with open(path, "w") as fh:
            fh.write("# Generated by workflow/scripts/py/make_sheets.py -- safe to edit.\n")
            fh.write("# Set use=false to drop a row without deleting it.\n")
            df.to_csv(fh, sep="\t", index=False)
        print(f"  wrote {path}  ({len(df)} rows)", file=sys.stderr)

    if not samples.empty:
        used = samples[samples["use"] == "true"]
        print(f"\n{len(used)} samples in {used['normgroup'].nunique()} normalisation groups; "
              f"{used['clone'].nunique()} clones, {used['mark'].nunique()} marks.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
