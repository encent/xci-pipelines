"""The sample / cooler universe, loaded once and shared by every rule module.

Everything a rule needs to know about *which* jobs exist comes from here, so
that adding a clone is a sheet edit and never a code edit.

Key concept: the **normgroup**, the unit of one csaw/edgeR-TMM run.
The original used one run per mark for WT, and one run *per clone-group* for the
degrons -- all writing into the same output namespace. There are 17 of them
(4 WT + 3 H3K27me3 + 3 H3K27ac + 2 CTCF + 1 Rad21 + 4 RNA-Seq).
"""

from __future__ import annotations

import os
from collections import defaultdict
from typing import Dict, Iterable, List, Optional

import pandas as pd

from . import errors, naming

SAMPLE_COLUMNS = [
    "sample", "mark", "clone", "condition", "allele", "replicate",
    "genotype", "snpsplit_tag", "normgroup", "input_type", "bam", "use",
]
COOLER_COLUMNS = [
    "cooler", "locus", "clone", "condition", "snpsplit_tag", "allele",
    "replicate", "path", "merge_loops", "merge_comps", "use",
]


class SampleSheet:
    def __init__(self, config: dict, samples_tsv: str, coolers_tsv: str):
        self.config = config
        self.samples = self._load_samples(samples_tsv)
        self.coolers = self._load_coolers(coolers_tsv)
        self._validate()

    # -- loading ----------------------------------------------------------
    @staticmethod
    def _load_samples(path: str) -> pd.DataFrame:
        if not os.path.exists(path):
            raise errors.sample_sheet_problem(
                f"{path} does not exist.",
                "Run  bash install.sh  to generate it from the metadata sheet.",
            )
        df = pd.read_csv(path, sep="\t", dtype=str, comment="#")
        missing = [c for c in SAMPLE_COLUMNS if c not in df.columns]
        if missing:
            raise errors.sample_sheet_problem(
                f"{path} is missing column(s): {', '.join(missing)}",
                f"Required columns are:\n    {'  '.join(SAMPLE_COLUMNS)}",
            )
        df["replicate"] = df["replicate"].astype(int)
        df["use"] = df["use"].astype(str).str.lower().isin(("true", "1", "yes"))
        return df[df["use"]].reset_index(drop=True)

    @staticmethod
    def _load_coolers(path: str) -> pd.DataFrame:
        if not os.path.exists(path):
            return pd.DataFrame(columns=COOLER_COLUMNS)
        df = pd.read_csv(path, sep="\t", dtype=str, comment="#").fillna("")
        df["use"] = df["use"].astype(str).str.lower().isin(("true", "1", "yes"))
        return df[df["use"]].reset_index(drop=True)

    # -- validation -------------------------------------------------------
    def _validate(self) -> None:
        # chrX must never be in the csaw restrict set (D-04)
        csaw = self.config["chromosomes"]["csaw"]
        if "chrX" in csaw:
            raise errors.chrx_in_csaw_set(csaw)

        # names must fit the frozen grammar
        for s in self.samples["sample"]:
            try:
                naming.parse_sample(s)
            except naming.NamingError as exc:
                raise errors.sample_sheet_problem(
                    str(exc),
                    "Sample ids must be  {mark}_{clone}_{condition}_{allele}_rep{n}\n"
                    "e.g.  H3K27me3_E6_WT_Xi_rep1",
                )

        # the G1/G2 -> Xa/Xi map is a property of the clone, not of the row
        for clone, grp in self.samples.groupby("clone"):
            allelic = grp[grp["allele"].isin(("Xa", "Xi"))]
            mapping = defaultdict(set)
            for _, r in allelic.iterrows():
                if r["snpsplit_tag"]:
                    mapping[r["snpsplit_tag"]].add(r["allele"])
            bad = {k: v for k, v in mapping.items() if len(v) > 1}
            if bad:
                rows = "\n".join(
                    f"    {r['sample']}\tsnpsplit_tag={r['snpsplit_tag']}\tallele={r['allele']}"
                    for _, r in allelic.iterrows()
                )
                raise errors.inconsistent_allele_map(clone, rows)

    # -- queries used by the rule modules ---------------------------------
    @property
    def normgroups(self) -> List[str]:
        return sorted(self.samples["normgroup"].unique())

    def normgroup_samples(self, normgroup: str, allele: Optional[str] = None) -> List[str]:
        df = self.samples[self.samples["normgroup"] == normgroup]
        if allele:
            df = df[df["allele"] == allele]
        return sorted(df["sample"])

    def normgroup_of(self, sample: str) -> str:
        hit = self.samples.loc[self.samples["sample"] == sample, "normgroup"]
        if hit.empty:
            raise KeyError(f"unknown sample {sample!r}")
        return hit.iloc[0]

    def normgroup_input_type(self, normgroup: str) -> str:
        df = self.samples[self.samples["normgroup"] == normgroup]
        return df["input_type"].iloc[0] if len(df) else "all"

    def normgroup_mark(self, normgroup: str) -> str:
        df = self.samples[self.samples["normgroup"] == normgroup]
        return df["mark"].iloc[0]

    @property
    def all_samples(self) -> List[str]:
        return sorted(self.samples["sample"])

    def samples_of(self, **kw) -> List[str]:
        df = self.samples
        for k, v in kw.items():
            df = df[df[k] == v] if not isinstance(v, (list, tuple)) else df[df[k].isin(v)]
        return sorted(df["sample"])

    @property
    def tracks(self) -> List[str]:
        """Every replicate-merged track id."""
        return sorted({naming.track_of(s) for s in self.all_samples})

    def tracks_of(self, mark: Optional[str] = None, allele: Optional[str] = None) -> List[str]:
        out = []
        for t in self.tracks:
            d = naming.parse_track(t)
            if mark and d["mark"] != mark:
                continue
            if allele and d["allele"] != allele:
                continue
            out.append(t)
        return out

    def reps_of_track(self, track: str) -> List[str]:
        d = naming.parse_track(track)
        return sorted(
            self.samples_of(
                mark=d["mark"], clone=d["clone"],
                condition=d["condition"], allele=d["allele"],
            )
        )

    def bam_of(self, sample: str) -> str:
        return self.samples.loc[self.samples["sample"] == sample, "bam"].iloc[0]

    @property
    def marks(self) -> List[str]:
        return sorted(self.samples["mark"].unique())

    @property
    def clones(self) -> List[str]:
        return sorted(self.samples["clone"].unique())

    # -- AcMe3 ------------------------------------------------------------
    def acme3_tracks(self) -> List[str]:
        """clone_condition_allele triples having BOTH H3K27ac and H3K27me3."""
        def key(t):
            d = naming.parse_track(t)
            return (d["clone"], d["condition"], d["allele"])

        ac = {key(t) for t in self.tracks_of(mark="H3K27ac")}
        me = {key(t) for t in self.tracks_of(mark="H3K27me3")}
        return sorted(f"AcMe3_{c}_{cond}_{a}" for c, cond, a in ac & me)

    # -- consensus --------------------------------------------------------
    def consensus_members(self, mark: str, group: str, allele: str) -> List[str]:
        conds = self.config["normalization"]["consensus"]["groups"][group]["conditions"]
        if mark == "AcMe3":
            return sorted(
                t for t in self.acme3_tracks()
                if naming.parse_track(t)["condition"] in conds
                and naming.parse_track(t)["allele"] == allele
            )
        return sorted(
            t for t in self.tracks_of(mark=mark, allele=allele)
            if naming.parse_track(t)["condition"] in conds
        )

    def consensus_names(self) -> List[str]:
        cfg = self.config["normalization"]["consensus"]
        out = []
        for mark in cfg["marks"]:
            for group in cfg["groups"]:
                for allele in ("Xa", "Xi"):
                    if self.consensus_members(mark, group, allele):
                        out.append(naming.consensus_id(mark, group, allele))
        return sorted(out)

    # -- coolers ----------------------------------------------------------
    @property
    def cooler_names(self) -> List[str]:
        return sorted(self.coolers["cooler"])

    @staticmethod
    def _merge_cell(value: str) -> List[str]:
        """A merge cell is comma-separated: a dTAG cooler belongs to both its
        specific compartment group and the pooled ``*_dTAG_*`` group."""
        return [v.strip() for v in str(value).split(",") if v.strip()]

    def merged_cooler_names(self, scope: str) -> List[str]:
        col = "merge_loops" if scope == "merged_loops" else "merge_comps"
        if col not in self.coolers.columns:
            return []
        names = set()
        for cell in self.coolers[col]:
            names.update(self._merge_cell(cell))
        return sorted(names)

    def merged_cooler_members(self, scope: str, name: str) -> List[str]:
        col = "merge_loops" if scope == "merged_loops" else "merge_comps"
        if col not in self.coolers.columns:
            return []
        hit = self.coolers[col].map(lambda c: name in self._merge_cell(c))
        return sorted(self.coolers.loc[hit, "cooler"])

    def cooler_path(self, name: str) -> str:
        """Source path of a cooler, or "" when it is produced by the pipeline
        (the replicate merges: `path` is empty, `merge_f3_replicates` builds it)."""
        hit = self.coolers.loc[self.coolers["cooler"] == name, "path"]
        return hit.iloc[0] if len(hit) else ""

    def replicate_members(self, name: str) -> List[str]:
        """Replicate coolers behind a replicate-merged cooler, or []."""
        pref = f"{name}_rep"
        return sorted(c for c in self.coolers["cooler"] if c.startswith(pref))

    def metaloci_datasets(self, run: str) -> List[str]:
        """Dataset list per METALoci run, matching the original's filters."""
        if run == "consensus":
            return sorted(
                n for n in self.merged_cooler_names("merged_comps")
                if "CTCF-dTAG" in n or "NodTAG-or-WT" in n
            )
        if run == "wt":
            return sorted(n for n in self.cooler_names if "_WT_" in n)
        if run == "degron":
            return sorted(
                n for n in self.cooler_names
                if "TAG" in n and not n.endswith(("_rep1", "_rep2"))
            )
        raise ValueError(f"unknown METALoci run {run!r}")

    def metaloci_signals_for(self, dataset: str) -> List[str]:
        """Per-dataset signal subset, from config.metaloci.sample_exclusions."""
        signals = list(self.config["metaloci"]["signals"])
        out = []
        for sig in signals:
            keep = True
            for rule in self.config["metaloci"].get("sample_exclusions", []):
                if "clone" in rule and rule["clone"] in dataset:
                    if sig not in rule.get("signals_except", []):
                        keep = False
                if rule.get("signal") == sig:
                    only = rule.get("only_clones")
                    excl = rule.get("exclude_clones")
                    if only and not any(c in dataset for c in only):
                        keep = False
                    if excl and any(c in dataset for c in excl):
                        keep = False
            if keep:
                out.append(sig)
        return out


def load(config: dict) -> SampleSheet:
    sheets = config["paths"]["sheets"]
    return SampleSheet(config, sheets["samples"], sheets["coolers"])
