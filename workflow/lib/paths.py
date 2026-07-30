"""Every path in the pipeline, derived from config. Nothing is hard-coded.

This module is the **frozen output contract** of the re-developed pipeline
(the design notes section 2.2). Panel rules and test comparators are written
against these functions, which is what lets the visualisation stage be built in
parallel with the feature stages.

Two invariants, both enforced by tests:
  1. every function here returns a path under ``config.paths.data_dir``;
  2. no rule outside ``90_visualise.smk`` may return a path under ``figures()``.
"""

from __future__ import annotations

import os
from typing import Optional


class Paths:
    """Resolves all pipeline paths from the config dict."""

    def __init__(self, config: dict):
        p = config["paths"]
        self.data = os.path.abspath(_expand(p["data_dir"], p))
        self.resources = os.path.abspath(_expand(p["resources_dir"], p))
        self.results = os.path.abspath(_expand(p["results_dir"], p))
        self.logs = os.path.abspath(_expand(p["log_dir"], p))
        self.tmp = os.path.abspath(_expand(p.get("tmpdir", "{data_dir}/tmp"), p))
        self.repo = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

    # -- generic ----------------------------------------------------------
    def work(self, *parts: str) -> str:
        return os.path.join(self.data, "work", *parts)

    def raw(self, *parts: str) -> str:
        return os.path.join(self.data, "raw", *parts)

    def res(self, *parts: str) -> str:
        return os.path.join(self.resources, *parts)

    def out(self, *parts: str) -> str:
        return os.path.join(self.results, *parts)

    def log(self, rule: str, name: str = "") -> str:
        return os.path.join(self.logs, rule, f"{name}.log" if name else "run.log")

    def benchmark(self, rule: str, name: str = "") -> str:
        return os.path.join(self.logs, "benchmarks", rule, f"{name}.tsv" if name else "run.tsv")

    def fixture(self, *parts: str) -> str:
        return os.path.join(self.repo, "resources", "fixtures", *parts)

    # -- 00 reference -----------------------------------------------------
    def genome(self, *parts: str) -> str:
        return self.res("genome", *parts)

    def chrom_sizes(self) -> str:
        return self.genome("mm10.chrom.sizes")

    def blacklist(self) -> str:
        return self.genome("mm10-blacklist.v2.bed")

    def fasta(self) -> str:
        return self.genome("mm10.fa")

    def gtf(self) -> str:
        return self.genome("Mus_musculus.GRCm38.102.chr.gtf")

    def gtf_noncoding(self) -> str:
        return self.genome("GRCm38.102_NC.gtf")

    def ctcf_motifs(self, chrx_only: bool = False) -> str:
        return self.genome("CTCF_mm10_X_only.bed" if chrx_only else "CTCF_mm10.bed")

    def tss(self) -> str:
        return self.genome("mm10_TSS.bed")

    def roi_table(self) -> str:
        return self.genome("regions_of_interest.tsv")

    def ice_blacklist(self, locus: str) -> str:
        return self.genome("ice_blacklists", f"{locus}.bed")

    def chrom_subset(self, scope: str) -> str:
        return self.genome(f"chroms_{scope}.txt")

    def allelic_gtf(self, clone: str, ext: str = "gtf") -> str:
        return self.genome("allelic_gtf", f"GRCm38.102_{clone}.{ext}")

    # -- 10 signal --------------------------------------------------------
    def bam(self, sample: str) -> str:
        return self.raw("bam", f"{sample}.bam")

    def scalefactors(self, normgroup: str, kind: str = "factors") -> str:
        """kind in {counts, totals, factors, allsamples}."""
        ext = "tsv.gz" if kind == "counts" else "tsv"
        return self.work("scalefactors", f"{normgroup}.{kind}.{ext}")

    def rep10(self, mark: str, sample: str) -> str:
        return self.work("tracks", "rep10", mark, f"{sample}.bw")

    def merged20(self, mark: str, track: str) -> str:
        return self.work("tracks", "merged20", mark, f"{track}.bw")

    def merged5k(self, mark: str, track: str) -> str:
        return self.work("tracks", "merged5k", mark, f"{track}.bw")

    def rep5k(self, mark: str, sample: str) -> str:
        return self.work("tracks", "rep5k", mark, f"{sample}.bw")

    def bedgraph5k(self, mark: str, track: str) -> str:
        return self.work("tracks", "bedgraph5k", mark, f"{track}.bed")

    def acme3_bw20(self, track: str) -> str:
        return self.work("tracks", "acme3_bw20", f"{track}.bw")

    def acme3_bg5k(self, track: str) -> str:
        return self.work("tracks", "acme3_bg5k", f"{track}.bed")

    def consensus5k(self, name: str) -> str:
        return self.work("tracks", "consensus5k", f"{name}.bed")

    # -- 20 contacts ------------------------------------------------------
    def cool_raw(self, name: str) -> str:
        return self.work("contacts", "raw", f"{name}.cool")

    def cool_corrected(self, name: str) -> str:
        return self.work("contacts", "corrected", f"{name}.cool")

    def cool_merged(self, scope: str, name: str) -> str:
        """scope in {merged_loops, merged_comps}."""
        return self.work("contacts", scope, f"{name}.cool")

    def mcool(self, scope: str, name: str) -> str:
        """scope in {individual, merged_loops, merged_comps}."""
        return self.work("contacts", "mcool", scope, f"{name}.mcool")

    def hic(self, scope: str, name: str) -> str:
        return self.work("contacts", "hic", scope, f"{name}.hic")

    # -- 30 valleys -------------------------------------------------------
    def valleys(self, chrom: str, track: str, ext: str = "bed") -> str:
        return self.work("features", "valleys", chrom, f"{track}_valleys.{ext}")

    def valley_states(self, chrom: str, track: str) -> str:
        return self.work("features", "valleys", chrom, f"{track}_states.tsv.gz")

    def valleys_fig4d_fixture(self, track: str) -> str:
        """The published valley BED for one of the four Fig 4d tracks (D-18).

        Read-only fixture; the pipeline's own call for the same track lives at
        `valleys()` and is always produced, so the two can be compared.
        """
        return self.fixture("valleys_fig4d", f"{track}_valleys.bed")

    def valleys_filtered(self, chrom: str, track: str) -> str:
        return self.work("features", "valleys_filtered", chrom, f"{track}_valleys.bed")

    def boundaries(self, variant: str, track: str) -> str:
        """variant in {all, motif_yes, motif_no}."""
        return self.work("features", "boundaries", variant, f"{track}_boundaries.bed")

    def boundary_ctcf(self, track: str) -> str:
        return self.work("features", "boundary_ctcf", f"{track}.tsv")

    def gene_content(self, track: str) -> str:
        return self.work("features", "gene_content", f"{track}.tsv")

    def valley_overlap(self, clone: str, scope: str) -> str:
        """scope in {all, escaping}."""
        return self.work("features", "valley_overlap", f"{clone}_{scope}.tsv")

    def stackup(self, track: str, signal: str, variant: str) -> str:
        return self.work("features", "stackups", f"{track}_{signal}_{variant}.npz")

    def boundary_profile(self, track: str, signal: str) -> str:
        """panel_boundary_analysis input: NO background subtraction (R-2)."""
        return self.work("features", "boundary_profile", f"{track}_{signal}.npz")

    def boundary_dtag(self, clone: str, signal: str) -> str:
        return self.work("features", "boundary_dtag", f"{clone}_{signal}.npz")

    def density(self, clone: str, allele: str, win: str, mask: str) -> str:
        return self.work("features", "density", f"{clone}_{allele}_{win}_{mask}.tsv")

    def valley_sizes(self) -> str:
        return self.work("features", "valley_sizes.tsv")

    def valley_xa_xi(self, track: str) -> str:
        return self.work("features", "valley_xa_xi", f"{track}.tsv")

    # -- 31 hic features --------------------------------------------------
    def loops_raw(self, name: str, ext: str = "tsv") -> str:
        return self.work("features", "loops", "raw", name, f"{name}.{ext}")

    def loops_refined(self, name: str, ext: str = "tsv") -> str:
        return self.work("features", "loops", "refined", name, f"{name}.{ext}")

    def loop_stats(self, comparison: str, roi: str, suffix: str = "") -> str:
        stem = f"loop_stats{suffix}.tsv"
        return self.work("features", "loops", "stats", comparison, roi, stem)

    def pileup(self, comparison: str, roi: str, exp: str, sample: str) -> str:
        return self.work("features", "pileups", comparison, roi, exp, f"pileup_{sample}.npz")

    def pileup_score(self, comparison: str, roi: str) -> str:
        return self.work("features", "pileups", comparison, roi, "scores.tsv")

    def loop_anchor_stackup(self, comparison: str, roi: str, signal: str, pair: str) -> str:
        return self.work("features", "loops", "stackups", comparison, roi, f"{signal}_{pair}.npz")

    def gc_track(self, roi: str, name: str) -> str:
        return self.work("features", "compartments", "gc", roi, f"{name}.tsv")

    def gc_track_canonical(self, binkey: str) -> str:
        """The GC track computed ONCE per distinct bin table.

        `bioframe.frac_gc` sees only the bin table and the FASTA -- never the
        cooler or the ROI -- so every cooler sharing a bin table has the same
        GC track. `binkey` identifies the bin table (chromosome set +
        resolution), deliberately NOT a cooler name: the identity is a property
        of the binning, and naming it after one arbitrary cooler would hide that.
        """
        return self.work("features", "compartments", "gc", "_canonical", f"{binkey}.tsv")

    def eigs(self, roi: str, name: str, which: Optional[str] = None) -> str:
        base = self.work("features", "compartments", "eigs", roi, name)
        if which is None:
            return os.path.join(base, f"Comp_{name}.tsv")
        return os.path.join(base, f"Comp_{which}_{name}.bw")

    def saddle(self, roi: str, name: str, eig: str) -> str:
        return self.work("features", "compartments", "saddle", roi, name, f"Saddle_{eig}_{name}.npz")

    def saddle_values(self, roi: str, name: str) -> str:
        return self.work("features", "compartments", "saddle", roi, name, f"Saddle_values_{name}.tsv")

    def comps_refined(self, scope: str, roi: str, name: str, fname: str) -> str:
        """scope in {refined, refined_merged}."""
        return self.work("features", "compartments", scope, roi, name, fname)

    def saddle_strength_selected(self, scope: str, roi: str) -> str:
        return self.work("features", "compartments", scope, roi, f"saddle_strength_selected_{roi}.tsv")

    def allelic_ratio_stats(self, locus: str, roi: str, degron: bool = False) -> str:
        suffix = "_degron" if degron else ""
        return self.work("features", "allelic_ratio", f"{locus}_{roi}_allelic_ratio_stats{suffix}.csv")

    # -- 32 metaloci ------------------------------------------------------
    def ml_dir(self, run: str, dataset: str) -> str:
        return self.work("features", "metaloci", run, dataset)

    def ml_signals(self, run: str, dataset: str) -> str:
        return os.path.join(self.ml_dir(run, dataset), f"{dataset}.signals")

    def ml_mlo(self, run: str, dataset: str, start: int, end: int) -> str:
        return os.path.join(self.ml_dir(run, dataset), "chrX", "objects", f"chrX_{start}_{end}_0.mlo")

    def ml_moran(self, run: str, dataset: str) -> str:
        return os.path.join(self.ml_dir(run, dataset), "moran_info.txt")

    def ml_compartmentalization(self, run: str, roi: str) -> str:
        return self.work("features", "metaloci", run, f"compartmentalization_{roi}.tsv")

    def ml_compartmentalization_final(self) -> str:
        return self.work("features", "metaloci", "compartmentalization_full_final.tsv")

    # -- 40 staging -------------------------------------------------------
    def staged_bigwig(self, track: str) -> str:
        return self.work("staging", "bigwigs_for_coolbox", f"{track}.bw")

    def staged_boundary(self, track: str, ext: str = "bed") -> str:
        return self.work("staging", "boundaries_for_coolbox", f"{track}_valleys.{ext}")

    def fixtures_verified(self) -> str:
        return self.work("fixtures", ".verified")

    # -- 50 qc ------------------------------------------------------------
    def qc(self, *parts: str) -> str:
        return self.work("qc", *parts)

    # -- 90 visualise -----------------------------------------------------
    # NOTHING outside 90_visualise.smk may write under figures().
    def figures(self, group: str, name: str, ext: str = "svg") -> str:
        return self.out("figures", group, f"{name}.{ext}")

    def table(self, name: str) -> str:
        return self.out("tables", name)

    def manifest(self) -> str:
        return self.out("figure_manifest.tsv")

    def paper_numbers(self) -> str:
        return self.out("paper_numbers.tsv")

    def report(self) -> str:
        return self.out("REPORT.html")

    def provenance(self) -> str:
        return self.out("provenance.json")

    def visualise_done(self) -> str:
        return self.out(".visualise.done")


def _expand(template: str, paths_cfg: dict) -> str:
    """Expand ``{data_dir}`` inside a path template."""
    return template.format(data_dir=paths_cfg["data_dir"])
