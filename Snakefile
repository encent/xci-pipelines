# =====================================================================
#  XCI PIPELINE
#
#  One Snakemake workflow: allele-split BAMs + Capture Hi-C coolers
#  -> normalised tracks -> H3K27me3 valleys, Hi-C features, METALoci
#  -> ONE terminal visualisation stage that emits every paper panel.
#
#  Users should run  ./run.sh , not snakemake directly.
#  Design: the design notes
# =====================================================================
import os
import sys

sys.path.insert(0, os.path.join(workflow.basedir, "workflow"))

from lib import naming
from lib.paths import Paths
from lib import samples as _samples


configfile: "config/config.yaml"


# --- resolved once, shared by every rule module ------------------------
P = Paths(config)
SS = _samples.load(config)

DOWNSTREAM_CHROMS = config["chromosomes"]["downstream"]
VALLEY_CHROMS = DOWNSTREAM_CHROMS + config["chromosomes"].get("valley_control", [])
ROI_TYPES = list(naming.ROI_TYPES)
LOCI = [k for k in config["loci"] if k != "display_names"]

os.makedirs(P.tmp, exist_ok=True)


# --- wildcard constraints: an ambiguous match must fail at DAG build ---
wildcard_constraints:
    **naming.WC


# --- rule modules -------------------------------------------------------
include: "workflow/rules/00_reference.smk"
include: "workflow/rules/10_signal.smk"
include: "workflow/rules/20_contacts.smk"
include: "workflow/rules/30_valleys.smk"
include: "workflow/rules/31_hic_features.smk"
include: "workflow/rules/32_metaloci.smk"
include: "workflow/rules/40_staging.smk"
include: "workflow/rules/50_qc.smk"
include: "workflow/rules/90_visualise.smk"


# --- default target -----------------------------------------------------
rule all:
    input:
        P.visualise_done(),


# --- convenience targets ------------------------------------------------
# Section ergonomics without section barriers: these are aliases with input
# lists, not sub-workflows, so a partial target still pulls its true
# prerequisites and nothing more.


rule reference:
    input:
        P.chrom_sizes(),
        P.blacklist(),
        P.roi_table(),
        expand(P.genome("ice_blacklists/{locus}.bed"), locus=LOCI),


rule signal:
    input:
        expand("{p}", p=[P.merged20(naming.parse_track(t)["mark"], t) for t in SS.tracks]),
        expand("{p}", p=[P.merged5k(naming.parse_track(t)["mark"], t) for t in SS.tracks]),
        expand("{p}", p=[P.bedgraph5k(naming.parse_track(t)["mark"], t) for t in SS.tracks]),
        expand("{p}", p=[P.acme3_bg5k(t) for t in SS.acme3_tracks()]),
        expand("{p}", p=[P.consensus5k(n) for n in SS.consensus_names()]),


rule contacts:
    input:
        expand("{p}", p=[P.mcool("individual", n) for n in SS.cooler_names]),
        expand("{p}", p=[P.mcool("merged_loops", n) for n in SS.merged_cooler_names("merged_loops")]),
        expand("{p}", p=[P.mcool("merged_comps", n) for n in SS.merged_cooler_names("merged_comps")]),


rule valleys:
    input:
        # Assembled in 30_valleys.smk, next to the rules and the eligibility
        # logic. Covers the whole module: calling, the gene filter, boundaries
        # and their motif split, gene content, CTCF status, both overlap
        # analyses, the size table, BOTH boundary-metaplot paths (R-2) and the
        # degron density tables.
        VALLEY_TARGETS,


rule hic:
    input:
        rules.contacts.input,
        # LOOP_SETS, not every merged-loop cooler: only the four
        # {locus}_{allele} sets were ever hand-curated, so only they have a
        # refined fixture. See the LOOP_SETS comment in 31_hic_features.smk.
        expand("{p}", p=[P.loops_refined(n) for n in LOOP_SETS]),
        expand("{p}", p=[P.ml_compartmentalization(r, "full") for r in config["metaloci"]["runs"]]),


rule figures:
    input:
        P.visualise_done(),
