#!/usr/bin/env python3
"""Run one panel script against the stubs, without Snakemake.

Snakemake injects a magic ``snakemake`` object into ``script:`` files. This
builds a minimal stand-in so a panel can be exercised end to end -- producing a
real SVG on disk -- from a shell, in seconds, with no DAG, no conda envs and no
heavy data.

That matters more than it sounds: a panel bug that only appears under a full
run costs a 30-hour round trip to find. This costs two seconds.

Usage::

    python tests/stubs/make_stubs.py --out /tmp/stub
    python tests/stubs/run_panel.py --data /tmp/stub --panel valley_gene_content_wt
    python tests/stubs/run_panel.py --data /tmp/stub --all

``--all`` runs every panel whose stub inputs exist and prints a pass/fail table.
"""

from __future__ import annotations

import argparse
import os
import sys
import traceback

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, "workflow"))

import yaml                                                    # noqa: E402

from lib import panels as _panels                              # noqa: E402
from lib.paths import Paths                                    # noqa: E402


class Fake(object):
    """The subset of the ``snakemake`` object a panel script touches."""

    def __init__(self, inputs, output, params, log, config, threads=1):
        self.input = _Input(inputs)
        self.output = _Output([output])
        self.params = _Params(params)
        self.log = [log]
        self.config = config
        self.threads = threads
        self.wildcards = _Params({})
        self.resources = _Params({})


class _Input(list):
    def __init__(self, items):
        super().__init__(str(i) for i in items)


class _Output(list):
    pass


class _Params(object):
    def __init__(self, mapping):
        self._m = dict(mapping)

    def __getattr__(self, name):
        try:
            return self._m[name]
        except KeyError:
            raise AttributeError(name)

    def __getitem__(self, key):
        return self._m[key]

    def get(self, key, default=None):
        return self._m.get(key, default)


#: Panel producer -> the script that draws it.
SCRIPTS = {
    "panel_valley_gene_content": "valley_gene_content.py",
    "panel_boundary_ctcf_pies": "boundary_ctcf_pies.py",
    "panel_boundary_analysis": "boundary_analysis.py",
    "panel_boundary_dtag": "boundary_dtag.py",
    "panel_loops_comps_grid": "loops_comps_grid.py",
    "panel_scatter_delta": "scatter_delta.py",
    "panel_violin_loops": "violin_loops.py",
    "panel_violin_saddle": "violin_saddle.py",
    "panel_density_ma": "density_ma.py",
    "panel_valley_density": "density_ma.py",
    "panel_valley_venn": "valley_venn.py",
    "panel_valley_locus": "valley_locus.py",
    "panel_valley_sizes": "valley_sizes.py",
    "panel_valley_xa_xi_venn": "valley_xa_xi_venn.py",
    "panel_metaloci_composite": "metaloci_composite.py",
    "panel_metaloci_violin": "metaloci_panels.py",
    "panel_metaloci_heatmap": "metaloci_panels.py",
    "panel_metaloci_scatter": "metaloci_panels.py",
    "panel_stackup_composite": "stackup_composite.py",
    "panel_pileup": "hic_panels.py",
    "panel_saddle": "hic_panels.py",
    "panel_eigenvector": "hic_panels.py",
    "panel_loop_overlay": "hic_panels.py",
    "panel_cooler_qc": "hic_panels.py",
    "panel_ma_plot": "qc_panels.py",
    "panel_frip": "qc_panels.py",
    "panel_peak_venn": "qc_panels.py",
    "panel_correlation": "qc_panels.py",
    "panel_fragment_sizes": "qc_panels.py",
    "panel_signal_enrichment": "qc_panels.py",
}

#: Inputs each producer needs, resolved against `Paths`. Kept here rather than
#: imported from 90_visualise.smk because that file is Snakemake syntax; the
#: duplication is small and this harness is a test, not the pipeline.
def inputs_for(P, panel, config):
    w = panel.wildcards
    producer = panel.producer
    wt = ["H3K27me3_{}_WT_Xi".format(c)
          for c in ("E6", "C5", "B1", "JTG", "CL30")]
    degron_mark = {"E6A7": "CTCF", "F3": "CTCF", "B1621": "Rad21"}
    nod = ["H3K27me3_{}_{}-NodTAG_Xi".format(c, m)
           for c, m in degron_mark.items()]
    dta = ["H3K27me3_{}_{}-dTAG_Xi".format(c, m)
           for c, m in degron_mark.items()]
    signals = ["H3K27me3", "H3K27ac", "CTCF", "RNA-Seq"]

    def sigs(track):
        clone = track.split("_")[1]
        out = [s for s in signals if not (s == "CTCF" and clone == "B1621")]
        if clone == "B1621":
            out.append("Rad21")
        return out

    if producer == "panel_valley_gene_content":
        tracks = wt if w.get("cohort") == "WT" else nod + dta
        return [P.gene_content(t) for t in tracks]
    if producer == "panel_boundary_ctcf_pies":
        cohort = w.get("cohort")
        tracks = wt if cohort == "WT" else (nod if cohort == "degron"
                                            else wt + nod)
        return ([P.boundary_ctcf(t) for t in tracks]
                + [P.gene_content(t) for t in tracks])
    if producer == "panel_boundary_analysis":
        clone = w["clone"]
        track = next((t for t in (wt + nod) if t.split("_")[1] == clone), None)
        if track is None:
            return []
        return ([P.boundary_profile(track, s) for s in sigs(track)]
                + [P.boundary_ctcf(track), P.gene_content(track)])
    if producer == "panel_boundary_dtag":
        clone = w["clone"]
        track = next((t for t in nod if t.split("_")[1] == clone), None)
        if track is None:
            return []
        return ([P.boundary_dtag(clone, s) for s in sigs(track)]
                + [P.boundary_ctcf(track), P.gene_content(track)])
    if producer == "panel_valley_venn":
        return [P.valley_overlap(c, w.get("scope", "escaping"))
                for c in degron_mark]
    if producer == "panel_violin_loops":
        return [P.pileup_score(w.get("comparison", "Xa_vs_Xi"),
                               w.get("roi", "full"))]
    if producer == "panel_violin_saddle":
        return [P.saddle_strength_selected(w.get("scope", "refined_merged"),
                                           w.get("roi", "full"))]
    if producer == "panel_density_ma":
        return [P.density(c, w.get("allele", "Xi"), "100kb", m)
                for c in degron_mark
                for m in ("antivalley", "allcoverage")]
    if producer == "panel_valley_density":
        return [P.density(w["clone"], w["allele"], "100kb", m)
                for m in ("antivalley", "allcoverage", "valley")]
    if producer == "panel_valley_locus":
        return [P.valley_states("chrX", w["track"])]
    if producer == "panel_valley_sizes":
        return [P.valley_sizes()]
    if producer == "panel_scatter_delta":
        roi = w.get("roi", "full")
        return ([P.pileup_score("Xa_vs_Xi", roi),
                 P.saddle_strength_selected("refined_merged", roi),
                 P.ml_compartmentalization_final()]
                + [P.allelic_ratio_stats(loc, roi)
                   for loc in ("Mecp2", "Jarid")])
    if producer == "panel_loops_comps_grid":
        locus, roi = w["locus"], w["roi"]
        comparison, kind = w["comparison"], w["kind"]
        mark = w.get("mark")
        out = [P.allelic_ratio_stats(locus, roi)]
        if comparison == "Xa_vs_Xi":
            coolers = [("{}_{}_WT".format(locus, c),
                        "{}_{}_{}_WT_{}_{}".format(a, locus, c, t, al),
                        "{}_{}_WT_{}_{}".format(locus, c, t, al))
                       for c in ("E6", "B1", "C5", "CL30", "JTG")
                       for a in ("Xa", "Xi")
                       for al, t in (("Xa", "G1"), ("Xi", "G2"))]
        else:
            degron_clones = (("E6A7", "CTCF"), ("F3", "CTCF"),
                             ("C5C10", "CTCF"), ("B1621", "Rad21"))
            coolers = [("{}_{}_G1_Xa".format(locus, c),
                        "{}_{}_{}-{}_G1_Xa".format(locus, c, m, cond),
                        "{}_{}_{}-{}_G1_Xa".format(locus, c, m, cond))
                       for c, m in degron_clones
                       if not mark or m == mark
                       for cond in ("NodTAG", "dTAG")]
        if kind in ("loops", "both"):
            out.append(P.pileup_score(comparison, roi))
            out += [P.pileup(comparison, roi, exp, sample)
                    for exp, sample, _ in coolers]
        if kind in ("comps", "both"):
            scope = ("refined_merged" if comparison == "Xa_vs_Xi"
                     else "refined")
            out.append(P.saddle_strength_selected(scope, roi))
            out += [P.saddle(roi, name, "E1") for _, _, name in coolers]
        return out
    if producer in ("panel_metaloci_violin", "panel_metaloci_heatmap",
                    "panel_metaloci_scatter"):
        return [P.ml_compartmentalization("wt", "full")]
    if producer == "panel_metaloci_composite":
        import yaml as _yaml
        with open(os.path.join(REPO, "resources", "fixtures",
                               "fig6_layout.yaml")) as fh:
            layout = _yaml.safe_load(fh)
        locus, slug = w["locus"], w["locus_slug"]
        out = [os.path.join(REPO, "resources", "fixtures",
                            "fig6_layout.yaml"),
               P.ml_compartmentalization_final(),
               P.ml_compartmentalization("consensus", "full")]
        for signal in layout["rows"]:
            for column in layout["columns"].get(slug, []):
                if column == "violin":
                    out.append(P.figures(
                        "extras", "metaloci_violin_{}_{}".format(locus, signal),
                        "pdf"))
                elif column == "Xa_consensus":
                    out.append(P.figures(
                        "extras", "metaloci_gaudi_consensus_{}_NodTAG-or-WT_Xa"
                        "_{}".format(locus, signal), "pdf"))
                else:
                    out.append(P.figures(
                        "extras", "metaloci_gaudi_wt_{}_{}_WT_G1_Xi_{}".format(
                            locus, column, signal), "pdf"))
        return out
    if producer == "panel_stackup_composite":
        return [P.stackup(w["track"], w["signal"], w["variant"])]
    if producer == "panel_valley_xa_xi_venn":
        return [P.valley_xa_xi(w["track"])]
    if producer == "panel_eigenvector":
        return [P.eigs(w["roi"], w["name"])]
    if producer == "panel_frip":
        return [P.qc("frip.tsv")]
    if producer == "panel_peak_venn":
        return [P.qc("peak_overlap.tsv")]
    if producer == "panel_correlation":
        return [P.qc("correlation", "{}.npz".format(w["normgroup"]))]
    if producer == "panel_fragment_sizes":
        import glob as _glob
        return sorted(_glob.glob(P.qc("fragment_sizes", "*.tsv")))
    if producer == "panel_signal_enrichment":
        return [P.qc("enrichment", "{}_{}.npz".format(w["mark"], w["anchor"]))]
    if producer == "panel_ma_plot":
        return [P.scalefactors(w["normgroup"], "counts"),
                P.scalefactors(w["normgroup"], "factors")]
    if producer == "panel_loop_overlay":
        name = w["name"]
        return [P.loops_raw(name), P.loops_refined(name)]
    if producer == "panel_saddle":
        return [P.saddle(w["roi"], w["name"], w["eig"]),
                P.saddle_values(w["roi"], w["name"])]
    if producer == "panel_pileup":
        return [P.pileup(w["comparison"], w["roi"], w["exp"], w["sample"])]
    return []


def run_one(panel, P, config, out_dir, ext="svg"):
    scripts_dir = os.path.join(REPO, "workflow", "scripts", "py", "panels")
    script = SCRIPTS.get(panel.producer)
    if not script:
        return "skip", "no script mapping for " + panel.producer

    inputs = [str(p) for p in inputs_for(P, panel, config)]
    if not inputs:
        return "skip", "no stub inputs defined for " + panel.producer
    missing = [p for p in inputs if not os.path.exists(p)]
    if missing:
        return "skip", "{} stub input(s) absent, e.g. {}".format(
            len(missing), os.path.basename(missing[0]))

    out_path = os.path.join(out_dir, "figures", panel.group,
                            "{}.{}".format(panel.id, ext))
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    log_path = os.path.join(out_dir, "logs", "{}.log".format(panel.id))
    os.makedirs(os.path.dirname(log_path), exist_ok=True)

    spec = dict(panel.wildcards or {})
    spec.update({
        "panel_id": panel.id, "stem": panel.id, "series": panel.series,
        "paper_ref": list(panel.paper_ref), "group": panel.group,
        "producer": panel.producer, "confidence": panel.confidence,
        "candidate": None, "candidates": list(panel.candidates),
        "notes": panel.notes,
        "libdir": os.path.join(REPO, "workflow"), "dpi": 300,
        "colours": os.path.join(REPO, "resources", "fixtures",
                                "figure_colours.yaml"),
        "display_names": {"Jarid": "Kdm5c", "Mecp2": "Mecp2"},
        "legacy": {"stackup_violin_test": "mannwhitneyu"},
    })
    params = {"spec": spec,
              "windows": {"Mecp2": [73315000, 74475000],
                          "Jarid": [150400000, 152350000]},
              "saddle_cfg": {"n_bins": 38, "strength_extent": 8},
              "pileup_cfg": {"flank": 100000},
              "display": os.path.join(REPO, "resources", "fixtures",
                                      "coolbox_display.yaml")}

    fake = Fake(inputs, out_path, params, log_path, config)
    path = os.path.join(scripts_dir, script)
    scope = {"snakemake": fake, "__file__": path, "__name__": "__main__"}
    saved = list(sys.path)
    try:
        with open(path) as fh:
            code = compile(fh.read(), path, "exec")
        exec(code, scope)
    except SystemExit as exc:
        if exc.code not in (0, None):
            return "fail", "SystemExit({})".format(exc.code)
    except Exception:
        return "fail", traceback.format_exc().strip().splitlines()[-1]
    finally:
        sys.path[:] = saved

    if not os.path.exists(out_path):
        return "fail", "no output file"
    size = os.path.getsize(out_path)
    if size == 0:
        return "fail", "output is empty"
    return "ok", "{} ({:,} bytes)".format(out_path, size)


def run_summaries(P, config, registry):
    """Exercise figure_manifest, paper_numbers and report end to end.

    These three are the stage's deliverables and none of them is a panel, so
    `--all` would never touch them. `paper_numbers` is run non-strict here
    because the stubs are synthetic: a MISMATCH is expected wherever a stub
    value was not deliberately set to the published one, and the point of the
    run is to prove the comparison machinery fires, not to check the numbers.
    """
    scripts_dir = os.path.join(REPO, "workflow", "scripts", "py", "panels")
    figures = []
    for panel in registry:
        for stem in panel.stems():
            for ext in panel.formats:
                path = P.figures(panel.group, stem, ext)
                if os.path.exists(path):
                    figures.append(path)

    degron_mark = {"E6A7": "CTCF", "F3": "CTCF", "B1621": "Rad21"}
    wt = ["H3K27me3_{}_WT_Xi".format(c)
          for c in ("E6", "C5", "B1", "JTG", "CL30")]
    nod = ["H3K27me3_{}_{}-NodTAG_Xi".format(c, m)
           for c, m in degron_mark.items()]
    dta = ["H3K27me3_{}_{}-dTAG_Xi".format(c, m)
           for c, m in degron_mark.items()]

    jobs = [
        ("figure_manifest.py", {"tsv": P.manifest()},
         {"registry": registry.to_records(), "settings": registry.settings,
          "spec": {"libdir": os.path.join(REPO, "workflow"),
                   "results": P.results}}, figures),
        ("paper_numbers.py", {"tsv": P.paper_numbers()},
         {"registry": os.path.join(REPO, "config", "panels.yaml"),
          "spec": {
              "libdir": os.path.join(REPO, "workflow"),
              "strict": False,
              "wt_tracks": wt, "degron_nodtag_tracks": nod,
              "degron_dtag_tracks": dta,
              "dtag_clones": list(degron_mark),
              "gene_content": {t: P.gene_content(t) for t in wt + nod + dta},
              "boundary_ctcf": {t: P.boundary_ctcf(t) for t in wt + nod},
              "valleys": {t: P.valleys("chrX", t) for t in wt + nod + dta},
              "valleys_filtered": {t: P.valleys_filtered("chrX", t)
                                   for t in wt + nod + dta},
              "valley_overlap": {c: P.valley_overlap(c, "escaping")
                                 for c in degron_mark},
              "pileup_scores": {
                  "Xa_vs_Xi": P.pileup_score("Xa_vs_Xi", "full"),
                  "dTAG_vs_NodTAG": P.pileup_score("dTAG_vs_NodTAG", "full")},
              "loop_stats": {},
              "loops_raw": {n: P.loops_raw(n) for n in
                            ("Jarid_Xa", "Jarid_Xi", "Mecp2_Xa", "Mecp2_Xi")},
              "loops_refined": {n: P.loops_refined(n) for n in
                                ("Jarid_Xa", "Jarid_Xi", "Mecp2_Xa",
                                 "Mecp2_Xi")},
              "saddle_strength": {
                  "Xa_vs_Xi": P.saddle_strength_selected("refined_merged",
                                                         "full"),
                  "dTAG_vs_NodTAG": P.saddle_strength_selected("refined",
                                                               "full")},
              "compartmentalization": P.ml_compartmentalization_final(),
          }}, []),
        ("report.py", {"html": P.report()},
         {"registry": registry.to_records(), "settings": registry.settings,
          "spec": {"libdir": os.path.join(REPO, "workflow"),
                   "results": P.results, "pipeline_version": "stub",
                   "max_inline_bytes": 6000000,
                   "eigenvector_orientation": P.work(
                       "features", "compartments",
                       "eigenvector_orientation.tsv"),
                   "provenance": P.provenance()}},
         figures),
    ]

    results = []
    for script, outputs, params, inputs in jobs:
        log_path = os.path.join(P.results, "logs",
                                script.replace(".py", ".log"))
        os.makedirs(os.path.dirname(log_path), exist_ok=True)
        fake = Fake(inputs, list(outputs.values())[0], params, log_path, config)
        fake.output = _Output(list(outputs.values()))
        for key, value in outputs.items():
            setattr(fake.output, key, value)
        fake.input = _Input(inputs)
        for key, value in (("manifest", P.manifest()),
                           ("numbers", P.paper_numbers()),
                           ("panels", inputs)):
            setattr(fake.input, key, value)

        path = os.path.join(scripts_dir, script)
        scope = {"snakemake": fake, "__file__": path, "__name__": "__main__"}
        saved = list(sys.path)
        try:
            with open(path) as fh:
                exec(compile(fh.read(), path, "exec"), scope)
            out = list(outputs.values())[0]
            size = os.path.getsize(out) if os.path.exists(out) else 0
            status = "ok" if size else "fail"
            results.append((status, script, "{} ({:,} bytes)".format(out, size)))
        except SystemExit as exc:
            results.append(("fail" if exc.code else "ok", script, str(exc)))
        except Exception:
            results.append(("fail", script, traceback.format_exc()))
        finally:
            sys.path[:] = saved

    for status, script, detail in results:
        print("{:5s} {:44s} {}".format(status, script, detail))
    return results


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--data", required=True, help="the stub data_dir")
    ap.add_argument("--panel", default="", help="one panel id")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--summaries", action="store_true",
                    help="also run figure_manifest, paper_numbers and report")
    ap.add_argument("--ext", default="svg", choices=("svg", "pdf", "png"))
    args = ap.parse_args()

    with open(os.path.join(REPO, "config", "config.yaml")) as fh:
        config = yaml.safe_load(fh)
    config["paths"]["data_dir"] = os.path.abspath(args.data)
    for key in ("resources_dir", "results_dir", "log_dir", "tmpdir"):
        config["paths"][key] = "{data_dir}/" + key.split("_")[0]
    config["paths"]["results_dir"] = "{data_dir}/results"
    P = Paths(config)

    registry = _panels.load(os.path.join(REPO, "config", "panels.yaml"),
                            config,
                            selectors=_stub_selectors())

    targets = [p for p in registry
               if not args.panel or p.id == args.panel]
    if args.panel and not targets:
        raise SystemExit("no panel {!r} in the registry".format(args.panel))
    if not args.all and not args.panel and not args.summaries:
        raise SystemExit("pass --panel <id>, --all or --summaries")
    if args.summaries and not (args.all or args.panel):
        return 1 if any(s == "fail"
                        for s, _, _ in run_summaries(P, config, registry)) else 0

    results = []
    for panel in targets:
        status, detail = run_one(panel, P, config, P.results, ext=args.ext)
        results.append((status, panel.id, detail))
        if status != "skip" or args.panel:
            print("{:5s} {:44s} {}".format(status, panel.id, detail))

    counts = {}
    for status, _, _ in results:
        counts[status] = counts.get(status, 0) + 1
    print("\n" + ", ".join("{} {}".format(v, k)
                           for k, v in sorted(counts.items())))
    if args.summaries:
        print("")
        for status, _, _ in run_summaries(P, config, registry):
            if status == "fail":
                counts["fail"] = counts.get("fail", 0) + 1
    return 1 if counts.get("fail") else 0


def _stub_selectors():
    """Selectors that do not need a sample sheet."""
    degron = ["E6A7", "F3", "B1621"]
    wt_tracks = ["H3K27me3_{}_WT_Xi".format(c)
                 for c in ("E6", "C5", "B1", "JTG", "CL30")]
    return {
        "wt_clones": ["E6", "C5", "B1", "JTG", "CL30"],
        "degron_clones": degron,
        "dtag_pair_clones": degron,
        "degron_nodtag_clones": degron,
        "me3_tracks": wt_tracks,
        "me3_xi_tracks": wt_tracks,
        "loci": ["Mecp2", "Jarid"],
        "rois": ["full", "escape", "non_escape"],
        "normgroups": ["WT_H3K27me3"],
        "coolers": ["Mecp2_E6_WT_G1_Xa"],
        "merged_loop_coolers": ["Jarid_Xa", "Jarid_Xi", "Mecp2_Xa", "Mecp2_Xi"],
        "merged_comp_coolers": ["Mecp2_NodTAG-or-WT_Xa"],
        "enrichment_rows": [{"mark": "CTCF", "anchor": "motif",
                             "_paper_ref": ["Fig 1e"], "_group": "figure_01"}],
        "fig6_composites": [
            {"locus": "Jarid", "locus_slug": "kdm5c", "run": "wt"},
            {"locus": "Mecp2", "locus_slug": "mecp2", "run": "wt"}],
        "eig_rows": [{"roi": "full", "name": "Mecp2_E6_WT_G1_Xa"}],
        "saddle_rows": [{"roi": "full", "name": "Mecp2_E6_WT_G2_Xi",
                         "eig": "E1"}],
        "pileup_rows": [{"comparison": "Xa_vs_Xi", "roi": "full",
                         "exp": "Mecp2_E6_WT",
                         "sample": "Xa_Mecp2_E6_WT_G1_Xa"}],
        "stackup_rows": [{"track": wt_tracks[0], "signal": "H3K27me3",
                          "variant": "all"}],
        "density_rows": [{"clone": c, "allele": "Xi"} for c in degron],
        "metaloci_rows": [{"run": "wt", "dataset": "Mecp2_E6_WT_G1_Xi",
                           "signal": "H3K27me3"}],
        "metaloci_kk_rows": [{"run": "wt", "dataset": "Mecp2_E6_WT_G1_Xi"}],
        "metaloci_signal_rows": [{"locus": l, "signal": s}
                                 for l in ("Mecp2", "Jarid")
                                 for s in ("H3K27me3", "H3K27ac", "CTCF",
                                           "AcMe3", "RNA-Seq", "Rad21")],
    }


if __name__ == "__main__":
    sys.exit(main())
