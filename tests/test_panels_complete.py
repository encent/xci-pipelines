#!/usr/bin/env python3
"""Check ``config/panels.yaml`` against the 25 target paper panels.

The registry is the only place that says which figure belongs to which
published panel, so it is the only place a panel can silently go missing.
This test is the guard, and it is deliberately cheap: YAML plus a text scan of
``90_visualise.smk``, no Snakemake, no config, no data (plan section 9.7).

What "unique" means here
------------------------
A paper reference does NOT map one-to-one onto a file. Fig 6 is one published
panel produced as two files (Mecp2 + Kdm5c). EFig 2e/2f/2g are three published
panels produced as *one* figure, five times over (once per WT clone). Fig 4f,
Fig 4g and EFig 6b are likewise three published panels of one figure.

So the registry carries a ``series`` field, and the invariant is:

    every one of the 25 target references is claimed by exactly ONE series.

Two series claiming "Fig 3a" is the failure this catches -- that is how a panel
gets drawn twice from two different code paths and nobody notices which one
went into the manuscript.

Run::

    python tests/test_panels_complete.py
    pytest tests/test_panels_complete.py        # also works
"""

from __future__ import annotations

import os
import re
import sys

import yaml

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REGISTRY = os.path.join(REPO, "config", "panels.yaml")
VIZ = os.path.join(REPO, "workflow", "rules", "90_visualise.smk")

#: The 25 target panels, from the design notes section 6.4 / the archaeology notes F.1.
#: `figure.pdf` (2026-07-25) is authoritative for the lettering, NOT
#: scripts/figure_legends.txt, which uses an older scheme and inverts 3a/3b.
TARGET_PANELS = [
    "Fig 2g", "Fig 2h",
    "Fig 3a", "Fig 3b", "Fig 3c",
    "Fig 4d", "Fig 4f", "Fig 4g",
    "Fig 5b", "Fig 5c", "Fig 5g", "Fig 5h",
    "Fig 6",
    "EFig 2e", "EFig 2f", "EFig 2g",
    "EFig 3a", "EFig 3b",
    "EFig 6a", "EFig 6b", "EFig 6c",
    "EFig 8c", "EFig 8d",
    "EFig 10c", "EFig 10d",
]

#: Published but outside the 25-panel target list (the design review P-3).
EXTRA_PAPER_PANELS = ["Fig 4e"]

REQUIRED_FIELDS = ("id", "group", "producer", "formats", "confidence")


def load_registry():
    with open(REGISTRY) as fh:
        doc = yaml.safe_load(fh)
    if not isinstance(doc, dict) or "panels" not in doc:
        raise AssertionError("{} has no `panels:` list".format(REGISTRY))
    return doc


def _merged(doc):
    defaults = doc.get("defaults") or {}
    for entry in doc["panels"]:
        merged = dict(defaults)
        merged.update(entry)
        yield merged


def _series(entry):
    return entry.get("series") or entry["id"]


# ---------------------------------------------------------------------------
# checks
# ---------------------------------------------------------------------------
def test_registry_parses():
    doc = load_registry()
    assert doc["panels"], "the registry is empty"


def test_required_fields_present():
    problems = []
    for entry in _merged(load_registry()):
        for field in REQUIRED_FIELDS:
            if field not in entry or entry[field] in (None, ""):
                problems.append("{}: missing `{}`".format(entry.get("id"), field))
        if entry.get("formats") and not isinstance(entry["formats"], list):
            problems.append("{}: `formats` must be a list".format(entry["id"]))
        if entry.get("paper_ref") is not None and not isinstance(
            entry["paper_ref"], list
        ):
            problems.append("{}: `paper_ref` must be a list".format(entry["id"]))
    assert not problems, "registry field problems:\n  " + "\n  ".join(problems)


def test_template_ids_unique():
    seen = {}
    dupes = []
    for entry in _merged(load_registry()):
        if entry["id"] in seen:
            dupes.append(entry["id"])
        seen[entry["id"]] = entry
    assert not dupes, "duplicate panel ids: {}".format(sorted(set(dupes)))


def test_all_25_target_panels_present():
    claims = {}
    for entry in _merged(load_registry()):
        for ref in entry.get("paper_ref") or []:
            claims.setdefault(ref, set()).add(_series(entry))

    missing = [ref for ref in TARGET_PANELS if ref not in claims]
    assert not missing, (
        "{} of the 25 target panels have no producer in the registry: {}".format(
            len(missing), missing
        )
    )

    n = len([r for r in claims if r in TARGET_PANELS])
    assert n == 25, "expected 25 target refs, registry claims {}".format(n)


def test_each_paper_ref_has_exactly_one_series():
    claims = {}
    for entry in _merged(load_registry()):
        for ref in entry.get("paper_ref") or []:
            claims.setdefault(ref, set()).add(_series(entry))

    ambiguous = {ref: sorted(s) for ref, s in claims.items() if len(s) > 1}
    assert not ambiguous, (
        "paper references claimed by more than one series -- the same published "
        "panel would be drawn twice from two code paths:\n  {}".format(ambiguous)
    )


def test_fig_4e_is_a_paper_panel():
    """the design review P-3: Fig 4e is published; it must not be filed as an extra."""
    claims = set()
    for entry in _merged(load_registry()):
        claims.update(entry.get("paper_ref") or [])
    for ref in EXTRA_PAPER_PANELS:
        assert ref in claims, (
            "{} is a published panel (Figure 4) and must carry a real "
            "paper_ref so it appears in figure_manifest.tsv".format(ref)
        )


def test_no_unknown_paper_refs():
    known = set(TARGET_PANELS) | set(EXTRA_PAPER_PANELS)
    unknown = set()
    for entry in _merged(load_registry()):
        for ref in entry.get("paper_ref") or []:
            if ref not in known:
                unknown.add(ref)
    assert not unknown, (
        "paper references that are neither one of the 25 nor a known extra: "
        "{}. Either it is a typo or TARGET_PANELS needs updating.".format(
            sorted(unknown)
        )
    )


def test_medium_confidence_panels_list_candidates():
    problems = []
    for entry in _merged(load_registry()):
        if (
            entry.get("confidence") == "medium"
            and not entry.get("candidates")
            and not entry.get("candidates_from_expand")
        ):
            problems.append(entry["id"])
    assert not problems, (
        "medium-confidence panels must list their `candidates:` so the pipeline "
        "can render all of them for a human to pick from, or set "
        "`candidates_from_expand: true` when the expand axis already is the "
        "candidate axis: {}".format(problems)
    )


def test_every_producer_rule_exists():
    if not os.path.exists(VIZ):
        raise AssertionError("{} does not exist".format(VIZ))
    with open(VIZ) as fh:
        text = fh.read()
    defined = set(re.findall(r"^rule\s+([A-Za-z_]\w*)\s*:", text, re.M))

    wanted = {e["producer"] for e in _merged(load_registry())}
    missing = sorted(wanted - defined)
    assert not missing, (
        "panels.yaml names producer rules that 90_visualise.smk does not "
        "define: {}".format(missing)
    )


def test_producers_are_panel_rules():
    bad = [
        e["id"] for e in _merged(load_registry())
        if not e["producer"].startswith("panel_")
    ]
    assert not bad, "producers must be `panel_*` rules; offenders: {}".format(bad)


def test_gene_content_panels_document_the_filtered_set():
    """TRAP P-2 -- the single most likely silent divergence in the pipeline."""
    hits = [
        e for e in _merged(load_registry())
        if e["producer"] == "panel_valley_gene_content"
    ]
    assert hits, "no panel_valley_gene_content entry (Fig 2g / Fig 4d)"
    for entry in hits:
        note = (entry.get("notes") or "").lower()
        assert "filter" in note, (
            "{}: the notes must say the panel uses the GENE-FILTERED valley "
            "set. Fed the raw HMM calls it draws 377/304/261/302/324 instead "
            "of the published 375/301/258/300/321 and looks correct.".format(
                entry["id"]
            )
        )


def test_boundary_analysis_is_not_fed_from_stackups():
    """TRAP R-2 / P-1 -- same inputs, different arithmetic."""
    hits = [
        e for e in _merged(load_registry())
        if e["producer"] == "panel_boundary_analysis"
    ]
    assert hits, "no panel_boundary_analysis entry (EFig 2e / 2f / 2g)"
    for entry in hits:
        inputs = " ".join(entry.get("inputs") or [])
        assert "boundary_profile" in inputs, (
            "{}: must be fed from P.boundary_profile".format(entry["id"])
        )
        assert "stackups/" not in inputs, (
            "{}: P.stackup subtracts a 100-shift random background; routing "
            "EFig 2e/2f there yields plausible curves that do not "
            "reproduce".format(entry["id"])
        )


def test_fig6_does_not_hardcode_the_grid():
    """TRAP F-6 -- the Fig 6 column order is not derivable."""
    hits = [
        e for e in _merged(load_registry())
        if e["producer"] == "panel_metaloci_composite"
    ]
    assert hits, "no panel_metaloci_composite entry (Fig 6)"
    for entry in hits:
        blob = " ".join(entry.get("inputs") or []) + " " + (entry.get("notes") or "")
        assert "fig6_layout.yaml" in blob, (
            "{}: the row/column order must come from "
            "resources/fixtures/fig6_layout.yaml (fixture F-6), never from "
            "code".format(entry["id"])
        )


ALL_TESTS = [v for k, v in sorted(globals().items()) if k.startswith("test_")]


def main():
    failures = []
    for fn in ALL_TESTS:
        try:
            fn()
        except AssertionError as exc:
            failures.append((fn.__name__, str(exc)))
        except Exception as exc:                              # pragma: no cover
            failures.append((fn.__name__, "{}: {}".format(type(exc).__name__, exc)))

    for name, msg in failures:
        print("FAIL  {}\n      {}\n".format(name, msg.replace("\n", "\n      ")))

    print("{}/{} checks passed".format(len(ALL_TESTS) - len(failures), len(ALL_TESTS)))
    if failures:
        return 1

    doc = load_registry()
    series = {_series(e) for e in _merged(doc)}
    print("registry: {} entries, {} series, 25/25 target panels + {}".format(
        len(doc["panels"]), len(series), ", ".join(EXTRA_PAPER_PANELS)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
