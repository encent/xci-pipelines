"""Load ``config/panels.yaml`` and turn it into concrete output paths.

``90_visualise.smk`` reads the registry and nothing else to decide what to
draw. This module is the reader. It has no Snakemake dependency so it can be
imported by tests and by the manifest/report scripts.

Three things happen here, in order:

1. **Expansion.** ``id: "boundary_analysis_{clone}"`` plus
   ``expand: {clone: "@wt_clones"}`` becomes five panels. Selectors (``@name``)
   are resolved by callables supplied by the rule module, which is where the
   sample sheet lives — the registry stays free of clone lists.

2. **Gating.** A panel is dropped, with a reason recorded, when

   * ``figures.panels`` in config excludes it (``all`` | ``paper`` | a list of
     paper references);
   * a ``requires_config`` key is falsy (``qc.enabled``, ``coolbox.enabled``);
   * a ``requires_signals`` entry is missing from ``metaloci.signals`` — so
     shortening that list SKIPS the METALoci panels instead of failing the DAG
     half an hour in;
   * none of its inputs can be produced yet. Modules land at different times;
     a panel whose upstream rule does not exist and whose input is not on disk
     is *blocked*, not an error. As each module lands, its panels light up with
     no edit here. Set ``settings.skip_panels_without_producers: false`` to
     turn a missing upstream back into a hard DAG failure.

3. **Path construction.** ``results/figures/{group}/{id}.{fmt}``, via
   ``paths.Paths.figures`` — never string-formatted here.

Every dropped panel still gets a row in ``figure_manifest.tsv`` with its
reason, which is the point: "why is Fig 5c not in my results" must be
answerable from one file.
"""

from __future__ import annotations

import itertools
import os
from typing import Callable, Dict, Iterable, List, Optional

import yaml

#: Panels with these statuses are not requested from Snakemake.
NOT_BUILT = ("skipped", "blocked", "disabled")


class PanelError(ValueError):
    """Raised when the registry is malformed. Always names the panel."""


class Panel(object):
    """One concrete figure: one id, one group, one producer, N formats."""

    __slots__ = (
        "id", "series", "paper_ref", "group", "producer", "wildcards",
        "formats", "confidence", "candidates", "notes", "inputs",
        "requires_signals", "requires_config", "status", "reason",
    )

    def __init__(self, **kw):
        for slot in self.__slots__:
            setattr(self, slot, kw.get(slot))
        self.status = kw.get("status") or "enabled"
        self.reason = kw.get("reason") or ""

    # -- output paths ------------------------------------------------------
    def stems(self) -> List[str]:
        """File stems this panel produces: the primary plus any candidates."""
        out = [self.id]
        for cand in self.candidates or []:
            out.append("{}__{}".format(self.id, cand))
        return out

    def files(self, figures: Callable[[str, str, str], str]) -> List[str]:
        return [
            figures(self.group, stem, ext)
            for stem in self.stems()
            for ext in self.formats
        ]

    @property
    def built(self) -> bool:
        return self.status not in NOT_BUILT

    def __repr__(self):                                        # pragma: no cover
        return "<Panel {} [{}] {}>".format(self.id, self.status, self.producer)


class Registry(object):
    def __init__(self, doc: dict, panels: List[Panel]):
        self.doc = doc
        self.panels = panels
        self.settings = doc.get("settings") or {}

    # -- queries -----------------------------------------------------------
    def __iter__(self):
        return iter(self.panels)

    def __len__(self):
        return len(self.panels)

    @property
    def enabled(self) -> List[Panel]:
        return [p for p in self.panels if p.built]

    def by_producer(self, producer: str, only_enabled: bool = True) -> List[Panel]:
        src = self.enabled if only_enabled else self.panels
        return [p for p in src if p.producer == producer]

    def producers(self) -> List[str]:
        return sorted({p.producer for p in self.panels})

    def get(self, panel_id: str) -> Panel:
        for p in self.panels:
            if p.id == panel_id:
                return p
        raise KeyError("no panel {!r} in the registry".format(panel_id))

    def by_stem(self, stem: str) -> Panel:
        """Resolve an output stem, candidate suffix included, back to its panel."""
        for p in self.panels:
            if stem in p.stems():
                return p
        raise KeyError("no panel produces the stem {!r}".format(stem))

    def stem_alternation(self, producer: str) -> str:
        """A regex alternation of every stem one producer owns.

        This is what keeps ~40 rules that all write
        ``results/figures/{group}/{name}.{ext}`` unambiguous: each rule
        constrains ``{name}`` to exactly its own registry entries, so no path
        can match two rules and Snakemake never has to guess.
        """
        import re
        stems = sorted(
            {s for p in self.by_producer(producer, only_enabled=False)
             for s in p.stems()},
            key=len, reverse=True,
        )
        if not stems:
            # An alternation that cannot match anything. A rule with no panels
            # must be unreachable rather than match everything.
            return r"(?!x)x"
        return "|".join(re.escape(s) for s in stems)

    def files(self, figures) -> List[str]:
        out = []
        for panel in self.enabled:
            out.extend(panel.files(figures))
        return out

    def series_of(self, panel_id: str) -> str:
        p = self.get(panel_id)
        return p.series or p.id

    # -- serialisation -----------------------------------------------------
    def to_records(self) -> List[dict]:
        """Plain dicts, for handing the EXPANDED registry to a script.

        ``figure_manifest`` and ``report`` must describe the registry the DAG
        actually used -- including which panels were gated out and why. They
        cannot re-derive that by calling ``load()`` again: expansion needs the
        selectors, which live in the rule module, and the gating needs the rule
        graph, which no longer exists inside a running job. So the rule passes
        this through ``params`` and the script rebuilds from it.
        """
        return [
            {slot: getattr(panel, slot) for slot in Panel.__slots__}
            for panel in self.panels
        ]


def from_records(records, settings=None) -> Registry:
    """Rebuild a Registry from :meth:`Registry.to_records` output."""
    panels = [Panel(**dict(record)) for record in records]
    return Registry({"settings": dict(settings or {})}, panels)


# ---------------------------------------------------------------------------
# loading
# ---------------------------------------------------------------------------
def load(
    path: str,
    config: dict,
    selectors: Optional[Dict[str, Callable]] = None,
    is_available: Optional[Callable[[str], bool]] = None,
    input_paths: Optional[Callable[[Panel], Iterable[str]]] = None,
    panel_output_stem: Optional[Callable[[str], Optional[str]]] = None,
) -> Registry:
    """Read the registry, expand it, gate it.

    Parameters
    ----------
    selectors
        ``{"wt_clones": callable_or_list}`` — resolves ``"@wt_clones"``.
    is_available
        ``path -> bool``: is there a rule that can make this file, or is it
        already on disk. Used for the "module not built yet" gate.
    input_paths
        ``panel -> [paths]``: the panel's upstream inputs. Supplied by the rule
        module, since that is where ``paths.py`` and the sample sheet are.
    panel_output_stem
        ``path -> stem or None``: recognises a path that is ANOTHER PANEL's
        output. Fig 6 is a montage of already-drawn Gaudi PDFs, so the
        visualisation stage has one internal dependency, and it cannot be
        resolved by asking the rule graph: ``load()`` runs while this module is
        still being parsed, so its own panel rules are not registered yet. Such
        inputs are resolved against the registry instead, to a fixed point.
    """
    with open(path) as fh:
        doc = yaml.safe_load(fh) or {}
    if "panels" not in doc:
        raise PanelError("{}: no `panels:` list".format(path))

    selectors = selectors or {}
    defaults = doc.get("defaults") or {}
    settings = doc.get("settings") or {}

    panels: List[Panel] = []
    seen: Dict[str, str] = {}
    for raw in doc["panels"]:
        for panel in _expand(raw, defaults, selectors):
            if panel.id in seen:
                raise PanelError(
                    "duplicate panel id {!r} (from templates {!r} and {!r}). "
                    "Ids are file stems; two panels cannot share one."
                    .format(panel.id, seen[panel.id], raw["id"])
                )
            seen[panel.id] = raw["id"]
            panels.append(panel)

    _gate_config(panels, config)
    if is_available is not None and input_paths is not None:
        skip = settings.get("skip_panels_without_producers", True)
        _gate_inputs(panels, is_available, input_paths, skip, panel_output_stem)

    return Registry(doc, panels)


def _expand(raw: dict, defaults: dict, selectors: dict) -> List[Panel]:
    entry = dict(defaults)
    entry.update(raw)
    if "id" not in entry:
        raise PanelError("a registry entry has no `id`")
    for field in ("group", "producer", "formats", "confidence"):
        if not entry.get(field):
            raise PanelError("{}: missing `{}`".format(entry["id"], field))

    rows = _rows(entry, selectors)

    out = []
    for row in rows:
        wildcards = dict(entry.get("wildcards") or {})
        overrides = {}
        for key, value in row.items():
            if key.startswith("_"):
                overrides[key[1:]] = value
            else:
                wildcards[key] = value
        try:
            panel_id = entry["id"].format(**wildcards)
        except KeyError as exc:
            raise PanelError(
                "{}: the id template needs wildcard {} which no `expand` / "
                "`expand_rows` / `wildcards` provides".format(entry["id"], exc)
            )
        candidates = list(entry.get("candidates") or [])
        out.append(Panel(
            id=overrides.get("id", panel_id),
            series=overrides.get("series", entry.get("series") or panel_id),
            paper_ref=list(overrides.get("paper_ref", entry.get("paper_ref") or [])),
            group=overrides.get("group", entry["group"]),
            producer=entry["producer"],
            wildcards=wildcards,
            formats=list(entry["formats"]),
            confidence=entry.get("confidence", "high"),
            candidates=candidates,
            notes=_squash(entry.get("notes") or ""),
            inputs=list(entry.get("inputs") or []),
            requires_signals=list(entry.get("requires_signals") or []),
            requires_config=list(entry.get("requires_config") or []),
        ))
    return out


def _rows(entry: dict, selectors: dict) -> List[dict]:
    """Turn ``expand`` / ``expand_rows`` into a list of wildcard dicts."""
    if entry.get("expand_rows"):
        rows = _resolve(entry["expand_rows"], selectors, entry["id"])
        if rows and not isinstance(rows[0], dict):
            raise PanelError(
                "{}: `expand_rows` selector must return dicts".format(entry["id"])
            )
        return list(rows)

    spec = entry.get("expand")
    if not spec:
        return [{}]

    keys = list(spec)
    values = [_resolve(spec[k], selectors, entry["id"]) for k in keys]
    return [dict(zip(keys, combo)) for combo in itertools.product(*values)]


def _resolve(value, selectors: dict, panel_id: str):
    if isinstance(value, str) and value.startswith("@"):
        name = value[1:]
        if name not in selectors:
            raise PanelError(
                "{}: unknown selector @{}. Add it to SELECTORS in "
                "workflow/rules/90_visualise.smk.".format(panel_id, name)
            )
        got = selectors[name]
        return list(got() if callable(got) else got)
    if isinstance(value, (list, tuple)):
        return list(value)
    return [value]


def _squash(text: str) -> str:
    """One-line notes: TSV and HTML both dislike embedded newlines."""
    return " ".join(str(text).split())


# ---------------------------------------------------------------------------
# gates
# ---------------------------------------------------------------------------
def _dotted(config: dict, key: str, default=None):
    node = config
    for part in key.split("."):
        if not isinstance(node, dict) or part not in node:
            return default
        node = node[part]
    return node


def _gate_config(panels: List[Panel], config: dict) -> None:
    figures = config.get("figures") or {}
    wanted = figures.get("panels", "all")
    signals = set((config.get("metaloci") or {}).get("signals") or [])

    for panel in panels:
        for key in panel.requires_config:
            if not _dotted(config, key):
                panel.status = "disabled"
                panel.reason = "config.{} is off".format(key)
                break
        if panel.status != "enabled":
            continue

        missing = [s for s in panel.requires_signals if s not in signals]
        if missing:
            panel.status = "skipped"
            panel.reason = "metaloci.signals lacks {}".format(", ".join(missing))
            continue

        if wanted in (None, "all"):
            continue
        if wanted == "paper":
            if not panel.paper_ref:
                panel.status = "disabled"
                panel.reason = "figures.panels: paper"
            continue
        if isinstance(wanted, (list, tuple)):
            want = {str(w).strip() for w in wanted}
            hit = (
                set(panel.paper_ref) & want
                or panel.id in want
                or (panel.series or "") in want
            )
            if not hit:
                panel.status = "disabled"
                panel.reason = "not in figures.panels"


def _gate_inputs(panels, is_available, input_paths, skip_when_unavailable,
                 panel_output_stem=None):
    """Block panels whose upstream cannot be built, in two kinds of pass.

    Pass 1 resolves everything the rule graph knows about. Panel-to-panel
    inputs (Fig 6 montaging the Gaudi PDFs) cannot be settled there, because
    this runs while 90_visualise.smk is still being parsed and its own rules
    are not registered yet. They are resolved against the registry in pass 2,
    iterated to a fixed point so that a blocked Gaudi panel correctly blocks
    the composite that montages it rather than producing a DAG error later.
    """
    resolved = {}
    for panel in panels:
        if not panel.built:
            continue
        try:
            needed = list(input_paths(panel))
        except Exception as exc:
            panel.status = "blocked"
            panel.reason = "input resolution failed: {}: {}".format(
                type(exc).__name__, exc)
            continue

        internal, external = [], []
        for item in needed:
            stem = panel_output_stem(item) if panel_output_stem else None
            (internal if stem else external).append(stem or item)
        resolved[panel.id] = internal

        missing = [p for p in external if not is_available(p)]
        if missing and skip_when_unavailable:
            _block(panel, missing)

    if not panel_output_stem:
        return

    by_stem = {}
    for panel in panels:
        for stem in panel.stems():
            by_stem[stem] = panel

    for _ in range(len(panels) + 1):
        changed = False
        for panel in panels:
            if not panel.built:
                continue
            missing = []
            for stem in resolved.get(panel.id, ()):
                dependency = by_stem.get(stem)
                if dependency is None:
                    missing.append(stem + " (no such panel)")
                elif not dependency.built:
                    missing.append("{} ({})".format(stem, dependency.status))
            if missing and skip_when_unavailable:
                _block(panel, missing, verb="depends on unavailable panel")
                changed = True
        if not changed:
            break


def _block(panel, missing, verb="no rule yet produces"):
    shown = ", ".join(os.path.basename(str(m)) for m in missing[:3])
    if len(missing) > 3:
        shown += ", +{} more".format(len(missing) - 3)
    panel.status = "blocked"
    panel.reason = "{} {}".format(verb, shown)


# ---------------------------------------------------------------------------
# manifest support
# ---------------------------------------------------------------------------
MANIFEST_COLUMNS = [
    "panel_id", "paper_ref", "file", "status", "n_items", "confidence", "notes",
]


def manifest_rows(registry: Registry, figures, n_items=None):
    """Yield the ``figure_manifest.tsv`` rows for every panel, built or not."""
    n_items = n_items or {}
    for panel in registry.panels:
        stems = panel.stems()
        for stem in stems:
            for ext in panel.formats:
                path = figures(panel.group, stem, ext)
                if panel.built:
                    if os.path.exists(path) and os.path.getsize(path) > 0:
                        status = "ok"
                    elif os.path.exists(path):
                        status = "empty"
                    else:
                        status = "missing"
                else:
                    status = panel.status
                note = panel.notes
                if panel.reason:
                    note = "{} [{}]".format(panel.reason, note) if note \
                        else panel.reason
                if stem != panel.id:
                    note = "candidate {}; {}".format(
                        stem.split("__", 1)[1], note)
                yield {
                    "panel_id": stem,
                    "paper_ref": ";".join(panel.paper_ref) or "-",
                    "file": path,
                    "status": status,
                    "n_items": n_items.get(stem, n_items.get(panel.id, "")),
                    "confidence": panel.confidence,
                    "notes": note,
                }
