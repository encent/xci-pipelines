"""`figure_manifest.tsv` — panel to paper panel, one row per emitted file.

Deliverable D-07. The original had no such file: reconstructing which SVG
became which published panel took the archaeology weeks, and two panels came
from an interactive session whose notebook was never saved.

    panel_id  paper_ref  file  status  n_items  confidence  notes

EVERY registry entry gets a row, built or not, because the question this file
has to answer is not "what did we draw" but "why is Fig 5c not in my results".
`status` is one of:

    ok        the file exists and is non-empty
    empty     the file exists and is zero bytes -- a crashed panel
    missing   requested but absent; look at the panel's log
    skipped   `requires_signals` not satisfied (a shortened metaloci.signals)
    blocked   nothing can produce its inputs yet, with the reason in `notes`
    disabled  switched off by config (qc.enabled, figures.panels)

`n_items` is what the panel actually drew -- clones, boundaries, loops. Panels
drop a one-line `.n_items` sidecar next to their output and this rule collects
it, rather than re-reading every input to recount.
"""

import os
import sys

spec = dict(snakemake.params.spec)
sys.path.insert(0, spec["libdir"])
sys.path.insert(0, os.path.join(spec["libdir"], "scripts", "py", "panels"))

from lib import panels as _panels                              # noqa: E402
from lib.paths import Paths                                    # noqa: E402

COLUMNS = _panels.MANIFEST_COLUMNS


def _sidecars(results):
    """{stem: (n_items, {key: value})} from the `.n_items` files."""
    out = {}
    figures = os.path.join(results, "figures")
    for root, _dirs, files in os.walk(figures):
        for name in files:
            if not name.endswith(".n_items"):
                continue
            stem = name[: -len(".n_items")]
            stem = os.path.splitext(stem)[0]
            try:
                with open(os.path.join(root, name)) as fh:
                    lines = [l.rstrip("\n") for l in fh if l.strip()]
            except OSError:
                continue
            if not lines:
                continue
            extras = dict(
                l.split("\t", 1) for l in lines[1:] if "\t" in l)
            out[stem] = (lines[0], extras)
    return out


log_path = snakemake.log[0]
os.makedirs(os.path.dirname(log_path) or ".", exist_ok=True)
with open(log_path, "w") as log:
    def say(msg=""):
        log.write(str(msg) + "\n")
        log.flush()

    P = Paths(snakemake.config)
    registry = _panels.from_records(snakemake.params.registry,
                                    snakemake.params.get("settings"))
    sidecars = _sidecars(P.results)
    say("registry: {} panels, {} sidecar counts".format(
        len(registry), len(sidecars)))

    n_items = {stem: value for stem, (value, _extra) in sidecars.items()}
    rows = list(_panels.manifest_rows(registry, P.figures, n_items=n_items))

    out_path = str(snakemake.output.tsv)
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(out_path, "w") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for row in rows:
            fh.write("\t".join(str(row.get(c, "")) for c in COLUMNS) + "\n")

    counts = {}
    for row in rows:
        counts[row["status"]] = counts.get(row["status"], 0) + 1
    say("{} rows: {}".format(
        len(rows), ", ".join("{} {}".format(v, k)
                             for k, v in sorted(counts.items()))))

    paper = [r for r in rows if r["paper_ref"] != "-"]
    bad = [r for r in paper if r["status"] not in ("ok",)]
    say("{} paper-panel rows, {} not ok".format(len(paper), len(bad)))
    for row in bad:
        say("  {:10s} {:40s} {:12s} {}".format(
            row["status"], row["panel_id"], row["paper_ref"], row["notes"][:80]))
    say("wrote {}".format(out_path))
