"""`REPORT.html` — one self-contained page with every panel inlined.

No external assets. SVG and PNG go in as `data:` URIs and the CSS is inline, so
the file can be emailed, copied to a laptop or opened from a USB stick and
still render. That is the point: the paper's figures currently live across two
projects and a dozen result directories, and a reader who wants to check them
should not have to mount a filesystem.

Sections follow the paper -- Figure 2, Figure 3, ... Extended Data -- with the
extras and QC last, each panel captioned with its paper reference, status,
n-value and the note from the registry.

THE EIGENVECTOR ORIENTATION TABLE IS SHOWN, NOT LISTED
------------------------------------------------------
`work/features/compartments/eigenvector_orientation.tsv` records every
eigenvector sign flip the pipeline applied. It is rendered near the top of the
report, in full, because of a property RD-3 pinned down: a GLOBAL sign flip
leaves the saddle corner score `(AA + BB) / (AB + BA)` **unchanged** — the
corners swap in pairs — while mirroring the plot. No numeric test on
`saddle_strength_selected_*.tsv` can catch it. `compartments*.json` fixes which
eigenvector was used but never its orientation, so a silent flip would corrupt
fifteen of the twenty-five panels and every automatic check would pass.

Reading that table is a human step with no substitute, which is why it is
displayed rather than merely produced.
"""

import base64
import html
import os
import sys

spec = dict(snakemake.params.spec)
sys.path.insert(0, spec["libdir"])
sys.path.insert(0, os.path.join(spec["libdir"], "scripts", "py", "panels"))

import pandas as pd                                            # noqa: E402

from lib import panels as _panels                              # noqa: E402
from lib.paths import Paths                                    # noqa: E402

MIME = {".svg": "image/svg+xml", ".png": "image/png", ".pdf": "application/pdf"}

#: Section order and titles.
GROUP_TITLES = [
    ("figure_01", "Figure 1"),
    ("figure_02", "Figure 2"),
    ("figure_03", "Figure 3"),
    ("figure_04", "Figure 4"),
    ("figure_05", "Figure 5"),
    ("figure_06", "Figure 6"),
    ("extended_02", "Extended Data Figure 2"),
    ("extended_03", "Extended Data Figure 3"),
    ("extended_06", "Extended Data Figure 6"),
    ("extended_08", "Extended Data Figure 8"),
    ("extended_10", "Extended Data Figure 10"),
    ("qc", "Quality control"),
    ("extras", "Other pipeline figures (not in the paper)"),
]

CSS = """
:root { color-scheme: light dark; }
* { box-sizing: border-box; }
body { font-family: -apple-system, "Segoe UI", Helvetica, Arial, sans-serif;
       margin: 0 auto; max-width: 1200px; padding: 2rem 1.25rem 6rem;
       line-height: 1.5; color: #1a1a1a; background: #fff; }
h1 { font-size: 1.8rem; margin: 0 0 .25rem; }
h2 { font-size: 1.3rem; margin: 2.5rem 0 .75rem; padding-bottom: .3rem;
     border-bottom: 2px solid #e5e5e5; }
h3 { font-size: 1rem; margin: 1.5rem 0 .4rem; }
.sub { color: #666; margin: 0 0 1.5rem; }
.panel { border: 1px solid #e0e0e0; border-radius: 6px; padding: 1rem;
         margin: 0 0 1.25rem; overflow: hidden; }
.panel img, .panel object { max-width: 100%; height: auto; display: block;
                            margin: .5rem auto; }
.meta { font-size: .82rem; color: #555; }
.meta code { background: #f4f4f4; padding: .1rem .3rem; border-radius: 3px;
             font-size: .95em; word-break: break-all; }
.note { font-size: .82rem; color: #444; margin-top: .5rem;
        border-left: 3px solid #ddd; padding-left: .6rem; }
.badge { display: inline-block; font-size: .72rem; font-weight: 600;
         padding: .12rem .45rem; border-radius: 3px; margin-right: .35rem;
         vertical-align: middle; }
.b-ok { background: #d7f0d7; color: #14501c; }
.b-missing, .b-empty { background: #fde2e1; color: #7a1d18; }
.b-blocked, .b-skipped, .b-disabled { background: #fdf0d5; color: #7a5310; }
.b-ref { background: #dde7fb; color: #1c3d75; }
.b-med { background: #f3e3fb; color: #55206f; }
.wrap { overflow-x: auto; }
table { border-collapse: collapse; font-size: .8rem; width: 100%; }
th, td { text-align: left; padding: .3rem .55rem; border-bottom: 1px solid #eee;
         white-space: nowrap; }
th { background: #fafafa; position: sticky; top: 0; }
tr.bad td { background: #fdf3f2; }
.callout { border: 2px solid #d97706; border-radius: 6px; padding: 1rem;
           margin: 1.5rem 0; background: #fffbeb; }
.callout h3 { margin-top: 0; color: #92400e; }
details { margin: .5rem 0; }
summary { cursor: pointer; font-size: .85rem; color: #444; }
@media (prefers-color-scheme: dark) {
  body { background: #16181c; color: #e6e6e6; }
  h2 { border-color: #333; }
  .panel { border-color: #333; }
  .meta, .note { color: #aaa; }
  .meta code { background: #24262b; }
  th { background: #1d1f24; }
  td, th { border-color: #2a2c31; }
  .note { border-color: #333; }
  .callout { background: #2a2213; border-color: #b45309; }
  .callout h3 { color: #fbbf24; }
  tr.bad td { background: #2c1d1c; }
}
"""


def esc(text):
    return html.escape(str(text), quote=True)


def embed(path, max_bytes):
    """Inline a figure as a data URI, or explain why it is only linked."""
    ext = os.path.splitext(path)[1].lower()
    if not os.path.exists(path):
        return '<p class="meta">file not produced</p>'
    size = os.path.getsize(path)
    if size == 0:
        return '<p class="meta">file is empty (0 bytes)</p>'
    if ext == ".pdf":
        return ('<p class="meta">PDF, {:.1f} KB - not inlined; open '
                '<code>{}</code></p>'.format(size / 1024.0, esc(path)))
    if size > max_bytes:
        return ('<p class="meta">{:.1f} MB exceeds the inline limit; open '
                '<code>{}</code></p>'.format(size / 1e6, esc(path)))
    with open(path, "rb") as fh:
        payload = base64.b64encode(fh.read()).decode("ascii")
    return '<img src="data:{};base64,{}" alt="{}">'.format(
        MIME.get(ext, "application/octet-stream"), payload,
        esc(os.path.basename(path)))


def table_html(frame, bad_column=None, bad_values=()):
    if frame is None or not len(frame):
        return '<p class="meta">no rows</p>'
    head = "".join("<th>{}</th>".format(esc(c)) for c in frame.columns)
    body = []
    for _, row in frame.iterrows():
        css = ""
        if bad_column and str(row.get(bad_column)) in bad_values:
            css = ' class="bad"'
        cells = "".join("<td>{}</td>".format(esc(v)) for v in row)
        body.append("<tr{}>{}</tr>".format(css, cells))
    return ('<div class="wrap"><table><thead><tr>{}</tr></thead>'
            "<tbody>{}</tbody></table></div>".format(head, "".join(body)))


log_path = snakemake.log[0]
os.makedirs(os.path.dirname(log_path) or ".", exist_ok=True)
with open(log_path, "w") as log:
    def say(msg=""):
        log.write(str(msg) + "\n")
        log.flush()

    P = Paths(snakemake.config)
    registry = _panels.from_records(snakemake.params.registry,
                                    snakemake.params.get("settings"))
    max_bytes = int(spec.get("max_inline_bytes", 6000000))

    manifest = pd.read_csv(str(snakemake.input.manifest), sep="\t")
    try:
        numbers = pd.read_csv(str(snakemake.input.numbers), sep="\t")
    except Exception:
        numbers = None

    by_stem = {}
    for _, row in manifest.iterrows():
        by_stem.setdefault(str(row["panel_id"]), []).append(row)

    parts = []
    parts.append("<h1>XCI pipeline - figure report</h1>")
    parts.append('<p class="sub">pipeline {} &middot; {} panels, {} files '
                 "&middot; every panel below is inlined; this page has no "
                 "external assets.</p>".format(
                     esc(spec.get("pipeline_version", "?")),
                     len(registry), len(manifest)))

    # --- the eigenvector orientation callout -----------------------------
    orient_path = spec.get("eigenvector_orientation")
    orient = None
    if orient_path and os.path.exists(orient_path):
        try:
            orient = pd.read_csv(orient_path, sep="\t")
        except Exception:
            orient = None
    flipped = 0
    if orient is not None:
        for column in ("flipped", "flip", "sign_flipped"):
            if column in orient.columns:
                flipped = int(orient[column].astype(str).str.lower()
                              .isin(("true", "1", "yes", "-1")).sum())
                break
    parts.append('<div class="callout">')
    parts.append("<h3>Eigenvector orientation - read this before the saddles</h3>")
    parts.append(
        "<p>A <strong>global</strong> eigenvector sign flip leaves the saddle "
        "corner score (AA + BB) / (AB + BA) <strong>unchanged</strong> - the "
        "corners swap in pairs - while mirroring the plot. No numeric check on "
        "<code>saddle_strength_selected_*.tsv</code> can catch it, and "
        "<code>compartments*.json</code> records which eigenvector was used but "
        "never its orientation. A silent flip would mirror every saddle and "
        "quietly corrupt 15 of the 25 panels, so this table is a human step "
        "with no substitute.</p>")
    if orient is None:
        parts.append('<p class="meta">Not produced in this run: <code>{}</code>'
                     "</p>".format(esc(orient_path or "-")))
    else:
        parts.append('<p class="meta">{} rows, {} flip(s) applied - '
                     "<code>{}</code></p>".format(
                         len(orient), flipped, esc(orient_path)))
        parts.append(table_html(orient))
    parts.append("</div>")

    # --- the numbers ------------------------------------------------------
    if numbers is not None:
        mism = numbers[numbers["status"] == "MISMATCH"] \
            if "status" in numbers.columns else numbers.iloc[0:0]
        checked = numbers[numbers["status"].isin(("ok", "MISMATCH"))] \
            if "status" in numbers.columns else numbers.iloc[0:0]
        parts.append("<h2>Published numbers</h2>")
        parts.append('<p class="meta">{} quantities, {} checked against the '
                     "published values, {} mismatched.</p>".format(
                         len(numbers), len(checked), len(mism)))
        if len(checked):
            parts.append(table_html(checked, "status", {"MISMATCH"}))
        parts.append("<details><summary>All {} quantities</summary>{}"
                     "</details>".format(len(numbers), table_html(numbers)))

    # --- the panels, section by section ----------------------------------
    seen_groups = set()
    for group, title in GROUP_TITLES:
        panels = [p for p in registry if p.group == group]
        if not panels:
            continue
        seen_groups.add(group)
        parts.append("<h2>{}</h2>".format(esc(title)))
        for panel in sorted(panels, key=lambda p: (p.paper_ref and
                                                   p.paper_ref[0] or "zz",
                                                   p.id)):
            for stem in panel.stems():
                rows = by_stem.get(stem, [])
                preferred = None
                for row in rows:
                    ext = os.path.splitext(str(row["file"]))[1].lower()
                    if ext in (".svg", ".png") and row["status"] == "ok":
                        preferred = row
                        break
                # `preferred` is a pandas Series; `a or b` on one raises
                # rather than testing for None.
                if preferred is None:
                    row = rows[0] if rows else None
                else:
                    row = preferred
                status = str(row["status"]) if row is not None else "missing"
                path = str(row["file"]) if row is not None else ""
                n_items = row["n_items"] if row is not None else ""

                parts.append('<div class="panel">')
                badges = ['<span class="badge b-{}">{}</span>'.format(
                    esc(status), esc(status))]
                for ref in panel.paper_ref:
                    badges.append('<span class="badge b-ref">{}</span>'.format(
                        esc(ref)))
                if panel.confidence == "medium":
                    badges.append('<span class="badge b-med">confidence: '
                                  "medium</span>")
                parts.append("<h3>{}{}</h3>".format("".join(badges), esc(stem)))
                parts.append(embed(path, max_bytes))
                parts.append('<p class="meta"><code>{}</code>{}</p>'.format(
                    esc(path),
                    " &middot; n = {}".format(esc(n_items))
                    if str(n_items) not in ("", "nan") else ""))
                if panel.candidates and stem == panel.id:
                    parts.append('<p class="meta">This panel is medium '
                                 "confidence: the archaeology could not "
                                 "identify which rendering the paper used. "
                                 "Candidates also emitted: {}</p>".format(
                                     esc(", ".join(panel.candidates))))
                if panel.notes:
                    parts.append('<p class="note">{}</p>'.format(
                        esc(panel.notes)))
                if panel.reason:
                    parts.append('<p class="note">{}</p>'.format(
                        esc(panel.reason)))
                parts.append("</div>")

    leftover = sorted({p.group for p in registry} - seen_groups)
    if leftover:
        say("groups with no section in GROUP_TITLES: {}".format(leftover))

    # --- the manifest ------------------------------------------------------
    parts.append("<h2>Figure manifest</h2>")
    parts.append('<p class="meta">Also written as <code>{}</code>. Every '
                 "registry entry appears, built or not.</p>".format(
                     esc(P.manifest())))
    parts.append(table_html(manifest, "status",
                            {"missing", "empty", "blocked"}))

    html_out = (
        "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n"
        '<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
        "<title>XCI pipeline - figure report</title>\n"
        "<style>{}</style>\n</head>\n<body>\n{}\n</body>\n</html>\n".format(
            CSS, "\n".join(parts)))

    out_path = str(snakemake.output.html)
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(out_path, "w") as fh:
        fh.write(html_out)

    say("wrote {} ({:.1f} KB, {} manifest rows, {} panels)".format(
        out_path, len(html_out) / 1024.0, len(manifest), len(registry)))
