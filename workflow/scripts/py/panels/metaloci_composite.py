"""Fig 6 — the METALoci Gaudi montage, one file per locus.

A PyMuPDF montage of PDFs that other panels already drew, exactly as
`02_03_make_composite_metaloci_figures copy 2.ipynb` did: it opened the written
`_gtp.pdf` files with `fitz` and placed them on a new page.

THE TRAP (F-6)
--------------
The row and column order is hard-coded in that notebook, is NOT alphabetical,
and is NOT derivable from the data. It lives in
`resources/fixtures/fig6_layout.yaml`:

    rows     CTCF, H3K27me3, H3K27ac, AcMe3, RNA-Seq
    columns  Kdm5c: violin, Xa_consensus, E6, C5, B1
             Mecp2: violin, Xa_consensus, E6, B1, C5, CL30, JTG

Note Kdm5c is E6, C5, B1 while Mecp2 is E6, B1, C5 -- the two loci genuinely
differ in the third and fourth columns. Getting this wrong produces a Fig 6
that is visually plausible and silently mislabelled, which is why it is a
versioned fixture and not a literal in this file.

Placement, also verbatim: column-major; the label
`f"{label}: {value:.2f}"` at fontsize 30 where `value` is the METALoci
compartmentalization; column 0 (the violin) advances y by 150 before placing,
every other column by 20.

One thing that looks like a bug and is not: the Xi panels are normalised
per clone (METALOCI_NEW) while the Xa panel is consensus-normalised
(METALOCI_NEW_CONSENSUS), so the two are not strictly on the same signal
scale. That is what the original did and what the paper shows.
"""

import os
import re
import sys

sys.path.insert(0, os.path.join(snakemake.params.spec["libdir"],
                                "scripts", "py", "panels"))
sys.path.insert(0, snakemake.params.spec["libdir"])

import yaml                                                    # noqa: E402
import _common as C                                            # noqa: E402
from lib import plotting as pl                                 # noqa: E402

spec = C.bootstrap(snakemake)
out_path = C.outputs(snakemake)


def _compartmentalization(paths):
    """dataset -> signal -> value, from every compartmentalization table."""
    out = {}
    for path in paths:
        if "compartmentalization" not in os.path.basename(path):
            continue
        try:
            table = C.read_table(path)
        except Exception:
            continue
        if not {"dataset", "signal", "compartmentalization"} <= set(table.columns):
            continue
        for _, row in table.iterrows():
            out.setdefault(str(row["dataset"]), {})[str(row["signal"])] = \
                float(row["compartmentalization"])
    return out


def _lookup(values, dataset, signal):
    """METALoci datasets carry `_rep1` / `_G1` variants; match the stem."""
    if dataset in values and signal in values[dataset]:
        return values[dataset][signal]
    for key, sigs in values.items():
        if key.startswith(dataset) and signal in sigs:
            return sigs[signal]
    stem = re.sub(r"_(G1|G2)", "", str(dataset))
    for key, sigs in values.items():
        if re.sub(r"_(G1|G2)", "", key).startswith(stem) and signal in sigs:
            return sigs[signal]
    return None


with C.logging(snakemake) as say:
    locus = spec.get("locus", "Mecp2")
    slug = spec.get("locus_slug", locus.lower())
    display = C.display_locus(locus, spec.get("display_names"))
    paths = [str(p) for p in snakemake.input]

    layout_path = next(p for p in paths if p.endswith("fig6_layout.yaml"))
    with open(layout_path) as fh:
        layout = yaml.safe_load(fh)
    rows = list(layout["rows"])
    columns = list(layout["columns"].get(slug, []))
    label_fmt = layout.get("panel_label", "{label}: {value:.2f}")
    say("panel {} ({}): {} rows x {} columns".format(
        spec["panel_id"], display, len(rows), len(columns)))
    say("  rows    {}".format(", ".join(rows)))
    say("  columns {}".format(", ".join(columns)))

    values = _compartmentalization(paths)
    pdfs = {os.path.basename(p): p for p in paths if p.endswith(".pdf")}

    # Rebuild the grid of file paths and labels the same way the input
    # function did, so a missing panel is reported by name rather than by
    # a KeyError deep inside fitz.
    grid, labels, cell_values = [], [], []
    for signal in rows:
        row_paths, row_labels, row_values = [], [], []
        for column in columns:
            if column == "violin":
                base = "metaloci_violin_{}_{}.pdf".format(locus, signal)
                row_labels.append("")
                row_values.append(0.0)
            else:
                match = [name for name in pdfs
                         if name.startswith("metaloci_gaudi_")
                         and name.endswith("_{}.pdf".format(signal))
                         and ("_{}_".format(column) in name
                              or (column == "Xa_consensus"
                                  and "consensus" in name))]
                base = match[0] if match else None
                dataset = ""
                if base:
                    dataset = base[len("metaloci_gaudi_"):
                                   -len("_{}.pdf".format(signal))]
                    dataset = dataset.split("_", 1)[1] if "_" in dataset else dataset
                row_labels.append("{}_{}".format(signal, column))
                row_values.append(_lookup(values, dataset, signal))
            row_paths.append(pdfs.get(base) if base else None)
        grid.append(row_paths)
        labels.append(row_labels)
        cell_values.append(row_values)

    missing = [(rows[i], columns[j])
               for i, row in enumerate(grid)
               for j, path in enumerate(row) if not path]
    for signal, column in missing:
        say("  MISSING cell: row {} column {}".format(signal, column))

    try:
        import fitz                                            # PyMuPDF
    except ImportError:
        pl.empty_panel(out_path,
                       "PyMuPDF (fitz) is not installed;\n"
                       "the Fig 6 montage needs it -- see environment.yml")
        C.write_n_items(snakemake, 0)
        raise SystemExit(0)

    # --- pass 1: sizes, column-major, exactly as the original ------------
    column_widths, layout_data = [], []
    max_height, max_width = 0, 0
    n_rows = 0
    for col_idx in range(len(columns)):
        col_data, row_heights = [], []
        for row_idx in range(len(rows)):
            path = grid[row_idx][col_idx]
            if not path:
                continue
            doc = fitz.open(path)
            page = doc[0]
            width, height = page.rect.width, page.rect.height
            doc.close()
            col_data.append((path, width, height, row_idx, col_idx))
            max_width = max(max_width, width)
            row_heights.append(height)
        max_height = max(sum(row_heights), max_height)
        n_rows = max(n_rows, len(row_heights))
        layout_data.append(col_data)
        column_widths.append(max_width)

    if not any(layout_data):
        pl.empty_panel(out_path, "no Gaudi or violin PDFs for {}".format(display))
        C.write_n_items(snakemake, 0)
        raise SystemExit(0)

    canvas_width = sum(column_widths)
    canvas_height = max_height + 20 * (n_rows + 1)

    # --- pass 2: place ----------------------------------------------------
    output = fitz.open()
    page = output.new_page(width=canvas_width, height=canvas_height)
    x_cursor = 0
    placed = 0
    for col_idx, col in enumerate(layout_data):
        y_cursor = 20
        for (path, width, height, row_idx, _) in col:
            label = labels[row_idx][col_idx]
            value = cell_values[row_idx][col_idx]
            if value is None:
                text = "{}: n/a".format(label)
            else:
                text = label_fmt.format(label=label, value=value)
            page.insert_text((x_cursor + 10, y_cursor + 10), text,
                             fontsize=30, color=(0, 0, 0))
            y_cursor += 150 if col_idx == 0 else 20
            rect = fitz.Rect(x_cursor, y_cursor,
                             x_cursor + width, y_cursor + height)
            doc = fitz.open(path)
            page.show_pdf_page(rect, doc, 0)
            doc.close()
            y_cursor += height
            placed += 1
        x_cursor += column_widths[col_idx]

    C.ensure_parent(out_path)
    if out_path.lower().endswith(".pdf"):
        output.save(out_path)
    else:
        # fitz writes SVG one page at a time.
        svg = page.get_svg_image()
        with open(out_path, "w") as fh:
            fh.write(svg)
    output.close()

    C.write_n_items(snakemake, placed,
                    extra={"rows": len(rows), "columns": len(columns),
                           "missing_cells": len(missing)})
    say("placed {} panels, {} missing -> {}".format(
        placed, len(missing), out_path))
