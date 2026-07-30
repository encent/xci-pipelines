"""`paper_numbers.tsv` — every quantity the paper quotes, with its source.

The user asked for "essentially all the figures, numbers and outcomes that are
used in the paper". The figures are the manifest's job; this file is the
numbers.

    quantity  value  expected  status  unit  source_table  paper_ref  panel_id  notes

`expected` is filled in for the quantities the archaeology recovered from the
published figures, and `status` compares the two:

    ok         computed == expected
    MISMATCH   computed != expected -- a reproduction failure
    -          no published value to check against; the row is informational
    missing    the source table does not exist in this run

With `settings.strict_paper_numbers: true` in `config/panels.yaml` (the
default) a MISMATCH fails this rule. That is deliberate: these are the numbers
that say whether the pipeline reproduced the paper, and a run that quietly
writes 377 where the paper says 375 is worse than a run that stops.

THE ASSERTION THAT MATTERS MOST
-------------------------------
Fig 2g plots the GENE-FILTERED valley counts: E6 375, C5 301, B1 258, JTG 300,
CL30 321. The raw HMM calls are 377 / 304 / 261 / 302 / 324. Both make a
plausible bar chart. If `panel_valley_gene_content` is ever re-pointed at the
raw set, this file is what notices.

Both sets are computed and both are written: the raw counts are also the
exact-diff target for `call_valleys`, so having them side by side turns "the
bars are wrong" into "the filter did not run" in one glance.
"""

import os
import sys

spec = dict(snakemake.params.spec)
sys.path.insert(0, spec["libdir"])
sys.path.insert(0, os.path.join(spec["libdir"], "scripts", "py", "panels"))

import pandas as pd                                            # noqa: E402

COLUMNS = ["quantity", "value", "expected", "status", "unit", "source_table",
           "paper_ref", "panel_id", "notes"]

ESCAPING = "escaping-gene-valley"

#: Published, from the verification notes section 8 and the archaeology notes F.1.
EXPECTED_GENE_FILTERED = {"E6": 375, "C5": 301, "B1": 258, "JTG": 300,
                          "CL30": 321}
EXPECTED_RAW_XI = {
    "B1_WT": 261, "C5_WT": 304, "CL30_WT": 324, "E6_WT": 377, "JTG_WT": 302,
    "B1621_Rad21-NodTAG": 284, "B1621_Rad21-dTAG": 281,
    "E6A7_CTCF-NodTAG": 321, "E6A7_CTCF-dTAG": 281,
    "F3_CTCF-NodTAG": 351, "F3_CTCF-dTAG": 205,
}
#: EFig 6b: escaping boundaries and their CTCF split.
EXPECTED_DTAG_BOUNDARIES = {"E6A7": (110, 80, 30), "F3": (36, 23, 13)}
#: EFig 6c Venns: (dTAG-only, shared, NodTAG-only).
EXPECTED_VENN = {"E6A7": (1, 44, 11), "F3": (1, 12, 6)}
#: Refined loop counts -- byte-identical to the frozen fixture.
EXPECTED_LOOPS_RAW = {"Jarid_Xa": 55, "Jarid_Xi": 32,
                      "Mecp2_Xa": 41, "Mecp2_Xi": 26}
EXPECTED_LOOPS_REFINED = {"Jarid_Xa": 14, "Jarid_Xi": 13,
                          "Mecp2_Xa": 10, "Mecp2_Xi": 12}
#: Fig 2h pie n-values, in the published clone order.
EXPECTED_FIG2H = [110, 96, 74, 76, 76]


def _label(track):
    parts = str(track).split("_")
    return "_".join(parts[1:-1]) if len(parts) >= 4 else str(track)


def _clone(track):
    parts = str(track).split("_")
    return parts[1] if len(parts) > 1 else str(track)


rows = []


def add(quantity, value, expected=None, unit="count", source_table="",
        paper_ref="", panel_id="", notes=""):
    if value is None:
        status = "missing"
        value = ""
    elif expected is None:
        status = "-"
    else:
        try:
            status = "ok" if int(value) == int(expected) else "MISMATCH"
        except (TypeError, ValueError):
            status = "ok" if str(value) == str(expected) else "MISMATCH"
    rows.append({
        "quantity": quantity,
        "value": value,
        "expected": "" if expected is None else expected,
        "status": status,
        "unit": unit,
        "source_table": source_table,
        "paper_ref": paper_ref,
        "panel_id": panel_id,
        "notes": notes,
    })
    return status


def _read(path):
    if not path or not os.path.exists(path):
        return None
    try:
        return pd.read_csv(path, sep="\t")
    except Exception:
        return None


def _count_bed(path):
    if not path or not os.path.exists(path):
        return None
    n = 0
    with open(path) as fh:
        for line in fh:
            if line.strip() and not line.startswith(("#", "track", "browser")):
                n += 1
    return n


log_path = snakemake.log[0]
os.makedirs(os.path.dirname(log_path) or ".", exist_ok=True)
with open(log_path, "w") as log:
    def say(msg=""):
        log.write(str(msg) + "\n")
        log.flush()

    # ---------------------------------------------------------------------
    # Fig 2g / Fig 4d -- the gene-filtered valley counts. TRAP P-2.
    # ---------------------------------------------------------------------
    for track, path in sorted((spec.get("gene_content") or {}).items()):
        table = _read(path)
        clone = _clone(track)
        is_wt = "_WT_" in track
        expected = EXPECTED_GENE_FILTERED.get(clone) if is_wt else None
        total = None if table is None else int(len(table))
        status = add(
            "valleys_gene_filtered.{}".format(_label(track)),
            total, expected, "valleys", path,
            "Fig 2g" if is_wt else "Fig 4d",
            "valley_gene_content_wt" if is_wt else "valley_gene_content_degron",
            "GENE-FILTERED set (P-2). The raw HMM calls differ by 1-3 per "
            "clone and make an equally plausible bar chart.",
        )
        if expected is not None:
            say("Fig 2g {:6s} filtered {} (published {}) -> {}".format(
                clone, total, expected, status))
        if table is None:
            continue
        counts = table["category"].value_counts() if "category" in table else {}
        for category in ("no-gene", "no-expressed-gene-valley",
                         "silent-gene-valley", ESCAPING):
            add("valley_class.{}.{}".format(_label(track), category),
                int(counts.get(category, 0)), None, "valleys", path,
                "Fig 2g" if is_wt else "Fig 4d",
                "valley_gene_content_wt" if is_wt else
                "valley_gene_content_degron")

    # ---------------------------------------------------------------------
    # Raw and filtered valley counts -- the exact-diff target for call_valleys
    # ---------------------------------------------------------------------
    for track, path in sorted((spec.get("valleys") or {}).items()):
        key = "_".join(str(track).split("_")[1:3])
        add("valleys_raw.{}".format(_label(track)), _count_bed(path),
            EXPECTED_RAW_XI.get(key), "valleys", path, "", "",
            "raw HMM calls, before the 3-gene filter; includes the "
            "chrX:0-3,285,000 assembly-gap valley, which the paper counts too")
    for track, path in sorted((spec.get("valleys_filtered") or {}).items()):
        add("valleys_filtered_bed.{}".format(_label(track)), _count_bed(path),
            None, "valleys", path, "", "",
            "after removing valleys overlapping Mid1, Tmem29, Firre")

    # ---------------------------------------------------------------------
    # Fig 2h / EFig 2g -- boundaries, CTCF status, Fisher p
    # ---------------------------------------------------------------------
    fig2h_order = [t for t in (spec.get("wt_tracks") or [])]
    fig2h_values = []
    for track, path in sorted((spec.get("boundary_ctcf") or {}).items()):
        table = _read(path)
        is_wt = "_WT_" in track
        if table is None:
            add("boundaries.{}".format(_label(track)), None, None,
                "boundaries", path, "Fig 2h" if is_wt else "", "")
            continue
        escaping = table[table["valley_class"] == ESCAPING] \
            if "valley_class" in table else table.iloc[0:0]
        n_esc = int(len(escaping))
        n_ctcf = int(escaping["has_ctcf"].astype(bool).sum()) \
            if "has_ctcf" in escaping else 0
        add("boundaries_total.{}".format(_label(track)), int(len(table)),
            None, "boundaries", path, "", "",
            "boundaries/all is 2 x the valley count")
        expected_pie = None
        if track in fig2h_order:
            index = fig2h_order.index(track)
            if index < len(EXPECTED_FIG2H):
                expected_pie = EXPECTED_FIG2H[index]
            fig2h_values.append((_label(track), n_esc))
        add("boundaries_escaping.{}".format(_label(track)), n_esc,
            expected_pie, "boundaries", path,
            "Fig 2h" if is_wt else "", "boundary_ctcf_pies_wt" if is_wt else "",
            "the pie n-value")
        add("boundaries_escaping_ctcf.{}".format(_label(track)), n_ctcf,
            None, "boundaries", path, "Fig 2h" if is_wt else "", "",
            "CTCF peak within 50 kb inward / 10 kb outward")
        if "fisher_p" in table and len(table):
            add("fisher_p_escapee_x_ctcf.{}".format(_label(track)),
                float(table["fisher_p"].iloc[0]), None, "p-value", path,
                "EFig 2g" if is_wt else "",
                "boundary_analysis_{}".format(_clone(track)) if is_wt else "",
                "two-sided Fisher exact on the 2x2 Escapee x CTCF table")

    if fig2h_values:
        say("Fig 2h escaping-boundary n per clone: {}  (published {})".format(
            ", ".join("{}={}".format(k, v) for k, v in fig2h_values),
            EXPECTED_FIG2H))

    # ---------------------------------------------------------------------
    # Fig 4f/4g + EFig 6b -- the degron boundary counts
    # ---------------------------------------------------------------------
    for clone, (n_all, n_pos, n_neg) in EXPECTED_DTAG_BOUNDARIES.items():
        track = next((t for t in (spec.get("degron_nodtag_tracks") or [])
                      if _clone(t) == clone), None)
        path = (spec.get("boundary_ctcf") or {}).get(track)
        table = _read(path)
        if table is None:
            add("dtag_escaping_boundaries.{}".format(clone), None, n_all,
                "boundaries", path or "", "EFig 6b", "boundary_dtag_" + clone)
            continue
        escaping = table[table["valley_class"] == ESCAPING]
        got_all = int(len(escaping))
        got_pos = int(escaping["has_ctcf"].astype(bool).sum())
        s1 = add("dtag_escaping_boundaries.{}".format(clone), got_all, n_all,
                 "boundaries", path, "EFig 6b", "boundary_dtag_" + clone)
        s2 = add("dtag_escaping_boundaries_ctcf.{}".format(clone), got_pos,
                 n_pos, "boundaries", path, "Fig 4f;Fig 4g",
                 "boundary_dtag_" + clone,
                 "Fig 4f/4g quote the CTCF+ subset")
        s3 = add("dtag_escaping_boundaries_no_ctcf.{}".format(clone),
                 got_all - got_pos, n_neg, "boundaries", path, "EFig 6b",
                 "boundary_dtag_" + clone)
        say("EFig 6b {:6s} {}/{}/{} (published {}/{}/{}) -> {} {} {}".format(
            clone, got_all, got_pos, got_all - got_pos, n_all, n_pos, n_neg,
            s1, s2, s3))

    # ---------------------------------------------------------------------
    # EFig 6c -- the Venn cardinalities
    # ---------------------------------------------------------------------
    for clone, path in sorted((spec.get("valley_overlap") or {}).items()):
        table = _read(path)
        expected = EXPECTED_VENN.get(clone)
        if table is None or table.empty:
            add("venn_escaping.{}".format(clone), None,
                None if expected is None else "/".join(map(str, expected)),
                "valleys", path or "", "EFig 6c", "valley_venn_dtag")
            continue
        row = table.iloc[0]
        got = (int(row["dtag_only"]), int(row["nodtag_overlap"]),
               int(row["nodtag_only"]))
        status = add(
            "venn_escaping.{}".format(clone), "/".join(map(str, got)),
            None if expected is None else "/".join(map(str, expected)),
            "dTAG-only/shared/NodTAG-only", path, "EFig 6c",
            "valley_venn_dtag",
            "the two overlap counts differ when valleys merge or split; "
            "dtag_overlap = {}".format(int(row["dtag_overlap"])))
        say("EFig 6c {:6s} {} (published {}) -> {}".format(
            clone, got, expected, status))

    # ---------------------------------------------------------------------
    # Loops -- raw and refined
    # ---------------------------------------------------------------------
    for name, path in sorted((spec.get("loops_raw") or {}).items()):
        table = _read(path)
        add("loops_raw.{}".format(name),
            None if table is None else int(len(table)),
            EXPECTED_LOOPS_RAW.get(name), "loops", path or "", "", "",
            "chromosight calls before manual refinement")
    for name, path in sorted((spec.get("loops_refined") or {}).items()):
        table = _read(path)
        status = add("loops_refined.{}".format(name),
                     None if table is None else int(len(table)),
                     EXPECTED_LOOPS_REFINED.get(name), "loops", path or "",
                     "", "loop_overlay_" + name,
                     "frozen fixture; must be byte-identical")
        say("loops_refined {:12s} -> {}".format(name, status))

    # ---------------------------------------------------------------------
    # Loop and compartment strength, and the violin statistics
    # ---------------------------------------------------------------------
    for comparison, path in sorted((spec.get("pileup_scores") or {}).items()):
        table = _read(path)
        if table is None:
            continue
        for _, row in table.iterrows():
            add("pileup_score.{}.{}".format(comparison, row["file"]),
                float(row["mean"]), None, "obs/exp", path,
                "EFig 3a" if comparison == "Xa_vs_Xi" else "",
                "violin_loops_xa_vs_xi" if comparison == "Xa_vs_Xi" else "",
                "what the paper's loop-strength violins plot; 01_03's "
                "loop_stats writes are commented out")

    for comparison, path in sorted((spec.get("saddle_strength") or {}).items()):
        table = _read(path)
        if table is None:
            continue
        cols = list(table.columns)
        for _, row in table.iterrows():
            add("saddle_strength.{}.{}".format(comparison, row[cols[0]]),
                float(row[cols[1]]), None, "(AA+BB)/(AB+BA)", path,
                "EFig 3b" if comparison == "Xa_vs_Xi" else "",
                "violin_saddle_xa_vs_xi" if comparison == "Xa_vs_Xi" else "",
                "a GLOBAL eigenvector sign flip leaves this score unchanged "
                "while mirroring the plot -- check "
                "work/features/compartments/eigenvector_orientation.tsv")

    comp = _read(spec.get("compartmentalization"))
    if comp is not None and "compartmentalization" in comp.columns:
        for _, row in comp.iterrows():
            add("metaloci_compartmentalization.{}.{}".format(
                row.get("dataset", "?"), row.get("signal", "?")),
                float(row["compartmentalization"]), None, "percent",
                spec.get("compartmentalization", ""), "Fig 6",
                "metaloci_composite",
                "(sq1 + sq3) / (q1+q2+q3+q4) * 100; the Fig 6 panel labels")

    # ---------------------------------------------------------------------
    # write
    # ---------------------------------------------------------------------
    out_path = str(snakemake.output.tsv)
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    frame = pd.DataFrame(rows, columns=COLUMNS)
    frame.to_csv(out_path, sep="\t", index=False)

    counts = frame["status"].value_counts().to_dict()
    say("")
    say("{} quantities: {}".format(
        len(frame), ", ".join("{} {}".format(v, k)
                              for k, v in sorted(counts.items()))))

    mismatches = frame[frame["status"] == "MISMATCH"]
    if len(mismatches):
        say("")
        say("MISMATCHES against the published values:")
        for _, row in mismatches.iterrows():
            say("  {:52s} computed {}  published {}".format(
                row["quantity"], row["value"], row["expected"]))

    say("wrote {}".format(out_path))

    if len(mismatches) and spec.get("strict", True):
        raise SystemExit(
            "paper_numbers: {} quantity/quantities do not match the published "
            "values (see {} and {}). Set `settings.strict_paper_numbers: "
            "false` in config/panels.yaml to record them and continue."
            .format(len(mismatches), out_path, log_path))
