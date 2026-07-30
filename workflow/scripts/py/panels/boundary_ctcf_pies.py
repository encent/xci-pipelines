"""Fig 2h — CTCF status of escapee valley boundaries, one pie per clone.

A boundary counts as CTCF-positive if a consensus CTCF peak lies within 50 kb
on the valley-interior side or 10 kb on the exterior side (`boundary_ctcf_status`
computes it; this panel only draws it). Published n per clone: 110, 96, 74, 76,
76 — those are ESCAPEE boundaries, not all boundaries.

A note on provenance, because the archaeology is genuinely ambiguous here.
the archaeology notes F.1 attributes Fig 2h to `ctcf_both_sides_proportion.svg`, the last
cell of `xci_valleys_check_LAST.ipynb`. That file is not a pie chart: it is the
proportion of valleys with a peak on BOTH sides, with Wald CIs. The published
pies with n = 110/96/74/76/76 are panel [0,2] of
`ctcf_boundary_analysis/*_boundary_analysis.svg`. This rule renders both:

    cohort=WT | degron     the per-clone escapee-boundary CTCF pies (Fig 2h)
    cohort=both_sides      the both-sides proportion plot

so whichever reading is right, the figure exists and the manifest says which is
which.
"""

import os
import sys

sys.path.insert(0, os.path.join(snakemake.params.spec["libdir"],
                                "scripts", "py", "panels"))
sys.path.insert(0, snakemake.params.spec["libdir"])

import numpy as np                                             # noqa: E402
import _common as C                                            # noqa: E402
from lib import plotting as pl                                 # noqa: E402

spec = C.bootstrap(snakemake)
out_path = C.outputs(snakemake)
ESCAPING = "escaping-gene-valley"


def _load(paths):
    """One record per clone, from the boundary_ctcf tables."""
    records = []
    for path in sorted(C.group_inputs(paths)):
        if os.sep + "boundary_ctcf" + os.sep not in path:
            continue
        track = os.path.basename(path)[: -len(".tsv")]
        table = C.read_table(path)
        cls = C.need(table, "valley_class", path)
        has_ctcf = C.need(table, "has_ctcf", path).astype(bool)
        escaping = cls == ESCAPING
        records.append({
            "track": track,
            "label": C.track_label(track),
            "n_escaping": int(escaping.sum()),
            "n_escaping_ctcf": int((escaping & has_ctcf).sum()),
            "n_total": int(len(table)),
            "table": table,
            "path": path,
        })
    return records


def _both_sides_row(record):
    """Per-valley L/R co-occurrence: FF, FT, TF, TT.

    A valley contributes one row; `side` says which end each boundary is.
    The published question is: among valleys with a peak on at least one side,
    how often is there one on both.
    """
    table = record["table"]
    side = C.need(table, "side", record["path"])
    has = C.need(table, "has_ctcf", record["path"]).astype(bool)
    start = C.need(table, "start", record["path"])
    chrom = C.need(table, "chrom", record["path"])

    left, right = {}, {}
    for ch, st, sd, hc in zip(chrom, start, side, has):
        (left if str(sd).upper().startswith("L") else right)[(ch, int(st))] = bool(hc)

    # Boundaries are the valley edges, so pair them up in genomic order.
    lefts = sorted(left)
    rights = sorted(right)
    pairs = list(zip(lefts, rights))
    ff = ft = tf = tt = 0
    for lkey, rkey in pairs:
        l, r = left[lkey], right[rkey]
        if l and r:
            tt += 1
        elif l:
            tf += 1
        elif r:
            ft += 1
        else:
            ff += 1
    return {"label": record["label"], "FF": ff, "FT": ft, "TF": tf, "TT": tt,
            "n_pairs": len(pairs)}


with C.logging(snakemake) as say:
    pl.style()
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    cohort = spec.get("cohort", "WT")
    records = _load(snakemake.input)
    say("panel {} (cohort {}): {} clones".format(
        spec["panel_id"], cohort, len(records)))

    if not records:
        pl.empty_panel(out_path, "no boundary_ctcf tables for cohort " + cohort)
        C.write_n_items(snakemake, 0)
        raise SystemExit(0)

    if cohort == "both_sides":
        # --- the proportion plot -----------------------------------------
        rows = [_both_sides_row(r) for r in records]
        from scipy.stats import binomtest

        labels, props, los, his, sigs = [], [], [], [], []
        for row in rows:
            non_ff = row["FT"] + row["TF"] + row["TT"]
            if not non_ff:
                continue
            prop = row["TT"] / non_ff
            se = np.sqrt(prop * (1 - prop) / non_ff)
            p = binomtest(row["TT"], non_ff, 0.5).pvalue
            labels.append(row["label"])
            props.append(prop)
            los.append(max(0.0, prop - 1.96 * se))
            his.append(min(1.0, prop + 1.96 * se))
            sigs.append(p < 0.05)
            say("  {:20s} FF={:4d} FT={:3d} TF={:3d} TT={:3d}  "
                "TT/(non-FF)={:.3f}  p={:.3g}".format(
                    row["label"], row["FF"], row["FT"], row["TF"], row["TT"],
                    prop, p))

        fig, ax = plt.subplots(figsize=(7, 4))
        x = np.arange(len(labels))
        ax.axhline(0.5, ls="--", color="grey", lw=0.8)
        ax.errorbar(x, props,
                    yerr=[np.array(props) - np.array(los),
                          np.array(his) - np.array(props)],
                    fmt="none", ecolor="grey", capsize=3, lw=1)
        ax.scatter(x, props,
                   c=["#2A9D8F" if s else "#E76F51" for s in sigs],
                   s=55, zorder=3)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=45, ha="right")
        ax.set_ylim(0, 1)
        ax.yaxis.set_major_formatter(
            plt.FuncFormatter(lambda v, _: "{:.0f}%".format(v * 100)))
        ax.set_ylabel("Proportion of valleys with a peak on both sides\n"
                      "(among those with a peak on >= 1 side)")
        ax.set_xlabel("Clone")
        ax.legend(handles=[Patch(color="#2A9D8F", label="p < 0.05"),
                           Patch(color="#E76F51", label="n.s.")],
                  title="binomial vs 50:50 (FF removed)", frameon=False,
                  loc="upper left")
        pl.despine(ax)
        ax.set_title("CTCF peak on both vs one side of a valley")
        fig.tight_layout()
        pl.save(fig, out_path, dpi=spec.get("dpi", 300))
        C.write_n_items(snakemake, len(labels))
        say("wrote {}".format(out_path))
        raise SystemExit(0)

    # --- the published pies ----------------------------------------------
    n = len(records)
    fig, axes = plt.subplots(1, n, figsize=(2.6 * n, 3.4))
    if n == 1:
        axes = [axes]

    for ax, rec in zip(axes, records):
        pos = rec["n_escaping_ctcf"]
        neg = rec["n_escaping"] - pos
        say("  {:20s} escaping={:4d}  CTCF+={:4d}  CTCF-={:4d}".format(
            rec["label"], rec["n_escaping"], pos, neg))
        if rec["n_escaping"] == 0:
            ax.text(0.5, 0.5, "no escaping\nboundaries", ha="center",
                    va="center", transform=ax.transAxes, fontsize=9)
            ax.axis("off")
        else:
            wedges, _, texts = ax.pie(
                [pos, neg],
                colors=[pl.CTCF_PIE_COLOR, pl.NO_CTCF_PIE_COLOR],
                autopct=pl.make_autopct(rec["n_escaping"]),
                startangle=90, counterclock=False, pctdistance=0.72,
                wedgeprops={"linewidth": 0.5, "edgecolor": "white"},
            )
            for t in texts:
                t.set_fontsize(8)
        ax.set_title("{}\nn = {}".format(rec["label"], rec["n_escaping"]),
                     fontsize=9)

    fig.legend(handles=[Patch(facecolor=pl.CTCF_PIE_COLOR,
                              label="CTCF <= 50 kb"),
                        Patch(facecolor=pl.NO_CTCF_PIE_COLOR,
                              label="No CTCF <= 50 kb")],
               loc="lower center", ncol=2, frameon=False, fontsize=8)
    fig.suptitle("CTCF status of escapee valley boundaries", fontsize=10)
    fig.tight_layout(rect=(0, 0.08, 1, 0.95))

    pl.save(fig, out_path, dpi=spec.get("dpi", 300))
    C.write_n_items(snakemake, len(records),
                    extra={r["label"]: r["n_escaping"] for r in records})
    say("wrote {}".format(out_path))
