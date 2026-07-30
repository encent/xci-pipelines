"""EFig 6c — escapee-valley Venns, dTAG vs Ctrl, one per degron clone.

Published cardinalities (dTAG-only / shared / NodTAG-only):
E6A7 1 / 44 / 11, F3 1 / 12 / 6.

Built from the GENE-FILTERED valley sets in both conditions, restricted to
valleys carrying an escaping gene (`scope: escaping`). The `scope: all`
variant is the non-published `valley_overlap_dTAG_vs_NodTAG.svg`.

The "shared" count is asymmetric in the source data -- `valley_overlap` records
`nodtag_overlap` and `dtag_overlap` separately, because one NodTAG valley can
overlap two dTAG valleys and vice versa. A Venn cannot show that, so the
intersection is drawn as `nodtag_overlap` and both numbers go into the log and
into paper_numbers.tsv. If they differ, that is real merging/splitting, not an
error.
"""

import os
import sys

sys.path.insert(0, os.path.join(snakemake.params.spec["libdir"],
                                "scripts", "py", "panels"))
sys.path.insert(0, snakemake.params.spec["libdir"])

import _common as C                                            # noqa: E402
from lib import plotting as pl                                 # noqa: E402

spec = C.bootstrap(snakemake)
out_path = C.outputs(snakemake)


with C.logging(snakemake) as say:
    pl.style()
    import matplotlib.pyplot as plt

    scope = spec.get("scope", "escaping")
    say("panel {} (scope {})".format(spec["panel_id"], scope))

    records = []
    for path in sorted(map(str, snakemake.input)):
        table = C.read_table(path)
        if table.empty:
            continue
        row = table.iloc[0]
        records.append({
            "clone": str(row.get("clone", os.path.basename(path))),
            "n_nodtag": int(row["n_nodtag"]),
            "n_dtag": int(row["n_dtag"]),
            "nodtag_overlap": int(row["nodtag_overlap"]),
            "dtag_overlap": int(row["dtag_overlap"]),
            "nodtag_only": int(row["nodtag_only"]),
            "dtag_only": int(row["dtag_only"]),
        })
        say("  {:8s} dTAG-only {:3d} | shared {:3d} (dTAG side {:3d}) | "
            "NodTAG-only {:3d}".format(
                records[-1]["clone"], records[-1]["dtag_only"],
                records[-1]["nodtag_overlap"], records[-1]["dtag_overlap"],
                records[-1]["nodtag_only"]))

    if not records:
        pl.empty_panel(out_path, "no valley_overlap tables (scope " + scope + ")")
        C.write_n_items(snakemake, 0)
        raise SystemExit(0)

    n = len(records)
    fig, axes = plt.subplots(1, n, figsize=(7 * n, 7))
    if n == 1:
        axes = [axes]

    for ax, rec in zip(axes, records):
        pl.venn2_counts(
            ax,
            only_a=rec["dtag_only"],
            shared=rec["nodtag_overlap"],
            only_b=rec["nodtag_only"],
            label_a="dTAG (n={})".format(rec["n_dtag"]),
            label_b="NodTAG (n={})".format(rec["n_nodtag"]),
            colors=(pl.DTAG_LINE_COLOR, pl.NODTAG_LINE_COLOR),
        )
        ax.set_title("{}\n{} valleys, dTAG vs Ctrl".format(rec["clone"], scope),
                     fontsize=11)

    fig.tight_layout()
    pl.save(fig, out_path, dpi=spec.get("dpi", 300))
    C.write_n_items(snakemake, n, extra={
        "{}_{}".format(r["clone"], k): r[k]
        for r in records
        for k in ("dtag_only", "nodtag_overlap", "nodtag_only")
    })
    say("wrote {}".format(out_path))
