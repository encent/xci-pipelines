"""Xa-vs-Xi valley overlap per clone — a secondary analysis, not in the paper.

DEVIATION D-4. The side vocabulary here is `L`/`R` everywhere; the original's
Xa/Xi overlap outputs used `left`/`right` and were not quantised to the 5 kb
grid. Normalise the vocabulary before diffing against ground truth.
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

    path = str(snakemake.input[0])
    table = C.read_table(path)
    track = spec.get("track", "?")
    say("panel {}: {}".format(spec["panel_id"], track))

    if table.empty or "category" not in table.columns:
        pl.empty_panel(out_path, "no Xa/Xi overlap categories for " + str(track))
        C.write_n_items(snakemake, 0)
        raise SystemExit(0)

    counts = table["category"].value_counts().to_dict()
    say("  " + "  ".join("{}={}".format(k, v) for k, v in sorted(counts.items())))

    def pick(*names):
        for name in names:
            if name in counts:
                return int(counts[name])
        return 0

    shared = pick("shared", "both", "common")
    xa_only = pick("Xa_only", "xa_only", "Xa-unique")
    xi_only = pick("Xi_only", "xi_only", "Xi-unique")
    if not (shared or xa_only or xi_only):
        keys = sorted(counts)
        xa_only = int(counts.get(keys[0], 0)) if keys else 0
        xi_only = int(counts.get(keys[1], 0)) if len(keys) > 1 else 0
        shared = int(counts.get(keys[2], 0)) if len(keys) > 2 else 0

    fig, ax = plt.subplots(figsize=(7, 7))
    pl.venn2_counts(ax, xa_only, shared, xi_only,
                    "Xa (n={})".format(xa_only + shared),
                    "Xi (n={})".format(xi_only + shared))
    ax.set_title("{}\nvalley overlap, Xa vs Xi".format(C.track_label(track)),
                 fontsize=11)
    fig.tight_layout()
    pl.save(fig, out_path, dpi=spec.get("dpi", 300))
    C.write_n_items(snakemake, xa_only + xi_only + shared,
                    extra={"Xa_only": xa_only, "shared": shared,
                           "Xi_only": xi_only})
    say("wrote {}".format(out_path))
