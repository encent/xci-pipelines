"""Valley size distributions — one histogram per clone.

100 bins, sizes below 500 kb (the original's cut for readability; the tail is
a handful of very long low-signal stretches). husl palette, one row per track.
Not a paper figure; it is the sanity check that the HMM has not collapsed into
one enormous valley, which is what a mis-scaled bigWig looks like.
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
MAX_SIZE = 500000

with C.logging(snakemake) as say:
    pl.style()
    import matplotlib.pyplot as plt

    table = C.read_table(str(snakemake.input[0]))
    path = str(snakemake.input[0])
    if table.empty:
        pl.empty_panel(out_path, "valley_sizes.tsv is empty")
        C.write_n_items(snakemake, 0)
        raise SystemExit(0)

    tracks = sorted(set(C.need(table, "track", path)))
    say("panel {}: {} tracks".format(spec["panel_id"], len(tracks)))

    try:
        import seaborn as sns
        palette = sns.color_palette("husl", len(tracks))
    except Exception:
        cmap = plt.get_cmap("hsv")
        palette = [cmap(i / max(1, len(tracks))) for i in range(len(tracks))]

    fig, axes = plt.subplots(len(tracks), 1, figsize=(10, 3 * len(tracks)),
                             squeeze=False)
    for ax_row, track, colour in zip(axes, tracks, palette):
        ax = ax_row[0]
        sizes = table.loc[table["track"] == track, "size"].astype(float)
        shown = sizes[sizes < MAX_SIZE]
        ax.hist(shown, bins=100, color=colour, alpha=0.7, edgecolor="black",
                linewidth=0.5)
        ax.set_xlabel("Valley size (bp)")
        ax.set_ylabel("Count")
        ax.set_title("{} (n={} valleys, {} shown below {} kb)".format(
            C.track_label(track), len(sizes), len(shown), MAX_SIZE // 1000))
        ax.grid(axis="y", alpha=0.3)
        say("  {:36s} n={:4d} median={:8.0f} max={:9.0f}".format(
            track, len(sizes),
            float(np.median(sizes)) if len(sizes) else 0.0,
            float(sizes.max()) if len(sizes) else 0.0))

    fig.suptitle("Distribution of valley sizes by clone", fontsize=14, y=1.0)
    fig.tight_layout()
    pl.save(fig, out_path, dpi=spec.get("dpi", 300))
    C.write_n_items(snakemake, len(tracks))
    say("wrote {}".format(out_path))
