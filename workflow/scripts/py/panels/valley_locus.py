"""Jarid / Mecp2 valley tracks — the valley caller's diagnostic figures.

Line plot of the 5 kb H3K27me3 signal across the locus window with the called
valley bins shaded red (`fill_between`, alpha 0.2). figsize (12, 4), dpi 300,
PNG. Windows come from `config.loci.*.view_window`:

    Jarid  chrX:150,400,000-152,350,000
    Mecp2  chrX: 73,315,000- 74,475,000

Not a paper figure, but the one that makes a moved valley boundary visible.
The per-bin state table (`P.valley_states`) exists for exactly this: the BED
alone says where the valleys ended up, the states say what the HMM saw.

DEVIATION D-2. The pipeline emits 36 degron figures where the ground truth has
34. The original wrapped `bioframe.merge` in a bare `try/except ... continue`,
so `H3K27me3_F3_CTCF-NodTAG_Gall` silently produced nothing. We do not swallow
it: an empty valley set gives a valid figure with no shading. Score the extras
as new artefacts, not failures.
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


with C.logging(snakemake) as say:
    pl.style()
    import matplotlib.pyplot as plt

    locus = spec.get("locus", "Mecp2")
    track = spec.get("track", "?")
    windows = dict(snakemake.params.windows)
    start, end = (int(v) for v in windows[locus])
    display = C.display_locus(locus, spec.get("display_names"))
    say("panel {}: {} at {} chrX:{:,}-{:,}".format(
        spec["panel_id"], track, display, start, end))

    states = C.read_table(str(snakemake.input[0]))
    path = str(snakemake.input[0])
    chrom = C.need(states, "chrom", path)
    bin_start = C.need(states, "start", path)
    bin_end = C.need(states, "end", path)
    value = C.need(states, "value", path).astype(float)
    state = C.need(states, "state", path)

    window = states[(chrom == "chrX") & (bin_end > start) & (bin_start < end)]
    say("  {} bins in the window".format(len(window)))

    fig, ax = plt.subplots(figsize=(12, 4))
    if len(window):
        x = (window["start"].values + window["end"].values) / 2.0
        y = window["value"].values.astype(float)
        ax.plot(x, y, lw=0.8, color="#333333")

        # The valley state is the LOW-mean state; call_valleys records it as
        # the argmin of the fitted means, and bins outside a training interval
        # keep -1 and can never be valleys.
        in_state = window["state"].values
        valley_state = None
        for candidate in sorted(set(int(s) for s in in_state if int(s) >= 0)):
            mean = np.nanmean(y[in_state == candidate])
            if valley_state is None or mean < valley_state[1]:
                valley_state = (candidate, mean)
        if valley_state is not None:
            mask = in_state == valley_state[0]
            ax.fill_between(x, 0, y, where=mask, color="red", alpha=0.2,
                            step=None, label="valley bins")
            say("  valley state {} ({} bins, mean {:.4f})".format(
                valley_state[0], int(mask.sum()), valley_state[1]))
            ax.legend(loc="upper right", fontsize=8)
    else:
        ax.text(0.5, 0.5, "no bins in the window", ha="center", va="center",
                transform=ax.transAxes)

    ax.set_xlim(start, end)
    ax.set_xlabel("chrX position (bp)")
    ax.set_ylabel("H3K27me3 (5 kb bins)")
    ax.set_title("{}  -  {}  chrX:{:,}-{:,}".format(
        C.track_label(track), display, start, end), fontsize=10)
    pl.despine(ax)
    fig.tight_layout()

    pl.save(fig, out_path, dpi=spec.get("dpi", 300))
    C.write_n_items(snakemake, len(window))
    say("wrote {}".format(out_path))
