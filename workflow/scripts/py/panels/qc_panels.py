"""Six QC / enrichment panels, dispatched on the producing rule.

    panel_ma_plot            csaw normalisation MA plots, per normgroup
    panel_frip               fraction of reads in peaks
    panel_peak_venn          replicate peak overlap
    panel_correlation        PCA + Spearman and Pearson heatmaps, per normgroup
    panel_fragment_sizes     fragment-size distributions
    panel_signal_enrichment  Fig 1e / Fig 2a / Fig 2d

One script rather than six because they share their whole structure and differ
only in which table they read; six near-identical files is how three of them
end up with different fonts.

`panel_signal_enrichment` is the odd one out: Fig 1e (CTCF), Fig 2a (H3K27ac)
and Fig 2d (H3K27me3) ARE published, which is why `signal_enrichment` in
50_qc.smk sits outside the `qc.enabled` gate while everything else here does
not.
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


def _ma_plot(say, paths, fig_maker):
    """MA plots from the csaw 10 kb bin-count matrix.

    A = log2 mean count across the pair, M = log2 ratio; one panel per sample
    pair against the reference sample, with the TMM scale factor drawn as a
    horizontal line. If the factors table is present the line is where the
    normalisation says M should sit.
    """
    import matplotlib.pyplot as plt

    counts_path = next((p for p in paths if ".counts." in os.path.basename(p)),
                       None)
    factors_path = next((p for p in paths if ".factors." in os.path.basename(p)),
                        None)
    if not counts_path:
        return None, "no csaw counts matrix"

    counts = C.read_table(counts_path)
    numeric = counts.select_dtypes(include=[np.number])
    samples = [c for c in numeric.columns
               if c.lower() not in ("start", "end", "width", "bin")]
    if len(samples) < 2:
        return None, "only {} sample column(s) in the count matrix".format(
            len(samples))

    factors = {}
    if factors_path and os.path.exists(factors_path):
        table = C.read_table(factors_path)
        cols = list(table.columns)
        if len(cols) >= 2:
            factors = {str(a): float(b)
                       for a, b in zip(table[cols[0]], table[cols[1]])
                       if np.isfinite(b)}

    reference = samples[0]
    others = samples[1:]
    ncol = min(3, len(others))
    nrow = int(np.ceil(len(others) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.6 * nrow),
                             squeeze=False)
    ref = numeric[reference].to_numpy(dtype=float)
    for i, sample in enumerate(others):
        ax = axes[i // ncol][i % ncol]
        y = numeric[sample].to_numpy(dtype=float)
        keep = (ref > 0) & (y > 0)
        A = np.log2((ref[keep] + y[keep]) / 2.0)
        M = np.log2(y[keep] / ref[keep])
        ax.scatter(A, M, s=1.5, alpha=0.3, linewidths=0, rasterized=True,
                   color="#4477aa")
        ax.axhline(0, color="k", lw=0.7, ls="--", alpha=0.5)
        if np.isfinite(np.median(M) if M.size else np.nan):
            ax.axhline(float(np.median(M)), color="tab:red", lw=0.9)
        if sample in factors and reference in factors:
            expected = np.log2(factors[sample] / factors[reference])
            ax.axhline(expected, color="tab:green", lw=0.9, ls=":")
        ax.set_title("{} vs {}\nn={} bins".format(sample, reference,
                                                  int(keep.sum())), fontsize=8)
        ax.set_xlabel("A = log2 mean count", fontsize=7)
        ax.set_ylabel("M = log2 ratio", fontsize=7)
        pl.despine(ax)
        say("  {:36s} n={:7d} median M={:+.4f}".format(
            sample, int(keep.sum()), float(np.median(M)) if M.size else np.nan))
    for j in range(len(others), nrow * ncol):
        axes[j // ncol][j % ncol].axis("off")
    fig.suptitle("csaw 10 kb bin counts - MA plots ({})".format(
        spec.get("normgroup", "")), fontsize=11)
    fig.tight_layout()
    return fig, len(others)


def _bar_from_table(say, path, value_col, label_col, ylabel, title):
    import matplotlib.pyplot as plt

    table = C.read_table(path)
    if table.empty:
        return None, "empty table " + os.path.basename(path)
    cols = list(table.columns)
    label_col = label_col if label_col in cols else cols[0]
    if value_col not in cols:
        numeric = table.select_dtypes(include=[np.number]).columns
        if not len(numeric):
            return None, "no numeric column in " + os.path.basename(path)
        value_col = numeric[0]

    labels = [str(x) for x in table[label_col]]
    values = table[value_col].astype(float).to_numpy()
    fig, ax = plt.subplots(figsize=(max(6, 0.32 * len(labels)), 4))
    ax.bar(np.arange(len(labels)), values, color="#4477aa")
    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, rotation=90, fontsize=6)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=10)
    pl.despine(ax)
    fig.tight_layout()
    for label, value in zip(labels, values):
        say("  {:44s} {:.4f}".format(label, value))
    return fig, len(labels)


def _correlation(say, path):
    """PCA scatter plus Spearman and Pearson heatmaps from the deepTools npz."""
    import matplotlib.pyplot as plt
    from scipy.stats import spearmanr

    npz = C.read_npz(path)
    key = next((k for k in ("matrix", "data", "arr_0") if k in npz), None)
    if key is None:
        return None, "no matrix in " + os.path.basename(path)
    matrix = np.asarray(npz[key], dtype=float)
    labels = ([str(x) for x in npz["labels"]] if "labels" in npz
              else ["s{}".format(i) for i in range(matrix.shape[1])])
    if matrix.shape[0] < matrix.shape[1]:
        matrix = matrix.T
    keep = np.all(np.isfinite(matrix), axis=1)
    matrix = matrix[keep]

    pearson = np.corrcoef(matrix, rowvar=False)
    spear = spearmanr(matrix).correlation
    spear = np.atleast_2d(spear)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    centred = matrix - matrix.mean(axis=0, keepdims=True)
    try:
        _, s, vt = np.linalg.svd(centred, full_matrices=False)
        pcs = vt[:2].T * s[:2]
        axes[0].scatter(pcs[:, 0], pcs[:, 1], s=60, color="#4477aa",
                        edgecolor="black")
        for (x, y), label in zip(pcs, labels):
            axes[0].text(x, y, label, fontsize=6, ha="center", va="bottom")
        axes[0].set_xlabel("PC1")
        axes[0].set_ylabel("PC2")
    except np.linalg.LinAlgError:
        axes[0].text(0.5, 0.5, "PCA did not converge", ha="center",
                     va="center", transform=axes[0].transAxes)
    axes[0].set_title("PCA of {} bins".format(matrix.shape[0]), fontsize=9)
    pl.despine(axes[0])

    for ax, data, name in ((axes[1], spear, "Spearman"),
                           (axes[2], pearson, "Pearson")):
        im = ax.imshow(data, cmap="RdYlBu_r", vmin=-1, vmax=1)
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=90, fontsize=6)
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels, fontsize=6)
        ax.set_title(name, fontsize=9)
        fig.colorbar(im, ax=ax, fraction=0.046)
    fig.tight_layout()
    say("  {} samples, {} finite bins".format(len(labels), matrix.shape[0]))
    return fig, len(labels)


def _fragment_sizes(say, paths):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9, 5))
    n = 0
    for path in sorted(paths):
        table = C.read_table(path)
        cols = list(table.columns)
        if len(cols) < 2:
            continue
        x = table[cols[0]].astype(float)
        y = table[cols[1]].astype(float)
        total = y.sum()
        ax.plot(x, y / total if total else y, lw=0.7, alpha=0.7,
                label=os.path.basename(path)[: -len(".tsv")])
        n += 1
    ax.set_xlabel("Fragment size (bp)")
    ax.set_ylabel("Fraction of fragments")
    ax.set_title("Fragment-size distributions ({} samples)".format(n),
                 fontsize=10)
    if n <= 20:
        ax.legend(fontsize=5, ncol=2)
    pl.despine(ax)
    fig.tight_layout()
    say("  {} samples".format(n))
    return fig, n


def _signal_enrichment(say, path, mark, anchor):
    """Fig 1e / 2a / 2d — mean profile plus the per-region heatmap."""
    import matplotlib.pyplot as plt

    npz = C.read_npz(path)
    key = next((k for k in ("matrix", "data", "arr_0", "stackup") if k in npz),
               None)
    if key is None:
        return None, "no matrix in " + os.path.basename(path)
    matrix = np.asarray(npz[key], dtype=float)
    if matrix.ndim == 3:
        matrix = matrix.reshape(-1, matrix.shape[-1])
    labels = ([str(x) for x in npz["labels"]] if "labels" in npz else None)

    fig, axes = plt.subplots(
        2, 1, figsize=(6, 8), gridspec_kw={"height_ratios": [1, 3]})
    with np.errstate(invalid="ignore"):
        profile = np.nanmean(matrix, axis=0)
    axes[0].plot(profile, color=pl.MARK_COLOR.get(mark, "#333333"), lw=1.6)
    axes[0].set_ylabel("mean {}".format(mark), fontsize=8)
    axes[0].set_xlim(0, len(profile) - 1)
    pl.boundary_xaxis(axes[0], nbins=len(profile), flank=200000)
    pl.despine(axes[0])

    order = np.argsort(-np.nan_to_num(matrix).mean(axis=1))
    im = axes[1].imshow(matrix[order], aspect="auto", cmap="Blues",
                        interpolation="none")
    axes[1].set_ylabel("{} regions".format(matrix.shape[0]), fontsize=8)
    pl.boundary_xaxis(axes[1], nbins=matrix.shape[1], flank=200000)
    fig.colorbar(im, ax=axes[1], fraction=0.046)
    fig.suptitle("{} enrichment at {} +/- 200 kb".format(mark, anchor),
                 fontsize=11)
    fig.tight_layout()
    say("  {} regions x {} bins".format(*matrix.shape))
    return fig, int(matrix.shape[0])


with C.logging(snakemake) as say:
    pl.style()
    producer = spec.get("producer", "")
    paths = [str(p) for p in snakemake.input]
    say("panel {} ({})".format(spec["panel_id"], producer))

    if producer == "panel_ma_plot":
        fig, n = _ma_plot(say, paths, None)
    elif producer == "panel_frip":
        fig, n = _bar_from_table(say, paths[0], "frip", "sample",
                                 "Fraction of reads in peaks",
                                 "FRiP by sample")
    elif producer == "panel_peak_venn":
        fig, n = _bar_from_table(say, paths[0], "n_overlap", "pair",
                                 "Overlapping peaks",
                                 "Peak overlap between replicates")
    elif producer == "panel_correlation":
        fig, n = _correlation(say, paths[0])
    elif producer == "panel_fragment_sizes":
        fig, n = _fragment_sizes(say, paths)
    elif producer == "panel_signal_enrichment":
        fig, n = _signal_enrichment(say, paths[0], spec.get("mark", "?"),
                                    spec.get("anchor", "tss"))
    else:
        fig, n = None, "qc_panels.py has no branch for " + producer

    if fig is None:
        pl.empty_panel(out_path, str(n))
        C.write_n_items(snakemake, 0)
        say("empty panel: {}".format(n))
    else:
        pl.save(fig, out_path, dpi=spec.get("dpi", 300))
        C.write_n_items(snakemake, n)
        say("wrote {}".format(out_path))
