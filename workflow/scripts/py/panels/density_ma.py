"""EFig 6a — H3K27me3 coverage, Ctrl vs dTAG, as MA and density scatters.

3 x 3 grid: rows are degron clones, columns are chr7 / chrX / both. The paper
uses the MA set; the density set is the same data plotted as y-vs-x.

MEDIUM CONFIDENCE. Twelve candidate panels exist across
allele {Xi, Xa, Gall} x mask {antivalley, allcoverage} x kind {density, MA},
and the archaeology could not tell which one is EFig 6a. All four Xi
candidates are rendered as `{id}__{candidate}`; a human picks against
`paper_figures/20260722_Extended_Data_Figure_6.pdf`.

Annotations, verbatim from `xci_density_plots.ipynb`:
  density panel  Pearson r, Lin's concordance correlation coefficient, and the
                 median log2 fold change over bins positive in both conditions
  MA panel       A = log2((x+y)/2), M = log2(y/x); a dashed line at M = 0 and a
                 red line at the median LFC
  colouring      log10 of a Gaussian KDE, subsampled to 5000 points with
                 `default_rng(42)` when there are more, points drawn in
                 increasing density order so the dense core is on top
  chr7 is the control chromosome: H3K27me3 there should not move on dTAG.
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

COLUMNS = [("chr7", "chr7 (control)"), ("chrX", "chrX"), (None, "both")]


def lin_ccc(x, y):
    """Lin's concordance correlation coefficient."""
    mx, my = x.mean(), y.mean()
    sx2, sy2 = x.var(), y.var()
    sxy = np.cov(x, y, ddof=0)[0, 1]
    return 2 * sxy / (sx2 + sy2 + (mx - my) ** 2)


def _density(x, y):
    from scipy import stats
    xy = np.vstack([x, y])
    if len(x) > 5000:
        rng = np.random.default_rng(42)
        sidx = rng.choice(len(x), 5000, replace=False)
        kde = stats.gaussian_kde(xy[:, sidx])
    else:
        kde = stats.gaussian_kde(xy)
    return kde(xy)


def density_scatter(ax, x, y, title="", xlabel="NodTAG (RPKM)",
                    ylabel="dTAG (RPKM)", cmap="jet", lim=None):
    from scipy import stats
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    if len(x) < 3:
        ax.set_title(title, fontsize=9, pad=3)
        ax.text(0.5, 0.5, "No data", ha="center", va="center",
                transform=ax.transAxes)
        return None
    d = _density(x, y)
    order = d.argsort()
    ax.scatter(x[order], y[order], c=np.log10(d[order] + 1e-10), cmap=cmap,
               s=1.5, rasterized=True, linewidths=0)
    r = stats.pearsonr(x, y)[0]
    ccc = lin_ccc(x, y)
    pos = (x > 0) & (y > 0)
    lfc = ("med LFC = {:+.2f}".format(np.median(np.log2(y[pos] / x[pos])))
           if pos.sum() else "med LFC = n/a")
    ax.text(0.97, 0.04, "r = {:.2f}  CCC = {:.2f}\n{}".format(r, ccc, lfc),
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7.5)
    if lim is None:
        lim = max(np.nanpercentile(x, 99.5), np.nanpercentile(y, 99.5)) * 1.08
    ax.plot([0, lim], [0, lim], "k--", lw=0.6, alpha=0.35)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(title, fontsize=9, pad=3)
    ax.set_xlabel(xlabel, fontsize=7)
    ax.set_ylabel(ylabel, fontsize=7)
    ax.tick_params(labelsize=7)
    return {"r": r, "ccc": ccc, "n": int(len(x))}


def ma_scatter(ax, x, y, title="", xlabel="log2 mean (RPKM)",
               ylabel="log2(dTAG / NodTAG)", cmap="jet"):
    mask = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    x, y = x[mask], y[mask]
    if len(x) < 3:
        ax.set_title(title, fontsize=9, pad=3)
        ax.text(0.5, 0.5, "No data", ha="center", va="center",
                transform=ax.transAxes)
        return None
    A = np.log2((x + y) / 2)
    M = np.log2(y / x)
    d = _density(A, M)
    order = d.argsort()
    ax.scatter(A[order], M[order], c=np.log10(d[order] + 1e-10), cmap=cmap,
               s=1.5, rasterized=True, linewidths=0)
    median_lfc = float(np.median(M))
    ax.axhline(0, color="k", lw=0.7, ls="--", alpha=0.4)
    ax.axhline(median_lfc, color="tab:red", lw=0.9, ls="-", alpha=0.8)
    ax.text(0.97, 0.04, "med LFC = {:+.2f}".format(median_lfc),
            transform=ax.transAxes, ha="right", va="bottom", fontsize=8,
            color="tab:red")
    ax.set_xlim(np.nanpercentile(A, 0.5), np.nanpercentile(A, 99.5))
    abs_m = np.nanpercentile(np.abs(M), 99.5)
    ax.set_ylim(-abs_m * 1.08, abs_m * 1.08)
    ax.set_title(title, fontsize=9, pad=3)
    ax.set_xlabel(xlabel, fontsize=7)
    ax.set_ylabel(ylabel, fontsize=7)
    ax.tick_params(labelsize=7)
    return {"median_lfc": median_lfc, "n": int(len(x))}


with C.logging(snakemake) as say:
    pl.style()
    import matplotlib.pyplot as plt

    candidate = spec.get("candidate")
    # candidate names are {mask}_{kind}; the primary file is the first one.
    default_mask = spec.get("mask", "antivalley")
    default_kind = "ma"
    if candidate:
        parts = candidate.split("_")
        default_mask = parts[0]
        default_kind = parts[-1]
    mask_wanted = default_mask
    kind = "ma" if default_kind.lower() in ("ma", "maplot") else "density"
    say("panel {}: mask={} kind={}".format(spec["panel_id"], mask_wanted, kind))

    tables = {}
    for path in sorted(map(str, snakemake.input)):
        stem = os.path.basename(path)[: -len(".tsv")]
        parts = stem.split("_")
        if len(parts) < 4:
            continue
        clone, allele, win, mask = parts[0], parts[1], parts[2], parts[3]
        if mask != mask_wanted:
            continue
        tables[clone] = C.read_table(path)

    if not tables:
        pl.empty_panel(out_path,
                       "no density tables with mask={}".format(mask_wanted))
        C.write_n_items(snakemake, 0)
        raise SystemExit(0)

    clones = sorted(tables)
    say("  clones: {}".format(", ".join(clones)))

    # A shared axis limit across every panel, as the original did in two passes.
    shared_lim = None
    if kind == "density":
        pooled = []
        for table in tables.values():
            pooled.append(np.asarray(table["nodtag"], dtype=float))
            pooled.append(np.asarray(table["dtag"], dtype=float))
        pooled = np.concatenate(pooled) if pooled else np.array([np.nan])
        finite = pooled[np.isfinite(pooled)]
        if finite.size:
            shared_lim = float(np.nanpercentile(finite, 99.5) * 1.08)
            say("  shared axis limit {:.4f} RPKM".format(shared_lim))

    fig, axes = plt.subplots(len(clones), len(COLUMNS),
                             figsize=(4 * len(COLUMNS), 4 * len(clones)),
                             squeeze=False)
    fig.suptitle(
        "H3K27me3 signal ({}): NodTAG vs dTAG\n"
        "{} bins in {} regions".format(
            spec.get("allele", "Xi"), spec.get("win", "100kb"), mask_wanted),
        fontsize=11,
    )

    stats_rows = {}
    for row, clone in enumerate(clones):
        table = tables[clone]
        for col, (chrom, col_label) in enumerate(COLUMNS):
            ax = axes[row][col]
            sub = table if chrom is None else table[table["chrom"] == chrom]
            x = np.asarray(sub["nodtag"], dtype=float)
            y = np.asarray(sub["dtag"], dtype=float)
            title = col_label if row == 0 else ""
            ylabel_prefix = "{}\n".format(clone) if col == 0 else ""
            if kind == "ma":
                info = ma_scatter(ax, x, y, title=title,
                                  ylabel=ylabel_prefix + "log2(dTAG / NodTAG)")
            else:
                info = density_scatter(ax, x, y, title=title,
                                       ylabel=ylabel_prefix + "dTAG (RPKM)",
                                       lim=shared_lim)
            if info:
                stats_rows["{}_{}".format(clone, col_label.split()[0])] = info
                say("    {:8s} {:14s} n={:6d} {}".format(
                    clone, col_label, info["n"],
                    " ".join("{}={:.4f}".format(k, v)
                             for k, v in info.items() if k != "n")))

    fig.tight_layout()
    pl.save(fig, out_path, dpi=spec.get("dpi", 300))
    C.write_n_items(snakemake, len(clones), extra={
        k: v.get("median_lfc", v.get("r")) for k, v in stats_rows.items()})
    say("wrote {}".format(out_path))
