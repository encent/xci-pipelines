"""The `Stackups_*` composites — WITH random-shift background subtraction.

CORRECTION R-2, restated where it matters most. This panel reads
``P.stackup(track, signal, variant)``, whose matrices have a 100-iteration
random circular-shift background subtracted
(``numpy.random.default_rng(seed = clone_index + 42)``). The
``panel_boundary_analysis`` panels read ``P.boundary_profile`` instead, which
does no such thing. The two rules take the same boundary anchors and the same
20 bp bigWigs; only the arithmetic differs.

**No panel here is in the paper.** the archaeology notes section 5B: none of the
``Stackups_*`` products reached publication. The branch stays alive because it
also writes ``boundaries/{all,motif_yes,motif_no}``, which ARE ground truth.
If you are looking for the figure a reviewer asked about, it is
``boundary_analysis_{clone}``, not this.

Violin significance honours ``config.legacy.stackup_violin_test``. The original
``annotate_violin_mwu`` announced Mann-Whitney in the title and called
``scipy.stats.ttest_rel`` (deviation D-5). We do not port the bug: the default
is a real Mann-Whitney, and setting the option to ``ttest_rel`` reproduces the
original if a comparison needs it.
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

    track = spec.get("track", "?")
    signal = spec.get("signal", "?")
    variant = spec.get("variant", "all")
    test_name = str(spec.get("legacy", {}).get(
        "stackup_violin_test", "mannwhitneyu"))
    say("panel {}: {} / {} / {}  (violin test: {})".format(
        spec["panel_id"], track, signal, variant, test_name))

    npz = C.read_npz(str(snakemake.input[0]))
    matrix = np.asarray(npz["matrix_final"], dtype=float) \
        if "matrix_final" in npz else np.asarray(npz["matrix"], dtype=float)
    raw = np.asarray(npz["matrix"], dtype=float) if "matrix" in npz else matrix
    profile = np.asarray(npz["profile"], dtype=float) if "profile" in npz \
        else np.nanmean(matrix, axis=0)
    profile_raw = np.asarray(npz["profile_raw"], dtype=float) \
        if "profile_raw" in npz else np.nanmean(raw, axis=0)
    n = int(npz["n_used"]) if "n_used" in npz else matrix.shape[0]

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4))

    axes[0].plot(profile_raw, color=pl.OTHER_COLOR, lw=1.4,
                 label="raw (n={})".format(n))
    axes[0].plot(profile, color=pl.MARK_COLOR.get(signal, "#333333"), lw=2.0,
                 label="background-subtracted")
    pl.boundary_xaxis(axes[0], nbins=len(profile))
    axes[0].set_ylabel(signal)
    axes[0].set_title("mean profile", fontsize=9)
    axes[0].legend(fontsize=7)
    pl.despine(axes[0])

    order = np.argsort(-np.nan_to_num(matrix).mean(axis=1))
    im = axes[1].imshow(matrix[order], aspect="auto", cmap="Blues",
                        interpolation="none")
    pl.boundary_xaxis(axes[1], nbins=matrix.shape[1])
    axes[1].set_ylabel("{} boundaries".format(matrix.shape[0]), fontsize=8)
    axes[1].set_title("stackup ({})".format(variant), fontsize=9)
    fig.colorbar(im, ax=axes[1], fraction=0.046)

    # Flanks vs centre: the quantity the original's violins compared.
    half = matrix.shape[1] // 2
    edge = max(1, matrix.shape[1] // 10)
    with np.errstate(invalid="ignore"):
        centre = np.nanmean(matrix[:, half - edge:half + edge], axis=1)
        flanks = np.nanmean(
            np.hstack([matrix[:, :edge], matrix[:, -edge:]]), axis=1)
    good = np.isfinite(centre) & np.isfinite(flanks)
    pvalue = float("nan")
    if good.sum() >= 3:
        from scipy import stats
        if test_name == "ttest_rel":
            pvalue = float(stats.ttest_rel(centre[good], flanks[good]).pvalue)
        elif test_name == "wilcoxon":
            pvalue = float(stats.wilcoxon(centre[good], flanks[good]).pvalue)
        else:
            pvalue = float(stats.mannwhitneyu(
                centre[good], flanks[good], alternative="two-sided").pvalue)
    axes[2].violinplot([centre[good], flanks[good]], positions=[0, 1],
                       showmeans=True)
    axes[2].set_xticks([0, 1])
    axes[2].set_xticklabels(["boundary", "flanks"])
    axes[2].set_ylabel(signal)
    axes[2].set_title("{}  p = {:.3g}\n(n = {})".format(
        test_name, pvalue, int(good.sum())), fontsize=9)
    pl.despine(axes[2])

    fig.suptitle("{}  {}  ({})  -- background-subtracted, NOT a paper panel"
                 .format(C.track_label(track), signal, variant), fontsize=10)
    fig.tight_layout()
    pl.save(fig, out_path, dpi=spec.get("dpi", 300))
    C.write_n_items(snakemake, n, extra={"violin_test": test_name,
                                         "p_value": pvalue})
    say("  n = {}, {} p = {:.6g}".format(n, test_name, pvalue))
    say("wrote {}".format(out_path))
