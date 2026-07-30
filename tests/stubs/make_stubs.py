#!/usr/bin/env python3
"""Synthetic feature files matching the real schemas, for exercising panels.

The visualisation stage was built against the frozen path contract in
``workflow/lib/paths.py`` while the feature modules were still landing, so it
needed data of the right SHAPE long before there was data of the right VALUE.
That is what this generates: files at exactly the paths ``paths.py`` returns,
with exactly the columns and ``.npz`` keys the panels read, filled with
plausible noise.

**These stubs prove structure, never values.** A panel that renders here is
wired correctly; whether it reproduces the paper is a question only real data
answers. The counts are deliberately set to the PUBLISHED ones (E6 375, C5 301,
B1 258, JTG 300, CL30 321; E6A7 110 boundaries 80/30; F3 36 boundaries 23/13;
EFig 6c 1/44/11 and 1/12/6) so that ``paper_numbers.tsv`` can be exercised
end to end, assertions included.

Usage::

    python tests/stubs/make_stubs.py --out /path/to/stub-data
    python tests/stubs/make_stubs.py --out ... --only gene_content,boundary_ctcf

Then run a panel against it with ``tests/stubs/run_panel.py``.
"""

from __future__ import annotations

import argparse
import gzip
import os
import sys

import numpy as np
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, "workflow"))

from lib.paths import Paths                                    # noqa: E402

RNG = np.random.default_rng(20260728)

WT_CLONES = ["E6", "C5", "B1", "JTG", "CL30"]
#: Published Fig 2g totals -- the GENE-FILTERED counts, not the raw HMM calls.
PUBLISHED_FILTERED = {"E6": 375, "C5": 301, "B1": 258, "JTG": 300, "CL30": 321}
PUBLISHED_RAW = {"E6": 377, "C5": 304, "B1": 261, "JTG": 302, "CL30": 324}
#: Published Fig 2h pie n-values, in clone order E6, C5, B1, JTG, CL30.
PUBLISHED_ESCAPING = {"E6": 110, "C5": 96, "B1": 74, "JTG": 76, "CL30": 76}
#: EFig 6b: (escaping boundaries, CTCF+, CTCF-)
DEGRON_BOUNDARIES = {"E6A7": (110, 80, 30), "F3": (36, 23, 13)}
#: EFig 6c: (dTAG-only, shared, NodTAG-only)
DEGRON_VENN = {"E6A7": (1, 44, 11), "F3": (1, 12, 6)}

DEGRON = {"E6A7": "CTCF", "F3": "CTCF", "B1621": "Rad21"}
SIGNALS = ["H3K27me3", "H3K27ac", "CTCF", "RNA-Seq"]
CHRX_LEN = 171031299
NBINS = 100


def wt_track(clone):
    return "H3K27me3_{}_WT_Xi".format(clone)


def degron_track(clone, which):
    return "H3K27me3_{}_{}-{}_Xi".format(clone, DEGRON[clone], which)


def _write(path, frame):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    frame.to_csv(path, sep="\t", index=False)
    return path


def _valleys(n, seed):
    """`n` non-overlapping 5 kb-quantised intervals on chrX."""
    rng = np.random.default_rng(seed)
    starts = np.sort(rng.choice(np.arange(0, CHRX_LEN // 5000 - 40), n,
                                replace=False)) * 5000
    widths = rng.integers(1, 12, size=n) * 5000
    return starts, starts + widths


# ---------------------------------------------------------------------------
# generators
# ---------------------------------------------------------------------------
def gene_content(P):
    """Fig 2g / 4d. Row counts ARE the published gene-filtered totals."""
    written = []
    plan = [(wt_track(c), PUBLISHED_FILTERED[c],
             PUBLISHED_ESCAPING[c] // 2) for c in WT_CLONES]
    for clone, (n_bnd, n_pos, _neg) in DEGRON_BOUNDARIES.items():
        plan.append((degron_track(clone, "NodTAG"), 300, n_bnd // 2))
        plan.append((degron_track(clone, "dTAG"), 260, n_bnd // 3))
    plan.append((degron_track("B1621", "NodTAG"), 284, 40))
    plan.append((degron_track("B1621", "dTAG"), 281, 38))

    for track, total, n_escaping in plan:
        starts, ends = _valleys(total, abs(hash(track)) % 10000)
        category = ["no-gene"] * total
        for i in range(min(n_escaping, total)):
            category[i] = "escaping-gene-valley"
        for i in range(n_escaping, min(n_escaping + total // 5, total)):
            category[i] = "silent-gene-valley"
        for i in range(n_escaping + total // 5,
                       min(n_escaping + total // 3, total)):
            category[i] = "no-expressed-gene-valley"
        frame = pd.DataFrame({
            "chrom": "chrX", "start": starts, "end": ends,
            "size": ends - starts,
            "has_escapee": [c == "escaping-gene-valley" for c in category],
            "has_silent": [c == "silent-gene-valley" for c in category],
            "has_gene": [c != "no-gene" for c in category],
            "category": category,
            "size_class": "medium",
        })
        written.append(_write(P.gene_content(track), frame))
    return written


def valleys(P):
    """Raw and gene-filtered BEDs. Raw counts ARE the published ones."""
    written = []
    plan = [(wt_track(c), PUBLISHED_RAW[c], PUBLISHED_FILTERED[c])
            for c in WT_CLONES]
    plan += [(degron_track("E6A7", "NodTAG"), 321, 318),
             (degron_track("E6A7", "dTAG"), 281, 279),
             (degron_track("F3", "NodTAG"), 351, 349),
             (degron_track("F3", "dTAG"), 205, 203),
             (degron_track("B1621", "NodTAG"), 284, 282),
             (degron_track("B1621", "dTAG"), 281, 279)]
    for track, raw, filtered in plan:
        starts, ends = _valleys(raw, abs(hash(track)) % 10000)
        frame = pd.DataFrame({"chrom": "chrX", "start": starts, "end": ends})
        for path in (P.valleys("chrX", track),):
            os.makedirs(os.path.dirname(path), exist_ok=True)
            frame.to_csv(path, sep="\t", index=False, header=False)
            written.append(path)
        path = P.valleys_filtered("chrX", track)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        frame.iloc[:filtered].to_csv(path, sep="\t", index=False, header=False)
        written.append(path)
    return written


def boundary_ctcf(P):
    """Fig 2h / EFig 2g. Escaping-boundary counts ARE the published ones."""
    from scipy.stats import fisher_exact
    written = []
    plan = [(wt_track(c), PUBLISHED_FILTERED[c], PUBLISHED_ESCAPING[c],
             int(PUBLISHED_ESCAPING[c] * 0.6)) for c in WT_CLONES]
    for clone, (n_bnd, n_pos, _neg) in DEGRON_BOUNDARIES.items():
        plan.append((degron_track(clone, "NodTAG"), 300, n_bnd, n_pos))
    plan.append((degron_track("B1621", "NodTAG"), 284, 60, 35))

    for track, n_valleys, n_escaping, n_escaping_ctcf in plan:
        n_total = 2 * n_valleys
        starts, ends = _valleys(n_valleys, abs(hash(track)) % 10000)
        rows = []
        escaping_left = n_escaping
        for start, end in zip(starts, ends):
            is_escaping = escaping_left > 0
            escaping_left -= 2 if is_escaping else 0
            for side, pos in (("L", start), ("R", end)):
                rows.append({
                    "chrom": "chrX", "start": int(pos), "end": int(pos) + 1,
                    "side": side,
                    "valley_class": ("escaping-gene-valley" if is_escaping
                                     else "no-gene"),
                })
        frame = pd.DataFrame(rows).iloc[:n_total]
        frame["is_escaping"] = frame["valley_class"] == "escaping-gene-valley"

        has_ctcf = np.zeros(len(frame), dtype=bool)
        escaping_idx = np.flatnonzero(frame["is_escaping"].to_numpy())
        other_idx = np.flatnonzero(~frame["is_escaping"].to_numpy())
        has_ctcf[escaping_idx[:n_escaping_ctcf]] = True
        has_ctcf[other_idx[:len(other_idx) // 4]] = True
        frame["has_ctcf"] = has_ctcf

        a = int((frame["is_escaping"] & frame["has_ctcf"]).sum())
        b = int((~frame["is_escaping"] & frame["has_ctcf"]).sum())
        c = int((frame["is_escaping"] & ~frame["has_ctcf"]).sum())
        d = int((~frame["is_escaping"] & ~frame["has_ctcf"]).sum())
        odds, p = fisher_exact([[a, b], [c, d]])
        frame["fisher_p"] = p
        frame["fisher_odds_ratio"] = odds
        frame["n_boundaries"] = len(frame)
        frame["n_escaping"] = int(frame["is_escaping"].sum())
        written.append(_write(P.boundary_ctcf(track), frame))
    return written


def boundary_profile(P):
    """EFig 2e/2f — NO background subtraction (R-2)."""
    written = []
    tracks = [wt_track(c) for c in WT_CLONES]
    tracks += [degron_track(c, "NodTAG") for c in DEGRON]
    for track in tracks:
        clone = track.split("_")[1]
        signals = [s for s in SIGNALS
                   if not (s == "CTCF" and clone == "B1621")]
        if clone == "B1621":
            signals.append("Rad21")
        ctcf = pd.read_csv(P.boundary_ctcf(track), sep="\t")
        n = len(ctcf)
        for signal in signals:
            x = np.linspace(-1, 1, NBINS)
            base = {"H3K27me3": 0.9, "H3K27ac": 0.10, "CTCF": 0.09,
                    "Rad21": 0.09, "RNA-Seq": 0.4}[signal]
            shape = base * (1.0 - 0.55 * np.exp(-(x ** 2) / 0.02))
            matrix = np.clip(
                shape[None, :] + RNG.normal(0, base * 0.25, (n, NBINS)), 0, None)
            path = P.boundary_profile(track, signal)
            os.makedirs(os.path.dirname(path), exist_ok=True)
            np.savez_compressed(
                path, matrix=matrix, background=np.zeros(0),
                matrix_final=matrix, profile=np.nanmean(matrix, axis=0),
                profile_raw=np.nanmean(matrix, axis=0),
                chroms=np.asarray(["chrX"] * n, dtype=str),
                positions=ctcf["start"].to_numpy(dtype=np.int64),
                sides=ctcf["side"].to_numpy(dtype=str),
                n_used=n, n_dropped=0)
            written.append(path)
    return written


def boundary_dtag(P):
    """Fig 4f/4g + EFig 6b."""
    written = []
    for clone in DEGRON:
        signals = [s for s in SIGNALS
                   if not (s == "CTCF" and clone == "B1621")]
        if clone == "B1621":
            signals.append("Rad21")
        n_all, n_ctcf = DEGRON_BOUNDARIES.get(clone, (60, 35, 25))[:2]
        for signal in signals:
            payload = {}
            base = {"H3K27me3": 0.9, "H3K27ac": 0.10, "CTCF": 0.09,
                    "Rad21": 0.09, "RNA-Seq": 0.4}[signal]
            x = np.linspace(-1, 1, NBINS)
            for subset, n in (("escaping", n_all), ("escaping_ctcf", n_ctcf)):
                for condition, scale in (("nodtag", 1.0), ("dtag", 0.65)):
                    shape = base * scale * (
                        1.0 - 0.55 * np.exp(-(x ** 2) / 0.02))
                    matrix = np.clip(
                        shape[None, :] + RNG.normal(0, base * 0.2, (n, NBINS)),
                        0, None)
                    payload["profile_{}_{}".format(subset, condition)] = \
                        np.nanmean(matrix, axis=0)
                    payload["matrix_{}_{}".format(subset, condition)] = matrix
                payload["n_{}".format(subset)] = n
                payload["positions_{}".format(subset)] = \
                    np.arange(n, dtype=np.int64) * 100000
                payload["sides_{}".format(subset)] = np.asarray(
                    ["L", "R"] * (n // 2 + 1), dtype=str)[:n]
            path = P.boundary_dtag(clone, signal)
            os.makedirs(os.path.dirname(path), exist_ok=True)
            np.savez_compressed(path, **payload)
            written.append(path)
    return written


def valley_overlap(P):
    """EFig 6c — cardinalities ARE the published ones."""
    written = []
    for clone, (dtag_only, shared, nodtag_only) in DEGRON_VENN.items():
        for scope in ("escaping", "all"):
            scale = 1 if scope == "escaping" else 5
            frame = pd.DataFrame([{
                "clone": clone, "scope": scope,
                "n_nodtag": (shared + nodtag_only) * scale,
                "n_dtag": (shared + dtag_only) * scale,
                "nodtag_overlap": shared * scale,
                "nodtag_only": nodtag_only * scale,
                "dtag_overlap": shared * scale,
                "dtag_only": dtag_only * scale,
                "pct_of_nodtag": 100.0 * shared / (shared + nodtag_only),
                "pct_of_dtag": 100.0 * shared / (shared + dtag_only),
            }])
            written.append(_write(P.valley_overlap(clone, scope), frame))
    for scope in ("escaping", "all"):
        frame = pd.DataFrame([{
            "clone": "B1621", "scope": scope, "n_nodtag": 40, "n_dtag": 35,
            "nodtag_overlap": 30, "nodtag_only": 10, "dtag_overlap": 30,
            "dtag_only": 5, "pct_of_nodtag": 75.0, "pct_of_dtag": 85.7,
        }])
        written.append(_write(P.valley_overlap("B1621", scope), frame))
    return written


def density(P):
    """EFig 6a — chr7 (control) and chrX, NodTAG vs dTAG."""
    written = []
    for clone in DEGRON:
        for allele in ("Xi", "Xa", "Gall"):
            for mask in ("antivalley", "allcoverage", "valley"):
                frames = []
                for chrom, n in (("chr7", 1450), ("chrX", 1700)):
                    nodtag = np.abs(RNG.lognormal(-0.6, 0.7, n))
                    factor = 0.55 if (chrom == "chrX" and allele == "Xi") else 1.0
                    dtag = np.abs(nodtag * factor
                                  * RNG.lognormal(0, 0.25, n))
                    with np.errstate(divide="ignore", invalid="ignore"):
                        lfc = np.log2(dtag / nodtag)
                    frames.append(pd.DataFrame({
                        "chrom": chrom,
                        "start": np.arange(n) * 100000,
                        "end": (np.arange(n) + 1) * 100000,
                        "nodtag": nodtag, "dtag": dtag, "log2_ratio": lfc,
                        "pearson_r": float(np.corrcoef(nodtag, dtag)[0, 1]),
                        "ccc": 0.8, "median_lfc": float(np.median(lfc)),
                    }))
                written.append(_write(P.density(clone, allele, "100kb", mask),
                                      pd.concat(frames, ignore_index=True)))
    return written


def valley_states(P):
    """The per-bin HMM state table behind the locus figures."""
    written = []
    for track in [wt_track(c) for c in WT_CLONES] + \
                 [degron_track(c, "NodTAG") for c in DEGRON]:
        n = CHRX_LEN // 5000 + 1
        start = np.arange(n) * 5000
        value = np.abs(RNG.lognormal(-0.4, 0.6, n))
        state = (value < np.quantile(value, 0.25)).astype(int)
        frame = pd.DataFrame({"chrom": "chrX", "start": start,
                              "end": start + 5000, "value": value,
                              "state": state})
        path = P.valley_states("chrX", track)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with gzip.open(path, "wt") as fh:
            frame.to_csv(fh, sep="\t", index=False)
        written.append(path)
    return written


def valley_sizes(P):
    rows = []
    for track in [wt_track(c) for c in WT_CLONES]:
        starts, ends = _valleys(PUBLISHED_RAW[track.split("_")[1]],
                                abs(hash(track)) % 10000)
        rows.append(pd.DataFrame({
            "track": track, "mark": "H3K27me3", "clone": track.split("_")[1],
            "condition": "WT", "allele": "Xi", "chrom": "chrX",
            "start": starts, "end": ends, "size": ends - starts,
            "size_class": "medium"}))
    return [_write(P.valley_sizes(), pd.concat(rows, ignore_index=True))]


def pileup_scores(P):
    """`scores.tsv` — the quantity Fig 3c and EFig 3a actually plot."""
    written = []
    rows = []
    for locus in ("Mecp2", "Jarid"):
        clones = (["E6", "B1", "C5", "CL30", "JTG"] if locus == "Mecp2"
                  else ["E6", "C5", "B1"])
        for clone in clones:
            for called_on in ("Xa", "Xi"):
                for allele, tag, strength in (("Xa", "G1", 1.55),
                                              ("Xi", "G2", 0.95)):
                    cooler = "{}_{}_WT_{}_{}".format(locus, clone, tag, allele)
                    rows.append({
                        "file": "{}_{}".format(called_on, cooler),
                        "cooler": cooler, "n_loops": 12,
                        "median": strength - 0.05,
                        "mean": strength + RNG.normal(0, 0.08),
                        "std": 0.25})
    written.append(_write(P.pileup_score("Xa_vs_Xi", "full"),
                          pd.DataFrame(rows)))

    rows = []
    for locus in ("Mecp2", "Jarid"):
        for clone, mark in (("E6A7", "CTCF"), ("F3", "CTCF"),
                            ("C5C10", "CTCF"), ("B1621", "Rad21")):
            for condition, strength in (("{}-NodTAG".format(mark), 1.6),
                                        ("{}-dTAG".format(mark), 1.1)):
                cooler = "{}_{}_{}_G1_Xa".format(locus, clone, condition)
                rows.append({"file": cooler, "cooler": cooler, "n_loops": 11,
                             "median": strength - 0.05,
                             "mean": strength + RNG.normal(0, 0.06),
                             "std": 0.22})
    written.append(_write(P.pileup_score("dTAG_vs_NodTAG", "full"),
                          pd.DataFrame(rows)))
    return written


def saddle_strength(P):
    """EFig 3b — one `file, value` row per cooler."""
    written = []
    for scope, comparison in (("refined_merged", "Xa_vs_Xi"),
                              ("refined", "dTAG_vs_NodTAG")):
        rows = []
        for locus in ("Mecp2", "Jarid"):
            if comparison == "Xa_vs_Xi":
                for clone in ["E6", "B1", "C5", "CL30", "JTG"]:
                    for allele, tag, value in (("Xa", "G1", 2.0),
                                               ("Xi", "G2", 2.5)):
                        rows.append({
                            "file": "{}_{}_WT_{}_{}".format(locus, clone, tag,
                                                            allele),
                            "value": value + RNG.normal(0, 0.15)})
            else:
                for clone, mark in (("E6A7", "CTCF"), ("F3", "CTCF"),
                                    ("B1621", "Rad21")):
                    for condition, value in (("{}-NodTAG".format(mark), 2.4),
                                             ("{}-dTAG".format(mark), 1.9)):
                        rows.append({
                            "file": "{}_{}_{}_G1_Xa".format(locus, clone,
                                                            condition),
                            "value": value + RNG.normal(0, 0.12)})
        written.append(_write(P.saddle_strength_selected(scope, "full"),
                              pd.DataFrame(rows)))
    return written


def _pileup_npz(path, peak):
    grid = np.linspace(-1, 1, 41)
    xx, yy = np.meshgrid(grid, grid)
    data = 1.0 + (peak - 1.0) * np.exp(-(xx ** 2 + yy ** 2) / 0.05)
    data = data * RNG.lognormal(0, 0.05, data.shape)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    np.savez_compressed(path, data=data, n_loops=12, flank=100000, n=12)
    return path


def _saddle_npz(path, n=38):
    ramp = np.linspace(-1, 1, n)
    saddle = np.exp(np.outer(ramp, ramp) * 0.8)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tl = float(saddle[:8, :8].mean())
    br = float(saddle[-8:, -8:].mean())
    tr = float(saddle[:8, -8:].mean())
    bl = float(saddle[-8:, :8].mean())
    np.savez_compressed(
        path, interaction_sum=saddle * 100,
        interaction_count=np.full_like(saddle, 100.0), saddle=saddle,
        n_bins=n, qrange=(0.025, 0.975),
        saddle_strength_profile=np.linspace(0, 1.2, n + 2),
        top_left=tl, bottom_right=br, top_right=tr, bottom_left=bl,
        score=float((tl + br) / (tr + bl)))
    return path


def degron_pileups_and_saddles(P):
    """The dTAG_vs_NodTAG matrices behind Fig 5b/5c/5g/5h and EFig 8c/d, 10c/d."""
    written = []
    for locus in ("Mecp2", "Jarid"):
        for clone, mark in (("E6A7", "CTCF"), ("F3", "CTCF"),
                            ("C5C10", "CTCF"), ("B1621", "Rad21")):
            for condition, peak in (("{}-NodTAG".format(mark), 1.6),
                                    ("{}-dTAG".format(mark), 1.15)):
                cooler = "{}_{}_{}_G1_Xa".format(locus, clone, condition)
                exp = "{}_{}_G1_Xa".format(locus, clone)
                written.append(_pileup_npz(
                    P.pileup("dTAG_vs_NodTAG", "full", exp, cooler), peak))
                written.append(_saddle_npz(P.saddle("full", cooler, "E1")))
                path = P.saddle_values("full", cooler)
                _write(path, pd.DataFrame([{"file": cooler, "E1": 2.4,
                                            "E2": 2.0, "E3": 1.7}]))
                written.append(path)
    return written


def saddles_and_pileups(P):
    """The .npz matrices behind Fig 3a/3b and the extras."""
    written = []
    for locus in ("Mecp2", "Jarid"):
        clones = ["E6", "B1", "C5", "CL30", "JTG"]
        for clone in clones:
            for called_on in ("Xa", "Xi"):
                for allele, tag in (("Xa", "G1"), ("Xi", "G2")):
                    cooler = "{}_{}_WT_{}_{}".format(locus, clone, tag, allele)
                    exp = "{}_{}_WT".format(locus, clone)
                    sample = "{}_{}".format(called_on, cooler)
                    grid = np.linspace(-1, 1, 41)
                    xx, yy = np.meshgrid(grid, grid)
                    peak = 1.6 if allele == "Xa" else 1.1
                    data = 1.0 + (peak - 1.0) * np.exp(-(xx ** 2 + yy ** 2) / 0.05)
                    data *= RNG.lognormal(0, 0.05, data.shape)
                    path = P.pileup("Xa_vs_Xi", "full", exp, sample)
                    os.makedirs(os.path.dirname(path), exist_ok=True)
                    np.savez_compressed(path, data=data, n_loops=12,
                                        flank=100000, n=12)
                    written.append(path)

            for allele, tag in (("Xa", "G1"), ("Xi", "G2")):
                cooler = "{}_{}_WT_{}_{}".format(locus, clone, tag, allele)
                n = 38
                ramp = np.linspace(-1, 1, n)
                saddle = np.exp(np.outer(ramp, ramp) * 0.8)
                isum = saddle * 100
                icount = np.full_like(saddle, 100.0)
                path = P.saddle("full", cooler, "E1")
                os.makedirs(os.path.dirname(path), exist_ok=True)
                np.savez_compressed(
                    path, interaction_sum=isum, interaction_count=icount,
                    saddle=saddle, n_bins=n, qrange=(0.025, 0.975),
                    saddle_strength_profile=np.linspace(0, 1.2, n + 2),
                    top_left=float(saddle[:8, :8].mean()),
                    bottom_right=float(saddle[-8:, -8:].mean()),
                    top_right=float(saddle[:8, -8:].mean()),
                    bottom_left=float(saddle[-8:, :8].mean()),
                    score=float((saddle[:8, :8].mean()
                                 + saddle[-8:, -8:].mean())
                                / (saddle[:8, -8:].mean()
                                   + saddle[-8:, :8].mean())))
                written.append(path)
                written.append(_write(
                    P.saddle_values("full", cooler),
                    pd.DataFrame([{"file": cooler, "E1": 2.3, "E2": 2.0,
                                   "E3": 1.8}])))
    return written


def eigenvectors(P):
    """`Comp_{name}.tsv` plus the orientation table the report surfaces."""
    written = []
    names = []
    for locus in ("Mecp2", "Jarid"):
        for clone in ("E6", "B1", "C5", "CL30", "JTG"):
            for allele, tag in (("Xa", "G1"), ("Xi", "G2")):
                names.append("{}_{}_WT_{}_{}".format(locus, clone, tag, allele))
    for name in names[:4]:
        n = 240
        start = np.arange(n) * 5000 + 73315000
        frame = pd.DataFrame({
            "chrom": "chrX", "start": start, "end": start + 5000,
            "E1": np.sin(np.linspace(0, 6, n)),
            "E2": np.cos(np.linspace(0, 4, n)) * 0.5,
            "E3": np.sin(np.linspace(0, 9, n)) * 0.3})
        written.append(_write(P.eigs("full", name), frame))
    orientation = pd.DataFrame([
        {"roi": "full", "name": n, "eigenvector": "E1",
         "gc_correlation": 0.62, "flipped": bool(i % 3 == 0)}
        for i, n in enumerate(names[:4])])
    written.append(_write(
        P.work("features", "compartments", "eigenvector_orientation.tsv"),
        orientation))
    return written


def valley_xa_xi(P):
    written = []
    for clone in WT_CLONES:
        track = wt_track(clone)
        n = 120
        starts, ends = _valleys(n, abs(hash(track)) % 9999)
        category = (["shared"] * (n // 2) + ["Xa_only"] * (n // 4)
                    + ["Xi_only"] * (n - n // 2 - n // 4))
        written.append(_write(P.valley_xa_xi(track), pd.DataFrame({
            "chrom": "chrX", "start": starts, "end": ends,
            "side": ["L", "R"] * (n // 2), "category": category})))
    return written


def stackups(P):
    """The background-SUBTRACTED matrices. NOT a paper panel (R-2)."""
    written = []
    for clone in WT_CLONES:
        track = wt_track(clone)
        for signal in ("H3K27me3", "H3K27ac", "CTCF"):
            for variant in ("all", "motif_yes", "motif_no"):
                n = 200
                x = np.linspace(-1, 1, NBINS)
                base = {"H3K27me3": 0.9, "H3K27ac": 0.10,
                        "CTCF": 0.09}[signal]
                shape = base * (1.0 - 0.5 * np.exp(-(x ** 2) / 0.02))
                matrix = np.clip(shape[None, :]
                                 + RNG.normal(0, base * 0.2, (n, NBINS)),
                                 0, None)
                background = np.full(NBINS, base * 0.1)
                path = P.stackup(track, signal, variant)
                os.makedirs(os.path.dirname(path), exist_ok=True)
                np.savez_compressed(
                    path, matrix=matrix, background=background,
                    matrix_final=matrix - background,
                    profile=np.nanmean(matrix - background, axis=0),
                    profile_raw=np.nanmean(matrix, axis=0),
                    chroms=np.asarray(["chrX"] * n, dtype=str),
                    positions=np.arange(n, dtype=np.int64) * 500000,
                    sides=np.asarray(["L", "R"] * (n // 2), dtype=str),
                    n_used=n, n_dropped=0)
                written.append(path)
    return written


def qc(P):
    """FRiP, peak overlap, correlation, fragment sizes, enrichment, csaw."""
    written = []
    samples = ["H3K27me3_{}_WT_Gall_rep{}".format(c, r)
               for c in WT_CLONES for r in (1, 2)]
    written.append(_write(P.qc("frip.tsv"), pd.DataFrame({
        "sample": samples, "frip": RNG.uniform(0.05, 0.45, len(samples))})))
    written.append(_write(P.qc("peak_overlap.tsv"), pd.DataFrame({
        "pair": ["{}_vs_{}".format(a, b)
                 for a, b in zip(samples[::2], samples[1::2])],
        "n_overlap": RNG.integers(2000, 9000, len(samples) // 2)})))
    for sample in samples:
        written.append(_write(P.qc("fragment_sizes", "{}.tsv".format(sample)),
                              pd.DataFrame({
                                  "size": np.arange(50, 500),
                                  "count": np.exp(-((np.arange(50, 500) - 180)
                                                    ** 2) / 4000) * 1e5})))
    n_bins, n_samples = 4000, len(samples)
    matrix = RNG.lognormal(0, 0.6, (n_bins, n_samples))
    path = P.qc("correlation", "WT_H3K27me3.npz")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    np.savez_compressed(path, matrix=matrix,
                        labels=np.asarray(samples, dtype=str))
    written.append(path)

    for mark, anchor in (("CTCF", "motif"), ("H3K27ac", "tss"),
                         ("H3K27me3", "tss")):
        x = np.linspace(-1, 1, NBINS)
        profile = np.exp(-(x ** 2) / 0.05)
        enrich = np.clip(profile[None, :] * RNG.lognormal(0, 0.5, (600, 1))
                         + RNG.normal(0, 0.05, (600, NBINS)), 0, None)
        path = P.qc("enrichment", "{}_{}.npz".format(mark, anchor))
        os.makedirs(os.path.dirname(path), exist_ok=True)
        np.savez_compressed(path, matrix=enrich)
        written.append(path)

    counts = pd.DataFrame(
        {s: RNG.poisson(120, 5000) for s in samples[:6]})
    counts.insert(0, "start", np.arange(5000) * 10000)
    path = P.scalefactors("WT_H3K27me3", "counts")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    counts.to_csv(path, sep="\t", index=False)
    written.append(path)
    written.append(_write(P.scalefactors("WT_H3K27me3", "factors"),
                          pd.DataFrame({"sample": samples[:6],
                                        "scale_factor": RNG.uniform(
                                            0.03, 0.07, 6)})))
    return written


def metaloci_pdfs(P):
    """Gaudi `_gtp` and violin PDFs, so the Fig 6 montage can be exercised.

    They are figures, so they live under results/figures/ exactly where the
    panels that draw them would put them -- panel_metaloci_composite montages
    other panels' output, and that is the one internal dependency in the
    visualisation stage.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    written = []
    layout_path = os.path.join(REPO, "resources", "fixtures",
                               "fig6_layout.yaml")
    import yaml
    with open(layout_path) as fh:
        layout = yaml.safe_load(fh)

    def _pdf(path, label):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        fig, ax = plt.subplots(figsize=(2.2, 2.2))
        theta = np.linspace(0, 2 * np.pi, 220)
        ax.scatter(np.cos(theta) * RNG.uniform(0.4, 1, 220),
                   np.sin(theta) * RNG.uniform(0.4, 1, 220),
                   c=RNG.integers(0, 4, 220), cmap="Set1", s=14)
        ax.set_title(label, fontsize=6)
        ax.axis("off")
        fig.savefig(path, bbox_inches="tight")
        plt.close(fig)
        return path

    slug_locus = {"kdm5c": "Jarid", "mecp2": "Mecp2"}
    for slug, columns in layout["columns"].items():
        locus = slug_locus[slug]
        for signal in layout["rows"]:
            for column in columns:
                if column == "violin":
                    path = P.figures("extras", "metaloci_violin_{}_{}".format(
                        locus, signal), "pdf")
                    written.append(_pdf(path, "violin {}".format(signal)))
                elif column == "Xa_consensus":
                    dataset = "{}_NodTAG-or-WT_Xa".format(locus)
                    path = P.figures("extras",
                                     "metaloci_gaudi_consensus_{}_{}".format(
                                         dataset, signal), "pdf")
                    written.append(_pdf(path, "{} Xa".format(signal)))
                else:
                    dataset = "{}_{}_WT_G1_Xi".format(locus, column)
                    path = P.figures("extras",
                                     "metaloci_gaudi_wt_{}_{}".format(
                                         dataset, signal), "pdf")
                    written.append(_pdf(path,
                                        "{} {}".format(signal, column)))
    return written


def allelic_ratio(P):
    written = []
    for locus in ("Mecp2", "Jarid"):
        frame = pd.DataFrame(
            {"B1": [0.807, 0.0424769728421052, 0.03],
             "C5": [3.158, 0.1857936804117647, 0.12],
             "CL30.7": [0.617, 0.0325232014210526, 0.02],
             "E6": [3.849, 0.2025864485789474, 0.18],
             "JTG": [2.346, 0.1234743895789474, 0.10]},
            index=["sum_all", "mean_all", "median_all"])
        path = P.allelic_ratio_stats(locus, "full")
        os.makedirs(os.path.dirname(path), exist_ok=True)
        frame.to_csv(path)
        written.append(path)
    return written


def loops(P):
    written = []
    counts = {"Jarid_Xa": (55, 14), "Jarid_Xi": (32, 13),
              "Mecp2_Xa": (41, 10), "Mecp2_Xi": (26, 12)}
    for name, (raw, refined) in counts.items():
        for path, n in ((P.loops_raw(name), raw),
                        (P.loops_refined(name), refined)):
            starts = np.sort(RNG.choice(np.arange(200), n, replace=False))
            frame = pd.DataFrame({
                "chrom1": "chrX", "start1": starts * 5000,
                "end1": starts * 5000 + 5000,
                "chrom2": "chrX", "start2": (starts + 20) * 5000,
                "end2": (starts + 20) * 5000 + 5000,
                "bin1": starts, "bin2": starts + 20, "score": 0.5})
            written.append(_write(path, frame))
    return written


def compartmentalization(P):
    rows = []
    for locus in ("Mecp2", "Jarid"):
        clones = (["E6", "B1", "C5", "CL30", "JTG"] if locus == "Mecp2"
                  else ["E6", "C5", "B1"])
        for clone in clones:
            for allele, tag, base in (("Xa", "G1", 22.0), ("Xi", "G2", 31.0)):
                for signal in ["CTCF", "H3K27me3", "H3K27ac", "AcMe3",
                               "RNA-Seq", "Rad21"]:
                    rows.append({
                        "dataset": "{}_{}_WT_{}_{}".format(locus, clone, tag,
                                                           allele),
                        "signal": signal,
                        "compartmentalization": base + RNG.normal(0, 3),
                        "pearsonr": 0.5, "pearsonp": 1e-17})
        for signal in ["CTCF", "H3K27me3", "H3K27ac", "AcMe3", "RNA-Seq",
                       "Rad21"]:
            rows.append({"dataset": "{}_NodTAG-or-WT_Xa".format(locus),
                         "signal": signal,
                         "compartmentalization": 20.0 + RNG.normal(0, 2),
                         "pearsonr": 0.5, "pearsonp": 1e-17})
    frame = pd.DataFrame(rows)
    written = [_write(P.ml_compartmentalization_final(), frame)]
    for run in ("wt", "degron", "consensus"):
        written.append(_write(P.ml_compartmentalization(run, "full"), frame))
    return written


GENERATORS = {
    "gene_content": gene_content,
    "valleys": valleys,
    "boundary_ctcf": boundary_ctcf,
    "boundary_profile": boundary_profile,
    "boundary_dtag": boundary_dtag,
    "valley_overlap": valley_overlap,
    "density": density,
    "valley_states": valley_states,
    "valley_sizes": valley_sizes,
    "pileup_scores": pileup_scores,
    "saddle_strength": saddle_strength,
    "saddles_and_pileups": saddles_and_pileups,
    "degron_pileups_and_saddles": degron_pileups_and_saddles,
    "eigenvectors": eigenvectors,
    "valley_xa_xi": valley_xa_xi,
    "stackups": stackups,
    "qc": qc,
    "metaloci_pdfs": metaloci_pdfs,
    "allelic_ratio": allelic_ratio,
    "loops": loops,
    "compartmentalization": compartmentalization,
}


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", required=True,
                    help="stub data_dir; files land at the real paths under it")
    ap.add_argument("--only", default="",
                    help="comma-separated subset of: " + ", ".join(GENERATORS))
    args = ap.parse_args()

    config = {"paths": {
        "data_dir": os.path.abspath(args.out),
        "resources_dir": "{data_dir}/resources",
        "results_dir": "{data_dir}/results",
        "log_dir": "{data_dir}/logs",
        "tmpdir": "{data_dir}/tmp",
    }}
    P = Paths(config)

    wanted = ([g.strip() for g in args.only.split(",") if g.strip()]
              or list(GENERATORS))
    total = 0
    for name in wanted:
        if name not in GENERATORS:
            raise SystemExit("unknown generator {!r}; known: {}".format(
                name, ", ".join(GENERATORS)))
        written = GENERATORS[name](P)
        total += len(written)
        print("{:24s} {:4d} files".format(name, len(written)))
    print("\n{} stub files under {}".format(total, P.data))
    print("These prove STRUCTURE, not values. A panel that renders here is "
          "wired correctly;\nwhether it reproduces the paper is a question "
          "only real data answers.")


if __name__ == "__main__":
    main()
