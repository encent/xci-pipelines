"""Broadcast each Gall scale factor onto that sample's Xa and Xi splits.

Why this exists at all
----------------------
csaw/edgeR computes one normalisation factor per *library*. A library is
sequenced once and then split by SNP into Gall (unsplit), Xa and Xi. All three
files therefore share a single scale factor, computed on the Gall BAM.

The original bug (decision D-03)
--------------------------------
``12_normalization_*.r`` broadcast the factors positionally::

    bam_df_1$scaleFactor <- gall_factors[rep(1:nrow(gall_factors), times = 3)]

over an *alphabetically sorted* ``list.files()`` result. That is correct only
when the full file list is exactly three times the Gall list and sorts into
three aligned blocks. The WT H3K27ac ``selected`` run drops one replicate
(CL30 rep2), leaving 9 Gall BAMs and 27 files that do NOT align -- so 5 of the
27 bigWigs were written with another sample's scale factor. It propagates into
the merged tracks, the AcMe3 ratio and the METALoci consensus.

This was verified numerically on the ground truth: the minimum positive value
of a bigWig on an autosome is exactly the applied scale factor, and
``H3K27ac_CL30_WT_Xa_rep1.bw`` carries 0.0425394 -- E6_rep1's factor.

Two code paths
--------------
``legacy: false`` (default)
    Join on ``(mark, clone, condition, replicate)``. Correct, and order-independent.

``legacy: true``
    Reproduce the positional broadcast bit-for-bit, so the testing team can
    demonstrate equivalence with the published files.
"""

import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))
from workflow.lib import naming  # noqa: E402

JOIN_KEYS = ["mark", "clone", "condition", "replicate"]

factors_path = snakemake.input.factors
sheet_path = snakemake.input.sheet
out_path = snakemake.output[0]
log_path = snakemake.log[0]
legacy = bool(snakemake.params.legacy)
normgroup = snakemake.params.normgroup

os.makedirs(os.path.dirname(out_path), exist_ok=True)
os.makedirs(os.path.dirname(log_path), exist_ok=True)


def _decompose(sample_ids: pd.Series) -> pd.DataFrame:
    rows = [naming.parse_sample(s) for s in sample_ids]
    out = pd.DataFrame(rows, index=sample_ids.index)
    return out.rename(columns={"rep": "replicate"})


with open(log_path, "w") as log:
    # `comment="#"` is required, not cosmetic. This file has two producers with
    # different formats: `csaw_norm_factors.R` writes a bare table, but
    # `recover_scale_factors.py` (the D-09 fixture path) writes an 8-line `#`
    # provenance preamble. Without this the fixture path dies on the preamble
    # with "Expected 1 fields in line 9" -- and the fixture path is the ONLY way
    # to test anything below `bam_coverage` on chrX, so it took out all of
    # phase 4 while `compute` mode looked fine.
    factors = pd.read_csv(
        factors_path, sep="\t", dtype={"sample": str}, comment="#"
    )

    # The two producers also disagree on the column NAME:
    #   csaw_norm_factors.R      -> scaleFactor   (camelCase, plus LibSize/NormFactor)
    #   recover_scale_factors.py -> scale_factor  (snake_case; LibSize/NormFactor
    #                               are deliberately EMPTY, because only their
    #                               product is recoverable from a finished bigWig)
    # Accept either and normalise, rather than making one producer win — both
    # formats are correct for their own purpose.
    if "scaleFactor" not in factors.columns:
        if "scale_factor" in factors.columns:
            factors = factors.rename(columns={"scale_factor": "scaleFactor"})
        else:
            raise SystemExit(
                f"{factors_path} has no scale-factor column.\n"
                f"Found: {list(factors.columns)}\n"
                "Expected `scaleFactor` (from csaw_norm_factors.R) or "
                "`scale_factor` (from recover_scale_factors.py)."
            )
    if factors["scaleFactor"].isna().any():
        bad = factors.loc[factors["scaleFactor"].isna(), "sample"].tolist()
        raise SystemExit(
            f"{factors_path} has empty scale factors for: {bad}\n"
            "A fixture row with no recovered value cannot be used; re-run "
            "recover_scale_factors.py for those samples."
        )

    sheet = pd.read_csv(sheet_path, sep="\t", dtype=str, comment="#")
    sheet["use"] = sheet["use"].astype(str).str.lower().isin(("true", "1", "yes"))
    sheet = sheet[sheet["use"] & (sheet["normgroup"] == normgroup)].copy()
    if sheet.empty:
        raise SystemExit(f"no samples in sheet for normgroup {normgroup!r}")
    sheet["replicate"] = sheet["replicate"].astype(int)

    gall = factors.copy()
    gall_meta = _decompose(gall["sample"])
    gall = pd.concat([gall, gall_meta[["allele"] + JOIN_KEYS]], axis=1)

    # The two producers also differ in WHICH ROWS they contain:
    #   csaw_norm_factors.R      -> Gall only (csaw counts the unsplit BAMs)
    #   recover_scale_factors.py -> every allele, since it reads every bigWig
    # So the fixture has three rows per (mark, clone, condition, replicate) and
    # the join below would not be many-to-one. Collapse to the per-library
    # factor, preferring the Gall row as authoritative.
    #
    # Where a fixture's Xa/Xi factor DISAGREES with its Gall factor, that is not
    # noise -- it is the D-03 mis-pairing recorded in the ground truth itself.
    # Report it rather than silently discarding it.
    if (gall["allele"] != "Gall").any():
        disagreements = []
        for keys, grp in gall.groupby(JOIN_KEYS):
            ref = grp.loc[grp["allele"] == "Gall", "scaleFactor"]
            if ref.empty:
                continue
            ref = ref.iloc[0]
            for _, r in grp[grp["allele"] != "Gall"].iterrows():
                if r["scaleFactor"] != ref:
                    disagreements.append(
                        f"    {r['sample']}  {r['scaleFactor']}  "
                        f"(Gall sibling: {ref})"
                    )
        if disagreements:
            log.write(
                f"NOTE: {len(disagreements)} allelic row(s) in {factors_path} "
                "carry a scale factor differing from their Gall sibling:\n"
                + "\n".join(disagreements)
                + "\n"
                "This is a SYMPTOM of the D-03 mis-pairing, but it is NOT the\n"
                "mis-pairing count and must not be reported as one. Two of the\n"
                "five genuinely mis-assigned WT_H3K27ac files are themselves\n"
                "`Gall` files, so a sibling comparison cannot see a corrupt\n"
                "*reference* -- and it additionally flags correct allelic files\n"
                "whose reference is one of those corrupt Gall rows. On\n"
                "WT_H3K27ac that inflates 5 genuine to 7 reported.\n"
                "The authoritative count is derived from the original R source\n"
                "(`12_normalization_H3K27ac_WT_compos.r` lines 114/136) and is\n"
                "asserted by tests/test_d03_mispairing_count.py in the test\n"
                "repo. the ground-truth inventory.md 1.6 lists the five files.\n"
                "The Gall value is used here, which is the correct per-library\n"
                "factor. Set legacy.h3k27ac_positional_scalefactors: true to\n"
                "reproduce the original mis-pairing instead.\n"
            )

        has_gall = set(map(tuple, gall.loc[gall["allele"] == "Gall", JOIN_KEYS].values))
        keep = gall["allele"] == "Gall"
        for idx, row in gall[~keep].iterrows():
            if tuple(row[JOIN_KEYS]) not in has_gall:
                # no Gall row for this library; fall back to an allelic one
                keep.loc[idx] = True
                has_gall.add(tuple(row[JOIN_KEYS]))
                log.write(
                    f"WARNING: no Gall row for {tuple(row[JOIN_KEYS])}; using "
                    f"{row['sample']}'s factor as the library factor.\n"
                )
        gall = gall[keep].reset_index(drop=True)

    # BUG-B: the flag is named for H3K27ac but WAS applied to every normgroup.
    # The original's positional broadcast only ever mis-paired the WT H3K27ac
    # `selected` run -- that is the ONLY live run with an odd Gall count. Every
    # other normgroup ran the same R code with an EVEN count, where
    # pairwise-then-triple is the identity and no mis-pairing occurs.
    #
    # Applying it globally silently corrupted H3K27me3, CTCF, Rad21 and RNA-Seq
    # -- 6/30 on each WT non-H3K27ac group and 8/12 on each dTAG group against
    # 30/30 and 12/12 -- and therefore every valley, stackup and METALoci signal
    # downstream. Scoped rather than renamed, so the blast radius matches the
    # defect being reproduced.
    legacy_normgroups = set(
        getattr(snakemake.params, "legacy_normgroups", ["WT_H3K27ac"])
    )
    if legacy and normgroup not in legacy_normgroups:
        log.write(
            f"legacy flag is ON but {normgroup} is not in legacy_normgroups "
            f"({sorted(legacy_normgroups)}), so the CORRECT name-based join is "
            "used here. The original only mis-paired the WT H3K27ac `selected` "
            "run; every other normgroup had an even Gall count, where the "
            "positional broadcast is the identity.\n"
        )
        legacy = False

    # BUG-C: in `fixture` mode the Gall row was taken as authoritative -- but
    # TWO of the nine WT H3K27ac Gall rows are themselves mis-written in the
    # ground truth, so the corrected path inherited two wrong library factors
    # and broadcast each to three files: 21/27 instead of 27/27.
    #
    # A library's three allele rows (Gall/Xa/Xi) should all carry one factor, so
    # take the MAJORITY over the three rather than trusting Gall. Where the
    # original mis-paired exactly one of the three, the other two agree and win.
    # This also retires the earlier claim that a Gall-sibling comparison cannot
    # detect a corrupt reference: a majority vote can.
    def _library_factor(rows):
        vals = [v for v in rows["scaleFactor"].tolist() if pd.notna(v)]
        if not vals:
            return None
        counts = {}
        for v in vals:
            counts[v] = counts.get(v, 0) + 1
        best = max(counts.items(), key=lambda kv: (kv[1], -vals.index(kv[0])))
        if best[1] == 1 and len(vals) > 1:
            log.write(
                f"WARNING: no majority among {vals} for this library; the three "
                "allele rows disagree pairwise. Falling back to the Gall row.\n"
            )
            return None
        return best[0]

    if legacy:
        # ---- the original's positional broadcast, reproduced exactly --------
        gall_sorted = gall.sort_values("sample").reset_index(drop=True)
        all_sorted = sheet.sort_values("sample").reset_index(drop=True)

        n_gall = len(gall_sorted)
        n_all = len(all_sorted)

        # R, 12_normalization_H3K27ac_WT_compos.r:114 --
        #   split(bam_df, (seq_len(nrow(bam_df)) - 1) %/% 2)   consecutive PAIRS
        #   lapply(..., function(x) x[rep(1:nrow(x), times = 3), ])  triple EACH
        #
        # PAIRWISE-THEN-TRIPLE, not triple-the-whole-list. The two are the same
        # permutation only when the Gall count is EVEN -- and an ODD count is
        # precisely the condition that produced the original defect, because the
        # `selected` H3K27ac run drops CL30 rep2 and leaves 9.
        #
        # An earlier version of this branch did `list(range(n_gall)) * 3`, the
        # global broadcast. It had never been executed, and it reproduced 2 of
        # 27 instead of 27/27.
        index = []
        for lo in range(0, n_gall, 2):
            index += list(range(lo, min(lo + 2, n_gall))) * 3

        if n_all != len(index):
            log.write(
                f"WARNING: legacy broadcast is misaligned for {normgroup}: "
                f"{n_all} files vs {n_gall} Gall x 3 = {len(index)}. "
                "This is the original bug being reproduced on purpose (D-03).\n"
            )
        index = (index + [index[-1]] * n_all)[:n_all]

        all_sorted["scaleFactor"] = [
            gall_sorted["scaleFactor"].iloc[i] for i in index
        ]
        all_sorted["scaleFactor_source"] = [
            gall_sorted["sample"].iloc[i] for i in index
        ]
        merged = all_sorted

        wrong = merged[
            merged.apply(
                lambda r: naming.parse_sample(r["scaleFactor_source"])["clone"]
                != r["clone"],
                axis=1,
            )
        ]
        log.write(
            f"mode=legacy (D-03 bug reproduction)  normgroup={normgroup}\n"
            f"  {len(merged)} samples, {len(wrong)} carrying another sample's factor\n"
        )
        for _, r in wrong.iterrows():
            log.write(f"    {r['sample']}  <- {r['scaleFactor_source']}\n")

    else:
        # ---- correct: explicit join ----------------------------------------
        merged = sheet.merge(
            gall[JOIN_KEYS + ["scaleFactor", "sample"]].rename(
                columns={"sample": "scaleFactor_source"}
            ),
            on=JOIN_KEYS,
            how="left",
            validate="many_to_one",
        )

        missing = merged[merged["scaleFactor"].isna()]
        if not missing.empty:
            detail = "\n".join(
                f"    {r['sample']}  (mark={r['mark']} clone={r['clone']} "
                f"condition={r['condition']} rep={r['replicate']})"
                for _, r in missing.iterrows()
            )
            raise SystemExit(
                f"no Gall scale factor for {len(missing)} sample(s) in "
                f"normgroup {normgroup}:\n{detail}\n\n"
                "Every allelic split needs a Gall sibling with the same "
                "(mark, clone, condition, replicate). Check config/samples.tsv."
            )
        log.write(
            f"mode=join on {'+'.join(JOIN_KEYS)}  normgroup={normgroup}\n"
            f"  {len(merged)} samples resolved from {len(gall)} Gall libraries\n"
        )

    cols = ["sample", "mark", "clone", "condition", "allele", "replicate",
            "scaleFactor", "scaleFactor_source"]
    merged[cols].sort_values("sample").to_csv(out_path, sep="\t", index=False)
    log.write(f"wrote {out_path}\n")
