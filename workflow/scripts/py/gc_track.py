"""GC phasing track for one (cooler, ROI) -- reused from the canonical file.

`bioframe.frac_gc` depends only on the bin table and the FASTA, never on the
contact matrix and never on the ROI. So this rule does not recompute it: it
asserts that this cooler's bin table is identical to the canonical one, then
reuses the canonical GC track.

The per-(roi, name) output paths are kept because they are the original's
layout and `eigs_cis` depends on them.

**The assertion is the point.** Reusing the canonical file is valid only while
every cooler shares one bin table -- true for this dataset (chrX only, 5 kb),
not true in general. A cooler at a different resolution, or spanning more than
one chromosome, has a genuinely different GC track. Rather than leave that as a
latent assumption behind a 200x speedup, this script fails loudly with the
actual mismatch, so a future dataset gets an error instead of quietly wrong
eigenvectors.
"""

import os
import shutil
import sys

import cooler
import pandas as pd


def main():
    smk = snakemake  # noqa: F821
    os.makedirs(os.path.dirname(smk.output.tsv), exist_ok=True)
    os.makedirs(os.path.dirname(smk.log[0]), exist_ok=True)

    uri = f"{smk.input.mcool}::resolutions/{smk.params.resolution}"
    mine = cooler.Cooler(uri).bins()[:][["chrom", "start", "end"]].reset_index(drop=True)
    canon = pd.read_csv(smk.input.bins, sep="\t").reset_index(drop=True)

    if not mine.equals(canon):
        detail = []
        if len(mine) != len(canon):
            detail.append(f"bin count: {len(mine)} here vs {len(canon)} canonical")
        mc, cc = sorted(mine["chrom"].unique()), sorted(canon["chrom"].unique())
        if mc != cc:
            detail.append(f"chromosomes: {mc} here vs {cc} canonical")
        for col in ("start", "end"):
            a, b = mine[col].to_numpy(), canon[col].to_numpy()
            n = min(len(a), len(b))
            bad = (a[:n] != b[:n]).nonzero()[0]
            if bad.size:
                i = int(bad[0])
                detail.append(f"first differing {col}: row {i}, {a[i]} here vs {b[i]} canonical")
        raise SystemExit(
            "\n".join(
                ["", "=" * 76,
                 f"ERROR  Bin table mismatch for cooler {smk.wildcards.name}",
                 "=" * 76, "",
                 "  The GC phasing track is shared across coolers because it depends",
                 "  only on the bin table and the genome FASTA -- not on the contact",
                 "  matrix. That sharing is valid only while every cooler has the SAME",
                 "  bin table, which is true for the published data (chrX only, 5 kb).",
                 "", "  This cooler's bin table differs:", ""]
                + [f"    {d}" for d in detail]
                + ["", "  What to do:", "",
                   "    Your coolers do not all share one binning, so each needs its own",
                   "    GC track. In workflow/rules/31_hic_features.smk widen the key in",
                   "    `_gc_bin_key()` so it distinguishes them -- it is currently",
                   "    {chromosomes}_{resolution}.",
                   "",
                   "    Do NOT delete this check. Without it the eigenvectors would be",
                   "    phased against the wrong GC track and every compartment call",
                   "    would be quietly wrong.",
                   "", "=" * 76, ""]
            )
        )

    shutil.copyfile(smk.input.gc, smk.output.tsv)

    with open(smk.log[0], "w") as fh:
        fh.write(
            f"{smk.wildcards.name} / {smk.wildcards.roi}: bin table matches the "
            f"canonical table ({len(mine)} bins, "
            f"{', '.join(sorted(mine['chrom'].unique()))} @ {smk.params.resolution} bp).\n"
            f"reused {smk.input.gc} -- the 2.8 GB FASTA was not reloaded.\n"
            "NOTE: content is ROI-independent; the {roi} directory only reproduces "
            "the original's layout.\n"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
