"""GC fraction per bin -- computed ONCE per distinct bin table.

`bioframe.frac_gc(bins[['chrom','start','end']], bioframe.load_fasta(mm10.fa))`
sees only the bin table and the FASTA. It never sees the cooler's contact
matrix, and it never sees the ROI. So every cooler that shares a bin table has
*the same* GC track, byte for byte.

For this dataset every cooler is chrX-only at 5 kb, so there is exactly one
distinct bin table and therefore exactly one distinct GC track. Verified: all
204 outputs -- and all 204 ground-truth files -- share the md5
`4ed6660643f07d805e74f4308c8623a3`.

The original recomputed it 204 times, reloading a 2.8 GB FASTA on each. That is
the same defect class as `index_bam`: recomputing something that already
exists, invisible to a correctness test because the answer is right, just paid
for 204 times.

**Why this is not simply hardcoded.** The identity is a property of *this
dataset*, not of the rule. Coolers at a different resolution, or spanning more
than one chromosome, genuinely produce different GC tracks, and a blind
copy would then be silently wrong. So the work is done once per bin table and
`gc_track` **asserts** that each consumer's bin table matches the canonical
one, failing loudly if it does not. The saving is kept; the assumption is
explicit and self-checking rather than latent.

The representative cooler is only a source of bins. The OUTPUT is keyed on the
bin table (`{chroms}_{resolution}`), never on a cooler name, so nothing reads
as "the GC track of cooler X".
"""

import os
import sys

import bioframe
import cooler


def main():
    smk = snakemake  # noqa: F821
    os.makedirs(os.path.dirname(smk.output.tsv), exist_ok=True)
    os.makedirs(os.path.dirname(smk.log[0]), exist_ok=True)

    uri = f"{smk.input.mcool}::resolutions/{smk.params.resolution}"
    clr = cooler.Cooler(uri)
    bins = clr.bins()[:][["chrom", "start", "end"]]

    genome = bioframe.load_fasta(smk.input.fasta)
    gc = bioframe.frac_gc(bins, genome)
    gc.to_csv(smk.output.tsv, index=False, sep="\t")

    # The bin table is the contract every consumer is checked against.
    sidecar = smk.output.tsv + ".bins"
    bins.to_csv(sidecar, index=False, sep="\t")

    with open(smk.log[0], "w") as fh:
        fh.write(
            f"bin key      : {smk.wildcards.binkey}\n"
            f"bins from    : {uri}  (representative only -- a source of bins)\n"
            f"chromosomes  : {', '.join(sorted(bins['chrom'].unique()))}\n"
            f"resolution   : {smk.params.resolution}\n"
            f"bins         : {len(bins)}\n"
            f"non-null GC  : {int(gc['GC'].notnull().sum())}\n"
            f"fasta        : {smk.input.fasta}\n"
            "\n"
            "Computed ONCE for this bin table. Every consumer asserts its own\n"
            "bin table matches the sidecar before reusing this file.\n"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
