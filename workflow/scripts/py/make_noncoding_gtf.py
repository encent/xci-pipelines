"""Non-coding feature track for coolbox's GFFNC track (Fig 4e, EFig 6d).

Ports the FIRST BLOCK of `01_07_get_allelic_ratio_gtf.py` verbatim. That block
is the sole producer of `GRCm38.102_NC.gtf` in the original.

  gff3_features = ['lnc_RNA', 'ncRNA', 'scRNA', 'snoRNA', 'miRNA', 'snRNA',
                   'rRNA', 'pseudogene', 'ncRNA_gene']
  nc_elems = bioframe.read_table(gff_file, schema="gff") \
                     .query('feature in @gff3_features and chrom=="X"')
  ...
  nc_concat[list(nc_elems.columns)].to_csv(..., sep="\t", header=None,
      index=None, quoting=csv.QUOTE_NONE, escapechar="\\")

Three things here are load-bearing and were got wrong once already, so they are
called out explicitly:

1. **The source is the chrX GFF3, not the genome-wide GTF.**
   `Mus_musculus.GRCm38.102.chromosome.X.gff3`. An earlier implementation
   derived this file from `Mus_musculus.GRCm38.102.chr.gtf` by excluding
   `gene_biotype "protein_coding"`. That is a different file: different source,
   different feature vocabulary, genome-wide rather than chrX, and Ensembl
   chromosome naming. RD-3 caught the divergence while porting `allelic_gtf`;
   the original wins.

2. **The filter is a feature WHITELIST, not a biotype exclusion.** GFF3
   `feature` (column 3), not GTF `gene_biotype`. `pseudogene` and `ncRNA_gene`
   are in the whitelist even though a biotype-exclusion approach would treat
   them quite differently.

3. **The chromosome stays Ensembl-named `X`, not `chrX`.** The original filters
   `chrom == "X"` and never renames. coolbox is given this file as-is, so
   renaming here would silently empty the track.

Only the nine original GFF columns are written -- the parsed attributes are
computed (they force the dtype coercion below) and then dropped, exactly as the
original does by selecting `list(nc_elems.columns)`.
"""

import csv
import subprocess

import bioframe
import bioframe.sandbox.gtf_io

GFF3_FEATURES = [
    "lnc_RNA",
    "ncRNA",
    "scRNA",
    "snoRNA",
    "miRNA",
    "snRNA",
    "rRNA",
    "pseudogene",
    "ncRNA_gene",
]

gff_file = snakemake.input.gff3
out_gtf = snakemake.output.gtf
out_bgz = snakemake.output.bgz


def say(msg: str) -> None:
    print(msg, flush=True)


nc_elems = (
    bioframe.read_table(gff_file, schema="gff")
    .query("feature in @GFF3_FEATURES and chrom == 'X'")
    .reset_index(drop=True)
)
say(f"{len(nc_elems)} non-coding features on X from {gff_file}")

if nc_elems.empty:
    raise SystemExit(
        f"No non-coding features matched in {gff_file}.\n"
        "Expected an Ensembl GRCm38.102 chromosome-X GFF3 whose `chrom` column "
        "is 'X' (not 'chrX') and whose column 3 uses the GFF3 feature "
        "vocabulary. Check genome.gff3_x in config/config.yaml."
    )

# Parsing the attributes is what coerces start/end to int in the original.
nc_attr = bioframe.sandbox.gtf_io.parse_gtf_attributes(
    nc_elems["attributes"], kv_sep="=", item_sep=";"
)
import pandas as pd  # noqa: E402  (imported late, as the original does)

nc_concat = pd.concat([nc_elems, nc_attr], axis=1)
nc_concat["start"] = nc_concat["start"].astype(int)
nc_concat["end"] = nc_concat["end"].astype(int)

# Only the original GFF columns are written; the attribute columns are dropped.
nc_concat[list(nc_elems.columns)].to_csv(
    out_gtf,
    sep="\t",
    header=None,
    index=None,
    quoting=csv.QUOTE_NONE,
    escapechar="\\",
)
say(f"wrote {out_gtf}")

# bgzip + tabix for coolbox. The original did this outside the script; there is
# no record of the exact invocation, so this is the standard one.
with open(out_bgz, "wb") as fh:
    sort = subprocess.run(
        ["sort", "-k1,1", "-k4,4n", out_gtf], capture_output=True, check=True
    )
    bg = subprocess.run(["bgzip", "-c"], input=sort.stdout, capture_output=True, check=True)
    fh.write(bg.stdout)
subprocess.run(["tabix", "-p", "gff", "-f", out_bgz], check=True)
say(f"wrote {out_bgz} + .tbi")
