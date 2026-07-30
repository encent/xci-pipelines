"""Valley overlap between the dTAG and NodTAG condition of one clone -- EFig 6c.

From ``xci_valley_overlap.ipynb``. Two things about this comparison are easy to
get wrong and both change the Venn cardinalities:

1. It runs on the GENE-FILTERED valley sets, not the raw ones.
2. It counts at the REGION level, not the fragment level. ``bioframe.overlap``
   returns one row per overlapping pair, so a NodTAG valley touching three dTAG
   valleys contributes three rows. The original takes ``.nunique()`` on the
   left index::

       ov = bioframe.overlap(nodtag, dtag, how='left', return_index=True)
       n_overlap = ov.loc[ov['index_'].notna(), 'index'].nunique()

   which is NodTAG-centric and guarantees ``nodtag_only + overlap == n_nodtag``
   exactly. The dTAG side is counted the same way with the frames swapped, so
   the two "overlap" numbers need NOT be equal -- that asymmetry is real and the
   Venn draws the NodTAG-centric one.

``scope``:
    ``all``       every valley in both filtered sets.
    ``escaping``  only valleys overlapping an escapee gene. EFig 6c reports
                  E6A7 (1 / 44 / 11) and F3 (1 / 12 / 6), i.e. 55 escaping
                  NodTAG valleys for E6A7 and 18 for F3 -- which is also the
                  110 and 36 boundaries of EFig 6b.

                  Each condition is classified with ITS OWN escapee delivery
                  (`E6A7_NodTAG` vs `E6A7_4DdTAG`). Escape status is precisely
                  what the degron perturbs, so using one list for both sides
                  would compare a set against itself under a different name.
                  `xci_valley_overlap.ipynb` cell 6 uses `bed_stem_dtag` and
                  `bed_stem_nodtag` separately.

Emits a one-row summary table::

    clone scope n_nodtag n_dtag nodtag_overlap nodtag_only dtag_overlap
    dtag_only pct_of_nodtag pct_of_dtag
"""

import os
import sys

import bioframe
import pandas as pd

log_path = snakemake.log[0]
os.makedirs(os.path.dirname(log_path), exist_ok=True)
os.makedirs(os.path.dirname(snakemake.output.tsv), exist_ok=True)

clone = snakemake.wildcards.clone
scope = snakemake.params.scope

with open(log_path, "w") as log:

    def say(msg: str) -> None:
        log.write(msg + "\n")
        log.flush()
        print(msg, file=sys.stderr)

    nodtag = bioframe.read_table(snakemake.input.nodtag, schema="bed3")
    dtag = bioframe.read_table(snakemake.input.dtag, schema="bed3")
    nodtag = nodtag[["chrom", "start", "end"]].reset_index(drop=True)
    dtag = dtag[["chrom", "start", "end"]].reset_index(drop=True)
    say(f"{clone} {scope}: {len(nodtag)} NodTAG valleys, {len(dtag)} dTAG valleys")

    if scope == "escaping":

        def read_escapees(path: str) -> pd.DataFrame:
            return pd.read_csv(
                path,
                sep="\t",
                header=None,
                names=["chrom", "start", "end", "gene_name", "score", "strand"],
            )[["chrom", "start", "end"]]

        # One escapee delivery per CONDITION, never one for both.
        escapees_nodtag = read_escapees(snakemake.input.escapees)
        escapees_dtag = read_escapees(snakemake.input.escapees_dtag)
        say(
            f"escapee genes: {len(escapees_nodtag)} NodTAG "
            f"({os.path.basename(snakemake.input.escapees)}), "
            f"{len(escapees_dtag)} dTAG "
            f"({os.path.basename(snakemake.input.escapees_dtag)})"
        )

        def escaping_only(valleys: pd.DataFrame, escapees: pd.DataFrame) -> pd.DataFrame:
            if len(valleys) == 0 or len(escapees) == 0:
                return valleys.iloc[:0]
            hits = bioframe.overlap(valleys, escapees, how="left", return_index=True)
            keep = sorted(set(hits.loc[hits["index_"].notna(), "index"]))
            return valleys.loc[keep].reset_index(drop=True)

        nodtag = escaping_only(nodtag, escapees_nodtag)
        dtag = escaping_only(dtag, escapees_dtag)
        say(f"  -> {len(nodtag)} NodTAG, {len(dtag)} dTAG escaping valleys")

    def region_level_overlap(left: pd.DataFrame, right: pd.DataFrame) -> int:
        if len(left) == 0 or len(right) == 0:
            return 0
        ov = bioframe.overlap(left, right, how="left", return_index=True)
        return int(ov.loc[ov["index_"].notna(), "index"].nunique())

    nodtag_overlap = region_level_overlap(nodtag, dtag)
    dtag_overlap = region_level_overlap(dtag, nodtag)

    summary = pd.DataFrame(
        [
            {
                "clone": clone,
                "scope": scope,
                "n_nodtag": len(nodtag),
                "n_dtag": len(dtag),
                "nodtag_overlap": nodtag_overlap,
                "nodtag_only": len(nodtag) - nodtag_overlap,
                "dtag_overlap": dtag_overlap,
                "dtag_only": len(dtag) - dtag_overlap,
                "pct_of_nodtag": (
                    100.0 * nodtag_overlap / len(nodtag) if len(nodtag) else 0.0
                ),
                "pct_of_dtag": (
                    100.0 * dtag_overlap / len(dtag) if len(dtag) else 0.0
                ),
            }
        ]
    )
    say(summary.to_string(index=False))

    summary.to_csv(snakemake.output.tsv, sep="\t", index=False)
