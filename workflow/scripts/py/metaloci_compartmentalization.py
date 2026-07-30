"""Compartment-like strength from the METALoci `.mlo` objects.

Transcribed from `01_10_metaloci_compartments.ipynb` cells 1 and 4. One output
table per (run, region); the original was run 3 regions x 3 directories = 9
times by hand-editing two variables.

THE NUMBER THAT REACHES THE PAPER
---------------------------------
    compartmentalization (%) = (N_significant_q1 + N_significant_q3)
                               / N_total_bins * 100

with significance `LMI_pvalue < 0.05` -- **strictly less than**, whereas
`metaloci figure`'s own `moran_info.txt` counter uses `<=`. Both are reproduced
where they belong: `<` here, `<=` in `metaloci_lm.py`. They differ only for
bins whose stored p is exactly 0.05, which float16 storage makes possible.

THE MERGE
---------
    pd.merge(mlobject['lmi_info'][signal],
             mlobject['lmi_geometry'],
             on=["bin_index", "moran_index"], how="inner")

then wrapped as a GeoDataFrame on its `geometry` column. Identical to
`metaloci/tools/figure.py`.

PER-QUADRANT POLYGON STATISTICS (`stats_polygons`)
--------------------------------------------------
For each quadrant's significant bins: `unary_union` of the Voronoi cells,
`polygonize` the result, drop polygons with `area <= 1e-2`, then report
`num_pols`, `mean_distance` (mean pairwise centroid distance, upper triangle
only) and `mean_area`.

Three quirks of the original, all reproduced because they are in the published
numbers:

* the `< 2 polygons` branch returns `num_pols = 1` **even when there are
  zero polygons**, and `mean_distance = nan`;
* with zero polygons `np.mean([])` is `nan` (with a RuntimeWarning), which is
  why the ground-truth tables have blank `mean_area` cells;
* centroids are computed from the *unfiltered* polygon list in the original and
  then discarded -- the filtered list is what everything uses. Harmless, kept
  out.

`pearsonr` PRECISION -- DO NOT "FIX" THIS
-----------------------------------------
`ZSig` / `ZLag` are stored as `np.half` (float16) by METALoci
(`spatial_stats/lmi.py`). Under numpy 1.26 value-based casting, `scipy.stats.
pearsonr` computes in the input dtype, so **`r` comes back as a float16** while
`p` is float64. That is exactly the ground truth's signature: `pearsonr` values
carry 3-5 significant digits (`0.2383`, `0.11786`, `0.0002823`) while
`pearsonp` carries 17 (`1.9426774879526065e-06`). Verified empirically in
`workflow/envs/metaloci.yaml`'s pins (numpy 1.26.4 / scipy 1.16.1 / pandas
2.0.3): `pearsonr(float16, float16) -> (np.float16, np.float64)`.

Casting the inputs to float64 "for accuracy" would change every `pearsonr`
value in the table. Do not do it. (Under numpy >= 2 the NEP-50 promotion rules
would silently give float64 instead -- another reason the env pin is load-
bearing rather than decorative.)

EXCLUSIONS (cell 4, verbatim)
-----------------------------
    if ('C5C10' in ksignal and signal != 'RNA-Seq') \
       or (signal == "Rad21" and 'B1621' not in ksignal) \
       or (signal == "CTCF" and "B1621" in ksignal):
        continue

applied upstream by `SS.metaloci_signals_for()`, so this script never sees an
excluded pair.

ROW ORDER
---------
The original iterated `os.listdir()` (arbitrary) x a hard-coded signal list.
Here it is `sorted(datasets)` x the config signal order, which is deterministic
but will NOT match the ground truth's row order. Compare on the
`(dataset, signal, merge)` key.
"""

import os
import sys
import warnings

import geopandas as gpd
import numpy as np
import pandas as pd
from scipy.spatial import distance
from scipy.stats import pearsonr
from shapely.ops import polygonize, unary_union

# jobs: [(dataset, mlo_path, [(mark, metaloci_signal_name), ...]), ...]
jobs = list(snakemake.params.jobs)
pvalue = float(snakemake.params.pvalue)
min_area = float(snakemake.params.min_polygon_area)
roi = snakemake.params.roi
run = snakemake.params.run
out_tsv = snakemake.output.tsv

QUADRANTS = (1, 2, 3, 4)

os.makedirs(os.path.dirname(snakemake.log[0]), exist_ok=True)
log = open(snakemake.log[0], "w")


def say(msg):
    log.write(f"{msg}\n")


def fail(msg):
    say(f"ERROR {msg}")
    log.close()
    sys.exit(f"metaloci_compartmentalization[{run}/{roi}]: {msg}")


def stats_polygons(gdf):
    """`01_10` cell 1, verbatim -- including the `num_pols = 1` on an empty set."""
    all_signi = unary_union(gdf.geometry)
    polygons = [poly for poly in polygonize(all_signi)]
    polygons = [poly for poly in polygons if poly.area > min_area]

    if len(polygons) < 2:
        with warnings.catch_warnings():
            # np.mean([]) on the empty case: the original emitted this warning
            # and wrote nan. Keep the nan, drop the noise.
            warnings.simplefilter("ignore", category=RuntimeWarning)
            mean_area = np.mean([poly.area for poly in polygons])
        return 1, np.nan, mean_area

    distances = distance.pdist([poly.centroid.coords[0] for poly in polygons])
    distances = distance.squareform(distances)
    distances = distances[np.triu_indices(len(polygons), k=1)]
    return len(polygons), distances.mean(), np.mean([poly.area for poly in polygons])


records = []
for dataset, mlo_path, signal_names in jobs:
    if not os.path.exists(mlo_path):
        fail(f"{mlo_path} does not exist")
    mlobject = pd.read_pickle(mlo_path)
    lmi_info = mlobject.get("lmi_info") or {}
    if mlobject.get("lmi_geometry") is None:
        fail(f"{mlo_path} has no lmi_geometry -- metaloci lm did not run")

    # `signal` is the bare mark (the table's `signal` column); `name` -- the
    # original's `ksignal` -- is the dataset-qualified key METALoci stored
    # under, i.e. `{mark}_{dataset_adj}`.
    for signal, name in signal_names:
        if name not in lmi_info:
            fail(
                f"{dataset}: signal {name!r} absent from {mlo_path}. "
                f"present: {sorted(lmi_info)}"
            )

        merged = pd.merge(
            lmi_info[name],
            mlobject["lmi_geometry"],
            on=["bin_index", "moran_index"],
            how="inner",
        )
        merged = gpd.GeoDataFrame(merged, geometry=merged.geometry)

        significant = {}
        for quadrant in QUADRANTS:
            significant[quadrant] = merged[
                (merged["moran_quadrant"] == quadrant) & (merged["LMI_pvalue"] < pvalue)
            ]
        non_significant = merged[merged["LMI_pvalue"] >= pvalue]

        compartmentalization = (
            (len(significant[1]) + len(significant[3])) / len(merged) * 100
        )

        poly_stats = {q: stats_polygons(significant[q]) for q in QUADRANTS}

        try:
            # float16 in, float16 r out. See the module docstring.
            r, p = pearsonr(merged["ZSig"], merged["ZLag"])
        except Exception:
            r, p = np.nan, np.nan

        row = {
            "dataset": dataset,
            "signal": signal,
            # `01_10`: "selected" if 'selected' in ksignal else "all". The
            # `_selected` variants were an abandoned branch; the live path is
            # always "all". Kept so the column exists with the right values.
            "merge": "selected" if "selected" in name else "all",
            "compartmentalization": compartmentalization,
            "pearsonr": r,
            "pearsonp": p,
            "significant_1": len(significant[1]),
            "significant_2": len(significant[2]),
            "significant_3": len(significant[3]),
            "significant_4": len(significant[4]),
            "non_significant": len(non_significant),
            "total_bins": len(merged),
        }
        for q in QUADRANTS:
            num_pols, mean_distance, mean_area = poly_stats[q]
            row[f"num_pols_{q}"] = num_pols
            row[f"mean_distance_{q}"] = mean_distance
            row[f"mean_area_{q}"] = mean_area
        records.append(row)
        say(
            f"{dataset}\t{signal}\tcomp={compartmentalization:.4f}\t"
            f"sq1={len(significant[1])} sq3={len(significant[3])} "
            f"n={len(merged)} r={r}"
        )

COLUMNS = [
    "dataset", "signal", "merge", "compartmentalization", "pearsonr", "pearsonp",
    "significant_1", "significant_2", "significant_3", "significant_4",
    "non_significant", "total_bins",
    "num_pols_1", "mean_distance_1", "mean_area_1",
    "num_pols_2", "mean_distance_2", "mean_area_2",
    "num_pols_3", "mean_distance_3", "mean_area_3",
    "num_pols_4", "mean_distance_4", "mean_area_4",
]

results = pd.DataFrame(records, columns=COLUMNS)
os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
results.to_csv(out_tsv, index=False, sep="\t")
say(f"wrote {len(results)} rows to {out_tsv}")
log.close()
