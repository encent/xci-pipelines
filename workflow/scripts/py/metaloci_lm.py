"""Run `metaloci lm` for a dataset's three regions, then derive `moran_info.txt`.

    metaloci lm -w {WD}/{ds} -s {WD}/{ds}/{ds}.signals \
                -g chrX:{start}-{end}_0 -m -t {threads} -f

Defaults left unset, exactly as the original: 9999 permutations, p-value
threshold 0.05, quadrants [1,3], no `-b`, no `-po`, no `-a`, and **no `-i`**
(the original never wrote the per-signal `moran_data/` tree).

Regions are processed in the original's order -- `full`, `escape`,
`non_escape` -- because `moran_info.txt` is append-ordered.

`-f` IS SAFE HERE
-----------------
`metaloci lm -f` clears only `lmi_info` and `lmi_geometry`
(`metaloci/tools/ml.py`); it never touches the `kk_*` fields, so it does not
undo `metaloci_transplant_kk`. (`metaloci layout -f` DOES remove and re-create
the object, which is why the DAG orders layout strictly before the transplant.)

WHY THIS SCRIPT WRITES moran_info.txt ITSELF
============================================
`the design notes` section 2.2 lists `moran_info.txt` as the output of
`metaloci_lm`. In METALoci 1.3.2 it is not: `moran_info.txt` is appended by
**`metaloci figure`** (`tools/figure.py`, the block after the composite image is
saved), and `metaloci lm` writes only the `.mlo`. Since `metaloci figure` is a
figure step and belongs to `90_visualise.smk`, and since `moran_info.txt` is
pure tabular data that the testing team diffs, this script recomputes it from
the `.mlo` pickles using figure.py's arithmetic verbatim:

    r_value, p_value  = scipy.stats.linregress(lmi_info[sig].Sig,
                                               lmi_info[sig].Lag)[2:4]
    q[quadrant-1]    += 1                       for every row of lmi_info[sig]
    sq[quadrant-1]   += 1                       where LMI_pvalue <= 0.05
    columns: region symbol gene_id signal r_value p_value sq1 q1 .. sq4 q4

Two byte-level fidelity details, both reproduced deliberately:

* figure.py builds the region frame as
  `pd.DataFrame({"coords": [region], "symbol": ["symbol"], "id": ["id"]})` and
  then writes `row.name` -- the *pandas index* -- into the `symbol` column. So
  every ground-truth row has `symbol == 0` and `gene_id == "id"`. Not a bug we
  are free to fix; it is the column content.
* figure.py's header and row f-strings contain backslash line continuations
  whose leading indentation lands **inside** the string: 36 spaces before
  `\\tsq4` in the header, 24 spaces after `{r_value}` in every data row. Both
  are present in the ground-truth files (verified byte-for-byte with `cat -A`
  on `METALOCI_NEW_CONSENSUS/Mecp2_NodTAG-or-WT_Xi/moran_info.txt`). Removing
  them would make every line differ.

DETERMINISM (plan section 9.3)
------------------------------
`r_value`, `sq1..sq4` and `q1..q4` are deterministic. `p_value` here is the
linregress p of Sig vs Lag and is deterministic too. The stochastic quantity is
`LMI_pvalue` *inside* the `.mlo`, from the 9999-permutation Moran test -- which
moves `sq1..sq4` at the significance boundary. See the module docstring of
`32_metaloci.smk`.

A NOTE ON PRECISION, so nobody "fixes" it
-----------------------------------------
METALoci stores `Sig`, `Lag`, `ZSig`, `ZLag`, `LMI_score` and `LMI_pvalue` as
`np.half` (float16) to shrink the pickle (`spatial_stats/lmi.py`). `linregress`
promotes to float64, so `r_value`/`p_value` here are full precision -- matching
the 16-17 significant digits in the ground truth. (`pearsonr` in
`metaloci_compartmentalization.py` does NOT promote; see that file.)

GROUND-TRUTH COMPARISON CAVEAT
------------------------------
The original `moran_info.txt` files are **append-only across re-runs**. They
contain duplicate `(region, signal)` rows -- up to 4 copies -- and stale signal
names from abandoned variants (e.g. `H3K27ac_E6_WT_Xi_selected` in
`METALOCI_NEW/Mecp2_E6_WT_G1_Xi/moran_info.txt`, 48 data rows where 15 are
expected). This script writes the file fresh: exactly
`len(regions) x len(signals)` rows, no duplicates. Compare on the
`(region, signal)` key against the LAST occurrence in the ground truth, never
line by line.
"""

import os
import subprocess
import sys

import pandas as pd
from scipy.stats import linregress

work_dir = snakemake.params.work_dir
signals_file = snakemake.input.signals
regions = list(snakemake.params.regions)          # [(roi, start, end, mlo_path), ...]
chrom = snakemake.params.chrom
threads = int(snakemake.threads)
signipval = float(snakemake.params.pvalue)
permutations = int(snakemake.params.permutations)
quadrants = list(snakemake.params.quadrants)
multiprocess = bool(snakemake.params.multiprocess)
force = bool(snakemake.params.force)
dataset = snakemake.params.dataset
out_moran = snakemake.output.moran

os.makedirs(os.path.dirname(snakemake.log[0]), exist_ok=True)
log = open(snakemake.log[0], "w")

# figure.py's literal whitespace. See the docstring -- these are not cosmetic.
HEADER_PAD = " " * 36
ROW_PAD = " " * 24
HEADER = (
    "region\tsymbol\tgene_id\tsignal\tr_value\tp_value"
    "\tsq1\tq1\tsq2\tq2\tsq3\tq3" + HEADER_PAD + "\tsq4\tq4\n"
)


def say(msg):
    log.write(f"{msg}\n")
    log.flush()


def fail(msg):
    say(f"ERROR {msg}")
    log.close()
    sys.exit(f"metaloci_lm[{dataset}]: {msg}")


with open(signals_file) as fh:
    signal_names = [ln.strip() for ln in fh if ln.strip()]
if not signal_names:
    fail(f"{signals_file} is empty")
say(f"signals ({len(signal_names)}): {', '.join(signal_names)}")

# ---------------------------------------------------------------- run lm --
for roi, start, end, _mlo in regions:
    coords = f"{chrom}:{start}-{end}_0"
    cmd = [
        "metaloci", "lm",
        "-w", work_dir,
        "-s", signals_file,
        "-g", coords,
        "-p", str(permutations),
        "-v", str(signipval),
        "-t", str(threads),
    ]
    if multiprocess:
        cmd.append("-m")
    if force:
        cmd.append("-f")
    # `-q` takes nargs="+", so it goes last -- nothing may follow its values.
    # The original left it unset; 1 3 IS the METALoci default, so passing it
    # explicitly is documentation, not a behaviour change.
    cmd += ["-q", *[str(x) for x in quadrants]]
    say(f"[{roi}] {' '.join(cmd)}")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    log.write(proc.stdout)
    log.write(proc.stderr)
    log.flush()
    if proc.returncode != 0:
        fail(f"metaloci lm failed for {coords} (exit {proc.returncode})")

# --------------------------------------------------- derive moran_info.txt --
os.makedirs(os.path.dirname(out_moran), exist_ok=True)
rows = 0
with open(out_moran, "w") as out:
    out.write(HEADER)
    for roi, start, end, mlo_path in regions:
        coords = f"{chrom}:{start}-{end}_0"
        if not os.path.exists(mlo_path):
            fail(f"{mlo_path} missing after metaloci lm")
        mlobject = pd.read_pickle(mlo_path)
        lmi_info = mlobject.get("lmi_info") or {}
        if mlobject.get("lmi_geometry") is None:
            fail(f"{mlo_path} has no lmi_geometry -- metaloci lm did not run")

        missing = [s for s in signal_names if s not in lmi_info]
        if missing:
            # This is the "METALoci silently dropped a signal" failure mode:
            # a signal listed in .signals whose bedGraph basename did not match.
            fail(
                f"[{roi}] signals absent from lmi_info: {missing}. "
                f"present: {sorted(lmi_info)}. This means metaloci prep never "
                "saw a bedGraph with that exact basename."
            )

        for signal in signal_names:
            df = lmi_info[signal]
            _, _, r_value, p_value, _ = linregress(df["Sig"], df["Lag"])
            sq = [0] * 4
            q = [0] * 4
            for _, lmi_row in df.iterrows():
                q[int(lmi_row.moran_quadrant) - 1] += 1
                if lmi_row.LMI_pvalue <= signipval:
                    sq[int(lmi_row.moran_quadrant) - 1] += 1
            q_string = "\t".join(f"{sq[i]}\t{q[i]}" for i in range(4))
            out.write(
                f"{coords}\t0\tid\t{signal}\t{r_value}"
                + ROW_PAD
                + f"\t{p_value}\t{q_string}\n"
            )
            rows += 1
            say(
                f"[{roi}] {signal}: r={r_value:.6f} p={p_value:.3e} "
                f"sq={sq} q={q} n={len(df)}"
            )

say(f"wrote {rows} rows to {out_moran}")
log.close()
