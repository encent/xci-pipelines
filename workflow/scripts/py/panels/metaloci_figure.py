"""`metaloci figure` — the external tool that only emits images.

    metaloci figure -w {dir} -s {ds}.signals -g chrX:{start}-{end}_0 \
                    -m -t {threads} -e

Runs in its own conda env (`workflow/envs/metaloci.yaml`) because METALoci's
pinned stack cannot coexist with the main environment. It is one of exactly two
places in the pipeline where an external tool draws rather than computes; the
other is coolbox. Both are invoked FROM the visualisation stage, which is what
keeps "figures only come from 90_visualise.smk" true rather than aspirational.

`metaloci figure` writes a whole family per signal -- `_hic`, `_signal`, `_kk`,
`_lmi`, `_gsp`, `_gtp`. `_gtp` is the Gaudi TYPE plot (LMI quadrant colouring),
and that is the Fig 6 panel type: `panel_metaloci_composite` montages exactly
those files.

`moran_info.txt` is NOT an output here even though `metaloci figure` writes it.
Plan section 2.2 attributes it to `metaloci lm`; the archaeology notes section D.3 has it
right that `figure` produces it. Since it is DATA, `metaloci_lm` recomputes it
upstream from the `.mlo` with `figure.py`'s arithmetic verbatim, and this rule
merely consumes it. Two rules declaring one output would be a DAG error, and
the layer split -- lm is data, figure is a picture -- is the point.

`-e` exports intermediates; `-m` is multiprocess. Note the argparse trap that
bit the original elsewhere: `-p -l 15` is plot-then-persistence-length, while
`-pl 15` parses as `-p` plus `-l 15` and silently loses the KK and
mixed-matrix PDFs. That flag pair belongs to `metaloci layout`, not here, but
the same parser is in play.
"""

import glob
import os
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.join(snakemake.params.spec["libdir"],
                                "scripts", "py", "panels"))
sys.path.insert(0, snakemake.params.spec["libdir"])

import _common as C                                            # noqa: E402
from lib import plotting as pl                                 # noqa: E402

spec = C.bootstrap(snakemake)
out_path = C.outputs(snakemake)

#: Which plot of the family each producer wants.
WANTED = {"panel_metaloci_gaudi": "gtp", "panel_metaloci_kk": "kk"}


with C.logging(snakemake) as say:
    run = spec.get("run", "")
    dataset = spec.get("dataset", "")
    signal = spec.get("signal", "")
    suffix = WANTED.get(spec.get("producer", ""), "gtp")
    mldir = str(snakemake.params.mldir)
    start, end = (int(v) for v in snakemake.params.window)
    figure_cfg = dict(snakemake.params.figure_cfg or {})

    say("panel {}: run={} dataset={} signal={} want _{}".format(
        spec["panel_id"], run, dataset, signal, suffix))

    signals = os.path.join(mldir, "{}.signals".format(dataset))
    region = "chrX:{}-{}_0".format(start, end)
    cmd = ["metaloci", "figure", "-w", mldir, "-s", signals, "-g", region,
           "-m", "-t", str(snakemake.threads)]
    if figure_cfg.get("export_intermediates", True):
        cmd.append("-e")
    if figure_cfg.get("zscore", False):
        cmd.append("-z")
    say("  $ " + " ".join(cmd))

    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
        say(proc.stdout[-8000:])
        if proc.stderr:
            say("stderr:\n" + proc.stderr[-8000:])
        rc = proc.returncode
    except FileNotFoundError:
        rc = 127
        say("  metaloci is not on PATH in this environment")
    except subprocess.TimeoutExpired:
        rc = 124
        say("  metaloci figure timed out")

    # Collect the requested plot out of METALoci's own tree.
    hits = []
    for pattern in ("**/*_{}.pdf".format(suffix), "**/{}/*.pdf".format(suffix.upper())):
        hits.extend(glob.glob(os.path.join(mldir, pattern), recursive=True))
    if signal:
        preferred = [h for h in hits if signal in os.path.basename(h)]
        hits = preferred or hits
    hits = sorted(set(hits))
    say("  {} candidate _{} file(s)".format(len(hits), suffix))

    C.ensure_parent(out_path)
    if hits and out_path.lower().endswith(".pdf"):
        shutil.copyfile(hits[0], out_path)
        say("  copied {} -> {}".format(hits[0], out_path))
        C.write_n_items(snakemake, len(hits))
    elif hits:
        try:
            import fitz
            doc = fitz.open(hits[0])
            with open(out_path, "w") as fh:
                fh.write(doc[0].get_svg_image())
            doc.close()
            C.write_n_items(snakemake, len(hits))
        except Exception as exc:
            pl.empty_panel(out_path, "cannot convert {} to {}: {}".format(
                os.path.basename(hits[0]),
                os.path.splitext(out_path)[1], exc))
            C.write_n_items(snakemake, 0)
    else:
        pl.empty_panel(out_path,
                       "metaloci figure produced no _{} plot\n"
                       "(exit {}) for {} / {} / {}".format(
                           suffix, rc, run, dataset, signal))
        C.write_n_items(snakemake, 0)
        say("  no _{} plot found; wrote a placeholder".format(suffix))
