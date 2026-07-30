"""Fig 4e (and the EFig 6d / Fig 1c-d style views) — coolbox browser zoom-ins.

Runs in its own conda env (`workflow/envs/coolbox.yaml`). coolbox and the main
environment cannot coexist (plan section 5), so this and `metaloci figure` are
the only two panels carrying a `conda:` directive. Both are invoked FROM the
visualisation stage, which is what keeps "figures come only from
90_visualise.smk" true rather than aspirational.

THE TRACK STACK, transcribed from `01_08_coolbox_diff_hic_only.py`
------------------------------------------------------------------
Per condition (NodTAG then dTAG, or Xa then Xi):

    Cool(mcool)                        cmap YlOrRd, matrix, MinValue -9 / MaxValue -1
      + HiCPeaksCoverage(loops.bedpe)  color #2255ff, line_width 5, side upper
    Spacer(0.5) + BigWig(Comp_selected)         the refined compartment track
    Spacer(0.5) + BigWig(mark)                  one per staged mark
      and for H3K27me3 only:
    Spacer(0.5) + BED(valleys)         labels False, collapsed, TrackHeight(0.3)
    GTFAllelic(allelic_gtf)            TrackHeight 12 for Mecp2, 7 for Jarid
    GFFNC(GRCm38.102_NC.gtf)           TrackHeight 12 -- ALWAYS 12

then once:

    HiCDiff(cool1, cool2)              normalize expect, diff_method diff,
                                       cmap bwr, MinValue -6 / MaxValue 6
    Selfish(cool1, cool2)              norm log, matrix
    XAxis()
    frame *= Feature(depth_ratio=0.35)

The two TrackHeight(12)s are NOT the same constant. The allelic GTF track is
12 for Mecp2 and 7 for Jarid, and lives in the F-7 fixture as
`gtf_track_height`; the non-coding GFFNC track is a flat 12 at all four call
sites (01_08 lines 354, 367, 459, 472). Collapsing them into one number would
silently change one panel's proportions.

Everything numeric here comes from `resources/fixtures/coolbox_display.yaml`
(fixture F-7). Wrong values render a published panel at a different contrast
and nothing complains, which is why they are versioned rather than inline.
"""

import os
import sys

sys.path.insert(0, os.path.join(snakemake.params.spec["libdir"],
                                "scripts", "py", "panels"))
sys.path.insert(0, snakemake.params.spec["libdir"])

import yaml                                                    # noqa: E402
import _common as C                                            # noqa: E402
from lib import plotting as pl                                 # noqa: E402

spec = C.bootstrap(snakemake)
out_path = C.outputs(snakemake)

#: 01_08 lines 354, 367, 459, 472 -- a flat 12 for the non-coding track.
GFFNC_TRACK_HEIGHT = 12


with C.logging(snakemake) as say:
    locus = spec.get("locus", "Mecp2")
    roi = spec.get("roi", "full")
    comparison = spec.get("comparison", "dTAG_vs_NodTAG")
    with open(str(snakemake.params.display)) as fh:
        display = yaml.safe_load(fh)

    say("panel {}: {} / {} / {}".format(spec["panel_id"], locus, roi,
                                        comparison))
    say("  F-7 display fixture: {}".format(display))

    paths = [str(p) for p in snakemake.input]
    mcools = sorted(p for p in paths if p.endswith(".mcool"))
    bigwigs = sorted(p for p in paths if p.endswith(".bw"))
    beds = sorted(p for p in paths if p.endswith(".bed"))
    bedpes = sorted(p for p in paths if p.endswith(".bedpe"))
    gtfs = sorted(p for p in paths if p.endswith(".gtf"))
    nc_gtf = next((p for p in gtfs if os.path.basename(p).endswith("_NC.gtf")),
                  None)
    allelic = [p for p in gtfs if p != nc_gtf]
    say("  {} mcools, {} bigWigs, {} valley BEDs, {} bedpe, {} allelic GTFs, "
        "GFFNC {}".format(len(mcools), len(bigwigs), len(beds), len(bedpes),
                          len(allelic), "yes" if nc_gtf else "MISSING"))
    if not nc_gtf:
        say("  WARNING: no GRCm38.102_NC.gtf among the inputs. The published "
            "panel carries a GFFNC track at TrackHeight(12); this render "
            "would be missing it.")

    try:
        from coolbox.api import (BED, BigWig, Cool, Feature, Frame, HiCDiff,
                                 HiCPeaksCoverage, MaxValue, MinValue, Spacer,
                                 Title, TrackHeight, XAxis)
    except Exception as exc:
        pl.empty_panel(
            out_path,
            "coolbox is not importable here ({}: {}).\n"
            "This panel needs workflow/envs/coolbox.yaml; run with "
            "--use-conda.".format(type(exc).__name__, exc))
        C.write_n_items(snakemake, 0)
        say("  coolbox unavailable: {}".format(exc))
        raise SystemExit(0)

    try:
        from coolbox.api import GFFNC, GTFAllelic, Selfish
    except Exception:
        GFFNC = GTFAllelic = Selfish = None
        say("  note: this coolbox build has no GFFNC / GTFAllelic / Selfish; "
            "those tracks are skipped and the panel is NOT publication-faithful")

    hic = display.get("hic", {})
    diff = display.get("diff", {})
    loops_cfg = display.get("loops", {})
    gtf_height = display.get("gtf_track_height", {}).get(locus, 12)
    valley_height = display.get("valley_track_height", 0.3)
    depth_ratio = display.get("depth_ratio", 0.35)
    dpi = int(display.get("dpi", spec.get("dpi", 300)))

    frame = Frame()
    cools = []
    n_tracks = 0
    for mcool in mcools:
        cool = (Cool(mcool, cmap=hic.get("cmap", "YlOrRd"),
                     style=hic.get("style", "matrix"))
                + MinValue(hic.get("min_value", -9))
                + MaxValue(hic.get("max_value", -1)))
        cools.append(cool)
        frame += cool + Title(os.path.basename(mcool))
        n_tracks += 1
        for bedpe in bedpes:
            frame += HiCPeaksCoverage(bedpe,
                                      color=loops_cfg.get("color", "#2255ff"),
                                      line_width=loops_cfg.get("line_width", 5),
                                      side=loops_cfg.get("side", "upper"))

        for bw in bigwigs:
            frame += Spacer(0.5) + BigWig(bw) + Title(os.path.basename(bw))
            n_tracks += 1
            if os.path.basename(bw).startswith("H3K27me3"):
                stem = os.path.splitext(os.path.basename(bw))[0]
                bed = next((b for b in beds
                            if os.path.basename(b).startswith(stem)), None)
                if bed:
                    frame += (Spacer(0.5)
                              + BED(bed, labels=False, display="collapsed")
                              + TrackHeight(valley_height)
                              + Title(os.path.basename(bed)))
                    n_tracks += 1

        if GTFAllelic is not None and allelic:
            frame += (GTFAllelic(allelic[0])
                      + Title(os.path.basename(allelic[0]))
                      + TrackHeight(gtf_height))
            n_tracks += 1
        if GFFNC is not None and nc_gtf:
            frame += (GFFNC(nc_gtf) + Title(os.path.basename(nc_gtf))
                      + TrackHeight(GFFNC_TRACK_HEIGHT))
            n_tracks += 1

    if len(cools) >= 2:
        frame += (HiCDiff(cools[0], cools[1],
                          normalize=diff.get("normalize", "expect"),
                          diff_method=diff.get("diff_method", "diff"),
                          cmap=diff.get("cmap", "bwr"),
                          style=diff.get("style", "matrix"))
                  + MinValue(diff.get("min_value", -6))
                  + MaxValue(diff.get("max_value", 6)))
        n_tracks += 1
        if Selfish is not None:
            frame += Selfish(cools[0], cools[1], norm="log", style="matrix")
            n_tracks += 1

    frame += XAxis()
    frame *= Feature(depth_ratio=depth_ratio)

    start, end = (int(v) for v in
                  snakemake.config["loci"][locus]["regions"][roi])
    region = "chrX:{}-{}".format(start, end)
    say("  plotting {} with {} tracks".format(region, n_tracks))

    C.ensure_parent(out_path)
    fig = frame.plot(region)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    C.write_n_items(snakemake, n_tracks,
                    extra={"region": region, "gffnc": bool(nc_gtf),
                           "gtf_track_height": gtf_height})
    say("wrote {}".format(out_path))
