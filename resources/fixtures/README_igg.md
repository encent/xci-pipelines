# F-9 — the IgG decision

`03_call_peaks.sh` hard-codes `USE_CONTROL=false` and references
`/mnt/scratch/nikolai/Antonia/bams_IgG/pkg/WT`, which **does not exist**. The
metadata sheet has no IgG rows and no IgG BAMs remain on disk.

**Decision: the IgG branch is intentionally dead.** MACS2 peak calling in this
pipeline runs without a control and is QC-only (`qc.enabled`, default `false`).

Recorded here so this is a decision on the record rather than something
rediscovered during testing.

## Related: the CTCF peaks that actually matter

Every boundary/valley analysis reads the **external** delivery
`resources/fixtures/CTCFpeak_per_clone/*_consensusPeaks_withMeanRatioandDirection.bed`,
**not** MACS2 output. Never wire `03_call_peaks.sh`-equivalent output into the
boundary analysis. (Plan correction N-5.)

## Related: WT H3K27ac narrow vs broad (deviation D-12)

The saved `03_call_peaks.sh` emits narrowPeak for H3K27ac, but the surviving
on-disk peaks are 30 `.broadPeak` + 30 `.gappedPeak` — a commented-out `elif`
branch was live at run time. Two generations exist. QC-only; do not chase.
