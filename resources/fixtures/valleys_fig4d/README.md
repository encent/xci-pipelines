# Fig 4d valley fixture — four CTCF-degron Xi tracks

The published H3K27me3 valley calls for the four tracks behind **Fig 4d**,
copied out of the original pipeline's output (read-only; nothing was written
back):

```
/mnt/scratch/nikolai/Antonia/all_bams/norm_bw_H3K27me3_compos_all/merged_5000/valleys_annotation/
```

**This fixture exists so that Fig 4d can be drawn exactly as published. It is
not a claim that the pipeline agrees with it.** The pipeline computes its own
valleys for these tracks on every run, and the two differ. The deltas are below
precisely so nobody assumes otherwise.

## The four tracks, and how far our computation is from each

At the frozen `valleys.blas_threads: 2`:

| track | pipeline | published | delta |
|---|---|---|---|
| `H3K27me3_E6A7_CTCF-NodTAG_Xi` | 293 | **321** | **−28** |
| `H3K27me3_E6A7_CTCF-dTAG_Xi` | 291 | **281** | **+10** |
| `H3K27me3_F3_CTCF-NodTAG_Xi` | 343 | **351** | **−8** |
| `H3K27me3_F3_CTCF-dTAG_Xi` | 204 | **205** | **−1** |

The severity spans more than an order of magnitude, which is why the four are
named individually rather than treated as a block. `F3_CTCF-dTAG` misses by a
**single valley**; `E6A7_CTCF-NodTAG` misses by 28.

## What this fixture is NOT for

**Only these four tracks.** The other seven Xi tracks — `B1`, `C5`, `CL30`,
`E6`, `JTG`, and **both** `B1621` Rad21 conditions — reproduce count-exactly
from computation and must **never** be served from a fixture. A fixture there
would mask a real regression.

**The four-track list is PROVISIONAL and must not be frozen yet.** Comparing
the pipeline's valleys against the sweep on the *original* bigWigs, 30 of 32
chrX tracks agree but **two flip**. `H3K27me3_B1621_Rad21-dTAG_Xi` gives **307
through the pipeline** against 281 from the original bigWig — so "both B1621
Rad21 conditions reproduce count-exactly" holds for the sweep but **not through
the pipeline**.

That points at the signal chain, not at the valley caller: `bam_coverage`
differs on 10.2 % of chrX bases and `merge_reps_5k` is 0 of 33 byte-identical,
with a systematic one-sided drift whose mechanism is unknown. Until that is
understood, the membership of this fixture cannot be settled — the mechanism is
built, the list is not final.

## Why Fig 2g does not need one

**Fig 2g reproduces from computation.** Running the pipeline's own three scripts
unmodified at `blas_threads=2` — `call_valleys` → `filter_valleys_genes` →
`valley_gene_content` — gives **375 / 301 / 258 / 300 / 321**, the published
sequence exactly.

And the obvious worry there is **measured absent, not assumed absent**: `B1` and
`CL30` are count-exact but *not* byte-identical, and both still produce the
correct gene-content total. So a ≤4-bin edge shift does not move a valley
between gene categories.

So the honest position is: **Fig 2g is reproduced from computation; Fig 4d is
the only published valley panel requiring this fixture, and one of its four
tracks misses by a single valley.**

## Switching it

`config.valleys.fig4d_source`:

* `fixture` *(default)* — Fig 4d is drawn from these files. The run log and
  `provenance.json` say so explicitly.
* `computed` — Fig 4d is drawn from the pipeline's own calls. The panel will
  then show 293 / 291 / 343 / 204 instead of the published numbers.

The pipeline computes its own valleys either way, so the comparison is always
available to the testing team.

## Provenance

BED3, `chrX` only, tab-separated, no header — identical in format to
`call_valleys` output. Checksummed in `../fixtures.sha256` and verified by
`check_fixtures` on every run.
