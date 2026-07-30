# Troubleshooting

Every problem below was hit for real while building and testing this pipeline.
They are ordered by how likely you are to meet them.

---

## 1. `install.sh` seems to hang at "Solving environment"

**It is not hung — it is conda's classic solver, and it may never finish.**

This environment pins about forty exact package versions across conda-forge and
bioconda. On the machine this pipeline was developed on, conda's classic solver
ran for **over forty minutes** without finishing, printing only
`Solving environment: ...working...`. mamba solved the identical file in about
two minutes.

`install.sh` handles this for you: it prefers `mamba`, falls back to
`micromamba`, and otherwise downloads a standalone micromamba into `.mamba/`
inside the repo. It only falls back to plain conda as a last resort, and warns
you when it does.

**If you build the environment by hand, use mamba:**

```bash
mamba env create -f environment.yml
```

Installing [Miniforge](https://github.com/conda-forge/miniforge) rather than
Miniconda avoids this entirely — it ships mamba by default.

---

## 2. `KeyError: 'loci'` or `KeyError: 'sheets'` on startup

You almost certainly passed `--configfile` twice:

```bash
# WRONG — silently discards the first file
snakemake --configfile config/config.yaml --configfile my_overrides.yaml

# RIGHT — one flag, space-separated; later files override earlier ones
snakemake --configfile config/config.yaml my_overrides.yaml
```

Snakemake's `--configfile` takes a *list*, so a second flag **replaces** the
first rather than merging with it. The pipeline then starts with a config that
is missing whole sections, and the failure surfaces from the `Snakefile` as a
missing key — which looks like a corrupt config file rather than a command-line
mistake.

---

## 3. `LockException`, or "Directory cannot be locked"

Snakemake takes an **exclusive lock on the working directory**, so only one run
at a time can use it. If you want two runs in parallel, give each its own
directory:

```bash
snakemake --directory /path/to/run_a ...
snakemake --directory /path/to/run_b ...
```

If a previous run was killed, the lock can be left behind. Clear it with:

```bash
snakemake --unlock
```

**Note:** `--directory` changes the working directory, so any *relative* paths in
your config (`config/samples.tsv`, `config/coolers.tsv`, `config/panels.yaml`)
will stop resolving. Make them absolute in your overrides file when you use it.

---

## 4. Valley calls differ between machines, or between runs

They should not, and if they do it is worth reporting. But be aware of the
mechanism, because it is subtle.

`call_valleys` fits a 2-state Gaussian HMM whose emission means are seeded by
k-means. That seeding sums partial results in an order that depends on **how many
threads the numerical libraries use** — so the same input, same code and same
seed can land on a different local optimum and shift whole valleys.

The pipeline pins this for you: `valleys.blas_threads` (default **2**) is applied
around the fit with `threadpoolctl`, and the *observed* thread count is written
to the log and to `results/provenance.json`.

Two things to know if you change it:

* **Setting `threads:` on a Snakemake rule does nothing here.** That only tells
  the scheduler how to allocate slots; it never reaches the numerical libraries.
* **The limit must apply to every runtime, not just BLAS.** Capping BLAS while
  leaving OpenMP at its default still changes the result. This is why the code
  calls `threadpool_limits(limits=N)` with no `user_api` argument — do not add
  one.

If you see unexpected valley differences, compare the observed thread counts in
the two `results/provenance.json` files first.

---

## 5. Coverage tracks differ slightly, and I sliced my BAMs to save time

**Do not slice BAMs to a single chromosome before running.** It changes the
answer and it is not faster.

deepTools estimates the fragment-length distribution from the reads it can see.
A `samtools view -b sample.bam chrX` slice gives a different estimate — measured
here as `maxPairedFragmentLength` 688 instead of 696 — and **27,170 differing
bases** on chrX alone.

Run `bam_coverage` on the full BAM and restrict downstream instead, with
`chromosomes.downstream` in the config.

---

## 6. "no scaleFactor row for &lt;sample&gt;" / scale-factor errors

Every allelic split needs an unsplit (`Gall`) sibling with the same
`(mark, clone, condition, replicate)`, because one library is sequenced once and
then split by SNP — so all three files share a single scale factor.

Check `config/samples.tsv` for a missing or mistyped `Gall` row.

If you are supplying pre-computed factors, note the file has two possible
formats: a `scale_factor` column (the product `1e6 / (LibSize × NormFactor)`, the
only quantity recoverable from a finished bigWig) or a `scaleFactor` column with
separate `LibSize`/`NormFactor`. Both are accepted. Comment lines beginning `#`
are ignored.

---

## 7. Normalisation on chrX alone gives different numbers

This is expected, not a bug. The TMM normalisation factors are computed on
**autosomes only** — chrX is deliberately excluded, because one X is inactivated
and *which* parental X differs between clones, so including it would make the
factors reflect X-inactivation rather than library size.

A chrX-only input therefore cannot reproduce the published factors at all. Use
`chromosomes.csaw` to control which chromosomes the factor estimation sees,
independently of `chromosomes.downstream`.

---

## 8. `samtools index` fails with "Permission denied"

Your BAM directory is read-only. That is fine and the pipeline supports it — it
copies any existing `.bai` into its own workspace rather than writing next to
your data, and never treats your input directory as writable.

If you hit this anyway, check that `paths.data_dir` points somewhere writable
with enough space.

---

## 9. A rule fails and I want to know where to look

* One log per step under `<data_dir>/logs/<rule>/`.
* Nothing is lost. `./run.sh` resumes from where it stopped.
* For valley mismatches, `<track>_states.tsv.gz` holds the per-bin signal and
  HMM state, so a shifted boundary can be localised to a specific bin rather than
  diffed as two BED files.

---

## 10. Still stuck

Open an issue with:

* the exact command you ran,
* the failing rule name from the Snakemake output,
* the contents of that rule's log file,
* `results/provenance.json` if it exists.

The provenance file records package versions, thread counts and config values,
and it is usually enough to reproduce the problem.
