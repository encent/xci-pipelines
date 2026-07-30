"""Select the curated eigenvector for one cooler and stage its outputs.

Reimplements ``01_06_compartments_refine.py``. The original read
`compartments{,_merged}_{locus_type}.json`, copied the chosen eigenvector's
bigWig and saddle SVGs into a `compartments_refined*` tree, and appended one
row to `saddle_strength_selected_{locus_type}.tsv`.

Two changes, both deliberate:

1. **The selected bigWig gets a stable name**, `Comp_selected_{name}.bw`,
   rather than the ground truth's `Comp_E2_{name}.bw`. A Snakemake output path
   cannot embed a value that is only known once the job runs -- and under
   `eigenvector_selection: auto` the choice IS data-dependent. `selection.tsv`
   records which eigenvector it is, so nothing is lost.
2. **No SVGs are copied.** They are drawn from the `.npz` in 90_visualise.smk.

PLAN 9.5. The JSON records WHICH eigenvector, never its ORIENTATION. The sign
canonicalisation applied by `eigs_cis` is carried into `selection.tsv` here, so
the choice and the flip appear in the same row. Under `auto` the choice is
max |Spearman(E, GC)| and all three correlations are recorded.
"""
import json
import os
import shutil
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from lib import errors  # noqa: E402

CURATE = """
        1. run   ./run.sh compartment_panels
        2. look at
           results/figures/qc/eigenvector_{name}_{roi}.svg
           and decide which of E1/E2/E3 tracks the A/B pattern
        3. add it to
           resources/fixtures/compartments/{fixture}
           as    "{name}": "E2"
"""


def main():
    smk = snakemake  # noqa: F821
    name = smk.wildcards.name
    roi = smk.wildcards.roi
    eig_labels = list(smk.params.eigs)

    os.makedirs(os.path.dirname(smk.output.tsv), exist_ok=True)
    os.makedirs(os.path.dirname(smk.log[0]), exist_ok=True)

    orientation = pd.read_csv(smk.input.orientation, sep="\t")
    values = pd.read_csv(smk.input.values, sep="\t").iloc[0]
    rho = dict(zip(orientation["eigenvector"], orientation["spearman_gc_final"]))
    flipped = dict(zip(orientation["eigenvector"], orientation["sign_flipped"]))

    with open(smk.log[0], "w") as fh:
        if smk.params.selection == "fixture":
            fixture = json.load(open(smk.input.fixture))
            if name not in fixture:
                raise errors.missing_fixture(
                    name=f"eigenvector choice for {name!r} ({roi})",
                    path=smk.input.fixture,
                    what="Which of the three eigenvectors tracks the A/B "
                    "compartment pattern was decided by eye, per cooler and "
                    "per region, by looking at the contact maps.",
                    curate=CURATE.format(
                        name=name, roi=roi,
                        fixture=os.path.basename(smk.input.fixture)),
                )
            chosen = fixture[name]
            source = "fixture"
        else:
            chosen = max(eig_labels, key=lambda e: abs(rho.get(e, 0.0)))
            source = "auto"
        fh.write(f"{name} / {roi}: {source} -> {chosen}\n")
        for e in eig_labels:
            fh.write(f"  {e}: Spearman(E, GC) = {rho.get(e, float('nan')):+.4f}"
                     f"  sign_flipped={bool(flipped.get(e, False))}"
                     f"  saddle_score = {values[e]}\n")

        if chosen not in eig_labels:
            raise ValueError(
                f"{smk.input.fixture} chooses {chosen!r} for {name}, but only "
                f"{eig_labels} were computed (hic.compartments.n_eigs)."
            )

        src_bw = dict(smk.params.eig_bw)[chosen]
        shutil.copyfile(src_bw, smk.output.bw)
        shutil.copyfile(smk.input.eigs, smk.output.tsv)
        fh.write(f"  {src_bw} -> {smk.output.bw}\n")

        pd.DataFrame([{
            "file": name,
            "roi": roi,
            "scope": smk.wildcards.cscope,
            "eigenvector": chosen,
            "source": source,
            "value": values[chosen],
            "sign_flipped": bool(flipped.get(chosen, False)),
            "spearman_gc": rho.get(chosen, float("nan")),
            **{f"spearman_gc_{e}": rho.get(e, float("nan")) for e in eig_labels},
            **{f"saddle_score_{e}": values[e] for e in eig_labels},
        }]).to_csv(smk.output.selection, sep="\t", index=False)

        if bool(flipped.get(chosen, False)):
            fh.write("WARNING: the selected eigenvector was sign-flipped to "
                     "GC-rich-positive. The ground truth's orientation is "
                     "unrecorded, so its saddle plot may be mirrored relative "
                     "to ours (plan 9.5, risk R-3). The visualisation stage "
                     "must surface this.\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
