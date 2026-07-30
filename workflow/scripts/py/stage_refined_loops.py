"""Stage the manually curated loop set (D-05).

The published loop calls are the raw chromosight output minus false positives
deleted by hand after visual inspection. Nothing in the original code produces
them, so `paper` mode ships them as fixtures and this script only chooses
between the fixture and the raw calls -- loudly.
"""
import hashlib
import os
import shutil
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from lib import errors  # noqa: E402

CURATE = """
        1. run   ./run.sh loop_overlays
        2. open  results/figures/qc/loop_overlay_{name}_before.pdf
           and delete the false-positive rows from
           work/features/loops/raw/{name}/{name}.tsv
        3. save the result as
           resources/fixtures/loops_refined/{name}.tsv
"""

NOTE_PAPER = """\
# Loop refinement: {name}

Source: **fixture** `resources/fixtures/loops_refined/{name}.tsv`
SHA256: `{sha}`

| | loops |
|---|---|
| raw (chromosight) | {n_raw} |
| refined (published) | {n_ref} |

The refined set is the raw set minus false positives deleted by hand after
visual inspection of the contact map. That judgement is not derivable from
code, so it ships as a versioned fixture rather than being recomputed.
Published counts: Jarid_Xa 14, Jarid_Xi 13, Mecp2_Xa 10, Mecp2_Xi 12.
"""

NOTE_AUTO = """\
# MANUAL CURATION NEEDED: {name}

Source: **raw chromosight calls, UNREFINED** ({n_raw} loops)

`fixtures.mode: auto` is set, so no hand-curated loop set was used. Every panel
drawn from this loop set is tagged `auto` in `results/figure_manifest.tsv` and
footnoted in `results/REPORT.html`.

The published figures used a curated set: false positives were deleted by hand
after looking at the contact map. To do that here:

1. `./run.sh loop_overlays`
2. open `results/figures/qc/loop_overlay_{name}_before.pdf`
3. delete the false-positive rows from
   `work/features/loops/raw/{name}/{name}.tsv`
4. save the result as `resources/fixtures/loops_refined/{name}.tsv`
5. set `fixtures.mode: paper` in `config/config.yaml` and re-run

For reference, the published pipeline kept these many loops:
Jarid_Xa 14 of 55, Jarid_Xi 13 of 32, Mecp2_Xa 10 of 41, Mecp2_Xi 12 of 26.
"""


def n_rows(path):
    with open(path) as fh:
        return max(0, sum(1 for _ in fh) - 1)


def main():
    smk = snakemake  # noqa: F821
    name = smk.wildcards.name
    fixture = smk.params.fixture
    os.makedirs(os.path.dirname(smk.output.tsv), exist_ok=True)
    os.makedirs(os.path.dirname(smk.log[0]), exist_ok=True)

    n_raw = n_rows(smk.input.raw)

    with open(smk.log[0], "w") as fh:
        if smk.params.mode == "paper":
            if not os.path.exists(fixture):
                raise errors.missing_fixture(
                    name=f"loop refinement for cooler {name!r}",
                    path=fixture,
                    what="The published loop sets were curated by hand -- false "
                    "positives were deleted after visual inspection of the "
                    "contact map.",
                    curate=CURATE.format(name=name),
                )
            shutil.copyfile(fixture, smk.output.tsv)
            sha = hashlib.sha256(open(fixture, "rb").read()).hexdigest()
            n_ref = n_rows(smk.output.tsv)
            fh.write(f"{name}: fixture {fixture}\n  sha256 {sha}\n"
                     f"  raw {n_raw} -> refined {n_ref}\n")
            note = NOTE_PAPER.format(name=name, sha=sha, n_raw=n_raw, n_ref=n_ref)
        else:
            shutil.copyfile(smk.input.raw, smk.output.tsv)
            fh.write(f"{name}: fixtures.mode={smk.params.mode}; using the RAW "
                     f"chromosight calls ({n_raw} loops), unrefined.\n"
                     f"  wrote {smk.output.note}\n")
            note = NOTE_AUTO.format(name=name, n_raw=n_raw)

    with open(smk.output.note, "w") as fh:
        fh.write(note)
    return 0


if __name__ == "__main__":
    sys.exit(main())
