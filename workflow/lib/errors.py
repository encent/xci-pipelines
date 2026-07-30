"""User-facing errors.

The target reader of this pipeline has our Nature Genetics paper open and has
never used a terminal. Every failure they can plausibly hit gets a written
message, never a stack trace, and every message follows the same shape:

    what went wrong  ->  why it matters  ->  exactly what to type next
"""

from __future__ import annotations

import sys
import textwrap

_WIDTH = 88


class XCIError(Exception):
    """Base for every error we present to a user."""

    def __init__(self, title: str, why: str, fix: str, detail: str = ""):
        self.title, self.why, self.fix, self.detail = title, why, fix, detail
        super().__init__(self.render())

    def render(self) -> str:
        out = ["", "=" * _WIDTH, f"ERROR  {self.title}", "=" * _WIDTH, ""]
        if self.detail:
            out += [_wrap(self.detail), ""]
        out += ["  Why this matters:", _wrap(self.why, indent="    "), ""]
        out += ["  What to do:", _indent_block(self.fix), "", "=" * _WIDTH, ""]
        return "\n".join(out)


def _wrap(text: str, indent: str = "  ") -> str:
    return "\n".join(
        textwrap.fill(p, _WIDTH - len(indent), initial_indent=indent, subsequent_indent=indent)
        for p in text.strip().split("\n\n")
    )


def _indent_block(text: str) -> str:
    return "\n".join("    " + ln for ln in textwrap.dedent(text).strip().splitlines())


# --------------------------------------------------------------------------
# Concrete errors
# --------------------------------------------------------------------------
def missing_input_dir(key: str, path: str) -> XCIError:
    return XCIError(
        title=f"Input directory not found: {key}",
        detail=f"config/config.yaml -> paths.inputs.bam_dirs.{key}\n  points at: {path}\n"
        "That directory does not exist, or you do not have permission to read it.",
        why="The pipeline reads your aligned BAM files from this directory. It never "
        "writes there -- it only creates symbolic links into its own working area -- "
        "but it cannot start without knowing where the data is.",
        fix=f"""
        Open config/config.yaml and set the correct path:

            paths:
              inputs:
                bam_dirs:
                  {key}: /the/real/path/to/your/bams

        If you do not have this data type at all, delete the '{key}:' line
        entirely -- the pipeline will skip everything that depends on it.
        """,
    )


def missing_bam_index(bam: str) -> XCIError:
    return XCIError(
        title="BAM file has no index",
        detail=f"{bam}\n  has no companion .bai file.",
        why="Reading coverage from a BAM requires a random-access index. Without it "
        "the normalisation step cannot run.",
        fix=f"""
        conda activate xci-pipeline
        samtools index "{bam}"

        To index everything at once:

            for f in /path/to/bams/*.bam; do samtools index "$f"; done
        """,
    )


def missing_fixture(name: str, path: str, what: str, curate: str) -> XCIError:
    return XCIError(
        title=f"Missing fixture: {name}",
        detail=f"Expected at: {path}",
        why=f"{what}\n\nThat judgement is not derivable from code, so the pipeline "
        "ships it as a versioned fixture rather than pretending to recompute it. "
        "You are normally seeing this because you are running on your own data "
        "rather than reproducing the paper.",
        fix=f"""
        Choose one:

        (a) Use the automatic fallback -- edit config/config.yaml:

                fixtures:
                  mode: auto

            Affected panels are then tagged "auto" in results/figure_manifest.tsv
            and footnoted in results/REPORT.html.

        (b) Curate it by hand:
        {curate}
        """,
    )


def chrx_in_csaw_set(chroms: list) -> XCIError:
    return XCIError(
        title="chrX must not be in the csaw chromosome set",
        detail=f"config/config.yaml -> chromosomes.csaw = {chroms}",
        why="The library-size normalisation (csaw/edgeR TMM) is computed on autosomes "
        "ONLY, and this is deliberate. One X chromosome is inactivated, and which "
        "parental X is inactivated differs between clones -- so including chrX makes "
        "the TMM factors reflect X-inactivation rather than library size, and every "
        "downstream track becomes wrong in a way that looks plausible.",
        fix="""
        Remove chrX from chromosomes.csaw. The default set is the one the paper used:

            chromosomes:
              csaw: [chr2, chr3, chr4, chr5, chr7, chr8, chr9, chr10,
                     chr13, chr14, chr15, chr17, chr18, chr19]

        chrX belongs in chromosomes.downstream, where it already is.
        """,
    )


def sample_sheet_problem(detail: str, fix: str) -> XCIError:
    return XCIError(
        title="Problem in config/samples.tsv",
        detail=detail,
        why="The sample sheet defines every job the pipeline will run. An inconsistency "
        "here would silently mis-pair samples -- exactly the class of bug that put a "
        "wrong scale factor on five files in the original pipeline.",
        fix=fix,
    )


def inconsistent_allele_map(clone: str, rows: str) -> XCIError:
    return XCIError(
        title=f"Inconsistent G1/G2 -> Xa/Xi mapping for clone {clone}",
        detail=f"config/samples.tsv rows:\n{rows}",
        why="Which parental X is inactivated is a property of the clone, fixed for all "
        "its samples. If one row says G1=Xi and another says G1=Xa for the same clone, "
        "one of them is wrong and allele-specific tracks would be swapped.",
        fix=f"""
        Open config/samples.tsv, find the rows for clone {clone}, and make the
        snpsplit_tag / allele pairing consistent. For reference, in the paper:

            E6, JTG, CL30  ->  Xi is C57BL/6J (G1)
            B1, C5         ->  Xi is Cast/EiJ (G2)
        """,
    )


def not_enough_disk(need_gb: int, have_gb: int, path: str) -> XCIError:
    return XCIError(
        title="Not enough free disk space",
        detail=f"{path}\n  free: {have_gb} GB     needed: about {need_gb} GB",
        why="A full run writes roughly 500-600 GB of intermediate tracks and contact "
        "maps. Running out halfway through wastes the whole run.",
        fix=f"""
        Either free up space, or point the pipeline somewhere larger:

            paths:
              data_dir: /a/disk/with/room

        To try the pipeline without the full data volume first:

            ./run.sh demo      # about 2 hours, a few GB
        """,
    )


def conda_missing() -> XCIError:
    return XCIError(
        title="conda was not found",
        detail="The 'conda' command is not on your PATH.",
        why="The pipeline installs its software with conda so that you get the exact "
        "package versions the paper used.",
        fix="""
        Install Miniforge (free, no licence restrictions):

            https://github.com/conda-forge/miniforge

        Then close and reopen your terminal and run  bash install.sh  again.
        """,
    )


def die(err: XCIError) -> None:
    """Print a user-facing error and exit non-zero."""
    sys.stderr.write(err.render())
    sys.exit(1)
