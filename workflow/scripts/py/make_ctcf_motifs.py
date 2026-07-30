"""CTCF motif occurrences (JASPAR MA0139.1), genome-wide and chrX-only.

Two modes, chosen by ``config.fixtures.ctcf_motifs``:

``fixture`` (default)
    Use the committed ``CTCF_mm10_X_only.bed.gz``. The original produced this
    set inside an R/AnnotationHub session that cannot be reproduced offline, and
    the motif set is an *input* to the valley-boundary analysis (it decides
    which boundaries count as CTCF-proximal), so it is versioned rather than
    recomputed. See the archaeology notes fixture inventory §J.7.

``recompute``
    Scan the mm10 FASTA with the JASPAR PWM. This will NOT reproduce the
    fixture exactly — threshold and background conventions differ between
    Biostrings' ``matchPWM`` and the scanner here — so it is offered for users
    running on their own genome, not for reproducing the paper.

Note on the genome-wide output in ``fixture`` mode: the fixture covers chrX
only, because chrX is all the paper analysis uses. The genome-wide file is
written as a copy of it and the restriction is logged loudly. Only the opt-in
QC enrichment rule reads the genome-wide file, and it is off by default.
"""

import gzip
import os
import shutil
import sys

log_path = snakemake.log[0]
os.makedirs(os.path.dirname(log_path), exist_ok=True)
os.makedirs(os.path.dirname(snakemake.output.chrx), exist_ok=True)

mode = snakemake.params.mode
fixture = snakemake.params.fixture
chrx_name = snakemake.params.chrx


def _log(handle, msg: str) -> None:
    handle.write(msg + "\n")
    print(msg, file=sys.stderr)


def _read_fixture(path: str) -> list[str]:
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        return [ln for ln in fh if ln.strip() and not ln.startswith(("#", "track"))]


def _scan_fasta(fa_path: str, jaspar_id: str, log) -> list[str]:
    """Fallback PWM scan. Deliberately simple and clearly not the fixture."""
    try:
        import numpy as np
        import pyjaspar  # noqa: F401
    except ImportError as exc:  # pragma: no cover - depends on user's env
        raise SystemExit(
            f"recompute mode needs pyjaspar and numpy ({exc}).\n"
            "Either install them, or set  fixtures.ctcf_motifs: fixture  in "
            "config/config.yaml to use the committed motif set."
        )

    from pyjaspar import jaspardb

    motif = jaspardb(release="JASPAR2020").fetch_motif_by_id(jaspar_id)
    pssm = motif.pssm
    threshold = pssm.max * 0.8
    _log(log, f"scanning with {jaspar_id}, length {len(motif)}, threshold {threshold:.3f}")

    from Bio import SeqIO

    rows = []
    for record in SeqIO.parse(fa_path, "fasta"):
        seq = record.seq.upper()
        for position, score in pssm.search(seq, threshold=threshold, both=True):
            start = position if position >= 0 else len(seq) + position
            strand = "+" if position >= 0 else "-"
            rows.append(
                f"{record.id}\t{start}\t{start + len(motif)}\t{jaspar_id}\t"
                f"{score:.4f}\t{strand}\n"
            )
    return rows


with open(log_path, "w") as log:
    if mode == "fixture":
        if not os.path.exists(fixture):
            raise SystemExit(
                f"fixture {fixture} is missing.\n"
                "Either restore it, or set  fixtures.ctcf_motifs: recompute ."
            )
        rows = _read_fixture(fixture)
        _log(log, f"mode=fixture  source={fixture}  intervals={len(rows)}")

        with open(snakemake.output.chrx, "w") as fh:
            fh.writelines(rows)

        # The fixture is chrX-only by construction.
        shutil.copyfile(snakemake.output.chrx, snakemake.output.genome_wide)
        _log(
            log,
            "NOTE: the committed motif fixture covers "
            f"{chrx_name} only, so the 'genome-wide' motif file is a copy of "
            "the chrX one. Every paper panel is chrX, so this affects nothing "
            "downstream; only the opt-in QC enrichment rule reads the "
            "genome-wide file. Set fixtures.ctcf_motifs: recompute for a true "
            "genome-wide scan.",
        )

    elif mode == "recompute":
        rows = _scan_fasta(snakemake.input.fa, snakemake.params.jaspar, log)
        rows.sort(key=lambda r: (r.split("\t")[0], int(r.split("\t")[1])))
        with open(snakemake.output.genome_wide, "w") as fh:
            fh.writelines(rows)
        chrx_rows = [r for r in rows if r.split("\t")[0] == chrx_name]
        with open(snakemake.output.chrx, "w") as fh:
            fh.writelines(chrx_rows)
        _log(
            log,
            f"mode=recompute  genome_wide={len(rows)}  {chrx_name}={len(chrx_rows)}\n"
            "WARNING: a recomputed motif set will NOT match the committed "
            "fixture, so valley-boundary CTCF status will differ from the "
            "published figures.",
        )

    else:
        raise SystemExit(
            f"unknown fixtures.ctcf_motifs mode {mode!r}; expected "
            "'fixture' or 'recompute'"
        )
