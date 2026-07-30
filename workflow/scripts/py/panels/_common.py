"""Shared plumbing for every ``panel_*`` script.

Panel scripts are the one place in the pipeline that reads *other modules'*
output schemas, so they are the one place most exposed to an upstream column
rename. The readers below are therefore deliberately tolerant about column
names and strict about meaning: a missing column raises with the file path, the
columns that were found, and the column that was wanted, because "KeyError:
'category'" at hour 30 of a run tells you nothing.

The schemas panels depend on, as of the modules that produced them:

``P.gene_content(track)``        chrom start end size has_escapee has_silent
                                has_gene category size_class
``P.boundary_ctcf(track)``       chrom start end side valley_class is_escaping
                                has_ctcf fisher_p fisher_odds_ratio
                                n_boundaries n_escaping
``P.boundary_profile(t, sig)``   npz: matrix background matrix_final profile
                                profile_raw chroms positions sides n_used
``P.boundary_dtag(clone, sig)``  npz: profile_{subset}_{nodtag,dtag},
                                matrix_*, n_{subset}, positions_*, sides_*
                                subset in {escaping, escaping_ctcf}
``P.valley_overlap(clone, s)``   clone scope n_nodtag n_dtag nodtag_overlap
                                nodtag_only dtag_overlap dtag_only pct_*
``P.density(...)``               chrom start end nodtag dtag log2_ratio
                                pearson_r ccc median_lfc
``P.valley_states(chrom, t)``    chrom start end value state
``P.valley_sizes()``             track mark clone condition allele chrom start
                                end size size_class
``P.loop_stats(cmp, roi, sfx)``  file experiment loop_set n_loops median mean std
``P.saddle(roi, name, eig)``     npz: saddle interaction_sum interaction_count
                                saddle_strength_profile score top_left ...
``P.saddle_values(roi, name)``   file E1 E2 E3          (the CORNER score)
``P.saddle_strength_selected``   file value
``P.pileup(...)``                npz: data n_loops flank n
``P.eigs(roi, name)``            chrom start end E1 E2 E3

Keep this docstring in step with reality; it is the contract the other
re-developers are coding against from the other side.
"""

from __future__ import annotations

import os
import sys
from contextlib import contextmanager

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# bootstrap
# ---------------------------------------------------------------------------
def bootstrap(smk):
    """Put ``workflow/`` on the path so ``from lib import plotting`` works.

    The rule passes ``params.spec["libdir"]`` rather than this script guessing
    from ``__file__``, because Snakemake may copy or exec the script from a
    different directory depending on how it was invoked.
    """
    spec = dict(smk.params.spec)
    libdir = spec.get("libdir")
    if libdir and libdir not in sys.path:
        sys.path.insert(0, libdir)
    return spec


@contextmanager
def logging(smk):
    """Tee stdout and stderr into the rule's log file.

    Snakemake redirects a ``shell:`` block's output but not a ``script:``
    block's, so without this a panel's diagnostics vanish.
    """
    path = smk.log[0] if smk.log else None
    if not path:
        yield lambda msg: None
        return
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w") as fh:
        def say(msg=""):
            fh.write(str(msg) + "\n")
            fh.flush()

        try:
            yield say
        except Exception as exc:
            fh.write("\nFAILED: {}: {}\n".format(type(exc).__name__, exc))
            import traceback
            traceback.print_exc(file=fh)
            raise


def outputs(smk):
    return str(smk.output[0])


def ensure_parent(path):
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    return path


# ---------------------------------------------------------------------------
# tolerant readers
# ---------------------------------------------------------------------------
def read_table(path, **kw):
    """Read a TSV, raising with context rather than a bare pandas error."""
    try:
        return pd.read_csv(path, sep="\t", **kw)
    except Exception as exc:
        raise RuntimeError(
            "cannot read {}: {}: {}".format(path, type(exc).__name__, exc))


def need(df, column, path, alternatives=()):
    """Return ``df[column]``, or the first present alternative, or explain."""
    for name in (column,) + tuple(alternatives):
        if name in df.columns:
            return df[name]
    raise KeyError(
        "{} has no column {!r}{}; it has {}. The panel and its producing rule "
        "have drifted apart -- see workflow/scripts/py/panels/_common.py for "
        "the schema panels expect.".format(
            path, column,
            " (nor {})".format(", ".join(map(repr, alternatives)))
            if alternatives else "",
            ", ".join(map(repr, df.columns)),
        )
    )


def read_npz(path):
    try:
        return np.load(path, allow_pickle=False)
    except Exception as exc:
        raise RuntimeError(
            "cannot read {}: {}: {}".format(path, type(exc).__name__, exc))


def track_label(track):
    """``H3K27me3_E6_WT_Xi`` -> ``E6_WT``, the label the original figures used."""
    parts = str(track).split("_")
    if len(parts) >= 4:
        return "_".join(parts[1:-1])
    return str(track)


def clone_of(track):
    parts = str(track).split("_")
    return parts[1] if len(parts) > 1 else str(track)


def condition_of(track):
    parts = str(track).split("_")
    return parts[2] if len(parts) > 2 else ""


def display_locus(locus, display_names):
    """`Jarid` is the code name; the paper says Kdm5c."""
    return (display_names or {}).get(locus, locus)


def group_inputs(paths, *tokens):
    """Filter an input list down to the paths containing every token."""
    out = []
    for path in map(str, paths):
        base = os.path.basename(path)
        if all(t in base for t in tokens):
            out.append(path)
    return out


def find_input(paths, *tokens):
    hit = group_inputs(paths, *tokens)
    return hit[0] if hit else None


# ---------------------------------------------------------------------------
# n_items, for the manifest
# ---------------------------------------------------------------------------
def write_n_items(smk, n_items, extra=None):
    """Record how many things a panel drew, next to its output.

    ``figure_manifest.tsv`` needs `n_items` per panel and cannot recompute it
    without re-reading every input, so each panel drops a one-line sidecar and
    the manifest collects them. The sidecar is not a Snakemake output: it is
    advisory, and a panel that crashes before writing it still produces a valid
    manifest row.
    """
    path = str(smk.output[0]) + ".n_items"
    try:
        ensure_parent(path)
        with open(path, "w") as fh:
            fh.write("{}\n".format(n_items))
            if extra:
                for key, value in dict(extra).items():
                    fh.write("{}\t{}\n".format(key, value))
    except Exception:
        pass
    return n_items
