"""BLAS thread control, and why it is not applied uniformly.

Defect V-C: `call_valleys` output depends on the ambient BLAS thread count.
hmmlearn 0.3.3 initialises `GaussianHMM` emission means with
`sklearn.cluster.KMeans` (`hmmlearn/hmm.py:311-315`), whose chunked parallel
reduction sums partial results in a thread-count-dependent order. A different
initialisation lands EM in a different local optimum, and a 2-state HMM turns
that into whole valleys appearing or vanishing:

    threads=1  378 valleys                    threads=2  377  == ground truth
    threads=4  377  == ground truth           threads=8  320, means REVERSED

Not the `n_iter=15` cap -- every fit reports `converged=True`.

**Neither `threads:` on a Snakemake rule nor env vars are the right tool.**
`threads:` only tells the scheduler how to allocate slots; it never reaches
OpenBLAS. Env vars work but are global, are read once at numpy import, and
cannot be scoped to one computation. `threadpoolctl.threadpool_limits` is the
correct instrument.

**Never pass `user_api="blas"`.** The OpenMP runtime is load-bearing too, and
capping BLAS alone fails silently:

    E6 Xi   BLAS=2 / OpenMP=48  ->  320 valleys
    E6 Xi   BLAS=2 / OpenMP=2   ->  377, byte-identical to ground truth

A test harness did exactly that and produced a plausible but wrong sweep table.
A bare `threadpool_limits(limits=N)` caps every runtime, which is what is
wanted. The config key is called `blas_threads` for continuity with the frozen
decision record; read it as "numeric threads".

Three controls must hold together, and any one alone reopens V-C: the thread
pin, scikit-learn's version pin (KMeans is the mechanism), and capping **both**
runtimes.

**Why `observe()` reports the effective count rather than the requested one.**
This is not defensive decoration. Logging the observed per-runtime count is what
caught the `user_api="blas"` defect above -- within hours of the logging being
introduced, because the log showed `openmp/openmp: 48` next to a request of 2.
Had it logged the request, the harness would have looked correct and its wrong
sweep table would have been believed. Do not simplify this to echo the request.

Why this is NOT applied uniformly across the pipeline
-----------------------------------------------------
It is tempting to pin every rule "for safety". That would be wrong here.

The Hi-C branch currently reproduces ground truth **byte-identically at the
ambient thread count** -- `gc_track` 27/27 exact, `Saddle_values` corner scores
81/81 at worst relative difference 0, eigenvector reconstruction to 1.91e-14.
Pinning those rules to a different count could *change* results that are
currently exact. A safety measure that risks breaking a verified reproduction is
not a safety measure.

So the policy is asymmetric, and deliberately so:

* `call_valleys` -- **pinned**, because it is demonstrably sensitive and
  unpinned it is machine-dependent.
* Hi-C and METALoci numerics -- **observed and recorded, not pinned**
  (`hic.blas_threads: null`, `metaloci.blas_threads: null`). They are known
  exact at the ambient count. If a future run diverges, the recorded thread
  count makes it diagnosable in one step instead of being a mystery.

The distinction rests on evidence, not on which feels safer, and it should be
revisited only with a measurement -- not with an intuition.
"""

from __future__ import annotations

import contextlib
import os
from typing import Dict, Optional, Tuple

try:
    import threadpoolctl
except ImportError:  # pragma: no cover
    threadpoolctl = None


def observe() -> Tuple[Dict[str, int], Dict[str, str]]:
    """Thread counts and library paths AS REPORTED, not as requested.

    Requested and effective counts can differ, so callers must log this rather
    than log what they asked for.

    Note: never use `nproc` to detect the ambient default -- it honours
    `OMP_NUM_THREADS` and will happily report the value someone already
    exported rather than the machine's core count.
    """
    if threadpoolctl is None:
        return {}, {}
    info = threadpoolctl.threadpool_info()
    counts = {f"{d['user_api']}/{d['internal_api']}": d["num_threads"] for d in info}
    libs = {d["internal_api"]: os.path.basename(d["filepath"] or "?") for d in info}
    return counts, libs


@contextlib.contextmanager
def limits(n: Optional[int], say=print, what: str = "computation"):
    """Limit BLAS/OpenMP threads to `n` for the enclosed block.

    `n=None` means DO NOT PIN -- observe and record only. That is the correct
    setting for a computation already verified exact at the ambient count.
    """
    counts, libs = observe()
    if n is None:
        say(f"BLAS [{what}]: NOT pinned (observe-only). observed={counts} libs={libs}")
        yield
        return

    if threadpoolctl is None:
        say(
            f"BLAS [{what}]: threadpoolctl is unavailable, so the requested "
            f"limit of {n} CANNOT be enforced. Output may not be reproducible "
            "across machines. Install threadpoolctl."
        )
        yield
        return

    with threadpoolctl.threadpool_limits(limits=n):
        counts, libs = observe()
        say(f"BLAS [{what}]: requested={n} observed={counts} libs={libs}")
        for api, got in counts.items():
            if got != n:
                say(
                    f"WARNING: {api} reports {got} threads, not the requested {n}. "
                    "Results will not be comparable to a run at the requested "
                    "count. Investigate before trusting this output."
                )
        yield
