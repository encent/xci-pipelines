"""Stage a dataset's signal bedGraphs under METALoci-legal basenames.

THE TRAP THIS RULE EXISTS TO CLOSE
----------------------------------
``metaloci prep`` takes the *name of a signal* from the **file basename**::

    "For single signal files, if header is omitted, the signal name will be
     the name of the file."                       -- metaloci prep --help

Nothing validates it. If the basename does not match what ``metaloci lm``
and ``metaloci figure`` later look up (``{mark}_{dataset_adj}``), the signal
is **silently dropped**: no error, no warning, just a Gaudi grid with one row
missing and a compartment-like strength computed over fewer signals. That is
the single easiest way to produce a plausible, wrong Fig 6.

So this script does two things and refuses to guess:

1. copies (symlinks) each source bedGraph to ``signal_input/{name}.bed``
   inside the dataset's METALoci working directory, where ``{name}`` is
   ``naming.metaloci_signal_name(mark, cooler)`` -- resolved by the rule, not
   here;
2. writes ``{dataset}.signals``, newline-separated, in the *same order*, which
   is what ``metaloci lm -s`` reads.

It hard-errors if the number of sources and names disagree, if a name repeats,
or if a name is not a legal METALoci signal token. A wrong wiring becomes a
crash at DAG time instead of a missing row at figure time.

Symlinks, not copies: the sources are 12-16 MB bedGraphs x ~200 dataset-signal
pairs. `readlink -f` resolution means the staged file's realpath is the true
source, which is how the testing team verifies this rule.

WHERE THE SOURCES COME FROM (section_reports/hic_metaloci.md section 7)
----------------------------------------------------------------------
wt         H3K27me3 / H3K27ac / CTCF / RNA-Seq  <- bedgraph5k (csaw-normalised)
           AcMe3                                <- acme3_bg5k
degron     same, degron tracks; Rad21 exists for B1621 only
consensus  all five                             <- consensus5k
           There is NO Rad21 consensus (pooling Rad21-dTAG with non-Rad21
           clones would be meaningless), so consensus datasets carry 5 signals.

Per-dataset subsetting is `SS.metaloci_signals_for()` in the rule file, which
reproduces the `01_10` cell-4 exclusions (C5C10 CUT&RUN dropped except RNA-Seq;
Rad21 only for B1621; CTCF never for B1621).
"""

import os
import sys

sources = list(snakemake.input.beds)
names = list(snakemake.params.signal_names)
signals_file = snakemake.output.signals
staged_dir = snakemake.output.staged
dataset = snakemake.params.dataset

os.makedirs(os.path.dirname(snakemake.log[0]), exist_ok=True)
log = open(snakemake.log[0], "w")


def fail(msg):
    log.write(f"ERROR {msg}\n")
    log.close()
    sys.exit(f"metaloci_signals[{dataset}]: {msg}")


if len(sources) != len(names):
    fail(f"{len(sources)} source bedGraphs but {len(names)} signal names")
if len(set(names)) != len(names):
    fail(f"duplicate signal names: {names}")
if not names:
    fail("no signals resolved -- every signal was excluded for this dataset")

os.makedirs(staged_dir, exist_ok=True)

for src, name in zip(sources, names):
    if "/" in name or name.endswith(".bed"):
        fail(f"illegal signal name {name!r}")
    dst = os.path.join(staged_dir, f"{name}.bed")
    real = os.path.realpath(src)
    if not os.path.exists(real):
        fail(f"source bedGraph missing: {src}")
    if os.path.islink(dst) or os.path.exists(dst):
        os.remove(dst)
    os.symlink(real, dst)
    log.write(f"staged {dst} -> {real}\n")
    # The basename is the signal name METALoci will use. Say it out loud in the
    # log so a mismatch is greppable rather than invisible.
    log.write(f"  signal name = {name}\n")

os.makedirs(os.path.dirname(signals_file), exist_ok=True)
with open(signals_file, "w") as fh:
    for name in names:
        fh.write(f"{name}\n")

log.write(f"wrote {len(names)} signal names to {signals_file}\n")
log.close()
