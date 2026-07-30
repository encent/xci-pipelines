#!/usr/bin/env bash
# =====================================================================
#  XCI PIPELINE — the run command.
#
#      ./run.sh                 everything, ending in all paper figures
#      ./run.sh demo            small test run, about 2 hours
#      ./run.sh panel Fig3a     just one paper panel
#      ./run.sh --help          the full list
# =====================================================================
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO"

CONFIG="config/config.yaml"
PROFILE="profiles/desktop"
EXTRA=()

# Snakemake 7.32 spells conda deployment --use-conda.
SM_CONDA="--use-conda"

usage() {
  cat <<'EOF'
XCI pipeline

  ./run.sh                    Run everything and produce all paper figures.
  ./run.sh figures            Same as above (the default target).
  ./run.sh demo               Small chrX-only run on 2 clones, about 2 hours.
                              Use this first to check your setup works.

  ./run.sh signal             Just the normalised signal tracks.
  ./run.sh contacts           Just the Hi-C contact maps.
  ./run.sh valleys            Just the H3K27me3 valleys (and what they need).
  ./run.sh hic                Loops, compartments and METALoci.

  ./run.sh panel Fig3a        Just one paper panel, named the way the paper
                              names it. Try:  Fig2g  Fig3a  Fig6  EFig3a
  ./run.sh report             Rebuild the HTML report from existing results.
  ./run.sh clean-figures      Delete and redraw the figures. Cheap: this never
                              touches the expensive intermediate files.

  ./run.sh dryrun             Show what would run, without running it.
  ./run.sh --help             This message.

Anything after a "--" is passed straight through to snakemake, e.g.
  ./run.sh -- --cores 8
EOF
}

die() { printf '\n%s\n\n' "$*" >&2; exit 1; }

need_env() {
  if ! command -v snakemake >/dev/null 2>&1; then
    cat >&2 <<EOF

========================================================================
ERROR  The xci-pipeline environment is not active.
========================================================================

  Why this matters:
    The pipeline runs inside a conda environment that pins the exact
    software versions used for the paper. Without it, results would not
    be reproducible even if they looked plausible.

  What to do:

    conda activate xci-pipeline
    ./run.sh

  If that fails because the environment does not exist yet:

    bash install.sh

========================================================================

EOF
    exit 1
  fi
}

resolve_panel() {
  # Map a paper reference (Fig3a) to the panel output file, via panels.yaml.
  local ref="$1"
  python3 - "$ref" <<'PY'
import sys, os, re
sys.path.insert(0, "workflow")
try:
    import yaml
except ImportError:
    sys.exit("PANEL_LOOKUP_UNAVAILABLE")

ref = sys.argv[1]
norm = lambda s: re.sub(r"[^a-z0-9]", "", s.lower())

if not os.path.exists("config/panels.yaml"):
    sys.exit("PANEL_REGISTRY_MISSING")

reg = yaml.safe_load(open("config/panels.yaml")) or []
if isinstance(reg, dict):
    reg = reg.get("panels", [])

hits = [p for p in reg
        if norm(ref) in {norm(r) for r in (p.get("paper_ref") or [])}
        or norm(ref) == norm(p.get("id", ""))]
if not hits:
    known = sorted({r for p in reg for r in (p.get("paper_ref") or [])})
    sys.stderr.write(
        "\n========================================================================\n"
        f"ERROR  Unknown figure panel: {ref}\n"
        "========================================================================\n\n"
        "  What to do:\n"
        "    Use one of the panel names below, spelled as the paper spells it.\n\n"
        "    " + "  ".join(known) + "\n\n"
        "    The full map from panel to output file is in results/figure_manifest.tsv\n"
        "    and, after a run, in results/figure_manifest.tsv\n\n"
        "========================================================================\n\n")
    sys.exit(2)

from lib.paths import Paths
cfg = yaml.safe_load(open("config/config.yaml"))
P = Paths(cfg)
for p in hits:
    for fmt in (p.get("formats") or ["svg"]):
        print(P.figures(p.get("group", "misc"), p["id"], fmt))
PY
}

# ---------------------------------------------------------------- args --
CMD="${1:-figures}"
[ $# -gt 0 ] && shift || true

TARGETS=()
case "$CMD" in
  -h|--help|help) usage; exit 0 ;;
  figures|"")     TARGETS=(figures) ;;
  signal|contacts|valleys|hic|reference)
                  TARGETS=("$CMD") ;;
  demo)           CONFIG="config/config.demo.yaml"; TARGETS=(figures) ;;
  dryrun)         EXTRA+=(--dry-run --quiet); TARGETS=(figures) ;;
  report)         TARGETS=(report_only) ;;
  clean-figures)
      need_env
      RESULTS="$(python3 -c "
import sys,yaml; sys.path.insert(0,'workflow')
from lib.paths import Paths
print(Paths(yaml.safe_load(open('config/config.yaml'))).out('figures'))")"
      printf 'Removing %s and redrawing...\n' "$RESULTS"
      rm -rf "$RESULTS"
      TARGETS=(figures)
      ;;
  panel)
      [ $# -ge 1 ] || die "Usage: ./run.sh panel Fig3a"
      mapfile -t TARGETS < <(resolve_panel "$1") || exit $?
      shift
      ;;
  --) TARGETS=(figures) ;;
  *)  die "Unknown command: $CMD

Run  ./run.sh --help  to see what is available." ;;
esac

# pass-through after --
while [ $# -gt 0 ]; do
  if [ "$1" = "--" ]; then shift; EXTRA+=("$@"); break; fi
  EXTRA+=("$1"); shift
done

need_env

# --------------------------------------------------------------- run ----
export PYTHONHASHSEED=0          # determinism (plan section 1.2)

# V-C: do NOT export OMP_NUM_THREADS here. Thread control is per-rule, via
# threadpoolctl in the scripts that need it, because pinning globally would
# also pin the Hi-C rules that are verified exact at the ambient count.
# The effective count is observed and logged per job, and recorded in
# results/provenance.json.

mkdir -p "$(dirname "$REPO/results")" 2>/dev/null || true
DATA_RESULTS="$(python3 -c "
import sys,yaml; sys.path.insert(0,'workflow')
from lib.paths import Paths
print(Paths(yaml.safe_load(open('$CONFIG'))).results)" 2>/dev/null || true)"
if [ -n "${DATA_RESULTS:-}" ]; then
  mkdir -p "$DATA_RESULTS"
  # convenience symlink so a new user finds their figures next to the code
  [ -e "$REPO/results" ] || ln -s "$DATA_RESULTS" "$REPO/results" 2>/dev/null || true
fi

printf '\nRunning: %s\n' "${TARGETS[*]}"
printf 'Config : %s\n' "$CONFIG"
printf 'Results: %s\n\n' "${DATA_RESULTS:-<unresolved>}"

set +e
snakemake \
  --profile "$PROFILE" \
  --configfile "$CONFIG" \
  $SM_CONDA \
  --rerun-triggers mtime \
  --keep-going \
  "${EXTRA[@]}" \
  "${TARGETS[@]}"
RC=$?
set -e

# ------------------------------------------------------------ summary ---
echo
if [ $RC -eq 0 ]; then
  echo "========================================================================"
  echo "Done."
  echo "========================================================================"
  [ -n "${DATA_RESULTS:-}" ] && {
    echo "  Figures : $DATA_RESULTS/figures/"
    echo "  Report  : $DATA_RESULTS/REPORT.html      <- open this in a browser"
    echo "  Which file is which paper panel:"
    echo "            $DATA_RESULTS/figure_manifest.tsv"
    echo "  Every number quoted in the paper:"
    echo "            $DATA_RESULTS/paper_numbers.tsv"
    if [ -f "$DATA_RESULTS/figure_manifest.tsv" ]; then
      SKIPPED=$(awk -F'\t' 'NR>1 && $4!="ok"' "$DATA_RESULTS/figure_manifest.tsv" | wc -l)
      [ "$SKIPPED" -gt 0 ] && {
        echo
        echo "  $SKIPPED panel(s) were skipped or are low-confidence. See the"
        echo "  status and notes columns of figure_manifest.tsv."
      }
    fi
  }
  echo
else
  echo "========================================================================"
  echo "The run did not finish."
  echo "========================================================================"
  echo "  The error above says which step failed. Logs for every step are in:"
  echo "    ${DATA_RESULTS%/results}/logs/"
  echo
  echo "  Common causes and fixes:  docs/TROUBLESHOOTING.md"
  echo "  Nothing is lost -- re-running ./run.sh resumes where it stopped."
  echo
fi
exit $RC
