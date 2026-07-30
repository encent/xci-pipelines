#!/usr/bin/env bash
# =====================================================================
#  XCI PIPELINE — install.
#
#  Run this once:   bash install.sh
#  Takes about 10 minutes and needs an internet connection.
#  Safe to re-run; it skips whatever is already done.
#
#  Uses mamba rather than conda to resolve dependencies -- see the
#  "Choosing a package installer" step for why that is not a preference.
# =====================================================================
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO"

ENV_NAME="xci-pipeline"

say()  { printf '\n\033[1m%s\033[0m\n' "$*"; }
ok()   { printf '  [ok] %s\n' "$*"; }
warn() { printf '  [!]  %s\n' "$*"; }

# ------------------------------------------------------------- conda ----
say "1/6  Checking conda"
if ! command -v conda >/dev/null 2>&1; then
  cat >&2 <<'EOF'

========================================================================
ERROR  conda was not found
========================================================================

  The 'conda' command is not on your PATH.

  Why this matters:
    The pipeline installs its software with conda so that you get the
    exact package versions used for the paper. Without conda we cannot
    guarantee you are running the same code we did.

  What to do:
    Install Miniforge (free, no licence restrictions, and it already
    includes the fast installer this script wants):

        https://github.com/conda-forge/miniforge

    Then CLOSE AND REOPEN your terminal and run:

        bash install.sh

========================================================================

EOF
  exit 1
fi
ok "conda $(conda --version 2>/dev/null | awk '{print $2}')"

# ----------------------------------------------------------- solver ----
# Which program actually resolves the dependencies.
#
# This matters far more than it looks. This environment pins ~40 exact
# versions across conda-forge and bioconda, and conda's built-in "classic"
# solver does not cope: on the machine this pipeline was developed on it ran
# for over 40 minutes and never finished, printing only
# "Solving environment: ...working..." the whole time -- no error, no
# progress, nothing to distinguish it from a hang. mamba solved the identical
# file in about two minutes.
#
# So: use mamba if it exists, else micromamba, else fetch a standalone
# micromamba (a single static binary, dropped inside this repo, touching
# nothing else on your system). Falling back to plain conda is the last
# resort and warns you first.
say "2/6  Choosing a package installer"
MM_LOCAL="$REPO/.mamba/bin/micromamba"

# BUG-E: Snakemake shells out to a LITERAL `mamba env create`, and micromamba
# has no `env create` subcommand -- it aborts on the prefix directory Snakemake
# pre-creates. So --conda-frontend mamba is not sufficient on a micromamba-only
# box. If the run will need --use-conda (it will: metaloci and coolbox have
# their own envs), REAL mamba is required.
if command -v mamba >/dev/null 2>&1; then
  SOLVER="mamba"; SOLVER_KIND="mamba"
  ok "using mamba $(mamba --version 2>/dev/null | head -1 | awk '{print $2}')"
  MAMBA_REAL=1
elif command -v micromamba >/dev/null 2>&1; then
  SOLVER="micromamba"; SOLVER_KIND="micromamba"
  ok "using micromamba $(micromamba --version 2>/dev/null)"
elif [ -x "$MM_LOCAL" ]; then
  SOLVER="$MM_LOCAL"; SOLVER_KIND="micromamba"
  ok "using the micromamba downloaded by a previous run"
else
  echo "  Neither mamba nor micromamba found. Downloading micromamba"
  echo "  (one static binary, ~5 MB, installed into $REPO/.mamba --"
  echo "   it does not modify your conda installation)."
  MM_ARCH="linux-64"
  case "$(uname -s)-$(uname -m)" in
    Darwin-arm64) MM_ARCH="osx-arm64" ;;
    Darwin-x86_64) MM_ARCH="osx-64" ;;
    Linux-aarch64) MM_ARCH="linux-aarch64" ;;
  esac
  if mkdir -p "$REPO/.mamba" && \
     curl -fsSL "https://micro.mamba.pm/api/micromamba/${MM_ARCH}/latest" \
       | tar -xvj -C "$REPO/.mamba" bin/micromamba >/dev/null 2>&1 && \
     [ -x "$MM_LOCAL" ]; then
    SOLVER="$MM_LOCAL"; SOLVER_KIND="micromamba"
    ok "micromamba ready"
  else
    SOLVER="conda"; SOLVER_KIND="conda"
    warn "could not fetch micromamba -- falling back to conda."
    echo "       EXPECT THIS TO BE SLOW: possibly 30-60 minutes, and it may"
    echo "       appear to hang on 'Solving environment'. That is normal for"
    echo "       conda's classic solver on an environment this pinned."
    echo "       To avoid it, install Miniforge (which bundles mamba) and"
    echo "       re-run this script."
  fi
fi

# --------------------------------------------------------- main env ----
say "3/6  Creating the '$ENV_NAME' environment"
ENV_PREFIX="$(conda info --base)/envs/$ENV_NAME"
if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME" || [ -d "$ENV_PREFIX" ]; then
  ok "'$ENV_NAME' already exists -- skipping (delete it with"
  echo "       conda env remove -n $ENV_NAME    if you want a clean rebuild)"
else
  echo "  This is the slow part: roughly 1 GB downloaded."
  echo "  With mamba it takes a few minutes."
  case "$SOLVER_KIND" in
    micromamba)
      "$SOLVER" create -y -p "$ENV_PREFIX" -f environment.yml -r "$(conda info --base)"
      ;;
    mamba)
      "$SOLVER" env create -f environment.yml
      ;;
    conda)
      conda env create -f environment.yml
      ;;
  esac
  ok "created"
fi

# shellcheck disable=SC1091
CONDA_BASE="$(conda info --base)"
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"

# ----------------------------------------------------- satellite envs ---
if [ "${MAMBA_REAL:-0}" != "1" ]; then
  cat >&2 <<'EOF'

  WARNING  Snakemake needs REAL mamba to build the helper environments.

    It shells out to a literal `mamba env create`. micromamba has no
    `env create` subcommand, so it aborts on the prefix directory Snakemake
    pre-creates -- and --conda-frontend mamba does not avoid that.

    The main xci-pipeline environment is fine. The metaloci and coolbox
    environments are NOT, which means Fig 6 and Fig 4e cannot be produced.

    To fix:   conda install -n base -c conda-forge mamba
    then re-run: bash install.sh

EOF
fi

say "4/6  Building the two helper environments"
echo "  These are managed by the pipeline; you never need to name them."
# Snakemake defaults to conda for its own envs, which would reintroduce the
# slow solver we just avoided. Point it at mamba when we have one.
CONDA_FRONTEND="conda"
if [ "$SOLVER_KIND" = "mamba" ]; then
  CONDA_FRONTEND="mamba"
elif [ "$SOLVER_KIND" = "micromamba" ]; then
  # Snakemake 7 accepts `mamba` as a frontend name; expose micromamba as one.
  mkdir -p "$REPO/.mamba/shim"
  printf '#!/usr/bin/env bash\nexec "%s" "$@"\n' "$SOLVER" > "$REPO/.mamba/shim/mamba"
  chmod +x "$REPO/.mamba/shim/mamba"
  export PATH="$REPO/.mamba/shim:$PATH"
  CONDA_FRONTEND="mamba"
fi
if snakemake --use-conda --conda-frontend "$CONDA_FRONTEND" \
     --conda-create-envs-only --cores 1 figures >/dev/null 2>&1; then
  ok "metaloci and coolbox environments ready"
else
  warn "could not pre-build them now. That is not fatal -- they will be built"
  echo "       on first use. If it fails again then, see docs/TROUBLESHOOTING.md"
fi

# --------------------------------------------------------- sheets ------
say "5/6  Building the sample sheets"
if [ -s config/samples.tsv ] && [ -s config/coolers.tsv ]; then
  ok "config/samples.tsv and config/coolers.tsv already exist -- keeping them"
  echo "       (re-make them with: python3 workflow/scripts/py/make_sheets.py --force)"
else
  python3 workflow/scripts/py/make_sheets.py --config config/config.yaml
fi

# ------------------------------------------------------- preflight -----
say "6/6  Checking your setup"
set +e
python3 tests/preflight.py
PF=$?
set -e

# --------------------------------------------------------- next --------
echo
echo "========================================================================"
if [ $PF -eq 0 ]; then
  echo "Installed."
  echo "========================================================================"
  cat <<'EOF'

  Next:

    1.  conda activate xci-pipeline

    2.  bash configure.sh          Tell the pipeline where your data is.
                                   (Skip this if you are on the machine the
                                   paper was analysed on -- the defaults are
                                   already correct.)

    3.  ./run.sh demo              A small ~2 hour run, to check it works.
        ./run.sh                   The real thing.

  A full run takes about 2-3 days on a 48-core desktop and writes roughly
  500 GB. ./run.sh demo first is strongly recommended.

EOF
else
  echo "Installed, but the setup check found problems."
  echo "========================================================================"
  cat <<'EOF'

  The messages above say what is wrong and how to fix it. Most often it is
  a path in config/config.yaml pointing at data that has moved.

  Fix those, then re-run:   bash install.sh

EOF
fi
