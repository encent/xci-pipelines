#!/usr/bin/env bash
# XCI PIPELINE — edit the configuration, then check it makes sense.
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$REPO"
CONFIG="config/config.yaml"

cat <<'EOF'

Opening config/config.yaml.

  You only need the blocks marked  [EDIT ME] :

    paths.inputs   where your BAM files and Hi-C coolers are
    paths.data_dir where the pipeline should write (needs ~600 GB)
    genome         leave as 'auto' to download the mm10 reference
    loci           the two chrX regions; leave alone to reproduce the paper
    figures.panels 'all', 'paper', or a list like [Fig3a, Fig6]

  Save and close the editor when you are done.

EOF
read -r -p "Press Enter to open the editor. " _ || true
"${EDITOR:-nano}" "$CONFIG"

echo
echo "Checking your configuration..."
python3 - <<'PY'
import sys, os
sys.path.insert(0, "workflow")
import yaml
from lib.paths import Paths
from lib import samples as S

cfg = yaml.safe_load(open("config/config.yaml"))
P = Paths(cfg)

try:
    ss = S.load(cfg)
except Exception as exc:
    print(exc); sys.exit(1)

n_panels = 0
if os.path.exists(cfg["paths"].get("panels", "config/panels.yaml")):
    reg = yaml.safe_load(open(cfg["paths"]["panels"])) or []
    if isinstance(reg, dict):
        reg = reg.get("panels", [])
    n_panels = sum(1 for p in reg if p.get("paper_ref"))

missing = [f"{k}: {v}" for k, v in cfg["paths"]["inputs"].get("bam_dirs", {}).items()
           if not os.path.isdir(v)]

print()
print("=" * 72)
print("OK — configuration is valid.")
print("=" * 72)
print(f"  {len(ss.all_samples)} samples in {len(ss.normgroups)} normalisation groups")
print(f"  {len(ss.clones)} clones, {len(ss.marks)} marks, {len(ss.tracks)} merged tracks")
print(f"  {len(ss.cooler_names)} Hi-C contact maps")
print(f"  chromosomes: {', '.join(cfg['chromosomes']['downstream'])} "
      f"(normalisation on {len(cfg['chromosomes']['csaw'])} autosomes)")
if n_panels:
    print(f"  {n_panels} paper panels will be produced")
print(f"  results will be written to: {P.results}")
if missing:
    print()
    print("  WARNING — these input directories do not exist yet:")
    for m in missing:
        print(f"    {m}")
    print("  The pipeline will skip whatever depends on them.")
print()
print("  Estimated for a full run: about 2-3 days, about 500 GB.")
print("  Try  ./run.sh demo  first (about 2 hours).")
print()
PY
