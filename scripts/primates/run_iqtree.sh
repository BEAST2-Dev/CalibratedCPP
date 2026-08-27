#!/usr/bin/env bash
# Maximum-likelihood tree for the nogapN alignment under the model used in the
# BEAST analyses (HKY + discrete gamma, 4 categories; see
# calibratedcpp-lphy/examples/primates_conditionedMRCAPrior.lphy).
#
# Results land in iqtree/; iqtree/<alignment>.treefile is the ML tree,
# rooted on Mus musculus so it can be compared with compare_rf.py.
#
# Usage: ./run_iqtree.sh [ALIGNMENT] [extra iqtree args...]
#   ALIGNMENT   basename in calibratedcpp-lphy/examples/data (default primates_nogapN)

set -euo pipefail

IQTREE="${IQTREE:-$HOME/WorkSpace/iqtree-3.1.3-macOS/bin/iqtree3}"
MODEL="${MODEL:-HKY+G4}"

cd "$(dirname "$0")"
NAME="${1:-primates_nogapN}"
[ $# -gt 0 ] && shift
ALN="../../calibratedcpp-lphy/examples/data/$NAME.nex"

if [ ! -f "$ALN" ]; then
    echo "alignment not found: $ALN" >&2
    exit 1
fi

if [ ! -x "$IQTREE" ]; then
    echo "iqtree not found at $IQTREE (override with \$IQTREE)" >&2
    exit 1
fi

mkdir -p iqtree
"$IQTREE" -s "$ALN" -m "$MODEL" -o Mus_musculus -B 1000 -T AUTO \
          --prefix "iqtree/$NAME" -redo "$@"
