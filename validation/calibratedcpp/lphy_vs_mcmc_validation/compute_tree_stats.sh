#!/bin/bash
# Compute tree statistics for all trees files using BEASTAppLauncher TreeStat2.
#
# Output (.treestreestats.log) is written alongside each input file.
#
# Usage: bash compute_tree_stats.sh

set -e

DIR="$(cd "$(dirname "$0")" && pwd)"

CTRL="$(mktemp /tmp/treestat_ctrl.XXXXXX)"
cat > "$CTRL" << 'EOF'
treestat2.statistics.TreeLength
treestat2.statistics.TreeHeight
treestat2.statistics.GammaStatistic
treestat2.statistics.CollessIndex
treestat2.statistics.B1Statistic
treestat2.statistics.CherryStatistic
EOF

run_treestat() {
    local trees="$1"
    applauncher TreeStat "$CTRL" "$trees" 2>&1 \
        | grep -v "^Loading\|^Package\|^Unexpected\|behavior may follow" \
        || true
}

# LPhy trees have a 'begin taxa;...end;' block that causes duplicate-taxa
# warnings in BEAST's TreeParser — strip it before computing stats.
LPHY_FILES=(fixStemLPhy fixRootLPhy fix4LeafStemLPhy fix4LeafRootLPhy)
for base in "${LPHY_FILES[@]}"; do
    echo "=== $base ==="
    TMP="$DIR/${base}_notaxa.trees"
    sed '/^begin taxa/,/^end;/d' "$DIR/${base}.trees" > "$TMP"
    run_treestat "$TMP"
    mv "${TMP%.trees}.treestreestats.log" "$DIR/${base}.treestreestats.log"
    rm -f "$TMP"
done

MCMC_FILES=(
    fixStemMCMC fixStemMCMC_not_cond
    fixRootMCMC fixRootMCMC_not_cond
    fix4LeafStemMCMC fix4LeafStemMCMC_not_cond
    fix4LeafRootMCMC fix4LeafRootMCMC_not_cond
)
for base in "${MCMC_FILES[@]}"; do
    echo "=== $base ==="
    run_treestat "$DIR/${base}.trees"
done

rm -f "$CTRL"

echo ""
echo "Done. Log files:"
ls -1 "$DIR"/*.treestreestats.log 2>/dev/null || echo "  (none found)"
