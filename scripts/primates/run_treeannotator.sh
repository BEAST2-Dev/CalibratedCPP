#!/usr/bin/env bash
# Summarise every .trees file in a folder with TreeAnnotator (CA heights, CCD0 topology).
# The summary lands next to its input as <stem>_summary.tree.
#
# Usage: ./run_treeannotator.sh [-f] [-b BURNIN] [DIR ...]
#   -f          re-annotate even if the summary is newer than the .trees file
#   -b BURNIN   burn-in percentage (default 10)
#   DIR         folders to scan (default: data-wait)

set -euo pipefail

TREEANNOTATOR="${TREEANNOTATOR:-$HOME/WorkSpace/beast3/bin/treeannotator}"
BURNIN=10
FORCE=0

while getopts "fb:h" opt; do
    case "$opt" in
        f) FORCE=1 ;;
        b) BURNIN="$OPTARG" ;;
        h) sed -n '2,10p' "$0"; exit 0 ;;
        *) exit 2 ;;
    esac
done
shift $((OPTIND - 1))

if [ ! -x "$TREEANNOTATOR" ]; then
    echo "treeannotator not found at $TREEANNOTATOR (override with \$TREEANNOTATOR)" >&2
    exit 1
fi

cd "$(dirname "$0")"
DIRS=("${@:-data-wait}")

for dir in "${DIRS[@]}"; do
    [ -d "$dir" ] || { echo "skipping $dir: not a directory" >&2; continue; }
    shopt -s nullglob
    for trees in "$dir"/*.trees; do
        out="${trees%.trees}_summary.tree"
        if [ "$FORCE" -eq 0 ] && [ -s "$out" ] && [ "$out" -nt "$trees" ]; then
            echo "up to date: $out"
            continue
        fi
        echo "annotating: $trees -> $out"
        "$TREEANNOTATOR" -height CA -topology CCD0 -burnin "$BURNIN" "$trees" "$out"
    done
    shopt -u nullglob
done
