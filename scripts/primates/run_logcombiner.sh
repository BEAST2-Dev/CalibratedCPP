#!/usr/bin/env bash
# Combine replicate BEAST runs (-rep1/-rep2/-rep3) into one trace and tree log
# per analysis, using the BEAST 3 LogCombiner.
#
# Usage: ./run_logcombiner.sh [RUN_DIR] [BURNIN_PERCENT]
#   RUN_DIR         directory holding <name>-rep<N>.log/.trees (default grid/out)
#   BURNIN_PERCENT  burnin discarded from each replicate (default 10)
#
# Output goes to $RUN_DIR/combined/<name>-combined.{log,trees}; runs that already
# have a complete output are skipped, so the script can be re-run to resume.

set -uo pipefail

LOGCOMBINER="${LOGCOMBINER:-$HOME/WorkSpace/beast3/bin/logcombiner}"
REPS="${REPS:-1 2 3}"
# LogCombiner sometimes writes its output but fails to exit; give up on a job
# once it reports "Wrote ... lines", or after this many seconds.
TIMEOUT="${TIMEOUT:-1200}"

cd "$(dirname "$0")"
RUN_DIR="${1:-grid/out}"
BURNIN="${2:-10}"

if [ ! -x "$LOGCOMBINER" ]; then
    echo "logcombiner not found at $LOGCOMBINER (override with \$LOGCOMBINER)" >&2
    exit 1
fi
if [ ! -d "$RUN_DIR" ]; then
    echo "run directory not found: $RUN_DIR" >&2
    exit 1
fi

cd "$RUN_DIR"
mkdir -p combined
status=0

for rep1 in *-rep1.log; do
    [ -e "$rep1" ] || { echo "no *-rep1.log files in $RUN_DIR" >&2; exit 1; }
    name="${rep1%-rep1.log}"

    for ext in log trees; do
        out="combined/$name-combined.$ext"
        err="combined/$name-$ext.stderr"

        inputs=()
        missing=""
        for r in $REPS; do
            f="$name-rep$r.$ext"
            [ -f "$f" ] && inputs+=(-log "$f") || missing="$missing $f"
        done
        if [ -n "$missing" ]; then
            echo "SKIP $name.$ext (missing:$missing)" >&2
            status=1
            continue
        fi

        if [ -s "$out" ] && grep -q '^Wrote' "$err" 2>/dev/null; then
            echo "skip $out (already complete)"
            continue
        fi

        echo "combining $name.$ext -> $out"
        "$LOGCOMBINER" -b "$BURNIN" "${inputs[@]}" -o "$out" >/dev/null 2>"$err" &
        pid=$!
        waited=0
        while [ "$waited" -lt "$TIMEOUT" ]; do
            grep -q '^Wrote' "$err" 2>/dev/null && break
            kill -0 "$pid" 2>/dev/null || break
            sleep 2
            waited=$((waited + 2))
        done
        sleep 2
        pkill -P "$pid" 2>/dev/null
        kill "$pid" 2>/dev/null
        wait "$pid" 2>/dev/null

        if grep -q '^Wrote' "$err" 2>/dev/null; then
            grep '^Wrote' "$err"
        else
            echo "FAILED $name.$ext (see $RUN_DIR/$err)" >&2
            status=1
        fi
    done
done

exit $status
