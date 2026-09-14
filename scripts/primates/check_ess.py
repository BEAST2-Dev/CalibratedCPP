#!/usr/bin/env python3
"""Effective sample sizes for the combined traces, using Tracer's estimator.

Reports the minimum ESS over the columns that matter (posterior, likelihood,
prior, tree height and every mrca.age) for each combined run, and separately for
each replicate if the -rep?.log files are still present, since a low per-replicate
ESS is hidden once chains are pooled.

Usage:
    python3 check_ess.py [RUN_DIR]        # default: grid/data
"""
import glob
import os
import sys

import numpy as np

THRESHOLD = 200
MAX_LAG = 2000


def ess(x):
    """Effective sample size, following beast.core.util.ESS."""
    x = np.asarray(x, dtype=float)
    n = len(x)
    if n < 10:
        return float("nan")
    d = x - x.mean()
    max_lag = min(n - 1, MAX_LAG)
    gamma = np.empty(max_lag)
    for lag in range(max_lag):
        gamma[lag] = np.dot(d[: n - lag], d[lag:]) / (n - lag)
    if gamma[0] == 0:
        return float("nan")
    var_stat = gamma[0]
    for lag in range(2, max_lag, 2):
        pair = gamma[lag - 1] + gamma[lag]
        if pair <= 0:
            break
        var_stat += 2.0 * pair
    return n * gamma[0] / var_stat


def columns(path, burnin=0.0):
    """{column: samples} for the trace columns worth checking."""
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
    if not header or header[0] != "Sample":
        return {}
    keep = [c for c in header
            if c in ("posterior", "likelihood", "prior", "tree.height")
            or c.startswith("mrca.age(")]
    idx = {c: header.index(c) for c in keep}
    data = {c: [] for c in keep}
    with open(path) as fh:
        fh.readline()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < len(header):
                continue
            for c in keep:
                try:
                    data[c].append(float(p[idx[c]]))
                except ValueError:
                    pass
    n = len(next(iter(data.values()), []))
    cut = int(n * burnin)
    return {c: np.array(v[cut:]) for c, v in data.items() if len(v[cut:]) > 10}


def report(paths, label, burnin):
    print("\n%s\n%s" % (label, "-" * len(label)))
    bad = []
    for path in paths:
        cols = columns(path, burnin)
        if not cols:
            continue
        scores = {c: ess(v) for c, v in cols.items()}
        worst = min(scores, key=lambda c: scores[c])
        n = len(next(iter(cols.values())))
        flag = "" if scores[worst] >= THRESHOLD else "   <-- below %d" % THRESHOLD
        print("  %-52s n=%5d  min ESS %7.0f  (%s)%s"
              % (os.path.basename(path), n, scores[worst], worst, flag))
        if scores[worst] < THRESHOLD:
            bad.append(os.path.basename(path))
    return bad


def main():
    run_dir = sys.argv[1] if len(sys.argv) > 1 else "grid/data"
    base = os.path.dirname(os.path.abspath(__file__))
    run_dir = run_dir if os.path.isabs(run_dir) else os.path.join(base, run_dir)

    combined = sorted(glob.glob(os.path.join(run_dir, "*-combined.log")))
    bad = report(combined, "combined runs (burn-in already removed)", 0.0)

    reps = sorted(glob.glob(os.path.join(run_dir, "*-rep?.log")))
    bad_reps = report(reps, "individual replicates (10%% burn-in)", 0.10)

    print("\n%d of %d combined runs below ESS %d%s"
          % (len(bad), len(combined), THRESHOLD,
             ("" if not bad else ":\n  " + "\n  ".join(bad))))
    print("%d of %d replicates below ESS %d%s"
          % (len(bad_reps), len(reps), THRESHOLD,
             ("" if not bad_reps else ":\n  " + "\n  ".join(bad_reps))))


if __name__ == "__main__":
    main()
