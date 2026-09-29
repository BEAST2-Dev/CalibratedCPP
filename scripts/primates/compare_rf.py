#!/usr/bin/env python3
"""Pairwise Robinson-Foulds distances among one dataset's six summary trees.

Prints the RF matrix, then for every pair that differs the clades found in only one
of the two trees. Writes the matrix to figures/<dataset>_rf.csv.

Usage:
    python3 compare_rf.py [--dataset=grid-codon-deepcal] [--unrooted]
"""

import csv
import itertools
import sys

import primates_common as pc

SHORT = {"calibrationPrior": "joint", "suggestedPrior": "deVries", "uniformPrior": "uniform"}


def splits(clades, taxa):
    """Each clade as the smaller side of its split, so root placement is ignored."""
    out = set()
    for c in clades:
        other = frozenset(taxa - c)
        if len(c) > 1 and len(other) > 1:
            out.add(min(c, other, key=lambda s: (len(s), sorted(s))))
    return out


def main():
    unrooted = "--unrooted" in sys.argv[1:]
    for a in sys.argv[1:]:
        if a.startswith("--dataset="):
            pc.use(a.split("=", 1)[1])
    trees = {}
    for model in pc.MODELS:
        for cond in ("true", "false"):
            name = "%s-%s" % (SHORT[model], "cond" if cond == "true" else "regular")
            clades = {t for t, *_ in pc.parse_summary_tree(pc.summary_tree_path(model, cond))}
            trees[name] = {c for c in clades if len(c) > 1}

    taxa = max((c for cl in trees.values() for c in cl), key=len)
    if unrooted:
        trees = {n: splits(c, taxa) for n, c in trees.items()}
    names = list(trees)
    dist = {(a, b): len(trees[a] ^ trees[b]) for a in names for b in names}

    width = max(len(n) for n in names)
    print("Robinson-Foulds distances (%s, %d taxa, max %d)\n" % (
        "unrooted splits" if unrooted else "rooted clades", len(taxa), 2 * (len(taxa) - 2)))
    print(" " * width + "".join("%4d" % i for i in range(1, len(names) + 1)))
    for i, a in enumerate(names, 1):
        print("%*s" % (width, a) + "".join("%4d" % dist[a, b] for b in names) + "   (%d)" % i)

    for a, b in itertools.combinations(names, 2):
        if not dist[a, b]:
            continue
        print("\n%s vs %s: RF = %d" % (a, b, dist[a, b]))
        for label, clades in ((a + " only", trees[a] - trees[b]), (b + " only", trees[b] - trees[a])):
            for c in sorted(clades, key=len):
                print("  %-*s %s" % (width + 5, label, ", ".join(sorted(c))))

    out = pc.fig_path("%s_rf.csv" % pc.DATASET)
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([""] + names)
        for a in names:
            w.writerow([a] + [dist[a, b] for b in names])
    print("\nwrote", out)


if __name__ == "__main__":
    main()
