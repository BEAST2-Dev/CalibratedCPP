#!/usr/bin/env python3
"""Read every grid .trees file once and cache what plot_grid_topology.py needs.

Writes grid/topology_support.json: the posterior probability of each resolution of
the two contested nodes, per analysis. The alternative resolutions are not in the
MCC tree, so their support cannot come from the summary tree.

Pure Python: runs anywhere the .trees files are. Reads the LogCombiner output
(three replicates, 10% burn-in each already removed).

Usage:
    python3 compute_grid_stats.py
"""

import json
import os

import grid_common as gc

# The two nodes on which the 18 analyses disagree (found by comparing the summary
# trees' RF distances), each given as the three possible resolutions of a three-way split.
GP = None  # Galeopterus + Primates, filled once the taxon list is known
CONFLICTS = {
    "euarchontoglires": [
        ("T1", "Tupaia+GaleoPrimates", "(Mus, (Tupaia, Galeopterus+Primates))"),
        ("T2", "Mus+Tupaia", "((Mus, Tupaia), Galeopterus+Primates)"),
        ("T3", "Mus+GaleoPrimates", "(Tupaia, (Mus, Galeopterus+Primates))"),
    ],
    "platyrrhini": [
        ("P1", "Aotus+Callithrix", "((Aotus, Callithrix), Cebus+Saimiri)"),
        ("P2", "Callithrix+CebusSaimiri", "(Aotus, (Callithrix, Cebus+Saimiri))"),
        ("P3", "Aotus+CebusSaimiri", "(Callithrix, (Aotus, Cebus+Saimiri))"),
    ],
}


def clade_sets(all_taxa):
    prim = set(all_taxa) - gc.OUTGROUP
    gp = prim | {"Galeopterus_variegatus"}
    return {
        "Tupaia+GaleoPrimates": gp | {"Tupaia_chinensis"},
        "Mus+Tupaia": {"Mus_musculus", "Tupaia_chinensis"},
        "Mus+GaleoPrimates": gp | {"Mus_musculus"},
        "Aotus+Callithrix": {"Aotus_nancymaae", "Callithrix_jacchus"},
        "Callithrix+CebusSaimiri": {"Callithrix_jacchus", "Cebus_capucinus_imitator", "Saimiri_boliviensis"},
        "Aotus+CebusSaimiri": {"Aotus_nancymaae", "Cebus_capucinus_imitator", "Saimiri_boliviensis"},
    }


def scan(path):
    """One pass over a .trees file: how many trees contain each target clade."""
    translate, trees = gc.read_trees(path, burnin=0.0)
    all_taxa = set(translate.values())
    targets = {k: frozenset(v) for k, v in clade_sets(all_taxa).items()}
    counts = {k: 0 for k in targets}
    for nwk in trees:
        seen = set()
        gc.walk_newick(nwk, translate, lambda taxa, height: seen.add(frozenset(taxa)))
        for k, t in targets.items():
            if t in seen:
                counts[k] += 1
    return counts, len(trees)


def main():
    topo = {}
    for st in gc.all_stems():
        path = os.path.join(gc.GRID_DATA, "%s-combined.trees" % st)
        if not os.path.exists(path):
            print("missing", path)
            continue
        counts, n = scan(path)
        topo[st] = {k: round(v / n, 4) for k, v in counts.items()}
        topo[st]["n"] = n
        print(st, "posterior trees:", n)

    with open(gc.TOPO_JSON, "w") as fh:
        json.dump({"conflicts": CONFLICTS, "support": topo}, fh, indent=1)
    print("wrote", gc.TOPO_JSON)


if __name__ == "__main__":
    main()
