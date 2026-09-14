#!/usr/bin/env python3
"""SI figure: the two nodes the 18 analyses disagree on, with the posterior probability
of every resolution in every analysis.

Top: the three resolutions of each conflict as small cladograms. Bottom: one row per
analysis (grouped by calibration scheme), one column per resolution, cell = posterior
probability, three replicates pooled.

Needs grid/topology_support.json from compute_grid_stats.py.

Usage:
    python3 plot_grid_topology.py
"""

import numpy as np

import grid_common as gc

CLADOGRAMS = {  # nested tuples of tip labels for each resolution code
    "T1": ((("Tupaia", "Gal.+Prim."), "Mus")),
    "T2": ((("Mus", "Tupaia"), "Gal.+Prim.")),
    "T3": ((("Mus", "Gal.+Prim."), "Tupaia")),
    "P1": ((("Aotus", "Callithrix"), ("Cebus", "Saimiri"))),
    "P2": ((("Callithrix", ("Cebus", "Saimiri")), "Aotus")),
    "P3": ((("Aotus", ("Cebus", "Saimiri")), "Callithrix")),
}
TITLES = {
    "euarchontoglires": "Conflict 1: Mus and Tupaia within Euarchontoglires",
    "platyrrhini": "Conflict 2: Aotus and Callithrix within Callitrichidae + Cebidae",
}
# Reference trees that have each resolution (see grid/rf_topologies.txt). The two
# published topologies are de Vries & Beck (2023) and the ASTRAL tree of Vanderpool et
# al. (2020); de Vries & Beck (2023) leave Conflict 2 as a polytomy.
REF_SUPPORT = {
    "T1": ["IQ-TREE: all 3 alignments", "de Vries & Beck (2023)", "Vanderpool et al. (2020)"],
    "T2": ["none"],
    "T3": ["none"],
    "P1": ["IQ-TREE: codon,", "    unpartition"],
    "P2": ["IQ-TREE:", "    codon-noGapN"],
    "P3": ["Vanderpool et al. (2020)"],
}
def depth(node):
    return 0 if isinstance(node, str) else 1 + max(depth(ch) for ch in node)


def draw(ax, node, x, x_tip, y_next, lw=0.9):
    """Draw a tip-aligned cladogram recursively; returns the y of this node."""
    if isinstance(node, str):
        y = y_next[0]
        y_next[0] -= 1.0
        ax.plot([x, x_tip], [y, y], "k", lw=lw)
        ax.text(x_tip + 0.2, y, node, va="center", fontsize=6.5)
        return y
    x_child = x_tip - 2.2 * (depth(node) - 1)
    ys = [draw(ax, ch, x_child, x_tip, y_next, lw) for ch in node]
    for cy in ys:
        ax.plot([x, x_child], [cy, cy], "k", lw=lw)
    ax.plot([x, x], [min(ys), max(ys)], "k", lw=lw)
    return (min(ys) + max(ys)) / 2


def cladogram(ax, code):
    ax.set_xlim(0, 13)
    ax.set_ylim(-3.2, 4.6)
    ax.axis("off")
    tree = CLADOGRAMS[code]
    draw(ax, tree, 1, 1 + 2.2 * depth(tree), [4.0])
    ax.set_title(code, fontsize=6.8, pad=0)
    ax.text(0, -0.2, "Reference trees:", fontsize=5.3, color="#444444", va="top", fontweight="bold")
    for i, line in enumerate(REF_SUPPORT[code]):
        ax.text(0, -1.05 - 0.85 * i, line, fontsize=5.3, color="#444444", va="top")


def main():
    plt = gc.figure_style()
    from matplotlib.colors import LinearSegmentedColormap

    data = gc.load_json(gc.TOPO_JSON)
    conflicts, support = data["conflicts"], data["support"]
    order = ["euarchontoglires", "platyrrhini"]

    rows = []  # (label, stem)
    for sch, _, _ in gc.SCHEMES:
        for aln, _ in gc.ALIGNMENTS:
            for cond, _, _ in gc.TREE_PRIORS:
                rows.append(("%s, %s" % (aln.replace("_", "-"), "calibrated" if cond == "" else "regular"),
                             gc.stem(aln, sch, cond)))
    per_scheme = len(gc.ALIGNMENTS) * len(gc.TREE_PRIORS)

    fig = plt.figure(figsize=(7.0, 5.9))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.55, 5.5], hspace=0.04, wspace=0.45)
    top = fig.add_subfigure(gs[0, :])
    tax = top.subplots(1, 6, gridspec_kw={"wspace": 0.15})
    k = 0
    for ci, cname in enumerate(order):
        for code, _, _ in conflicts[cname]:
            cladogram(tax[k], code)
            k += 1
        top.text(0.26 + 0.5 * ci, 1.02, TITLES[cname], ha="center", fontsize=7.5, fontweight="bold")

    cmap = LinearSegmentedColormap.from_list("b", ["#FFFFFF", "#0072B2"])
    for ci, cname in enumerate(order):
        codes = [c[0] for c in conflicts[cname]]
        keys = [c[1] for c in conflicts[cname]]
        M = np.array([[support[st][kk] for kk in keys] for _, st in rows])
        ax = fig.add_subplot(gs[1, ci])
        ax.imshow(M, cmap=cmap, vmin=0, vmax=1, aspect="auto")
        for r in range(M.shape[0]):
            for c in range(M.shape[1]):
                ax.text(c, r, "%.2f" % M[r, c], ha="center", va="center", fontsize=6.5,
                        color="white" if M[r, c] > 0.6 else "black")
        ax.set_xticks(range(len(codes)))
        ax.set_xticklabels(codes)
        ax.xaxis.tick_top()
        ax.set_yticks(range(len(rows)))
        ax.set_yticklabels([lab for lab, _ in rows] if ci == 0 else [])
        ax.tick_params(length=0)
        for b in range(per_scheme, len(rows), per_scheme):
            ax.axhline(b - 0.5, color="black", lw=0.8)
        for s in ax.spines.values():
            s.set_visible(False)
        if ci == 0:
            for si, (sch, _, col) in enumerate(gc.SCHEMES):
                ax.text(-1.02, si * per_scheme + (per_scheme - 1) / 2, gc.SCHEME_SHORT[sch], rotation=90,
                        ha="center", va="center", fontsize=7, color=col, fontweight="bold",
                        transform=ax.get_yaxis_transform())
    gc.save(fig, "grid_topology_conflicts")


if __name__ == "__main__":
    main()
