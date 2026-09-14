#!/usr/bin/env python3
"""The three grid alignments' per-clade age comparisons merged into one PDF.

Same content as the three per-alignment scatter figures: x = median node age
under the regular tree prior, y = under the calibrated tree prior, error bars
= 95% HPD. Rows are alignments, columns are calibration priors, and every panel
shares one set of axis limits so the rows can be compared directly.

Usage:
    python3 plot_clade_merged.py
"""
import os

import matplotlib.pyplot as plt

import plot_clade_condCal as pcc
import primates_common as pc

CONDITIONS = [
    ("grid-codon", "Codon-partitioned, full alignment"),
    ("grid-codon-nogapN", "Codon-partitioned, gap/N codons removed"),
    ("grid-unpartition", "Unpartitioned, gap/N columns removed"),
]


def main():
    nrow, ncol = len(CONDITIONS), len(pc.MODELS)
    fig, axes = plt.subplots(nrow, ncol, squeeze=False,
                             figsize=(4.8 * ncol, 5.0 * nrow))
    lo, hi = [], []
    for r, (dataset, subtitle) in enumerate(CONDITIONS):
        pc.use(dataset)
        names = pc.clade_names()
        for c, model in enumerate(pc.MODELS):
            ax = axes[r][c]
            pcc.panel(ax, model, names)
            if (r, c) != (0, 0) and ax.get_legend():
                ax.get_legend().remove()
            ax.set_title(pc.MODEL_LABELS[model] if r == 0 else "", fontsize=10)
            ax.set_xlabel("median age, regular tree prior (Ma)"
                          if r == nrow - 1 else "")
            if ax.get_xscale() == "log":
                lo.append(ax.get_xlim()[0])
                hi.append(ax.get_xlim()[1])
        axes[r][0].set_ylabel("median age, calibrated tree prior (Ma)", fontsize=9.5)

    for row in axes:
        for ax in row:
            if ax.get_xscale() == "log":
                ax.set_xlim(min(lo), max(hi))
                ax.set_ylim(min(lo), max(hi))

    fig.tight_layout(rect=(0, 0, 1, 0.955), h_pad=4.0)
    for r, (_, subtitle) in enumerate(CONDITIONS):
        # row 0 also carries the column titles, so its subtitle needs more clearance
        top = axes[r][0].get_position().y1 + (0.030 if r == 0 else 0.010)
        fig.text(0.5, min(top, 0.995), subtitle, ha="center", va="bottom",
                 fontsize=13, fontweight="bold")
    out = os.path.join(pc.BASE, "grid_clade_ages.pdf")
    fig.savefig(out)
    print("wrote", out)


if __name__ == "__main__":
    main()
