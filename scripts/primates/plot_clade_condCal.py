#!/usr/bin/env python3
"""Per-clade node ages under conditionOnCalibrations=false vs true.

One point per clade shared by the two summary trees: x = posterior median age
with conditioning off, y = with conditioning on, error bars = 95% HPD. Calibrated
clades are drawn separately from the non-calibrated ones, which are the clades the
question is really about.

Summary trees come from "treeannotator -height CA -topology CCD0".

Usage:
    python3 plot_clade_condCal.py
"""

import os
import sys

import matplotlib.pyplot as plt

import primates_common as pc



def panel(ax, model, names):
    t_false = pc.summary_tree_path(model, "false")
    t_true = pc.summary_tree_path(model, "true")
    if not (os.path.exists(t_false) and os.path.exists(t_true)):
        ax.text(0.5, 0.5, "%s\nno condFalse run yet" % pc.MODEL_LABELS[model],
                ha="center", va="center", transform=ax.transAxes, color="grey")
        ax.set_xticks([])
        ax.set_yticks([])
        return

    false = {c[0]: c[1:] for c in pc.parse_summary_tree(t_false)}
    true = {c[0]: c[1:] for c in pc.parse_summary_tree(t_true)}
    calib = set(pc.calibrated_clades(model))
    shared = sorted(set(false) & set(true), key=lambda c: -false[c][0])

    for is_cal, colour, marker, label in [
            (False, "#3b6ea5", "o", "other clades"),
            (True, "#c0392b", "s", "calibrated")]:
        sel = [c for c in shared if (c in calib) == is_cal]
        if not sel:
            continue
        x = [false[c][0] for c in sel]
        y = [true[c][0] for c in sel]
        xerr = [[x[i] - false[c][1] for i, c in enumerate(sel)],
                [false[c][2] - x[i] for i, c in enumerate(sel)]]
        yerr = [[y[i] - true[c][1] for i, c in enumerate(sel)],
                [true[c][2] - y[i] for i, c in enumerate(sel)]]
        ax.errorbar(x, y, xerr=xerr, yerr=yerr, fmt=marker, ms=4, color=colour,
                    ecolor=colour, elinewidth=0.8, alpha=0.75, capsize=0, label=label)

    lo = min(min(false[c][1] for c in shared), min(true[c][1] for c in shared))
    hi = max(max(false[c][2] for c in shared), max(true[c][2] for c in shared))
    lo, hi = lo * 0.8, hi * 1.25
    ax.plot([lo, hi], [lo, hi], "k--", lw=0.8, zorder=0)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal")
    ax.set_title("%s  (%d clades)" % (pc.MODEL_LABELS[model], len(shared)), fontsize=10)
    ax.set_xlabel("median age, conditioning off (Ma)")
    ax.legend(fontsize=8, loc="upper left", frameon=False)

    for c in shared:
        if c in calib:
            ax.annotate(names.get(c, ""), (false[c][0], true[c][0]), fontsize=6.5,
                        xytext=(5, -3), textcoords="offset points", color="#c0392b")


def main():
    if len(sys.argv) > 1 and sys.argv[1].startswith("--dataset="):
        pc.use(sys.argv.pop(1).split("=", 1)[1])
    out = os.path.join(pc.BASE, "%s_cladeCondCal.png" % pc.DATASET)
    names = pc.clade_names()
    fig, axes = plt.subplots(1, len(pc.MODELS), figsize=(5 * len(pc.MODELS), 5.8))
    for ax, model in zip(axes, pc.MODELS):
        panel(ax, model, names)
    axes[0].set_ylabel("median age, conditioning on (Ma)")
    fig.suptitle("Effect of conditioning the tree prior on calibrations, per clade", fontsize=12)
    fig.tight_layout(rect=[0, 0.02, 1, 0.95])
    fig.savefig(out, dpi=200)
    print("wrote", out)


if __name__ == "__main__":
    main()
