#!/usr/bin/env python3
"""SI figure: posterior mean and 95% HPD of the eight calibrated nodes under the three
calibration schemes, against the fossil bounds; one panel per alignment.

Reads grid/mrca_age_comparison_grid.csv (from compare_mrca_ages.py grid/data) and the
bounds from the LPhy header of the grid XML.

Usage:
    python3 plot_grid_bounds.py [--cond=true|false]     # which tree prior to show
"""

import csv
import sys

import grid_common as gc


def main(cond):
    plt = gc.figure_style()
    import matplotlib.ticker as mt

    bounds = gc.calibration_bounds()
    rows = {}
    with open(gc.MRCA_CSV) as fh:
        for r in csv.DictReader(fh):
            rows[(r["clade"], r["run"])] = r

    n = len(gc.CALIBRATED)
    offsets = (0.24, 0.0, -0.24)
    fig, axes = plt.subplots(1, 3, figsize=(7.0, 3.3), sharey=True, gridspec_kw={"wspace": 0.08})
    for ax, (aln, atitle) in zip(axes, gc.ALIGNMENTS):
        for i, (clade, label) in enumerate(gc.CALIBRATED):
            y = n - 1 - i
            lo, hi = bounds[clade]
            ax.add_patch(plt.Rectangle((lo, y - 0.4), hi - lo, 0.8, color="#E6E6E6", lw=0, zorder=0))
            for (sch, _, col), off in zip(gc.SCHEMES, offsets):
                r = rows[(clade, gc.stem(aln, sch, cond) + "-combined")]
                ax.plot([float(r["hpd95_lower"]), float(r["hpd95_upper"])], [y + off, y + off],
                        color=col, lw=1.3, solid_capstyle="butt", zorder=2)
                ax.plot(float(r["mean"]), y + off, "o", ms=3, color=col, mec="white", mew=0.4, zorder=3)
        ax.set_xscale("log")
        ax.set_xlim(3, 260)
        ticks = [5, 10, 20, 50, 100, 200]
        ax.xaxis.set_major_locator(mt.FixedLocator(ticks))
        ax.xaxis.set_major_formatter(mt.FixedFormatter([str(t) for t in ticks]))
        ax.xaxis.set_minor_formatter(mt.NullFormatter())
        ax.tick_params(which="minor", length=0)
        ax.set_ylim(-0.6, n - 0.4)
        ax.set_title(atitle.replace(", ", ",\n"), fontsize=7.5, pad=4)
        ax.tick_params(axis="y", length=0)
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
        ax.grid(axis="x", color="#EEEEEE", lw=0.5, zorder=0)
    axes[0].set_yticks(range(n))
    axes[0].set_yticklabels([lab for _, lab in gc.CALIBRATED][::-1])
    axes[1].set_xlabel("Node age (Ma), posterior mean and 95% HPD")
    handles = [plt.Line2D([], [], color=c, marker="o", ms=3, lw=1.3, mec="white") for _, _, c in gc.SCHEMES]
    handles.append(plt.Rectangle((0, 0), 1, 1, color="#E6E6E6"))
    fig.legend(handles, [s[1] for s in gc.SCHEMES] + ["Fossil bounds (min–max)"], loc="lower center",
               ncol=4, frameon=False, fontsize=6.8, bbox_to_anchor=(0.5, -0.09), columnspacing=1.4)
    gc.save(fig, "grid_bounds" + ("_condFalse" if cond else ""))


if __name__ == "__main__":
    cond = ""
    for a in sys.argv[1:]:
        if a.startswith("--cond=") and a.split("=", 1)[1].lower() == "false":
            cond = "_condFalse"
    main(cond)
