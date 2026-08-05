"""
Compare LPhy-simulated trees vs BEAST sample-from-prior trees.

Input: .treestreestats.log files produced by compute_tree_stats.sh
       (TreeStat2 tab-delimited output, one row per tree).

Statistics:
    tree_length  – sum of branch lengths          (continuous)
    root_age     – root height                    (continuous)
    gamma        – Pybus & Harvey gamma statistic (continuous)
    colless      – normalised Colless imbalance   (continuous)
    b1           – B1 imbalance index             (continuous)
    cherries     – number of cherries             (discrete → integer histogram)

Usage (from this directory):
    python compare_distributions.py

Requires: numpy, scipy, matplotlib
"""

import argparse
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from analytical_cherry_prob import cherryProb, cherryProb_root

# ── Column mapping (TreeStat2 header → internal name) ────────────────────────

TREESTAT_COLS = {
    "Tree Length":            "tree_length",
    "Tree Height":            "root_age",
    "Gamma":                  "gamma",
    "Colless tree-imbalance": "colless",
    "B1":                     "b1",
    "Cherry count":           "cherries",
}

STAT_NAMES  = ["tree_length", "root_age", "gamma", "colless", "b1", "cherries"]
CONTINUOUS  = ["tree_length", "root_age", "gamma", "colless", "b1"]
DISCRETE_S  = ["cherries"]

STAT_LABELS = {
    "tree_length": "Tree Length",
    "root_age":    "Tree Height",
    "gamma":       "Gamma Statistic",
    "colless":     "Colless Index",
    "b1":          "B1 Index",
    "cherries":    "Number of Cherries",
}
DISCRETE = {"cherries"}

COLORS = {"lphy": "#56B4E9", "mcmc": "#E69F00", "mcmc_nc": "#009E73"}


# ── Reader ────────────────────────────────────────────────────────────────────

def read_log(path, burnin=0.0, thin=1):
    """Read a TreeStat2 .treestreestats.log file; return dict of stat→np.array."""
    rows = []
    header = None
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if header is None:
                header = parts
                continue
            rows.append(parts)

    col_idx = {}
    for col_i, col_name in enumerate(header):
        if col_name in TREESTAT_COLS:
            col_idx[col_i] = TREESTAT_COLS[col_name]

    data = {s: [] for s in STAT_NAMES}
    for row in rows:
        for ci, sname in col_idx.items():
            try:
                data[sname].append(float(row[ci]))
            except (IndexError, ValueError):
                data[sname].append(float("nan"))

    n_burn = int(len(rows) * burnin)
    result = {}
    for s in STAT_NAMES:
        arr = np.array(data[s], dtype=float)[n_burn:][::thin]
        result[s] = arr
    return result


# ── Comparison printing ───────────────────────────────────────────────────────

def print_comparison(stats_a, stats_b, label_a, label_b, stat_subset=None):
    names = stat_subset if stat_subset is not None else STAT_NAMES
    w = 20
    print(f"\n{'Statistic':<{w}} {label_a+' mean':>14} {label_b+' mean':>14} {'t p':>10}")
    print("-" * 60)
    for s in names:
        a = stats_a[s][np.isfinite(stats_a[s])]
        b = stats_b[s][np.isfinite(stats_b[s])]
        if len(a) < 2 or len(b) < 2:
            print(f"{STAT_LABELS[s]:<{w}} {'n/a':>14} {'n/a':>14} {'n/a':>10}")
            continue
        t_p = stats.ttest_ind(a, b, equal_var=False).pvalue
        print(f"{STAT_LABELS[s]:<{w}} {a.mean():>14.4f} {b.mean():>14.4f} {t_p:>10.4f}")


# ── Axis plotting ─────────────────────────────────────────────────────────────

def plot_ax(ax, series, s, show_legend=False):
    """series: list of (values_array, label, color_key) tuples."""
    clean = [(v[np.isfinite(v)], lbl, col) for v, lbl, col in series]
    if any(len(v) < 2 for v, _, _ in clean):
        ax.set_title(STAT_LABELS[s], fontsize=8)
        return

    if s in DISCRETE:
        lo = int(min(v.min() for v, _, _ in clean))
        hi = int(max(v.max() for v, _, _ in clean))
        bins = np.arange(lo - 0.5, hi + 1.5, 1.0)
        for vals, lbl, col in clean:
            ax.hist(vals, bins=bins, density=True, alpha=0.5,
                    color=COLORS[col], label=lbl, edgecolor="none")
        ax.set_ylabel("Proportion", fontsize=7)
    else:
        for vals, lbl, col in clean:
            ax.hist(vals, bins=40, density=True, alpha=0.5,
                    color=COLORS[col], label=lbl, edgecolor="none")
        ax.set_ylabel("Density", fontsize=7)

    ax.set_title(STAT_LABELS[s], fontsize=8, fontweight="bold")
    ax.set_xlabel(STAT_LABELS[s], fontsize=7)
    ax.tick_params(labelsize=6)
    if s in DISCRETE:
        ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
    else:
        ax.ticklabel_format(useOffset=False, axis="x")
        ax.xaxis.set_major_locator(plt.MaxNLocator(4, prune="both"))
        ax.xaxis.set_major_formatter(plt.FormatStrFormatter("%.2f"))
    if show_legend:
        ax.legend(fontsize=6)


# ── Figure builders ───────────────────────────────────────────────────────────

def make_figure_100leaf(stats_a, stats_b, label_a, label_b, title, stats_nc=None):
    """Two rows of three: continuous top, continuous+discrete bottom."""
    top_row = ["tree_length", "root_age", "gamma"]
    bot_row = ["colless", "b1", "cherries"]
    fig, axes = plt.subplots(2, 3, figsize=(6.5, 4.5))
    for i, s in enumerate(top_row):
        ser = [(stats_a[s], label_a, "lphy"), (stats_b[s], label_b, "mcmc")]
        if stats_nc is not None:
            ser.append((stats_nc[s], "MCMC (not conditioned)", "mcmc_nc"))
        plot_ax(axes[0, i], ser, s, show_legend=False)
    for i, s in enumerate(bot_row):
        ser = [(stats_a[s], label_a, "lphy"), (stats_b[s], label_b, "mcmc")]
        if stats_nc is not None:
            ser.append((stats_nc[s], "MCMC (not conditioned)", "mcmc_nc"))
        plot_ax(axes[1, i], ser, s, show_legend=False)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=len(handles),
               fontsize=7, frameon=False, bbox_to_anchor=(0.5, 0.0))
    fig.subplots_adjust(bottom=0.14)
    plt.tight_layout(rect=[0, 0.10, 1, 1])
    return fig


def make_figure_4leaf_combined(stem_data, root_data, stat_subset,
                               stem_analytical_cherry_p=None,
                               root_analytical_cherry_p=None):
    """2×3 figure combining origin- and root-conditioned 4-leaf scenarios.

    stem_data / root_data: dicts with keys 'lphy', 'mcmc', 'nc' (optional),
                           each mapping to a stats dict.
    stem_analytical_cherry_p: if given, draw the analytical P on the cherry
                              panel of the Origin Conditioned row.
    root_analytical_cherry_p: if given, draw the analytical P on the cherry
                              panel of the Root Conditioned row.
    """
    n = len(stat_subset)
    fig, axes = plt.subplots(2, n, figsize=(2 * n, 5))
    fig.subplots_adjust(hspace=0.55, top=0.88, bottom=0.18,
                        left=0.12, right=0.97, wspace=0.55)

    row_info = [
        ("Origin Conditioned", stem_data),
        ("Root Conditioned",   root_data),
    ]

    legend_handles = None
    analytical_handle = None
    for row, (row_title, data) in enumerate(row_info):
        for col, s in enumerate(stat_subset):
            ax = axes[row, col]
            series = [(d[s], lbl, col_key)
                      for d, lbl, col_key in
                      [(data["lphy"], "LPhy", "lphy"),
                       (data["mcmc"], "MCMC", "mcmc")]
                      + ([(data["nc"], "MCMC (not conditioned)", "mcmc_nc")]
                         if data.get("nc") is not None else [])]
            plot_ax(ax, series, s, show_legend=False)
            ax.set_title("")
            if legend_handles is None:
                legend_handles, legend_labels = ax.get_legend_handles_labels()

            analytical_p = (stem_analytical_cherry_p if row == 0
                            else root_analytical_cherry_p)
            if s == "cherries" and analytical_p is not None:
                p2 = analytical_p
                p1 = 1.0 - p2
                line, = ax.plot([0.5, 1.5], [p1, p1], color="black",
                                linestyle="--", linewidth=1.2)
                ax.plot([1.5, 2.5], [p2, p2], color="black",
                        linestyle="--", linewidth=1.2)
                if analytical_handle is None:
                    analytical_handle = line

    # Shared legend below the figure
    all_handles = legend_handles or []
    all_labels  = legend_labels  if legend_handles else []
    if analytical_handle is not None:
        all_handles = all_handles + [analytical_handle]
        all_labels  = all_labels  + ["Analytical"]
    if all_handles:
        fig.legend(all_handles, all_labels,
                   loc="lower center", ncol=len(all_handles),
                   fontsize=7, frameon=False,
                   bbox_to_anchor=(0.5, 0.0))

    # Place row titles in figure coordinates
    for row, (row_title, _) in enumerate(row_info):
        pos0 = axes[row, 0].get_position()
        pos1 = axes[row, n - 1].get_position()
        x_center = (pos0.x0 + pos1.x1) / 2
        y_top    = pos0.y1

        # Row title: centred above the row
        fig.text(x_center, y_top + 0.01,
                 row_title,
                 fontsize=9, fontweight="bold",
                 va="bottom", ha="center")

    return fig


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stem-lphy", default="fixStemLPhy.treestreestats.log")
    parser.add_argument("--stem-mcmc", default="fixStemMCMC.treestreestats.log")
    parser.add_argument("--stem-nc",   default="fixStemMCMC_not_cond.treestreestats.log")
    parser.add_argument("--root-lphy", default="fixRootLPhy.treestreestats.log")
    parser.add_argument("--root-mcmc", default="fixRootMCMC.treestreestats.log")
    parser.add_argument("--root-nc",   default="fixRootMCMC_not_cond.treestreestats.log")
    parser.add_argument("--mcmc-burnin", type=float, default=0.10)
    parser.add_argument("--mcmc-thin",   type=int,   default=5)
    args = parser.parse_args()

    LEAF4_STATS = ["tree_length", "root_age", "cherries"]

    # ── 100-leaf scenarios ────────────────────────────────────────────────────
    leaf100_scenarios = [
        ("Fixed-stem  n=100  λ=2 μ=1 ρ=0.1 stemAge=3",
         args.stem_lphy, args.stem_mcmc, args.stem_nc, "fixStem_distributions.pdf"),
        ("Fixed-root  n=100  λ=2 μ=1 ρ=0.1 rootAge=3",
         args.root_lphy, args.root_mcmc, args.root_nc, "fixRoot_distributions.pdf"),
    ]

    for title, lphy_path, mcmc_path, nc_path, out_pdf in leaf100_scenarios:
        print(f"\n{'='*60}\nScenario: {title}\n{'='*60}")
        try:
            lphy_stats = read_log(lphy_path, burnin=0.0, thin=1)
            mcmc_stats = read_log(mcmc_path, burnin=args.mcmc_burnin,
                                  thin=args.mcmc_thin)
            nc_stats = None
            try:
                nc_stats = read_log(nc_path, burnin=args.mcmc_burnin,
                                    thin=args.mcmc_thin)
            except FileNotFoundError:
                pass
            print(f"  LPhy: {len(lphy_stats['tree_length'])}  "
                  f"MCMC: {len(mcmc_stats['tree_length'])}")
            print_comparison(lphy_stats, mcmc_stats, "LPhy", "MCMC")
            if nc_stats is not None:
                print_comparison(lphy_stats, nc_stats, "LPhy",
                                 "MCMC (not conditioned)")
            fig = make_figure_100leaf(lphy_stats, mcmc_stats, "LPhy", "MCMC",
                                      title, stats_nc=nc_stats)
            fig.savefig(out_pdf, bbox_inches="tight")
            plt.close(fig)
            print(f"  Saved: {out_pdf}")
        except FileNotFoundError as e:
            print(f"  File not found: {e} — skipping")

    # ── 4-leaf combined figure ────────────────────────────────────────────────
    def load_4leaf(lphy_path, mcmc_path, nc_path):
        data = {}
        try:
            data["lphy"] = read_log(lphy_path, burnin=0.0, thin=1)
            data["mcmc"] = read_log(mcmc_path, burnin=args.mcmc_burnin,
                                    thin=args.mcmc_thin)
        except FileNotFoundError as e:
            print(f"  Missing: {e}")
            return None
        try:
            data["nc"] = read_log(nc_path, burnin=args.mcmc_burnin,
                                  thin=args.mcmc_thin)
        except FileNotFoundError:
            data["nc"] = None
        return data

    print(f"\n{'='*60}\n4-leaf combined figure\n{'='*60}")
    stem_data = load_4leaf("fix4LeafStemLPhy.treestreestats.log",
                           "fix4LeafStemMCMC.treestreestats.log",
                           "fix4LeafStemMCMC_not_cond.treestreestats.log")
    root_data = load_4leaf("fix4LeafRootLPhy.treestreestats.log",
                           "fix4LeafRootMCMC.treestreestats.log",
                           "fix4LeafRootMCMC_not_cond.treestreestats.log")

    if stem_data and root_data:
        for label, data in [("Origin cond", stem_data), ("Root cond", root_data)]:
            print(f"\n  {label}:")
            print_comparison(data["lphy"], data["mcmc"], "LPhy", "MCMC",
                             stat_subset=LEAF4_STATS)
            if data["nc"]:
                print_comparison(data["lphy"], data["nc"], "LPhy",
                                 "MCMC (not conditioned)", stat_subset=LEAF4_STATS)
        analytical_stem_p = cherryProb(
            birth_rate=2.0, death_rate=1.0, rho=0.5, time=3.0, cherryTime=1.0
        )
        analytical_root_p = cherryProb_root(
            birth_rate=2.0, death_rate=1.0, rho=0.5, time=3.0, cherryTime=1.0
        )
        print(f"\n  Analytical P (stem conditioned) = {analytical_stem_p:.4f}")
        print(f"  Analytical P (root conditioned) = {analytical_root_p:.4f}")
        fig = make_figure_4leaf_combined(stem_data, root_data, LEAF4_STATS,
                                         stem_analytical_cherry_p=analytical_stem_p,
                                         root_analytical_cherry_p=analytical_root_p)
        fig.savefig("fix4Leaf_distributions.pdf", bbox_inches="tight")
        plt.close(fig)
        print("  Saved: fix4Leaf_distributions.pdf")


if __name__ == "__main__":
    main()
