"""
Plot JMH scaling benchmark results from scaling_benchmark.csv.

Usage:
    python plot_scaling_benchmark.py
"""

import csv
import os
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CSV_PATH   = os.path.join(SCRIPT_DIR, "scaling_benchmark.csv")

series: dict[str, dict[str, list]] = {}  # {name: {"x": [], "y": [], "err": []}}

with open(CSV_PATH, newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        # Strip package prefix to get short name: "bd", "erlang", "gammaWarm", "gammaCold"
        full_name = row["Benchmark"].strip('"')
        name = full_name.split(".")[-1]

        x   = int(row["Param: nLeaves"].strip('"'))
        y   = float(row["Score"].strip('"'))
        err = float(row["Score Error (99.9%)"].strip('"'))

        if name not in series:
            series[name] = {"x": [], "y": [], "err": []}
        series[name]["x"].append(x)
        series[name]["y"].append(y)
        series[name]["err"].append(err)

# Sort each series by x
for s in series.values():
    order = np.argsort(s["x"])
    s["x"]   = [s["x"][i]   for i in order]
    s["y"]   = [s["y"][i]   for i in order]
    s["err"] = [s["err"][i] for i in order]

# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------
LABELS = {
    "bd":        "Constant-rate BD",
    "erlang":    "Erlang k=3",
    "erlang10":    "Erlang k=10",
    "gammaCold": "Gamma VIDE",
}
COLORS = {
    "bd":        "#2196F3",   # blue
    "erlang":    "#4CAF50",   # green
    "erlang10":  "#8F00FF",  # purple
    "gammaCold": "#FF5722",   # deep orange
}
MARKERS = {
    "bd":        "o",
    "erlang":    "s",
    "erlang10":  "^",
    "gammaCold": "D",
}

# ---------------------------------------------------------------------------
# Figure 1 – log-log overview (all four series)
# ---------------------------------------------------------------------------
fig1, ax1 = plt.subplots(figsize=(7, 5))

for name in ["bd", "erlang", "erlang10", "gammaCold"]:
    if name not in series:
        continue
    d = series[name]
    ax1.errorbar(
        d["x"], d["y"],
        yerr=d["err"],
        label=LABELS[name],
        color=COLORS[name],
        marker=MARKERS[name],
        markersize=6,
        linewidth=1.8,
        capsize=4,
        capthick=1.2,
    )

ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.set_xlabel("Number of leaves", fontsize=12)
ax1.set_ylabel("Average time (µs / evaluation)", fontsize=12)
ax1.legend(fontsize=10)
ax1.xaxis.set_major_formatter(ticker.ScalarFormatter())
ax1.yaxis.set_major_formatter(ticker.LogFormatterSciNotation())
ax1.set_xticks([10, 25, 50, 100, 250, 500, 1000, 10000])
ax1.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.6)
fig1.tight_layout()
out1 = os.path.join(SCRIPT_DIR, "scaling_benchmark_loglog.pdf")
fig1.savefig(out1, dpi=150)
print(f"Saved: {out1}")

# ---------------------------------------------------------------------------
# Figure 2 – two-panel: BD+Erlang (linear) | Gamma cold/warm (linear)
# ---------------------------------------------------------------------------
fig2, (ax_fast, ax_slow) = plt.subplots(1, 2, figsize=(12, 5))

# Left panel: fast methods
for name in ["bd", "erlang", "erlang10"]:
    if name not in series:
        continue
    d = series[name]
    ax_fast.errorbar(
        d["x"], d["y"],
        yerr=d["err"],
        label=LABELS[name],
        color=COLORS[name],
        marker=MARKERS[name],
        markersize=6,
        linewidth=1.8,
        capsize=4,
        capthick=1.2,
    )
ax_fast.set_xlabel("Number of leaves", fontsize=12)
ax_fast.set_ylabel("Average time (µs / evaluation)", fontsize=12)
ax_fast.set_title("Closed-form methods", fontsize=12)
ax_fast.legend(fontsize=10)
ax_fast.set_xticks([0, 100, 200, 300, 400, 500, 1000, 10000])
ax_fast.ticklabel_format(style="sci", axis="y", scilimits=(0, 0), useMathText=True)
ax_fast.grid(True, linestyle="--", linewidth=0.5, alpha=0.6)

# Right panel: VIDE (gamma) methods
for name in ["gammaCold"]:
    if name not in series:
        continue
    d = series[name]
    ax_slow.errorbar(
        d["x"], d["y"],
        yerr=d["err"],
        label=LABELS[name],
        color=COLORS[name],
        marker=MARKERS[name],
        markersize=6,
        linewidth=1.8,
        capsize=4,
        capthick=1.2,
    )
ax_slow.set_xlabel("Number of leaves", fontsize=12)
ax_slow.set_ylabel("Average time (µs / evaluation)", fontsize=12)
ax_slow.set_title("Numerical VIDE (Gamma)", fontsize=12)
ax_slow.legend(fontsize=10)
ax_slow.set_xticks([0, 100, 200, 300, 400, 500, 1000, 10000])
ax_slow.ticklabel_format(style="sci", axis="y", scilimits=(0, 0), useMathText=True)
ax_slow.grid(True, linestyle="--", linewidth=0.5, alpha=0.6)

fig2.tight_layout()
out2 = os.path.join(SCRIPT_DIR, "scaling_benchmark_panels.pdf")
fig2.savefig(out2, dpi=150)
print(f"Saved: {out2}")

plt.show()
