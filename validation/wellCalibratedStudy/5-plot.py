#!/usr/bin/env python3
import os, re, glob, concurrent.futures
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, binom

prefix="fixStemInference"
input_csv   = f"./xmls_summary/{prefix}_diversification_turnover_results.csv"
logs_glob   = f"./xmls/{prefix}-*.log"
burnin_frac = 0.1
out_png     = f"./xmls_summary/{prefix}_diversification_turnover_calibration_from_logs.png"
out_pdf     = f"./xmls_summary/{prefix}_diversification_turnover_calibration_from_logs.pdf"

log_col_diversification = "diversification"
log_col_turnover     = "turnover"
csv_col_r           = "r"
csv_col_true_diversification = "true_diversification"
csv_col_true_turnover     = "true_turnover"

# ============================================================
# Helper functions
# ============================================================
def infer_r_from_filename(path):
    fname = os.path.basename(path)
    pats = [
        re.compile(r"[_-]r(\d+)(?!\d)"),
        re.compile(r"[_-]rep(?:licate)?(\d+)(?!\d)"),
        re.compile(r"(\d+)\.log$")
    ]
    for p in pats:
        m = p.search(fname)
        if m:
            return int(m.group(1))
    return None

def read_one_log(path):
    try:
        df = pd.read_csv(path, comment="#", sep=r"\s+",
                         usecols=[log_col_diversification, log_col_turnover])
    except Exception:
        return None
    r_id = infer_r_from_filename(path)
    if r_id is None:
        return None
    n = len(df)
    if n == 0:
        return None
    start = int(np.floor(burnin_frac * n))
    lam = df[log_col_diversification].iloc[start:].to_numpy()
    turnover  = df[log_col_turnover].iloc[start:].to_numpy()
    return (r_id, lam, turnover)

# ============================================================
# Parallel reading
# ============================================================
paths = glob.glob(logs_glob)
post_diversification, post_turnover = {}, {}

print(f"Reading {len(paths)} log files in parallel...")
with concurrent.futures.ThreadPoolExecutor(max_workers=os.cpu_count()) as ex:
    for result in ex.map(read_one_log, paths):
        if result is None:
            continue
        r_id, lam, turnover = result
        if lam.size == 0 or turnover.size == 0:
            continue
        post_diversification[r_id] = lam
        post_turnover[r_id] = turnover
print(f"Loaded {len(post_diversification)} valid logs.")

# ============================================================
# Load true values
# ============================================================
truth_df = pd.read_csv(input_csv)

# keep r values that match logs, including 0
valid_r = sorted(set(truth_df[csv_col_r]).intersection(post_diversification.keys(), post_turnover.keys()))

truth_df = truth_df[truth_df[csv_col_r].isin(valid_r)].copy()
truth_df.sort_values(csv_col_r, inplace=True)

R = truth_df[csv_col_r].to_numpy()

true_diversification = truth_df[csv_col_true_diversification].to_numpy()
true_turnover     = truth_df[csv_col_true_turnover].to_numpy()

# ============================================================
# Compute ranks, medians, coverage
# ============================================================
def _adjusted_gamma(counts_sim, N, p, prob):
    lower = binom.cdf(counts_sim, N, p[None, :])
    upper = binom.sf(counts_sim - 1, N, p[None, :])
    ptail = np.minimum(lower, upper)
    return np.quantile(ptail.min(axis=1), 1.0 - prob)

def ecdf_null_band(N, prob=0.95, n_grid=200, M=2000, seed=0):
    rng  = np.random.default_rng(seed)
    u    = np.linspace(1.0/(n_grid+1), n_grid/(n_grid+1), n_grid)
    sims = rng.random((M, N))
    counts = (sims[:, None, :] <= u[None, :, None]).sum(axis=2)
    g  = _adjusted_gamma(counts, N, u, prob)
    return u, binom.ppf(g, N, u) / N, binom.ppf(1.0 - g, N, u) / N

def coverage_null_band(N, alphas, prob=0.95, n_grid=200, M=2000, seed=0):
    u, loE, hiE = ecdf_null_band(N, prob=prob, n_grid=n_grid, M=M, seed=seed)
    a = np.asarray(alphas, float)
    return np.interp(a, u, loE), np.interp(a, u, hiE)

def compute_from_samples(true_vals, posts):
    n = len(true_vals)
    ranks   = np.empty(n)
    medians = np.empty(n)
    q025    = np.empty(n)
    q975    = np.empty(n)
    coverage_flags = np.empty(n, dtype=bool)
    alphas  = np.linspace(0, 1, 21)
    cover   = np.empty_like(alphas)

    for i, (r, t) in enumerate(zip(R, true_vals)):
        s = posts[int(r)]
        ranks[i]   = np.mean(s <= t)
        medians[i] = np.median(s)
        q025[i], q975[i] = np.quantile(s, [0.025, 0.975])
        coverage_flags[i] = (t >= q025[i]) and (t <= q975[i])

    all_samples = [posts[int(r)] for r in R]
    for j, a in enumerate(alphas):
        loq, hiq = (1-a)/2, 1-(1-a)/2
        covered = [(np.quantile(s, loq) <= t <= np.quantile(s, hiq))
                   for s, t in zip(all_samples, true_vals)]
        c = np.mean(covered)
        cover[j] = c

    return {
        "ranks": ranks,
        "medians": medians,
        "q025": q025,
        "q975": q975,
        "coverage_flags": coverage_flags,
        "alphas": alphas,
        "coverage": cover
    }

print("Computing calibration statistics...")
diversification_stats = compute_from_samples(true_diversification, post_diversification)
turnover_stats = compute_from_samples(true_turnover, post_turnover)
print("Calibration stats computed.")

# ============================================================
# Plot 3×2 figure
# ============================================================
fig, axes = plt.subplots(3, 2, figsize=(10, 10))
pairs = [
    ("diversification (λ - μ)", true_diversification, diversification_stats),
    ("turnover (µ/λ))", true_turnover,     turnover_stats),
]

for j, (label, truths, st) in enumerate(pairs):
    # Coverage
    ax = axes[0, j]
    Ncov = len(st["ranks"])
    lo, hi = coverage_null_band(Ncov, st["alphas"])
    ax.fill_between(st["alphas"], lo, hi, color="0.8", alpha=0.6,
                    label="95% calibrated region")
    ax.plot([0, 1], [0, 1], "r-", lw=1.5)
    ax.plot(st["alphas"], st["coverage"], "ko", ms=4)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("Credibility level α", fontsize=14)
    ax.set_ylabel("Coverage probability", fontsize=14)
    ax.set_title(label, fontsize=16)

    # ECDF
    ax = axes[1, j]
    r_sorted = np.sort(st["ranks"])
    ecdf_y = np.arange(1, len(r_sorted)+1) / len(r_sorted)
    ug, lo, hi = ecdf_null_band(len(r_sorted))
    ax.fill_between(ug, lo, hi, color="cyan", alpha=0.25,
                    label="95% calibrated region")
    ax.plot([0,1], [0,1], "r-")
    ax.step(r_sorted, ecdf_y, where="post", color="k")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("Posterior rank of true value", fontsize=14)
    ax.set_ylabel("ECDF", fontsize=14)

    # Median vs true
    ax = axes[2, j]
    inside = st["coverage_flags"]
    colors = np.where(inside, "cyan", "red")
    yerr = np.vstack([
        st["medians"] - st["q025"],
        st["q975"] - st["medians"]
    ])
    for xi, yi, ye, c in zip(truths, st["medians"], yerr.T, colors):
        ax.errorbar(xi, yi, yerr=[[ye[0]], [ye[1]]], fmt="o",
                    color="k", ecolor=c, elinewidth=2, alpha=0.9, markersize=5)

    lo = min(truths.min(), st["medians"].min())
    hi = max(truths.max(), st["medians"].max())
    ax.plot([lo, hi], [lo, hi], "r-")

    # Pearson R
    r_val, _ = pearsonr(truths, st["medians"])
    ax.text(0.05, 0.9, f"r={r_val:.3f}", transform=ax.transAxes)

    # Coverage percentage annotation
    pct = 100 * np.mean(st["coverage_flags"])
    ax.text(0.05, 0.82, f"Coverage = {pct:.1f}%", transform=ax.transAxes)

    ax.set_xlabel(f"True value ({label})", fontsize=14)
    ax.set_ylabel("Posterior median", fontsize=14)

plt.tight_layout()
os.makedirs(os.path.dirname(out_png), exist_ok=True)
plt.savefig(out_png, dpi=300, bbox_inches="tight")
plt.savefig(out_pdf, bbox_inches="tight")
plt.close(fig)

print(f"Saved combined figure:\n  {out_png}\n  {out_pdf}")
