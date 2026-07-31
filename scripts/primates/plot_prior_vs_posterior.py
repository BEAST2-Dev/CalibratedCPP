#!/usr/bin/env python3
"""
Overlay each calibrated node's marginal PRIOR distribution (from a sample-from-prior run,
where the alignment has been reduced to a single missing/'?' site so the likelihood
contributes nothing and MCMC explores the tree prior alone) against its POSTERIOR
distribution (from the matching real-data run), per node.

This is the diagnostic for whether a discrepancy between a posterior node age and the
dos Reis ground truth is coming from the DATA (posterior pulls away from an accurate
prior) or is already baked into the MODEL (the tree-prior process reshapes the specified
per-node calibration density -- e.g. Uniform(a,b) -- into a different marginal, because the
birth-death/CPP process and the OTHER calibrations interact). See Heled & Drummond (2012)
for why a tree prior's marginal on a calibrated node is not generally equal to the density
you typed into that node's calibration.

Both runs must come from the same XML family (same taxa, same calibrations, same
TaxonSetN identifiers) so that mrca.age(TaxonSetN) columns line up directly -- no tree
parsing or node-identity matching is needed. Self-contained (no imports from other
scripts in this folder).

Usage:
    python plot_prior_vs_posterior.py

Looks for (relative to this script's directory), trying data/ then the script's own
directory, and .txt then .log as the extension:
    sample-from-prior_allCalibrations(.txt|.log)                 vs  allCalibrations(.txt|.log)
    sample-from-prior_allCalibrations_suggestedPriors(.txt|.log) vs  allCalibrations_suggestedPriors(.txt|.log)

Any pair whose files aren't found yet is skipped with a message, so this can be re-run as
runs complete without editing anything.
"""
import os

import matplotlib.pyplot as plt
import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
BURNIN_FRACTION = 0.1


def read_log_columns(log_path, columns, burnin_fraction):
    with open(log_path) as f:
        header = f.readline().strip().split("\t")
    col_idx = {c: header.index(c) for c in columns if c in header}
    missing = [c for c in columns if c not in col_idx]
    if missing:
        raise ValueError(f"Columns not found in log file: {missing}")

    data = {c: [] for c in columns}
    with open(log_path) as f:
        f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(header):
                continue
            for c in columns:
                data[c].append(float(parts[col_idx[c]]))

    n = len(next(iter(data.values())))
    burnin = int(n * burnin_fraction)
    return {c: np.array(v[burnin:]) for c, v in data.items()}


SEARCH_DIRS = [os.path.join(BASE, "data"), BASE]

# TaxonSetN -> (fossil-paper node number, short clade label, lower bound, upper bound), in
# Ma. Matches the 12 MRCAPriors shared by allCalibrations.xml / allCalibrations_suggestedPriors.xml
# (de Vries & Beck Table 1; bounds are identical between the uniform and suggested-prior
# versions -- only the shape of the density between the bounds differs).
TAXONSET_INFO = {
    "TaxonSet1":  (1,  "Euarchontoglires",        65.79,   125.816),
    "TaxonSet2":  (3,  "Euarchonta",               65.79,   125.816),
    "TaxonSet3":  (5,  "Primates",                 55.935,  66.095),
    "TaxonSet5":  (13, "Cercopithecidae",          12.47,   25.235),
    "TaxonSet6":  (14, "Colobinae",                8.125,   15.0),
    "TaxonSet7":  (15, "Cercopithecinae",          6.5,     15.0),
    "TaxonSet8":  (16, "Papionini",                5.33,    12.51),
    "TaxonSet9":  (18, "Hominoidea",               13.4,    25.235),
    "TaxonSet10": (19, "Hominidae",                12.3,    25.235),
    "TaxonSet11": (20, "Homo+Pan",                 4.631,   15.0),
    "TaxonSet12": (23, "Callitrichidae+Cebidae",   13.183,  34.5),
    "TaxonSet13": (24, "Cebidae",                  13.032,  34.5),
}

RUN_PAIRS = [
    {
        "name": "allCalibrations (uniform priors)",
        "prior_base": "sample-from-prior_allCalibrations",
        "posterior_base": "allCalibrations",
        "out": os.path.join(BASE, "allCalibrations_prior_vs_posterior.pdf"),
    },
    {
        "name": "allCalibrations_suggestedPriors (offset-exponential + uniform)",
        "prior_base": "sample-from-prior_allCalibrations_suggestedPriors",
        "posterior_base": "allCalibrations_suggestedPriors",
        "out": os.path.join(BASE, "allCalibrations_suggestedPriors_prior_vs_posterior.pdf"),
    },
    {
        "name": "allCalibrations (13 calibrations with uniform prior)",
        "prior_base": "sample-from-prior_allCalibrations_13calibrations",
        "posterior_base": "allCalibrations_13calibrations",
        "out": os.path.join(BASE, "allCalibrations_13calibrations_prior_vs_posterior.pdf"),
    },
]


def find_log(base_name):
    for d in SEARCH_DIRS:
        for ext in (".txt", ".log"):
            path = os.path.join(d, base_name + ext)
            if os.path.exists(path):
                return path
    return None


def plot_pair(name, prior_log, posterior_log, out_pdf):
    print(f"=== {name} ===")
    print(f"  prior:     {prior_log}")
    print(f"  posterior: {posterior_log}")

    with open(prior_log) as f:
        prior_header = f.readline().strip().split("\t")
    with open(posterior_log) as f:
        post_header = f.readline().strip().split("\t")

    mrca_cols = [c for c in TAXONSET_INFO
                 if f"mrca.age({c})" in prior_header and f"mrca.age({c})" in post_header]
    mrca_cols = [f"mrca.age({c})" for c in mrca_cols]
    if not mrca_cols:
        print("  No shared mrca.age(...) columns found between prior and posterior logs; skipping.\n")
        return

    prior_data = read_log_columns(prior_log, mrca_cols, BURNIN_FRACTION)
    post_data = read_log_columns(posterior_log, mrca_cols, BURNIN_FRACTION)

    entries = []
    for col in mrca_cols:
        tsid = col[len("mrca.age("):-1]
        node_number, label, lo, hi = TAXONSET_INFO[tsid]
        entries.append((node_number, label, lo, hi, prior_data[col], post_data[col]))
    entries.sort(key=lambda e: e[0])

    n = len(entries)
    ncols = 3
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 3.5 * nrows), squeeze=False)
    axes = axes.flatten()

    for ax, (node_number, label, lo, hi, prior_samples, post_samples) in zip(axes, entries):
        lo_range = min(lo, prior_samples.min(), post_samples.min())
        hi_range = max(hi, prior_samples.max(), post_samples.max())
        bins = np.linspace(lo_range, hi_range, 60)

        ax.hist(prior_samples, bins=bins, density=True, color="#bab0ac", alpha=0.75,
                label="Prior (sample-from-prior)")
        ax.hist(post_samples, bins=bins, density=True, color="#4e79a7", alpha=0.6,
                label="Posterior")
        ax.axvline(lo, color="#e15759", linestyle=":", linewidth=1.2, alpha=0.8)
        ax.axvline(hi, color="#e15759", linestyle=":", linewidth=1.2, alpha=0.8,
                   label="Specified calibration bounds")
        ax.set_title(f"Node {node_number}: {label}", fontsize=9)
        ax.set_xlabel("Age (Ma)")
        ax.set_yticks([])
        ax.legend(fontsize=6.5, loc="upper right")

        prior_med, post_med = np.median(prior_samples), np.median(post_samples)
        print(f"  Node {node_number} ({label}): prior median={prior_med:.2f}, "
              f"posterior median={post_med:.2f}, bounds=[{lo},{hi}]")

    for ax in axes[len(entries):]:
        ax.axis("off")

    fig.suptitle(f"{name}: prior vs. posterior node ages", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_pdf)
    plt.close(fig)
    print(f"  Wrote {out_pdf}\n")


def main():
    for pair in RUN_PAIRS:
        prior_log = find_log(pair["prior_base"])
        posterior_log = find_log(pair["posterior_base"])
        if prior_log is None or posterior_log is None:
            missing = []
            if prior_log is None:
                missing.append(f"prior ({pair['prior_base']}.txt/.log)")
            if posterior_log is None:
                missing.append(f"posterior ({pair['posterior_base']}.txt/.log)")
            print(f"=== {pair['name']} ===\n  Missing: {', '.join(missing)} -- skipping.\n")
            continue
        plot_pair(pair["name"], prior_log, posterior_log, pair["out"])


if __name__ == "__main__":
    main()
