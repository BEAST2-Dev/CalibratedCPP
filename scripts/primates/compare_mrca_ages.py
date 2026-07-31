#!/usr/bin/env python3
"""
Build a comparison table of estimated MRCA node ages (mean + 95% HPD interval) across every
BEAST log file in data/ (and this script's own directory), so different runs -- different
calibration schemes, priors, alignments, starting trees -- can be compared side by side for
the same node.

Self-contained (no imports from other scripts in this folder).

Usage:
    python compare_mrca_ages.py

Output:
    mrca_age_comparison.csv          -- long format: one row per (run, node)
    printed to stdout                -- the same data, grouped by node for quick reading

Any file without at least one "Sample" + "mrca.age(...)" column is skipped (so this can be
pointed at a directory containing non-BEAST-log files without erroring).
"""
import csv
import glob
import os

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
SEARCH_DIRS = [os.path.join(BASE, "data"), BASE]
BURNIN_FRACTION = 0.1
OUT_CSV = os.path.join(BASE, "mrca_age_comparison.csv")

# TaxonSetN -> (fossil-paper node number, short clade label), for labelling only. Covers all
# 13 de Vries & Beck calibrations used across this project's XML family (including
# TaxonSet4/Lorisiformes). A TaxonSetN not in this map is still included in the table, just
# labelled by its raw id.
TAXONSET_INFO = {
    "TaxonSet1":  (1,  "Euarchontoglires"),
    "TaxonSet2":  (3,  "Euarchonta"),
    "TaxonSet3":  (5,  "Primates"),
    "TaxonSet4":  (6,  "Lorisiformes"),
    "TaxonSet5":  (13, "Cercopithecidae"),
    "TaxonSet6":  (14, "Colobinae"),
    "TaxonSet7":  (15, "Cercopithecinae"),
    "TaxonSet8":  (16, "Papionini"),
    "TaxonSet9":  (18, "Hominoidea"),
    "TaxonSet10": (19, "Hominidae"),
    "TaxonSet11": (20, "Homo+Pan"),
    "TaxonSet12": (23, "Callitrichidae+Cebidae"),
    "TaxonSet13": (24, "Cebidae"),
}


def label_for(taxonset_id):
    info = TAXONSET_INFO.get(taxonset_id)
    if info is None:
        return None, taxonset_id
    return info


def hpd_interval(samples, mass=0.95):
    """Narrowest interval containing `mass` fraction of the (sorted) samples."""
    s = np.sort(samples)
    n = len(s)
    interval_idx = int(np.floor(mass * n))
    if interval_idx >= n:
        return s[0], s[-1]
    widths = s[interval_idx:] - s[: n - interval_idx]
    best = np.argmin(widths)
    return s[best], s[best + interval_idx]


def discover_log_files():
    """{run_name: path}, deduplicated by basename (data/ takes priority over BASE)."""
    found = {}
    for d in SEARCH_DIRS:
        for path in sorted(glob.glob(os.path.join(d, "*.txt"))) + sorted(glob.glob(os.path.join(d, "*.log"))):
            name = os.path.splitext(os.path.basename(path))[0]
            found.setdefault(name, path)
    return found


def read_mrca_columns(path):
    """Return {TaxonSetId: np.array(post-burnin samples)} for every mrca.age(...) column,
    or None if this file isn't a BEAST log (no Sample/mrca.age columns)."""
    with open(path) as f:
        first_line = f.readline()
    if not first_line.strip():
        return None
    header = first_line.strip().split("\t")
    if not header or header[0] != "Sample":
        return None
    mrca_cols = [c for c in header if c.startswith("mrca.age(")]
    if not mrca_cols:
        return None

    col_idx = {c: header.index(c) for c in mrca_cols}
    data = {c: [] for c in mrca_cols}
    with open(path) as f:
        f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(header):
                continue
            for c in mrca_cols:
                try:
                    data[c].append(float(parts[col_idx[c]]))
                except ValueError:
                    pass

    n = len(next(iter(data.values()), []))
    if n == 0:
        return None
    burnin = int(n * BURNIN_FRACTION)

    out = {}
    for c in mrca_cols:
        taxonset_id = c[len("mrca.age("):-1]
        arr = np.array(data[c][burnin:])
        if len(arr) > 0:
            out[taxonset_id] = arr
    return out


def main():
    log_files = discover_log_files()
    print(f"Found {len(log_files)} candidate log file(s) in {', '.join(SEARCH_DIRS)}\n")

    rows = []  # (node_number, clade_label, taxonset_id, run_name, n, mean, hpd_lo, hpd_hi)
    for run_name, path in sorted(log_files.items()):
        cols = read_mrca_columns(path)
        if cols is None:
            continue
        for taxonset_id, samples in cols.items():
            node_number, label = label_for(taxonset_id)
            mean = float(np.mean(samples))
            hpd_lo, hpd_hi = hpd_interval(samples)
            rows.append((node_number, label, taxonset_id, run_name, len(samples), mean, hpd_lo, hpd_hi))

    if not rows:
        print("No BEAST log files with mrca.age(...) columns found.")
        return

    # Write long-format CSV.
    rows.sort(key=lambda r: (r[0] is None, r[0] if r[0] is not None else 0, r[2], r[3]))
    with open(OUT_CSV, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["node_number", "clade", "taxonset", "run", "n_samples", "mean", "hpd95_lower", "hpd95_upper"])
        for node_number, label, taxonset_id, run_name, n, mean, hpd_lo, hpd_hi in rows:
            w.writerow([node_number if node_number is not None else "", label, taxonset_id, run_name,
                        n, f"{mean:.4f}", f"{hpd_lo:.4f}", f"{hpd_hi:.4f}"])
    print(f"Wrote {OUT_CSV}\n")

    # Print grouped by node, one mini-table per node, for quick visual comparison across runs.
    current_key = None
    for node_number, label, taxonset_id, run_name, n, mean, hpd_lo, hpd_hi in rows:
        key = (node_number, taxonset_id)
        if key != current_key:
            if current_key is not None:
                print()
            tag = f"Node {node_number}" if node_number is not None else "Node ?"
            print(f"=== {tag}: {label} ({taxonset_id}) ===")
            print(f"  {'run':<45} {'n':>7} {'mean':>10} {'95% HPD':>22}")
            current_key = key
        print(f"  {run_name:<45} {n:>7} {mean:>10.2f}   [{hpd_lo:>8.2f}, {hpd_hi:>8.2f}]")
    print()


if __name__ == "__main__":
    main()
