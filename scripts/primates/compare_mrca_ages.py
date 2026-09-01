#!/usr/bin/env python3
"""
Build a comparison table of estimated MRCA node ages (mean + 95% HPD interval) across every
BEAST log file in data/, so different runs -- different
calibration schemes, priors, alignments, starting trees -- can be compared side by side for
the same node.

Clade labels are resolved from each run's own XML (via primates_common): the TaxonSetN
numbering is not shared between the nogapN and codon XML families, nor between models.

Usage:
    python compare_mrca_ages.py [DIR ...]     # default: data/ and this folder

Output:
    mrca_age_comparison.csv          -- long format: one row per (run, node)
    printed to stdout                -- the same data, grouped by node for quick reading

Any file without at least one "Sample" + "mrca.age(...)" column is skipped (so this can be
pointed at a directory containing non-BEAST-log files without erroring).
"""
import csv
import glob
import os
import sys

import numpy as np

import primates_common as pc

BASE = os.path.dirname(os.path.abspath(__file__))
SEARCH_DIRS = [os.path.join(BASE, "data")]
BURNIN_FRACTION = 0.1
OUT_CSV = os.path.join(BASE, "mrca_age_comparison.csv")

# Fossil-paper node numbers, keyed by clade name (stable, unlike the TaxonSetN ids).
NODE_NUMBERS = {
    "Euarchontoglires": 1, "Euarchonta": 3, "Primates": 5, "Lorisiformes": 6,
    "Cercopithecidae": 13, "Colobinae": 14, "Cercopithecinae": 15, "Papionini": 16,
    "Hominoidea": 18, "Hominidae": 19, "Homo_Pan": 20,
    "Callitrichidae_Cebidae": 23, "Cebidae": 24,
}

_LABELS = {}


def labels_for_run(run_name):
    """{TaxonSetId: (node_number, clade_name)} resolved from this run's own XML."""
    if run_name in _LABELS:
        return _LABELS[run_name]
    out = {}
    for dataset in pc.DATASETS:
        pc.use(dataset)
        for model in pc.MODELS:
            if any(pc.stem(model, c, p) == run_name
                   for c in ("true", "false") for p in (False, True)):
                names = pc.clade_names()
                for tid, taxa in pc.taxon_sets(model).items():
                    if names.get(taxa):
                        out[tid] = (NODE_NUMBERS.get(names[taxa]), names[taxa])
                break
        if out:
            break
    _LABELS[run_name] = out
    return out


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


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if argv:
        global SEARCH_DIRS, OUT_CSV
        SEARCH_DIRS = [d if os.path.isabs(d) else os.path.join(BASE, d) for d in argv]
        OUT_CSV = os.path.join(BASE, "mrca_age_comparison_%s.csv"
                               % "_".join(os.path.basename(d.rstrip("/")) for d in SEARCH_DIRS))
    log_files = discover_log_files()
    print(f"Found {len(log_files)} candidate log file(s) in {', '.join(SEARCH_DIRS)}\n")

    rows = []  # (node_number, clade_label, taxonset_id, run_name, n, mean, hpd_lo, hpd_hi)
    for run_name, path in sorted(log_files.items()):
        cols = read_mrca_columns(path)
        if cols is None:
            continue
        labels = labels_for_run(run_name)
        for taxonset_id, samples in cols.items():
            node_number, label = labels.get(taxonset_id, (None, taxonset_id))
            mean = float(np.mean(samples))
            hpd_lo, hpd_hi = hpd_interval(samples)
            rows.append((node_number, label, taxonset_id, run_name, len(samples), mean, hpd_lo, hpd_hi))

    if not rows:
        print("No BEAST log files with mrca.age(...) columns found.")
        return

    # Write long-format CSV.
    rows.sort(key=lambda r: (r[0] is None, r[0] if r[0] is not None else 0, r[1], r[3]))
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
        key = (node_number, label)
        if key != current_key:
            if current_key is not None:
                print()
            tag = f"Node {node_number}" if node_number is not None else "Node ?"
            print(f"=== {tag}: {label} ===")
            print(f"  {'run':<45} {'n':>7} {'mean':>10} {'95% HPD':>22}")
            current_key = key
        print(f"  {run_name:<45} {n:>7} {mean:>10.2f}   [{hpd_lo:>8.2f}, {hpd_hi:>8.2f}]")
    print()


if __name__ == "__main__":
    main()
