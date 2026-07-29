#!/usr/bin/env python3
"""
Extract a curated, moderate-sized subset of sites from a NEXUS DNA matrix: the
contiguous window of a given size (default ~100kb) with the lowest overall missing-data
rate (N/-/? characters, summed across all taxa), rather than an arbitrary prefix. This
alignment is a concatenation of many loci with very uneven missing-data coverage across
taxa, so picking a fixed-size window by minimum missing data is a much better way to get
a "clean" moderate-sized subset than just taking the first N% of sites.

Usage:
    python truncate_alignment.py [input.nex] [output.nex] [window_size]

Defaults:
    input.nex   = primates.nex
    output.nex  = primates_low_missing_100kb.nex
    window_size = 100000
"""
import re
import sys

import numpy as np

MISSING_BYTES = [b"N", b"-", b"?"]


def parse_nexus_matrix(text):
    """Return (ntax, nchar, format_line, taxa) where taxa is a list of (name, sequence)."""
    dims = re.search(r"Dimensions\s+ntax=(\d+)\s+nchar=(\d+)\s*;", text, re.IGNORECASE)
    if not dims:
        raise ValueError("Could not find a Dimensions line (ntax=.../nchar=...)")
    ntax, nchar = int(dims.group(1)), int(dims.group(2))

    fmt = re.search(r"(Format[^\n]*;)", text, re.IGNORECASE)
    format_line = fmt.group(1).strip() if fmt else "Format datatype=nucleotide gap=-;"

    matrix_start = re.search(r"Matrix\s*", text, re.IGNORECASE)
    if not matrix_start:
        raise ValueError("Could not find the Matrix keyword")

    lines = text[matrix_start.end():].splitlines()
    taxa = []
    for line in lines:
        line = line.strip()
        if not line or line in (";", "End;") or line.lower().startswith("end"):
            continue
        if line == ";":
            break
        name, seq = line.split(None, 1)
        taxa.append((name, seq.strip().rstrip(";").strip()))
        if len(taxa) == ntax:
            break

    if len(taxa) != ntax:
        raise ValueError(f"Expected {ntax} taxa, found {len(taxa)}")
    for name, seq in taxa:
        if len(seq) != nchar:
            raise ValueError(f"Taxon {name}: expected {nchar} characters, got {len(seq)}")

    return ntax, nchar, format_line, taxa


def find_best_window(taxa, nchar, window_size):
    """Return (start, size, total_missing) for the window_size-site contiguous window
    (out of nchar) with the fewest missing (N/-/?) characters summed across all taxa."""
    window_size = min(window_size, nchar)

    missing_per_site = np.zeros(nchar, dtype=np.int32)
    for _, seq in taxa:
        arr = np.frombuffer(seq.upper().encode("ascii"), dtype="S1")
        missing_per_site += np.isin(arr, MISSING_BYTES)

    # Prefix sums give O(1) window-sum lookups instead of O(window_size) per window.
    cumsum = np.concatenate(([0], np.cumsum(missing_per_site)))
    n_windows = nchar - window_size + 1
    window_sums = cumsum[window_size:] - cumsum[:n_windows]

    best_start = int(np.argmin(window_sums))
    return best_start, window_size, int(window_sums[best_start])


def write_nexus(path, ntax, nchar, format_line, taxa):
    name_width = max(len(name) for name, _ in taxa) + 3
    with open(path, "w", newline="\n") as fh:
        fh.write("#NEXUS\n\n")
        fh.write("BEGIN DATA;\n")
        fh.write(f"\tDimensions ntax={ntax} nchar={nchar};\n")
        fh.write(f"\t{format_line}\n")
        fh.write("\tMatrix\n")
        for name, seq in taxa:
            fh.write(f"\t{name.ljust(name_width)}{seq}\n")
        fh.write(";\n")
        fh.write("END;\n")


def main():
    input_path = sys.argv[1] if len(sys.argv) > 1 else "primates.nex"
    output_path = sys.argv[2] if len(sys.argv) > 2 else "primates_low_missing_100kb.nex"
    window_size = int(sys.argv[3]) if len(sys.argv) > 3 else 100_000

    with open(input_path) as fh:
        text = fh.read()

    ntax, nchar, format_line, taxa = parse_nexus_matrix(text)
    start, size, total_missing = find_best_window(taxa, nchar, window_size)

    subset = [(name, seq[start:start + size]) for name, seq in taxa]
    write_nexus(output_path, ntax, size, format_line, subset)

    total_chars = ntax * size
    print(f"Wrote {output_path}: {ntax} taxa, {size} sites "
          f"(best window [{start}, {start + size}) of {nchar}), "
          f"missing data = {total_missing}/{total_chars} ({total_missing / total_chars:.3%})")


if __name__ == "__main__":
    main()
