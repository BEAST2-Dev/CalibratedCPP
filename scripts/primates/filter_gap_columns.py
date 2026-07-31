#!/usr/bin/env python3
"""
Build a NEXUS alignment by filtering out COLUMNS (sites), rather than truncating to a
contiguous region like truncate_alignment.py does. A column survives only if every single
taxon has a real base call there -- if any one taxon has a missing character at that site,
the whole column is dropped. Surviving columns are concatenated in their original order
(so, unlike truncate_alignment.py, the result is not one unbroken genomic region -- it's
whatever sites pass the filter, scattered across the full alignment). No target output
length is enforced; the result is whatever survives.

Two filter modes ('?' is deliberately ignored -- it doesn't occur in this project's
alignments, and semantically it's neither a gap nor an ambiguous base call, just a generic
"unknown" symbol, so it isn't a clean fit for either tier):
    nogap   -- drop a column if any taxon has '-' (a true alignment gap) there
    nogapN  -- drop a column if any taxon has '-' or 'N' (gap or ambiguous base) there

Usage:
    python filter_gap_columns.py [input.nex] [output.nex] [mode]

Defaults:
    input.nex  = primates.nex
    output.nex = primates_nogap.nex (or primates_nogapN.nex, matching mode)
    mode       = nogap
"""
import re
import sys

import numpy as np

MODE_BYTES = {
    "nogap":  [b"-"],
    "nogapN": [b"-", b"N"],
}


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


def find_surviving_columns(taxa, nchar, missing_bytes):
    """Return a boolean mask (nchar,), True where NO taxon has a missing character."""
    any_missing = np.zeros(nchar, dtype=bool)
    for _, seq in taxa:
        arr = np.frombuffer(seq.upper().encode("ascii"), dtype="S1")
        any_missing |= np.isin(arr, missing_bytes)
    return ~any_missing


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
    mode = sys.argv[3] if len(sys.argv) > 3 else "nogap"
    if mode not in MODE_BYTES:
        raise ValueError(f"mode must be one of {list(MODE_BYTES)}, got {mode!r}")
    output_path = sys.argv[2] if len(sys.argv) > 2 else f"primates_{mode}.nex"

    with open(input_path) as fh:
        text = fh.read()

    ntax, nchar, format_line, taxa = parse_nexus_matrix(text)
    mask = find_surviving_columns(taxa, nchar, MODE_BYTES[mode])
    n_surviving = int(mask.sum())

    idx = np.where(mask)[0]
    subset = [(name, "".join(seq[i] for i in idx)) for name, seq in taxa]
    write_nexus(output_path, ntax, n_surviving, format_line, subset)

    print(f"Wrote {output_path}: {ntax} taxa, {n_surviving}/{nchar} sites survived "
          f"the '{mode}' filter ({n_surviving / nchar:.3%}) -- "
          f"dropped any column where a taxon had {MODE_BYTES[mode]}")


if __name__ == "__main__":
    main()
