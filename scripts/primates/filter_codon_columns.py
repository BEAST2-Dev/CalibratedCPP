#!/usr/bin/env python3
"""
Build a NEXUS alignment by filtering out whole CODONS, keeping the reading frame intact.

filter_gap_columns.py drops individual columns, which is fine for an unpartitioned
analysis but destroys the codon structure: after it runs, site i is no longer at a
predictable codon position, so a charset like "1-.\\3" no longer selects 1st positions.
This script works on triplets instead -- a codon survives only if all three of its
columns are clean in every taxon -- so the surviving sites stay in frame and the codon
positions of the output are, in order, 1, 2, 3, 1, 2, 3, ...

Two filter modes ('?' is deliberately ignored, matching filter_gap_columns.py):
    nogap   -- drop a codon if any taxon has '-' in it
    nogapN  -- drop a codon if any taxon has '-' or 'N' in it

Usage:
    python filter_codon_columns.py [input.nex] [output.nex] [mode] [frame_start]

Defaults:
    input.nex   = data/primates.nex
    output.nex  = data/primates_codon_nogapN.nex (matching mode)
    mode        = nogapN
    frame_start = 1   (1-based column where the first full codon begins)
"""
import re
import sys

import numpy as np

MODE_CHARS = {
    "nogap": ["-"],
    "nogapN": ["-", "N"],
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

    taxa = []
    for line in text[matrix_start.end():].splitlines():
        line = line.strip()
        if not line or line == ";" or line.lower().startswith("end"):
            continue
        name, seq = line.split(None, 1)
        taxa.append((name, seq.strip().rstrip(";").strip()))
        if len(taxa) == ntax:
            break

    if len(taxa) != ntax:
        raise ValueError("Expected %d taxa, parsed %d" % (ntax, len(taxa)))
    for name, seq in taxa:
        if len(seq) != nchar:
            raise ValueError("%s has %d sites, expected %d" % (name, len(seq), nchar))
    return ntax, nchar, format_line, taxa


def keep_codons(taxa, nchar, bad_chars, frame_start):
    """Boolean mask over codons; a codon is kept only if no taxon has a bad char in it."""
    matrix = np.array([list(seq.upper()) for _, seq in taxa])
    offset = frame_start - 1
    ncodon = (nchar - offset) // 3
    trimmed = nchar - offset - ncodon * 3
    if trimmed:
        print("warning: %d trailing site(s) are not a whole codon and were dropped" % trimmed)

    body = matrix[:, offset:offset + ncodon * 3]
    bad = np.isin(body, [c.upper() for c in bad_chars])
    # any bad character anywhere in a taxon x codon block condemns the whole codon
    bad_codon = bad.reshape(len(taxa), ncodon, 3).any(axis=2).any(axis=0)
    return body, ~bad_codon


def write_nexus(path, format_line, taxa, sequences):
    width = max(len(name) for name, _ in taxa) + 4
    with open(path, "w") as out:
        out.write("#NEXUS\n\nBEGIN DATA;\n")
        out.write("\tDimensions ntax=%d nchar=%d;\n" % (len(taxa), len(sequences[0])))
        out.write("\t%s\n\tMatrix\n" % format_line)
        for (name, _), seq in zip(taxa, sequences):
            out.write("\t%s%s\n" % (name.ljust(width), seq))
        out.write("\t;\nEND;\n")


def main():
    args = sys.argv[1:]
    in_path = args[0] if len(args) > 0 else "data/primates.nex"
    mode = args[2] if len(args) > 2 else "nogapN"
    if mode not in MODE_CHARS:
        sys.exit("unknown mode %r; use one of %s" % (mode, ", ".join(MODE_CHARS)))
    out_path = args[1] if len(args) > 1 else "data/primates_codon_%s.nex" % mode
    frame_start = int(args[3]) if len(args) > 3 else 1

    with open(in_path) as fh:
        text = fh.read()
    ntax, nchar, format_line, taxa = parse_nexus_matrix(text)

    body, keep = keep_codons(taxa, nchar, MODE_CHARS[mode], frame_start)
    if not keep.any():
        sys.exit("no codon survived the %s filter" % mode)

    kept_idx = np.repeat(keep, 3)
    sequences = ["".join(row[kept_idx]) for row in body]

    write_nexus(out_path, format_line, taxa, sequences)
    print("%s: %d taxa, %d sites (%d codons)" % (in_path, ntax, nchar, len(keep)))
    print("%s: %d sites (%d codons kept, %.1f%%), frame preserved as 1,2,3,..."
          % (out_path, len(sequences[0]), keep.sum(), 100.0 * keep.sum() / len(keep)))


if __name__ == "__main__":
    main()
