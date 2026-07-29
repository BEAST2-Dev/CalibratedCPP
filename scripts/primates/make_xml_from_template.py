#!/usr/bin/env python3
"""
Build a new BEAST XML from a template by swapping in a different alignment (NEXUS) and
starting tree (raw newick), while keeping everything else (calibrations, priors,
operators, etc.) from the template unchanged.

Usage:
    python make_xml_from_template.py <template.xml> <alignment.nex> <tree.nwk> <output.xml>
"""
import re
import sys


def parse_nexus_matrix(text):
    """Return {taxon_name: sequence} from a NEXUS alignment.

    Handles both the single-block "BEGIN DATA; Dimensions ntax=.. nchar=..; ... Matrix"
    style, and the two-block "begin taxa; dimensions ntax=..; ... begin characters;
    dimensions nchar=..; ... matrix" style (e.g. LPhy's *_true_D.nexus output) — ntax and
    nchar are searched for independently rather than requiring them on the same line.
    """
    ntax_m = re.search(r"ntax\s*=\s*(\d+)", text, re.IGNORECASE)
    nchar_m = re.search(r"nchar\s*=\s*(\d+)", text, re.IGNORECASE)
    if not ntax_m or not nchar_m:
        raise ValueError("Could not find ntax=... and nchar=... in the NEXUS file")
    ntax, nchar = int(ntax_m.group(1)), int(nchar_m.group(1))

    matrix_start = re.search(r"Matrix\s*", text, re.IGNORECASE)
    if not matrix_start:
        raise ValueError("Could not find the Matrix keyword")

    lines = text[matrix_start.end():].splitlines()
    taxa = {}
    for line in lines:
        line = line.strip()
        if not line or line in (";", "End;") or line.lower().startswith("end"):
            continue
        if line == ";":
            break
        name, seq = line.split(None, 1)
        taxa[name] = seq.strip().rstrip(";").strip()
        if len(taxa) == ntax:
            break

    if len(taxa) != ntax:
        raise ValueError(f"Expected {ntax} taxa, found {len(taxa)}")
    for name, seq in taxa.items():
        if len(seq) != nchar:
            raise ValueError(f"Taxon {name}: expected {nchar} characters, got {len(seq)}")

    return taxa


def replace_alignment(xml_text, taxon_sequences):
    """Replace each <sequence ... taxon="X" ... value="..."/> with the new sequence for X."""
    pattern = re.compile(r'(<sequence\b[^>]*\btaxon="([^"]+)"[^>]*\bvalue=")([^"]*)("[^>]*/>)')

    missing = []

    def repl(m):
        prefix, taxon, _old_value, suffix = m.groups()
        if taxon not in taxon_sequences:
            missing.append(taxon)
            return m.group(0)
        return f"{prefix}{taxon_sequences[taxon]}{suffix}"

    new_text, n = pattern.subn(repl, xml_text)
    if missing:
        raise ValueError(f"Taxa in template XML not found in new alignment: {missing}")
    if n == 0:
        raise ValueError("No <sequence ... taxon=... value=.../> elements found in template XML")
    return new_text, n


def replace_tree(xml_text, newick):
    """Replace the newick attribute of <stateNode id="tree" ...> with the new tree."""
    pattern = re.compile(r'(<stateNode id="tree"[^>]*\bnewick=")([^"]*)(")')
    new_text, n = pattern.subn(lambda m: f"{m.group(1)}{newick}{m.group(3)}", xml_text)
    if n == 0:
        raise ValueError('Could not find <stateNode id="tree" ... newick="..."> in template XML')
    return new_text, n


def main():
    if len(sys.argv) != 5:
        print(__doc__)
        sys.exit(1)

    template_path, alignment_path, tree_path, output_path = sys.argv[1:5]

    with open(template_path) as f:
        xml_text = f.read()
    with open(alignment_path) as f:
        alignment_text = f.read()
    with open(tree_path) as f:
        newick = f.read().strip()

    taxon_sequences = parse_nexus_matrix(alignment_text)
    xml_text, n_seq = replace_alignment(xml_text, taxon_sequences)
    xml_text, n_tree = replace_tree(xml_text, newick)

    with open(output_path, "w") as f:
        f.write(xml_text)

    print(f"Wrote {output_path}: replaced {n_seq} sequences and {n_tree} starting tree "
          f"from template {template_path}")


if __name__ == "__main__":
    main()
