#!/usr/bin/env python3
"""Robinson-Foulds distances between the primates summary trees and reference trees.

Reference topology: de Vries and Beck (2023), "Fossil calibrations for primate
divergences", Palaeontologia Electronica 26(2):a25
(https://palaeo-electronica.org/content/2023/3777-primate-fossil-calibrations),
figure 1, restricted to the 29 taxa of this dataset and encoded in
primates_deVries.newick.
"""

import argparse
import os
import re

import primates_common as pc

REF_TREES = {
    "deVries": os.path.join(pc.BASE, "refTrees", "primates_deVries.newick"),
    "ASTRAL": os.path.join(pc.BASE, "refTrees", "ASTRAL_Tree_AVGdates.tre"),
    "IQ-TREE-nogapN": os.path.join(pc.BASE, "iqtree", "primates_nogapN.treefile"),
    "IQ-TREE-full": os.path.join(pc.BASE, "iqtree", "primates.treefile"),
}

_META = re.compile(r"\[[^\]]*\]")


def _clades(newick, translate=None, outgroup=None):
    """{frozenset(taxa)} for every internal node, trivial ones excluded.

    A multifurcating root (as written by IQ-TREE, which leaves the root
    trifurcating even with -o) is resolved by grouping the children that do not
    hold `outgroup`, so the clade counts match a rooted tree.
    """
    out = set()
    root_children = []
    pos = 0

    def skip_label():
        nonlocal pos
        while pos < len(newick) and newick[pos] not in "(),;":
            pos += 1

    def node():
        nonlocal pos
        while pos < len(newick) and newick[pos] in " \t\n":
            pos += 1
        if newick[pos] == "(":
            pos += 1
            taxa = set()
            depth0 = pos == 1
            while True:
                child = node()
                taxa |= child
                if depth0:
                    root_children.append(child)
                if newick[pos] == ",":
                    pos += 1
                    continue
                break
            pos += 1  # ')'
            skip_label()
            out.add(frozenset(taxa))
            return taxa
        start = pos
        while pos < len(newick) and newick[pos] not in "(),:":
            pos += 1
        label = newick[start:pos].strip()
        skip_label()
        return {translate.get(label, label) if translate else label}

    node()
    if outgroup and len(root_children) > 2:
        rest = set()
        for child in root_children:
            if outgroup not in child:
                rest |= child
        out.add(frozenset(rest))
    return {c for c in out if len(c) > 1}


def read_tree(path, outgroup="Mus_musculus"):
    """Clade set of the first tree in a newick or NEXUS file."""
    with open(path) as fh:
        text = fh.read()
    translate = None
    if "Translate" in text:
        block = text.split("Translate")[1].split(";")[0]
        translate = dict(re.findall(r"(\w+)\s+([\w.]+)\s*,?", block))
    m = re.search(r"^\s*tree\s[^=]*=", text, re.M | re.I)
    newick = _META.sub("", text[m.end():] if m else text).strip()
    if ";" in newick:
        newick = newick[:newick.index(";")]
    return _clades(newick, translate, outgroup)


def rf(a, b, taxa=None):
    """RF distance: clades present in exactly one of the two trees.

    With `taxa` given the comparison is unrooted: each clade is replaced by the
    smaller side of the split it induces, so trees that differ only in where the
    root sits come out identical.
    """
    if taxa is not None:
        a, b = _splits(a, taxa), _splits(b, taxa)
    return len(a ^ b)


def _splits(clades, taxa):
    out = set()
    for c in clades:
        other = frozenset(taxa - c)
        if len(c) > 1 and len(other) > 1:
            out.add(min(c, other, key=lambda s: (len(s), sorted(s))))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dataset", default=pc.DATASET, choices=sorted(pc.DATASETS),
                    help="which set of summary trees to compare")
    ap.add_argument("--unrooted", action="store_true",
                    help="ignore root placement (compare splits, not clades)")
    ap.add_argument("--ref", action="append", metavar="NAME=FILE",
                    help="extra reference tree; repeatable")
    args = ap.parse_args()
    pc.use(args.dataset)

    refs = dict(REF_TREES)
    for spec in args.ref or []:
        name, _, path = spec.partition("=")
        refs[name] = path or name
    trees = {n: read_tree(p) for n, p in refs.items() if os.path.exists(p)}
    if not trees:
        raise SystemExit("no reference tree found")
    for model in pc.MODELS:
        for cond in ("true", "false"):
            path = pc.summary_tree_path(model, cond)
            if os.path.exists(path):
                trees[model + pc.COND[cond]] = read_tree(path)

    names = list(trees)
    taxa = {t for c in trees[names[0]] for t in c}
    for name in names[1:]:
        other = {t for c in trees[name] for t in c}
        if other != taxa:
            raise SystemExit("%s: taxon set differs from %s (%s)" % (
                name, names[0], sorted(other ^ taxa)))

    key = taxa if args.unrooted else None
    width = max(len(n) for n in names)
    print("Robinson-Foulds distances (%s, %d taxa)\n" % (
        "unrooted splits" if args.unrooted else "rooted clades", len(taxa)))
    print(" " * width + "".join("%4d" % i for i in range(1, len(names) + 1)))
    for i, a in enumerate(names, 1):
        print("%*s" % (width, a) + "".join(
            "%4d" % rf(trees[a], trees[b], key) for b in names) + "   (%d)" % i)

    for ref in refs:
        if ref not in trees:
            continue
        print("\nvs %s:" % ref)
        for name in names:
            if name == ref:
                continue
            d = rf(trees[ref], trees[name], key)
            print("  %-*s RF = %2d  (normalised %.3f)" % (
                width, name, d, d / (2.0 * (len(taxa) - 2))))
            for label, clades in (("missing", trees[ref] - trees[name]),
                                  ("extra", trees[name] - trees[ref])):
                for c in sorted(clades, key=len):
                    print("    %-7s %s" % (label, ", ".join(sorted(c))))


if __name__ == "__main__":
    main()
