#!/usr/bin/env python3
"""
Replace the initial tree in primates_uniform_prior.xml with the tree sampled from
primates_narrow.lphy (i.e., the first tree in primates_narrow_tree.trees).

Usage:
    python change_initial_tree.py

The script reads the first newick from primates_narrow_tree.trees (NEXUS format
produced by SLPhy) and writes it into the stateNode[@id="tree"] attribute of
primates_uniform_prior.xml in place.
"""
import os
import re
from lxml import etree as ET


def extract_newick_from_trees(trees_file):
    """Return the newick of the first tree in a NEXUS .trees file.

    SLPhy writes lines of the form:
        tree TREE_0= [&R] <newick>;
    The [&R] annotation is stripped so BEAST2 receives a plain newick string.
    """
    with open(trees_file) as fh:
        for line in fh:
            line = line.strip()
            if line.lower().startswith("tree "):
                # Remove optional [&R] / [&U] annotation then grab everything up to the trailing ;
                newick = re.sub(r"\[&[RU]\]\s*", "", line.split("=", 1)[1]).strip()
                if newick.endswith(";"):
                    newick = newick[:-1]
                return newick
    raise ValueError(f"No 'tree' line found in {trees_file}")


def replace_initial_tree(target_file, newick, output_file):
    """Overwrite the newick attribute in target_file and save to output_file."""
    parser = ET.XMLParser(remove_blank_text=True)
    tree = ET.parse(target_file, parser)
    root = tree.getroot()

    node = root.find('.//stateNode[@id="tree"]')
    if node is None:
        raise ValueError(f"No <stateNode id='tree'> found in {target_file}")

    old_newick = node.get("newick", "")
    node.set("newick", newick)

    tree.write(output_file, pretty_print=True, xml_declaration=True, encoding="UTF-8")
    print(f"Replaced initial tree in {output_file}")
    print(f"  Old newick prefix: {old_newick[:80]}...")
    print(f"  New newick prefix: {newick[:80]}...")


if __name__ == "__main__":
    base = os.path.dirname(os.path.abspath(__file__))
    trees_file = os.path.join(base, "primates_narrow_tree.trees")
    target_xml = os.path.join(base, "primates_uniform_prior.xml")

    newick = extract_newick_from_trees(trees_file)
    replace_initial_tree(target_xml, newick, target_xml)
