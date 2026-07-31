#!/usr/bin/env python3
"""
Generate all 8 primates XML variants from the two base XMLs.

  {consistent_prior, uniform_prior}
    x {conditioned, not-conditioned}
    x {with sequences, sample-from-prior (sequences replaced with "?")}

Naming convention:
  [sample-from-prior_]primates_{consistent,uniform}_prior[_not-conditioned].xml

The two base XMLs (with sequences, conditioned) are written as part of the set
so all 8 files share the same generation pipeline.

Usage:
    python make_sample_from_prior.py
"""
import copy
import os
from lxml import etree as ET

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

INPUTS = [
    "primates_consistent_prior.xml",
    "primates_uniform_prior.xml",
]


def set_not_conditioned(root):
    """Add conditionOnCalibrations="false" to CalibratedBirthDeathSkylineModel."""
    nodes = root.xpath(
        './/distribution[@spec="calibratedcpp.CalibratedBirthDeathSkylineModel"]'
    )
    if not nodes:
        print("  WARNING: CalibratedBirthDeathSkylineModel not found.")
    for node in nodes:
        node.set("conditionOnCalibrations", "false")
    return len(nodes)


def strip_sequence_data(root):
    sequences = root.findall('.//sequence[@value]')
    for seq in sequences:
        seq.set("value", "?")
    return len(sequences)


def normalize_log_filenames(root):
    """Replace hardcoded log file names with $(filebase) tokens."""
    for logger in root.findall('.//logger[@fileName]'):
        name = logger.get("fileName")
        if name.endswith(".log"):
            logger.set("fileName", "$(filebase).log")
        elif name.endswith(".trees"):
            logger.set("fileName", "$(filebase).trees")


def make_output_name(stem, conditioned: bool, with_sequences: bool) -> str:
    parts = []
    if not with_sequences:
        parts.append("sample-from-prior")
    parts.append(stem)
    if not conditioned:
        parts.append("not-conditioned")
    return "_".join(parts) + ".xml"


def process_all(input_file):
    parser = ET.XMLParser(remove_blank_text=True)
    base_tree = ET.parse(input_file, parser)
    stem = os.path.splitext(os.path.basename(input_file))[0]

    for conditioned in (True, False):
        for with_sequences in (True, False):
            root = copy.deepcopy(base_tree.getroot())

            if not conditioned:
                set_not_conditioned(root)
            if not with_sequences:
                strip_sequence_data(root)
            normalize_log_filenames(root)

            out_name = make_output_name(stem, conditioned, with_sequences)
            out_path = os.path.join(BASE_DIR, out_name)

            out_tree = ET.ElementTree(root)
            out_tree.write(
                out_path, pretty_print=True, xml_declaration=True, encoding="UTF-8"
            )
            cond_label = "conditioned" if conditioned else "not-conditioned"
            seq_label = "with sequences" if with_sequences else "sample-from-prior"
            print(f"  {cond_label:<17}  {seq_label:<17}  -> {out_name}")


if __name__ == "__main__":
    for filename in INPUTS:
        input_path = os.path.join(BASE_DIR, filename)
        print(f"\n{filename}")
        process_all(input_path)
