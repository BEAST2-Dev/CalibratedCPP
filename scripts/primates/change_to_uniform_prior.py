#!/usr/bin/env python3
"""
Convert primates_consistent_prior.xml -> primates_uniform_prior.xml by replacing
the CalibrationPrior distribution with individual uniform MRCAPrior distributions.

The CalibrationPrior block is removed and each CalibrationCladePrior calibration
becomes a MRCAPrior with a Uniform(lower, upper) distribution.
"""
import os
from lxml import etree as ET


def get_calibration_info(cal_node):
    """Return (taxonset_id, lower, upper) for a CalibrationCladePrior element."""
    # taxa="@TaxonSetX" attribute reference
    taxa_ref = cal_node.get("taxa")
    if taxa_ref and taxa_ref.startswith("@"):
        taxset_id = taxa_ref[1:]
    else:
        taxa_el = cal_node.find(".//taxa")
        taxset_id = taxa_el.get("id") if taxa_el is not None else None

    def get_age(tag):
        el = cal_node.find(f"./{tag}")
        if el is not None:
            return el.get("value") or (el.text.strip() if el.text else None)
        return None

    lower = get_age("lowerAge")
    upper = get_age("upperAge")
    return taxset_id, lower, upper


def make_mrca_prior(idx, taxset_id, lower, upper):
    """Build a <distribution MRCAPrior> element with an inline Uniform distr."""
    dist = ET.Element(
        "distribution",
        attrib={
            "id": f"MRCAPrior{idx}",
            "spec": "beast.base.spec.evolution.tree.MRCAPrior",
            "monophyletic": "true",
            "taxonset": f"@{taxset_id}",
            "tree": "@tree",
        },
    )
    ET.SubElement(
        dist,
        "distr",
        attrib={
            "id": f"Uniform{idx}",
            "spec": "beast.base.spec.inference.distribution.Uniform",
            "lower": str(lower),
            "upper": str(upper),
        },
    )
    return dist


def add_mrca_logs(root, mrca_ids):
    """Append <log idref="MRCAPriorX"/> entries to the file Logger."""
    file_logger = root.find('.//logger[@id="Logger"]')
    if file_logger is None:
        print("WARNING: Could not find <logger id='Logger'> – skipping log entries.")
        return
    for mrca_id in mrca_ids:
        ET.SubElement(file_logger, "log", attrib={"idref": mrca_id})
    print(f"Added {len(mrca_ids)} log entries to Logger.")


def convert(input_file, output_file):
    parser = ET.XMLParser(remove_blank_text=True)
    tree = ET.parse(input_file, parser)
    root = tree.getroot()

    cal_prior_nodes = root.xpath(
        './/distribution[@spec="calibrationprior.CalibrationPrior"]'
    )
    if not cal_prior_nodes:
        print("WARNING: No CalibrationPrior distribution found – nothing to convert.")
        return

    all_mrca_ids = []
    for cal_prior in cal_prior_nodes:
        parent = cal_prior.getparent()
        insert_idx = list(parent).index(cal_prior)

        calibrations = cal_prior.findall(
            './/calibration[@spec="calibrationprior.CalibrationCladePrior"]'
        )

        new_elems = []
        for i, cal in enumerate(calibrations):
            taxset_id, lower, upper = get_calibration_info(cal)
            if not all([taxset_id, lower, upper]):
                print(
                    f"  WARNING: skipping calibration id='{cal.get('id')}' "
                    f"– missing lower={lower}, upper={upper}, taxset={taxset_id}"
                )
                continue
            new_elems.append(make_mrca_prior(i, taxset_id, lower, upper))

        parent.remove(cal_prior)
        for i, elem in enumerate(new_elems):
            parent.insert(insert_idx + i, elem)

        all_mrca_ids.extend(e.get("id") for e in new_elems)
        print(
            f"Replaced CalibrationPrior id='{cal_prior.get('id')}' "
            f"with {len(new_elems)} MRCAPrior element(s)."
        )

    add_mrca_logs(root, all_mrca_ids)

    tree.write(output_file, pretty_print=True, xml_declaration=True, encoding="UTF-8")
    print(f"Written: {output_file}")


if __name__ == "__main__":
    base = os.path.dirname(os.path.abspath(__file__))
    input_file = os.path.join(base, "primates_consistent_prior.xml")
    output_file = os.path.join(base, "primates_uniform_prior.xml")
    convert(input_file, output_file)
