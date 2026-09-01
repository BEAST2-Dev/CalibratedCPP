#!/usr/bin/env python3
"""Shared helpers for the primates conditionOnCalibrations figures.

Handles the three model XMLs (calibrationPrior / suggestedPrior / uniformPrior),
their TaxonSet definitions, and TreeAnnotator "-height CA" summary trees.
"""

import os
import re
import xml.etree.ElementTree as ET

BASE = os.path.dirname(os.path.abspath(__file__))
XML_ROOT = os.path.join(BASE, "..", "..", "calibratedcpp-beast", "src", "test",
                        "resources", "calibratedcpp", "examples", "primates", "xmls")

MODELS = ["calibrationPrior", "suggestedPrior", "uniformPrior"]
MODEL_LABELS = {
    "calibrationPrior": "Joint (Beta-LogNormal) calibration prior",
    "suggestedPrior": "de Vries and Beck independent calibration priors",
    "uniformPrior": "Uniform independent calibration priors",
}

# Shared palette, so the same colour means the same thing in every figure.
BLUE, ORANGE = "#1f77b4", "#ff7f0e"

# One entry per alignment analysed. "suffix" is appended to the model name in every
# file stem, "cond" gives the conditionOnCalibrations=false marker, and "prior" is the
# folder of sample-from-prior runs (None when that dataset has none).
DATASETS = {
    "nogapN": {
        "post": "data",
        "prior": "data",
        "xmls": "",
        "suffix": "",
        "cond": {"true": "", "false": "-condFalse"},
    },
    "codon": {
        "post": "codons",
        "prior": "codons",
        "xmls": "codons",
        "suffix": "_codon",
        "cond": {"true": "", "false": "_condFalse"},
    },
    "codon-nogapN": {
        "post": "codons",
        "prior": "codons",
        "xmls": "noGapNCodon",
        "suffix": "_codon_nogapN",
        "cond": {"true": "", "false": "_condFalse"},
    },
}

DATASET = os.environ.get("PRIMATES_DATASET", "nogapN")


def use(name):
    """Switch the dataset every path below refers to."""
    global DATASET, COND, PRIOR_DATA, POST_DATA, XMLS
    if name not in DATASETS:
        raise SystemExit("unknown dataset %r; known: %s" % (name, ", ".join(DATASETS)))
    DATASET = name
    d = DATASETS[name]
    COND = d["cond"]
    POST_DATA = os.path.join(BASE, d["post"])
    PRIOR_DATA = os.path.join(BASE, d["prior"]) if d["prior"] else None
    XMLS = os.path.join(XML_ROOT, d["xmls"])


def has_prior_runs():
    return DATASETS[DATASET]["prior"] is not None


def stem(model, cond, prior_only=False):
    return (model + DATASETS[DATASET]["suffix"] + COND[cond]
            + ("-fromPrior" if prior_only else ""))


def model_xml(model, cond="true"):
    return os.path.join(XMLS, stem(model, cond) + ".xml")


def log_path(model, cond, prior_only):
    if prior_only and not has_prior_runs():
        return ""
    return os.path.join(PRIOR_DATA if prior_only else POST_DATA,
                        stem(model, cond, prior_only) + ".txt")


def trees_path(model, cond, prior_only=False):
    if prior_only and not has_prior_runs():
        return ""
    return os.path.join(PRIOR_DATA if prior_only else POST_DATA,
                        stem(model, cond, prior_only) + ".trees")


def summary_tree_path(model, cond):
    return os.path.join(POST_DATA, stem(model, cond) + "_summary.tree")


def taxon_sets(model):
    """{TaxonSetN: frozenset(taxon names)} for one model XML."""
    root = ET.parse(model_xml(model)).getroot()
    out = {}
    for el in root.iter():
        i = el.get("id", "")
        if i.startswith("TaxonSet"):
            names = {t.get("idref") or t.get("id") for t in el.iter() if t.tag == "taxon"}
            if names:
                out[i] = frozenset(names)
    return out


_CLADE_RE = re.compile(r"^\s*([A-Za-z_]+) = \[(.*?)\];", re.S | re.M)


def clade_names():
    """{frozenset(taxa): clade name} from the LPhy headers of the calibrationPrior XMLs.

    Names for uncalibrated clades (Primates, Colobinae, ...) only appear in the
    nogapN XML, so both it and the current dataset's XML are read.
    """
    names = {}
    paths = [os.path.join(XML_ROOT, "calibrationPrior.xml"), model_xml("calibrationPrior")]
    for path in dict.fromkeys(paths):
        if not os.path.exists(path):
            continue
        with open(path) as fh:
            head = fh.read(20000)
        comment = head[head.find("<!--"):head.find("-->")]
        for name, body in _CLADE_RE.findall(comment):
            taxa = frozenset(re.findall(r'"([^"]+)"', body))
            if taxa:
                names[taxa] = name
    # the root clade is defined via D.getTaxaNames() in LPhy, so it has no literal list
    allt = taxon_sets("calibrationPrior")["TaxonSet"]
    names.setdefault(allt, "Euarchontoglires")
    return names


def calibrated_clades(model):
    """{frozenset(taxa): TaxonSetN} for the sets this model actually calibrates.

    TaxonSet (the alignment's own set) and TaxonSet1 both hold all 29 taxa; only
    the ones referenced by a calibration/MRCAPrior are returned.
    """
    root = ET.parse(model_xml(model)).getroot()
    sets = taxon_sets(model)
    used = set()
    for el in root.iter():
        spec = el.get("spec", "")
        ref = el.get("taxonset") or el.get("taxa")
        if ref and ref.startswith("@TaxonSet") and (
                "MRCAPrior" in spec or "CalibrationCladePrior" in spec):
            used.add(ref[1:])
    return {sets[i]: i for i in used if i in sets}


# ---------------------------------------------------------------- summary trees

_META = re.compile(r"\[&([^\]]*)\]")


def _split_top(s):
    depth = 0
    part = []
    for ch in s:
        if ch == "," and depth == 0:
            yield "".join(part)
            part = []
            continue
        if ch in "([{":
            depth += 1
        elif ch in ")]}":
            depth -= 1
        part.append(ch)
    yield "".join(part)


def _parse_meta(text):
    out = {}
    for key, val in re.findall(r"([\w%.]+)=(\{[^}]*\}|[^,]+)", text):
        if val.startswith("{"):
            out[key] = [float(x) for x in val[1:-1].split(",")]
        else:
            try:
                out[key] = float(val)
            except ValueError:
                out[key] = val
    return out


def parse_summary_tree(path):
    """[(frozenset(taxa), median_age, hpd_low, hpd_high)] for every internal node.

    Uses the CAheight_* annotations written by "treeannotator -height CA".
    """
    with open(path) as fh:
        text = fh.read()
    # the block runs up to the first ";"; the last entry has no trailing comma
    translate = dict(re.findall(r"(\d+)\s+([\w.]+)",
                                text.split("Translate")[1].split(";")[0]))
    newick = text[text.index("tree TREE"):]
    newick = newick[newick.index("=") + 1:].strip().rstrip(";")

    clades = []

    def walk(node):
        node = node.strip()
        if not node.startswith("("):
            label = node.split("[")[0].split(":")[0].strip()
            return frozenset([translate.get(label, label)])
        close = node.rindex(")")
        m = _META.search(node[close:])
        meta = m.group(1) if m else ""
        body = node[1:close]
        taxa = set()
        for child in _split_top(body):
            taxa |= walk(child)
        taxa = frozenset(taxa)
        d = _parse_meta(meta)
        if "CAheight_median" in d:
            hpd = d.get("CAheight_95%_HPD", [float("nan")] * 2)
            clades.append((taxa, d["CAheight_median"], hpd[0], hpd[1]))
        return taxa

    walk(newick)
    return clades


# ---------------------------------------------------------------------- logs

def read_log(path, burnin=0.1):
    """{column: [values]} from a BEAST trace file, post burn-in."""
    header, rows = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if header is None:
                header = fields
                continue
            try:
                rows.append([float(x) for x in fields])
            except ValueError:
                continue
    rows = rows[int(len(rows) * burnin):]
    return {name: [r[i] for r in rows] for i, name in enumerate(header) if rows and i < len(rows[0])}


# ------------------------------------------------- clade ages from a tree file

_STRIP_META = re.compile(r"\[[^\]]*\]")
_TOKEN = re.compile(r"[(),;]|[^(),;:]+(?::[-\d.eE+]+)?")


def _mrca_height(newick, targets):
    """Height of the MRCA of `targets` in an ultrametric newick string.

    Post-order walk carrying (number of targets below, node height). The first
    node holding all of them is the MRCA.
    """
    pos = 0
    n = len(newick)
    found = [None]

    def node():
        nonlocal pos
        if newick[pos] == "(":
            pos += 1
            count, height = 0, 0.0
            while True:
                c, h, blen = node()
                count += c
                height = max(height, h + blen)
                if newick[pos] == ",":
                    pos += 1
                    continue
                break
            pos += 1  # ')'
        else:
            start = pos
            while pos < n and newick[pos] not in "(),:":
                pos += 1
            label = newick[start:pos].strip()
            count = 1 if label in targets else 0
            height = 0.0
        blen = 0.0
        if pos < n and newick[pos] == ":":
            pos += 1
            start = pos
            while pos < n and newick[pos] not in "(),;":
                pos += 1
            blen = float(newick[start:pos])
        if found[0] is None and count == len(targets):
            found[0] = height
        return count, height, blen

    node()
    return found[0]


def clade_age_trace(trees_path, taxa, burnin=0.1, cache=True):
    """Posterior sample of the MRCA age of `taxa`, read straight from a .trees file.

    Needed because BEAST only logs mrca.age() for clades that carry a calibration,
    and the clades of interest here (e.g. crown Primates) are not always calibrated.
    Results are cached beside the tree file.
    """
    import hashlib
    import numpy as np

    stem = os.path.basename(trees_path)[:-len(".trees")]
    digest = hashlib.md5("\n".join(sorted(taxa)).encode()).hexdigest()[:8]
    key = "%s.%s" % (stem, digest)
    cache_file = os.path.join(os.path.dirname(trees_path), ".ages_%s.npy" % key)
    if cache and os.path.exists(cache_file) and os.path.getmtime(cache_file) > os.path.getmtime(trees_path):
        ages = np.load(cache_file)
    else:
        translate, ages = {}, []
        targets = None
        with open(trees_path) as fh:
            in_translate = False
            for line in fh:
                s = line.strip()
                if s.startswith("Translate"):
                    in_translate = True
                    continue
                if in_translate:
                    if s.startswith("tree ") or s == ";":
                        in_translate = False
                    else:
                        num, name = s.rstrip(",;").split()
                        translate[num] = name
                if s.startswith("tree "):
                    if targets is None:
                        inv = {v: k for k, v in translate.items()}
                        targets = {inv.get(t, t) for t in taxa}
                    nwk = _STRIP_META.sub("", s[s.index("=") + 1:]).strip().rstrip(";")
                    ages.append(_mrca_height(nwk, targets))
        ages = np.asarray(ages, dtype=float)
        if cache:
            np.save(cache_file, ages)
    return ages[int(len(ages) * burnin):]


use(DATASET)


def set_shared_xlim(axes, samples, tail=0.025, pad=0.03):
    """Give every panel the same x range, trimmed to the bulk of `samples`.

    The default keeps every 95% HPD in the figure fully visible and clips only
    what lies outside it; without trimming, one long tail (the suggested MRCA
    prior reaches ~260 Ma) squeezes all six panels into the left third.
    """
    import numpy as np

    if not len(samples):
        return
    lo = min(np.quantile(s, tail) for s in samples)
    hi = max(np.quantile(s, 1 - tail) for s in samples)
    span = hi - lo
    for row in axes:
        for ax in row:
            ax.set_xlim(lo - pad * span, hi + pad * span)
