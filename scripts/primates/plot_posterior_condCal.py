#!/usr/bin/env python3
"""Age of one clade under conditionOnCalibrations=true vs false.

2 rows (sample from prior / with data) x 3 columns (joint, de Vries and Beck, uniform
calibration priors). Each panel overlays the calibrated and the regular tree prior as
filled KDEs, with the K-Pg boundary marked.

Ages are read from the .trees files rather than the traces: BEAST only logs mrca.age()
for calibrated clades, and crown Primates is not calibrated in these runs.

Self-contained: everything it needs is in this file.

Usage:
    python3 plot_posterior_condCal.py                     # nogapN, Primates
    python3 plot_posterior_condCal.py Colobinae
    python3 plot_posterior_condCal.py --dataset=codon Primates
"""

import hashlib
import os
import re
import sys
import xml.etree.ElementTree as ET

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde

BASE = os.path.dirname(os.path.abspath(__file__))
XML_ROOT = os.path.join(BASE, "..", "..", "calibratedcpp-beast", "src", "test",
                        "resources", "calibratedcpp", "examples", "primates", "xmls")

MODELS = ["calibrationPrior", "suggestedPrior", "uniformPrior"]
MODEL_LABELS = {
    "calibrationPrior": "Joint (Beta-LogNormal) calibration prior",
    "suggestedPrior": "de Vries and Beck independent calibration priors",
    "uniformPrior": "Uniform independent calibration priors",
}
BLUE, ORANGE = "#1f77b4", "#ff7f0e"
BOUNDARY = 66.0
COLOURS = {"true": BLUE, "false": ORANGE}
COND_LABELS = {"true": "Calibrated tree prior", "false": "Regular tree prior"}
ROWS = [("sample from prior", True), ("with data", False)]

# One entry per alignment analysed. "suffix" is appended to the model name in every
# file stem, "cond" gives the conditionOnCalibrations=false marker, and "prior" is the
# folder of sample-from-prior runs (None when that dataset has none).
DATASETS = {
    "nogapN": {"post": "data", "prior": "data", "xmls": "", "suffix": "",
               "cond": {"true": "", "false": "-condFalse"}},
    "codon": {"post": "codons", "prior": "codons", "xmls": "codons", "suffix": "_codon",
              "cond": {"true": "", "false": "_condFalse"}},
    "codon-nogapN": {"post": "codons-nogapN", "prior": "codons-nogapN",
                     "xmls": "noGapNCodon", "suffix": "_codon_nogapN",
                     "cond": {"true": "", "false": "_condFalse"}},
}
DATASET = os.environ.get("PRIMATES_DATASET", "nogapN")


def use(name):
    """Switch the dataset every path below refers to."""
    global DATASET, D, POST_DATA, PRIOR_DATA, XMLS
    if name not in DATASETS:
        raise SystemExit("unknown dataset %r; known: %s" % (name, ", ".join(DATASETS)))
    DATASET = name
    D = DATASETS[name]
    POST_DATA = os.path.join(BASE, D["post"])
    PRIOR_DATA = os.path.join(BASE, D["prior"]) if D["prior"] else None
    XMLS = os.path.join(XML_ROOT, D["xmls"])


def has_prior_runs():
    return D["prior"] is not None


def stem(model, cond, prior_only=False):
    return (model + D["suffix"] + D["cond"][cond] + ("-fromPrior" if prior_only else ""))


def trees_path(model, cond, prior_only=False):
    if prior_only and not has_prior_runs():
        return ""
    return os.path.join(PRIOR_DATA if prior_only else POST_DATA,
                        stem(model, cond, prior_only) + ".trees")


# ------------------------------------------------------------- clade names

_CLADE_RE = re.compile(r"^\s*([A-Za-z_]+) = \[(.*?)\];", re.S | re.M)


def clade_names():
    """{frozenset(taxa): clade name} from the LPhy headers of the calibrationPrior XMLs.

    Names for uncalibrated clades (Primates, Colobinae, ...) only appear in the
    nogapN XML, so both it and the current dataset's XML are read.
    """
    names = {}
    own = os.path.join(XMLS, stem("calibrationPrior", "true") + ".xml")
    for path in dict.fromkeys([os.path.join(XML_ROOT, "calibrationPrior.xml"), own]):
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
    root = ET.parse(own).getroot()
    for el in root.iter():
        if el.get("id") == "TaxonSet":
            allt = frozenset(t.get("idref") or t.get("id")
                             for t in el.iter() if t.tag == "taxon")
            names.setdefault(allt, "Euarchontoglires")
    return names


# ------------------------------------------- clade ages from a .trees file

_STRIP_META = re.compile(r"\[[^\]]*\]")


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


def clade_age_trace(path, taxa, burnin=0.1, cache=True):
    """Posterior sample of the MRCA age of `taxa`, cached beside the tree file."""
    name = os.path.basename(path)[:-len(".trees")]
    digest = hashlib.md5("\n".join(sorted(taxa)).encode()).hexdigest()[:8]
    cache_file = os.path.join(os.path.dirname(path), ".ages_%s.%s.npy" % (name, digest))
    if cache and os.path.exists(cache_file) and os.path.getmtime(cache_file) > os.path.getmtime(path):
        ages = np.load(cache_file)
    else:
        translate, ages, targets = {}, [], None
        with open(path) as fh:
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
                        num, taxon = s.rstrip(",;").split()
                        translate[num] = taxon
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


def set_shared_xlim(axes, samples, tail=0.025, pad=0.03):
    """Give every panel the same x range, trimmed to the bulk of `samples`.

    The default keeps every 95% HPD in the figure fully visible and clips only
    what lies outside it; without trimming, one long tail (the suggested MRCA
    prior reaches ~260 Ma) squeezes all six panels into the left third.
    """
    if not len(samples):
        return
    lo = min(np.quantile(s, tail) for s in samples)
    hi = max(np.quantile(s, 1 - tail) for s in samples)
    span = hi - lo
    for row in axes:
        for ax in row:
            ax.set_xlim(lo - pad * span, hi + pad * span)


# ------------------------------------------------------------------ figure

def panel(ax, model, prior_only, taxa, clade):
    """Draw one panel; return the age samples it plotted, for shared axis limits."""
    drawn = []
    for cond in ("true", "false"):
        path = trees_path(model, cond, prior_only)
        if not os.path.exists(path):
            continue
        ages = clade_age_trace(path, taxa)
        if len(ages) < 20 or ages.std() < 1e-9:
            continue
        grid = np.linspace(ages.min(), ages.max(), 500)
        dens = gaussian_kde(ages)(grid)
        ax.plot(grid, dens, color=COLOURS[cond], lw=1.6, label=COND_LABELS[cond])
        ax.fill_between(grid, dens, color=COLOURS[cond], alpha=0.25)
        drawn.append(ages)

    if not drawn:
        ax.text(0.5, 0.5, "no runs yet", ha="center", va="center",
                transform=ax.transAxes, color="grey")
        return drawn

    ax.axvline(BOUNDARY, color="k", ls="--", lw=1.0,
               label="K-Pg boundary (%g Ma)" % BOUNDARY)
    ax.legend(fontsize=8, loc="upper right", framealpha=0.9)
    ax.set_xlabel("%s age (Ma)" % clade)
    return drawn


def main():
    if len(sys.argv) > 1 and sys.argv[1].startswith("--dataset="):
        use(sys.argv.pop(1).split("=", 1)[1])
    clade = sys.argv[1] if len(sys.argv) > 1 else "Primates"
    names = clade_names()
    match = [t for t, n in names.items() if n == clade]
    if not match:
        sys.exit("unknown clade %r; known: %s" % (clade, ", ".join(sorted(names.values()))))
    taxa = match[0]

    rows = [r for r in ROWS if has_prior_runs() or not r[1]]
    fig, axes = plt.subplots(len(rows), len(MODELS), squeeze=False, sharex=True,
                             figsize=(6.0 * len(MODELS), 4.6 * len(rows)))
    drawn = []
    for r, (row_label, prior_only) in enumerate(rows):
        for c, model in enumerate(MODELS):
            ax = axes[r][c]
            drawn += panel(ax, model, prior_only, taxa, clade)
            ax.set_title("%s - %s" % (MODEL_LABELS[model], row_label), fontsize=12)
        axes[r][0].set_ylabel("posterior density" if not prior_only else "prior density")

    set_shared_xlim(axes, drawn)

    out = os.path.join(BASE, "%s_condCal_%s.png" % (DATASET, clade))
    fig.tight_layout()
    fig.savefig(out, dpi=180)
    print("wrote", out)


use(DATASET)

if __name__ == "__main__":
    main()
