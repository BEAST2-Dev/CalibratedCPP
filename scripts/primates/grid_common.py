"""Shared pieces for the grid figures (compute_grid_stats.py and plot_grid_*.py).

The grid is the 3 alignments x 3 calibration schemes x 2 tree priors design:

    primates_<alignment>[_<scheme>][_condFalse]-{combined,fromPrior}

Everything here is pure Python so that compute_grid_stats.py can run on a machine
without matplotlib (e.g. the cluster where the .trees files live).
"""

import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
GRID = os.path.join(HERE, "grid")
GRID_DATA = os.path.join(GRID, "data")        # LogCombiner output + summary trees
FIG_DIR = os.path.join(GRID, "figures")
XML_DIR = os.path.join(GRID, "xmls")   # copy of the grid XMLs; falls back to the BEAST repo
if not os.path.isdir(XML_DIR):
    XML_DIR = os.path.join(HERE, "..", "..", "calibratedcpp-beast", "src", "test", "resources",
                           "calibratedcpp", "examples", "primates", "xmls", "primatesGrid")

TOPO_JSON = os.path.join(GRID, "topology_support.json")
MRCA_CSV = os.path.join(GRID, "mrca_age_comparison_grid.csv")

# ------------------------------------------------------------------ the design
ALIGNMENTS = [
    ("codon", "Codon-partitioned, full alignment"),
    ("codon_noGapN", "Codon-partitioned, gap/N codons removed"),
    ("unpartition", "Unpartitioned, gap/N columns removed"),
]
SCHEMES = [                      # (stem suffix, label, colour)
    ("", "Joint calibration prior", "#0072B2"),
    ("_suggested", "Independent priors (de Vries & Beck)", "#D55E00"),
    ("_uniform", "Independent priors (uniform)", "#009E73"),
]
SCHEME_SHORT = {"": "Joint prior", "_suggested": "de Vries & Beck", "_uniform": "Uniform"}
TREE_PRIORS = [                  # (stem suffix, label, colour)
    ("", "Tree prior conditioned on calibrations", "#0072B2"),
    ("_condFalse", "Unconditioned tree prior", "#D55E00"),
]
DEFAULT_ALIGNMENT = "codon_noGapN"
BURNIN = 0.1

OUTGROUP = {"Mus_musculus", "Tupaia_chinensis", "Galeopterus_variegatus"}

CALIBRATED = [   # display order, top to bottom
    ("Euarchontoglires", "Euarchontoglires"),
    ("Callitrichidae_Cebidae", "Callitrichidae + Cebidae"),
    ("Hominoidea", "Hominoidea"),
    ("Hominidae", "Hominidae"),
    ("Cercopithecidae", "Cercopithecidae"),
    ("Cercopithecinae", "Cercopithecinae"),
    ("Papionini", "Papionini"),
    ("Homo_Pan", "Homo–Pan"),
]


def stem(alignment, scheme="", cond=""):
    return "primates_%s%s%s" % (alignment, scheme, cond)


def all_stems():
    for aln, _ in ALIGNMENTS:
        for sch, _, _ in SCHEMES:
            for cond, _, _ in TREE_PRIORS:
                yield stem(aln, sch, cond)



def calibration_bounds(xml_path=None):
    """{clade: (lower, upper)} from the LPhy header of a grid XML."""
    xml_path = xml_path or os.path.join(XML_DIR, stem(DEFAULT_ALIGNMENT) + "-rep1.xml")
    with open(xml_path) as fh:
        head = fh.read(20000)
    out = {}
    for name, hi, lo in re.findall(r"calibration\(taxa=(\w+), upper=([\d.]+), lower=([\d.]+)\)", head):
        out[name] = (float(lo), float(hi))
    return out


# ------------------------------------------------------------ reading .trees
def read_trees(path, burnin=BURNIN):
    """(translate, [newick]) from a BEAST .trees file, burn-in removed."""
    translate, trees, in_tr = {}, [], False
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if s.lower().startswith("translate"):
                in_tr = True
                continue
            if in_tr:
                if s == ";":
                    in_tr = False
                    continue
                num, name = s.rstrip(",;").split()[:2]
                translate[num] = name
                continue
            if s.startswith("tree "):
                trees.append(s[s.index("("):].rstrip(";"))
    return translate, trees[int(len(trees) * burnin):]


def _skip_meta(s, i):
    if i < len(s) and s[i] == "[":
        return s.index("]", i) + 1
    return i


_NUM = re.compile(r"[0-9.eE+-]+")


def walk_newick(s, translate, on_clade):
    """Depth-first walk; calls on_clade(taxa_set, height) for every internal node.

    Heights are measured from the tips (tree assumed ultrametric, as BEAST writes).
    """
    i = 0

    def node():
        nonlocal i
        if s[i] == "(":
            i += 1
            taxa, height = set(), 0.0
            while True:
                t, h, bl = node()
                taxa |= t
                height = max(height, h + bl)
                if s[i] == ",":
                    i += 1
                    continue
                i += 1  # ')'
                break
            while i < len(s) and s[i] not in ":,);":
                i = _skip_meta(s, i) if s[i] == "[" else i + 1
            bl = 0.0
            if i < len(s) and s[i] == ":":
                i = _skip_meta(s, i + 1)
                m = _NUM.match(s, i)
                bl = float(m.group())
                i = m.end()
            on_clade(taxa, height)
            return taxa, height, bl
        m = re.compile(r"[^:,)\[]+").match(s, i)
        name = m.group()
        i = _skip_meta(s, m.end())
        bl = 0.0
        if s[i] == ":":
            i = _skip_meta(s, i + 1)
            m = _NUM.match(s, i)
            bl = float(m.group())
            i = m.end()
        return {translate.get(name, name)}, 0.0, bl

    node()


def load_json(path):
    with open(path) as fh:
        return json.load(fh)


def figure_style():
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans", "axes.linewidth": 0.6,
                         "xtick.major.width": 0.6, "ytick.major.width": 0.6,
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    os.makedirs(FIG_DIR, exist_ok=True)
    return plt


def save(fig, name):
    out = os.path.join(FIG_DIR, name + ".pdf")
    fig.savefig(out, bbox_inches="tight")
    print("wrote", out)
