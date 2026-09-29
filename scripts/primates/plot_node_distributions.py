#!/usr/bin/env python3
"""Posterior age of every internal node under the three calibration schemes.

One small panel per clade: the posterior MRCA-age density under each calibration
scheme (joint, de Vries and Beck, uniform), all with the same tree prior. Beneath
them, the calibration prior suggested for that node in de Vries and Beck (Table 1):
uniform, offset exponential (5% beyond the soft maximum) or minimum bound only;
filled if the calibration is used in these runs, outlined if not. Clades are those of
the summary tree, in taxonomic order.

Usage:
    python3 plot_node_distributions.py [--dataset=grid-codon-deepcal] [--cond=true]
        [--format=pdf]
"""

import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.stats import gaussian_kde

import primates_common as pc

SCHEMES = [("calibrationPrior", "Joint calibration prior", "#0072B2"),
           ("suggestedPrior", "de Vries & Beck priors", "#D55E00"),
           ("uniformPrior", "Uniform priors", "#009E73")]
FOSSIL = "#7a7a7a"
STRIP = (-0.40, 0.26)   # (baseline, height) of the calibration-prior strip, x peak height

# de Vries and Beck (2023) Table 1 for the nodes sampled here:
# clade -> (node number, shape, minimum, maximum). Shapes: "uniform", "offexp", "min".
# Their node 23 (Callitrichidae + Cebidae) is our sampled crown Platyrrhini.
PAPER_PRIORS = {
    "Euarchontoglires": (1, "offexp", 65.79, 125.816),
    "Euarchonta": (3, "offexp", 65.79, 125.816),
    "Primatomorpha": (4, "min", 55.935, None),
    "Primates": (5, "offexp", 55.935, 66.095),
    "Strepsirrhini": (6, "uniform", 36.573, 55.8),
    "Haplorhini": (10, "min", 41.0, None),
    "Anthropoidea": (11, "uniform", 33.9, 56.035),
    "Catarrhini": (12, "uniform", 25.193, 35.102),
    "Cercopithecidae": (13, "offexp", 12.47, 25.235),
    "Colobinae": (14, "uniform", 8.125, 15.0),
    "Cercopithecinae": (15, "uniform", 6.5, 15.0),
    "Papionini": (16, "uniform", 5.33, 12.51),
    "Hominoidea": (18, "offexp", 13.4, 25.235),
    "Hominidae": (19, "offexp", 12.3, 25.235),
    "Hominini": (20, "uniform", 4.631, 15.0),
    "Platyrrhini": (23, "uniform", 13.183, 34.5),
    "Cebidae": (24, "uniform", 13.032, 34.5),
}

_COLOBINI = {"Colobus_angolensis_palliatus", "Piliocolobus_tephrosceles"}
_RHINO = {"Rhinopithecus_bieti", "Rhinopithecus_roxellana"}
_MACACA = {"Macaca_mulatta", "Macaca_fascicularis", "Macaca_nemestrina"}
_PAPIO = {"Papio_anubis", "Theropithecus_gelada"}
_MANGABEY = {"Cercocebus_atys", "Mandrillus_leucophaeus"}
_PAN = {"Pan_paniscus", "Pan_troglodytes"}
_HOMININAE = _PAN | {"Homo_sapiens", "Gorilla_gorilla"}
_HOMINOIDEA = _HOMININAE | {"Pongo_abelii", "Nomascus_leucogenys"}
_CERCOPITHECIDAE = _COLOBINI | _RHINO | _MACACA | _PAPIO | _MANGABEY | {"Chlorocebus_sabaeus"}
_PLATYRRHINI = {"Cebus_capucinus_imitator", "Saimiri_boliviensis", "Aotus_nancymaae",
                "Callithrix_jacchus"}
_ANTHROPOIDEA = _HOMINOIDEA | _CERCOPITHECIDAE | _PLATYRRHINI
_STREPSIRRHINI = {"Otolemur_garnettii", "Microcebus_murinus", "Propithecus_coquereli"}
_PRIMATES = _ANTHROPOIDEA | _STREPSIRRHINI | {"Carlito_syrichta"}

# (clade, common name, taxa) in display order.
CLADES = [
    ("Euarchontoglires", "root (primates, treeshrew, colugo + mouse)",
     _PRIMATES | {"Galeopterus_variegatus", "Tupaia_chinensis", "Mus_musculus"}),
    ("Euarchonta", "primates, colugo + treeshrew",
     _PRIMATES | {"Galeopterus_variegatus", "Tupaia_chinensis"}),
    ("Primatomorpha", "primates + colugo", _PRIMATES | {"Galeopterus_variegatus"}),
    ("Primates", "crown primates", _PRIMATES),
    ("Strepsirrhini", "lemurs + galago", _STREPSIRRHINI),
    ("Lemuroidea", "mouse lemur + sifaka", {"Microcebus_murinus", "Propithecus_coquereli"}),
    ("Haplorhini", "tarsier + monkeys and apes", _ANTHROPOIDEA | {"Carlito_syrichta"}),
    ("Anthropoidea", "monkeys and apes", _ANTHROPOIDEA),
    ("Platyrrhini", "New World monkeys", _PLATYRRHINI),
    ("Cebidae", "capuchin + squirrel monkey", {"Cebus_capucinus_imitator", "Saimiri_boliviensis"}),
    ("Aotus + Callithrix", "night monkey + marmoset", {"Aotus_nancymaae", "Callithrix_jacchus"}),
    ("Catarrhini", "Old World monkeys and apes", _HOMINOIDEA | _CERCOPITHECIDAE),
    ("Hominoidea", "apes", _HOMINOIDEA),
    ("Hominidae", "great apes", _HOMININAE | {"Pongo_abelii"}),
    ("Homininae", "African great apes", _HOMININAE),
    ("Hominini", "human + chimpanzees", _PAN | {"Homo_sapiens"}),
    ("Pan", "chimpanzee + bonobo", _PAN),
    ("Cercopithecidae", "Old World monkeys", _CERCOPITHECIDAE),
    ("Colobinae", "leaf monkeys", _COLOBINI | _RHINO),
    ("Colobini", "African colobus monkeys", _COLOBINI),
    ("Rhinopithecus", "snub-nosed monkeys", _RHINO),
    ("Cercopithecinae", "cheek-pouched monkeys", _CERCOPITHECIDAE - _COLOBINI - _RHINO),
    ("Papionini", "macaques, baboons + kin", _MACACA | _PAPIO | _MANGABEY),
    ("Macaca", "macaques", _MACACA),
    ("M. mulatta + M. fascicularis", "rhesus + long-tailed macaque",
     {"Macaca_mulatta", "Macaca_fascicularis"}),
    ("Papionina", "baboons, gelada, mandrills, mangabeys", _PAPIO | _MANGABEY),
    ("Papio + Theropithecus", "baboon + gelada", _PAPIO),
    ("Cercocebus + Mandrillus", "mangabey + drill", _MANGABEY),
]
NCOLS = 4


def density(ax, ages, colour):
    """Draw the density; return (0.1%, 99.9% quantiles, peak height)."""
    grid = np.linspace(*np.quantile(ages, [0.001, 0.999]), 300)
    d = gaussian_kde(ages)(grid)
    ax.fill_between(grid, 0, d, color=colour, alpha=0.12, lw=0)
    ax.plot(grid, d, color=colour, lw=1.2)
    return grid[0], grid[-1], d.max()


def paper_prior(ax, prior, used, top, xmax):
    """Draw a de Vries and Beck calibration prior in the strip under the densities."""
    _, shape, lo, hi = prior
    y0, h = STRIP[0] * top, STRIP[1] * top
    style = dict(facecolor=FOSSIL if used else "white", edgecolor=FOSSIL, lw=1,
                 alpha=0.75 if used else 1)
    if shape == "uniform":
        ax.fill_between([lo, hi], y0, y0 + h, **style)
    elif shape == "offexp":
        x = np.linspace(lo, xmax, 200)
        mean = (hi - lo) / np.log(20)   # 5% beyond the soft maximum
        ax.fill_between(x, y0, y0 + h * np.exp(-(x - lo) / mean), **style)
        ax.plot([hi, hi], [y0, y0 + h * 0.8], color="#333333", lw=1, ls=(0, (2, 1.5)))
    else:
        ax.plot([lo, lo], [y0, y0 + h], color=FOSSIL, lw=1.5)
        ax.annotate("", xy=(xmax, y0 + h / 2), xytext=(lo, y0 + h / 2),
                    arrowprops=dict(arrowstyle="-|>", color=FOSSIL, lw=1.2,
                                    mutation_scale=8, shrinkA=0, shrinkB=0))


def main():
    opts = {"dataset": "grid-codon-deepcal", "cond": "true", "format": "pdf"}
    for a in sys.argv[1:]:
        k, v = a.lstrip("-").split("=", 1)
        opts[k] = v
    pc.use(opts["dataset"])
    cond = opts["cond"]

    named = {frozenset(t) for *_, t in CLADES}
    for model, *_ in SCHEMES:
        unnamed = {t for t, *_ in pc.parse_summary_tree(pc.summary_tree_path(model, cond))} - named
        if unnamed:
            sys.exit("%s summary tree has unnamed clades: %s"
                     % (model, [sorted(t) for t in unnamed]))

    used = set(pc.calibrated_clades("calibrationPrior"))

    plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans", "axes.linewidth": 0.6,
                         "pdf.fonttype": 42})
    nrows = -(-len(CLADES) // NCOLS)
    fig, axes = plt.subplots(nrows, NCOLS, figsize=(10.5, 1.9 * nrows + 0.7))

    for ax, (clade, common, taxa) in zip(axes.flat, CLADES):
        taxa = frozenset(taxa)
        spans, top = [], 0.0
        for model, _, colour in SCHEMES:
            post = pc.clade_age_trace(pc.trees_path(model, cond), taxa, burnin=0.0)
            lo, hi, peak = density(ax, post, colour)
            spans += [lo, hi]
            top = max(top, peak)
        prior = PAPER_PRIORS.get(clade)
        if prior:
            spans += [prior[2]] + ([prior[3]] if prior[3] else [])
        lo, hi = min(spans), max(spans)
        pad = 0.06 * (hi - lo)
        xlim = (max(0, lo - pad), hi + pad)
        if prior:
            paper_prior(ax, prior, taxa in used, top, xlim[1])
        ax.set_xlim(*xlim)
        ax.set_ylim((STRIP[0] - 0.05) * top, 1.08 * top)
        ax.set_title(clade + (" (#%d)" % prior[0] if prior else ""), fontsize=9,
                     fontweight="bold", loc="left", pad=11)
        ax.text(0, 1.02, common, transform=ax.transAxes, fontsize=7, color="#666666",
                va="bottom")
        ax.set_yticks([])
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        ax.tick_params(axis="x", labelsize=7, length=2)

    for ax in axes.flat[len(CLADES):]:
        ax.axis("off")

    handles = [Line2D([], [], color=c, lw=1.5, label="Posterior, " + lab[0].lower() + lab[1:]
                      if lab.startswith(("Joint", "Uniform")) else "Posterior, " + lab)
               for _, lab, c in SCHEMES]
    handles += [Patch(facecolor=FOSSIL, alpha=0.75, edgecolor=FOSSIL),
                Patch(facecolor="white", edgecolor=FOSSIL),
                Line2D([], [], color="#333333", lw=1, ls=(0, (2, 1.5)))]
    labels = [h.get_label() for h in handles[:3]] + [
        "Calibration prior (de Vries & Beck), used", "Same, not used in these runs",
        "Soft maximum (offset exponential)"]
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False, fontsize=8.5,
               bbox_to_anchor=(0.5, 1.0))
    fig.supxlabel("Node age (Ma)", fontsize=9)
    fig.tight_layout(rect=(0, 0, 0.965, 0.95), h_pad=1.6, w_pad=2.5)

    out = pc.fig_path("%s_node_distributions_%s.%s"
                      % (pc.DATASET, "cond" if cond == "true" else "regular", opts["format"]))
    fig.savefig(out, dpi=250)
    print("wrote", out)


if __name__ == "__main__":
    main()
