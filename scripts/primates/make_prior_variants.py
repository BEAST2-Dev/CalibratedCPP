#!/usr/bin/env python3
"""
Expand each calibrationPrior XML into the full grid of primates analyses.

One XML per alignment goes in; six come out, covering the three calibration schemes and
both settings of conditionOnCalibrations:

    <stem>.xml                  <stem>_condFalse.xml
    <stem>_uniform.xml          <stem>_uniform_condFalse.xml
    <stem>_suggested.xml        <stem>_suggested_condFalse.xml

--from-prior additionally writes a <name>-fromPrior.xml for each, with
sampleFromPrior="true" on the MCMC element, to check what the calibrations imply before the
data is allowed to speak.

--replicates N writes each of the six as N identical copies named <name>-rep1..N.xml, for
running the same analysis as independent chains (BEAST's -seed differentiates them). A
"-repN" suffix already on the input is likewise kept at the end of every output name.

Everything except the calibration block is copied verbatim, so the alignment, substitution
model, clock and tree prior are identical across all six -- the only thing that varies is
what is being compared.

The schemes:

    calibrationPrior   one joint calibrationprior.CalibrationPrior over all clades
    uniform            one MRCAPrior(Uniform(lower, upper)) per clade
    suggested          as uniform, except the clades de Vries and Beck give an
                       offset-exponential shape to, which become
                       MRCAPrior(OffsetReal(offset, Exponential(mean)))

The exponential mean follows the paper's convention of a 5% probability of the divergence
being older than the soft maximum, so for age = offset + Exp(mean):

    mean = (upper - offset) / -ln(0.05)

matching compute_offset_exponential_mean.py. Which clades take that shape is decided by
their (lower, upper) bounds, which identify the paper's nodes uniquely in this dataset.

The starting tree is checked against the calibration bounds before anything is written. The
joint CalibrationPrior is a smooth density and tolerates a clade outside its interval, but
MRCAPrior(Uniform(lower, upper)) returns -Infinity there, so BEAST would fail to initialise
with "Could not find a proper state". The LPhy simulator does emit such trees, so convert
with another seed if this reports a violation.

Log file names are rewritten to BEAST's $(filebase) so every output writes its traces and
trees under its own name rather than the name it was copied from.

A partitioned XML also gets clock.rate wired back into its branch-rate model. The LPhy
scripts all pass mu=mRate to PhyloCTMC, but LPhyBEAST drops it on the multi-partition path,
leaving UCRelaxedClockModel on its default rate of 1.0 and mRate an orphan: sampled and
operated on, but read by nothing. branchRates then has to carry the absolute rate while its
own prior holds it at mean 1, so the clock is fought from both ends. Remove this once the
conversion is fixed upstream.

Usage:
    python make_prior_variants.py IN.xml [IN2.xml ...] [-o OUTDIR] [--replicates N]
                                  [--from-prior]
"""

import argparse
import math
import os
import re

TAIL_PROBABILITY = 0.05  # P(age > upper), per de Vries and Beck

# (lower, upper) of the calibrations the paper gives an offset-exponential shape to.
# Keyed by bounds because the XML carries no clade names -- these pairs are unique here.
OFFSET_EXPONENTIAL_BOUNDS = {
    (65.79, 125.816),   # Euarchontoglires / Euarchonta
    (55.935, 66.095),   # Primates
    (12.47, 25.235),    # Cercopithecidae
    (13.4, 25.235),     # Hominoidea
    (12.3, 25.235),     # Hominidae
}

CALIBRATION_PRIOR = re.compile(
    r'(?P<indent>[ \t]*)<distribution id="CalibrationPrior"[^>]*>.*?</distribution>\n',
    re.S,
)
CLADE = re.compile(
    r'<calibration\b[^>]*taxa="@(?P<taxonset>\w+)".*?'
    r'<upperAge\b[^>]*value="(?P<upper>[\d.eE+-]+)".*?'
    r'<lowerAge\b[^>]*value="(?P<lower>[\d.eE+-]+)"',
    re.S,
)
MODEL = re.compile(r'<distribution id="[^"]*" spec="calibratedcpp.CalibratedBirthDeathSkylineModel"')


def offset_exponential_mean(offset, upper, tail_probability=TAIL_PROBABILITY):
    """The Exp mean putting `tail_probability` of the mass above `upper`."""
    return (upper - offset) / -math.log(tail_probability)


def _num(x):
    """Render a float without a trailing '.0' clutter, as the generated XMLs do."""
    return repr(round(x, 6)).rstrip("0").rstrip(".") if isinstance(x, float) else str(x)


def _mrca_prior(idx, taxonset, lower, upper, scheme, indent):
    """One MRCAPrior element, uniform or offset-exponential depending on the bounds."""
    sfx = "" if idx == 0 else str(idx)
    pad = indent + "    "
    head = (f'{indent}<distribution id="MRCAPrior{sfx}" '
            f'spec="beast.base.spec.evolution.tree.MRCAPrior" monophyletic="true" '
            f'taxonset="@{taxonset}" tree="@tree">\n')

    if scheme == "suggested" and (lower, upper) in OFFSET_EXPONENTIAL_BOUNDS:
        mean = offset_exponential_mean(lower, upper)
        body = (
            f'{pad}<distr id="OffsetReal{sfx}" spec="beast.base.spec.inference.distribution.OffsetReal">\n'
            f'{pad}    <distribution id="Exponential{sfx}" spec="beast.base.spec.inference.distribution.Exponential">\n'
            f'{pad}        <mean id="ExpMean{sfx}" spec="beast.base.spec.inference.parameter.RealScalarParam" '
            f'domain="PositiveReal" value="{_num(mean)}"/>\n'
            f'{pad}    </distribution>\n'
            f'{pad}    <offset id="ExpOffset{sfx}" spec="beast.base.spec.inference.parameter.RealScalarParam" '
            f'domain="NonNegativeReal" value="{_num(lower)}"/>\n'
            f'{pad}</distr>\n'
        )
    else:
        body = (
            f'{pad}<distr id="Uniform{sfx}" spec="beast.base.spec.inference.distribution.Uniform">\n'
            f'{pad}    <lower id="UniformLower{sfx}" spec="beast.base.spec.inference.parameter.RealScalarParam" '
            f'domain="NonNegativeReal" value="{_num(lower)}"/>\n'
            f'{pad}    <upper id="UniformUpper{sfx}" spec="beast.base.spec.inference.parameter.RealScalarParam" '
            f'domain="NonNegativeReal" value="{_num(upper)}"/>\n'
            f'{pad}</distr>\n'
        )
    return head + body + f'{indent}</distribution>\n'


def _clade_ages(newick):
    """{frozenset(taxa): age} for every node of a newick string with branch lengths."""
    pos = 0

    def node():
        nonlocal pos
        if newick[pos] == "(":
            pos += 1
            kids = []
            while True:
                kids.append(node())
                if newick[pos] == ",":
                    pos += 1
                else:
                    break
            pos += 1                      # closing bracket
        else:
            name = re.match(r"[^:,()]+", newick[pos:]).group(0)
            pos += len(name)
            kids = []
        m = re.match(r":([\d.eE+-]+)", newick[pos:])
        length = float(m.group(1)) if m else 0.0
        if m:
            pos += len(m.group(0))
        taxa = frozenset([name]) if not kids else frozenset().union(*(k[0] for k in kids))
        age = 0.0 if not kids else max(k[1] + k[2] for k in kids)
        ages[taxa] = age
        return taxa, age, length

    ages = {}
    node()
    return ages


def check_starting_tree(xml, path):
    """Warn if the starting tree breaks a calibration, which would give -Infinity later."""
    tree = re.search(r'newick="([^"]+)"', xml)
    if not tree:
        return
    ages = _clade_ages(tree.group(1).rstrip(";"))
    for m in re.finditer(r'<calibration\b[^>]*taxa="@(\w+)".*?'
                         r'<upperAge\b[^>]*value="([\d.eE+-]+)".*?'
                         r'<lowerAge\b[^>]*value="([\d.eE+-]+)"', xml, re.S):
        taxonset, upper, lower = m.group(1), float(m.group(2)), float(m.group(3))
        block = re.search(r'<calibrations id="%s"[^>]*>(.*?)</calibrations>' % taxonset,
                          xml, re.S)
        if not block:
            continue
        clade = frozenset(re.findall(r'idref="([^"]+)"', block.group(1)))
        age = ages.get(clade)
        if age is None:
            print(f"  WARNING {os.path.basename(path)}: {taxonset} is not monophyletic in the "
                  f"starting tree; MRCAPrior(monophyletic) will reject it")
        elif not lower <= age <= upper:
            print(f"  WARNING {os.path.basename(path)}: {taxonset} starts at {age:.3f}, outside "
                  f"[{lower}, {upper}]; the swapped MRCAPriors will give -Infinity. Reconvert "
                  f"the LPhy script with another seed.")


def swap_calibrations(xml, scheme):
    """Replace the joint CalibrationPrior with one MRCAPrior per clade."""
    m = CALIBRATION_PRIOR.search(xml)
    if not m:
        raise SystemExit("no <distribution id=\"CalibrationPrior\"> block found")

    clades = [(c.group("taxonset"), float(c.group("lower")), float(c.group("upper")))
              for c in CLADE.finditer(m.group(0))]
    if not clades:
        raise SystemExit("CalibrationPrior block contains no calibrations")

    indent = m.group("indent")
    priors = "".join(_mrca_prior(i, ts, lo, up, scheme, indent)
                     for i, (ts, lo, up) in enumerate(clades))
    xml = xml[:m.start()] + priors + xml[m.end():]

    # the tracelog logged the joint prior; log each MRCAPrior instead
    refs = "".join(f'{indent}    <log idref="MRCAPrior{"" if i == 0 else i}"/>\n'
                   for i in range(len(clades)))
    xml, n = re.subn(r'[ \t]*<log idref="CalibrationPrior"/>\n', refs, xml)
    if n != 1:
        raise SystemExit(f'expected one <log idref="CalibrationPrior"/>, found {n}')
    return xml


def set_condition(xml, condition):
    """Set conditionOnCalibrations on the tree prior (the attribute defaults to true)."""
    xml = re.sub(r'\s*conditionOnCalibrations="[a-z]*"', "", xml, count=1)
    if condition:
        return xml
    m = MODEL.search(xml)
    if not m:
        raise SystemExit("CalibratedBirthDeathSkylineModel element not found")
    return xml[:m.end()] + ' conditionOnCalibrations="false"' + xml[m.end():]


BRANCH_RATE_MODEL = re.compile(
    r'<branchRateModel\b(?![^>]*\bclock\.rate=)(?P<attrs>[^>]*'
    r'spec="beast\.base\.spec\.evolution\.branchratemodel\.UCRelaxedClockModel"[^>]*)>')


def wire_clock_rate(xml, rate="mRate"):
    """Give the relaxed clock its overall rate back, if the conversion left it out."""
    if f'id="{rate}"' not in xml:
        return xml
    return BRANCH_RATE_MODEL.sub(
        lambda m: f'<branchRateModel clock.rate="@{rate}"{m.group("attrs")}>', xml, count=1)


def sample_from_prior(xml):
    """Ignore the likelihood, so the chain samples what the priors alone imply."""
    xml = re.sub(r'\s*sampleFromPrior="[a-z]*"', "", xml, count=1)
    return re.sub(r'(<run\b[^>]*\bspec="MCMC")', r'\1 sampleFromPrior="true"', xml, count=1)


def use_filebase(xml):
    """Point every logger at $(filebase), so each copy writes under its own name."""
    return re.sub(r'fileName="[^"]*\.(log|trees)"', r'fileName="$(filebase).\1"', xml)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+", help="calibrationPrior XML, one per alignment")
    ap.add_argument("-o", "--outdir", help="where to write (default: beside each input)")
    ap.add_argument("--replicates", type=int, default=0, metavar="N",
                    help="also write each output as N copies named -rep1..N")
    ap.add_argument("--from-prior", action="store_true",
                    help="also write a -fromPrior.xml of each, sampling the prior only")
    args = ap.parse_args()

    for path in args.inputs:
        source = wire_clock_rate(use_filebase(open(path).read()))
        check_starting_tree(source, path)
        outdir = args.outdir or os.path.dirname(path) or "."
        os.makedirs(outdir, exist_ok=True)
        stem = os.path.splitext(os.path.basename(path))[0]
        # a replicate marker stays last, after the scheme and conditioning suffixes
        stem, rep = re.match(r"(.*?)(-rep\d+)?$", stem).groups()
        rep = rep or ""

        for scheme in ("calibration", "uniform", "suggested"):
            xml = source if scheme == "calibration" else swap_calibrations(source, scheme)
            for condition in (True, False):
                name = stem
                if scheme != "calibration":
                    name += "_" + scheme
                if not condition:
                    name += "_condFalse"
                body = set_condition(xml, condition)
                suffixes = ([rep] if not args.replicates
                            else [f"-rep{k}" for k in range(1, args.replicates + 1)])
                for suffix in suffixes:
                    out = os.path.join(outdir, name + suffix + ".xml")
                    with open(out, "w") as fh:
                        fh.write(body)
                    print("wrote", out)
                if args.from_prior:
                    # one per configuration: the prior does not need replicate chains
                    out = os.path.join(outdir, name + rep + "-fromPrior.xml")
                    with open(out, "w") as fh:
                        fh.write(sample_from_prior(body))
                    print("wrote", out)


if __name__ == "__main__":
    main()
