#!/usr/bin/env python3
"""
Rewrite a joint-calibration-prior LPhy script into the independent-prior schemes.

The three primates schemes differ only in how the eight clade calibrations enter the
model; everything else (data, substitution model, clock, tree prior) is identical:

    calibrationPrior  cX = calibration(taxa=X, upper=U, lower=L)
                      calibrations ~ ConditionedMRCAPrior(calibrations=[...])

    uniformPrior      cX ~ UniformMRCA(taxa=X, upper=U, lower=L)
                      calibrations = toArray(calibrations=[...])

    suggestedPrior    as uniformPrior, except the clades de Vries and Beck give an
                      offset-exponential shape to:
                      cX ~ OffsetExponentialMRCA(taxa=X, offset=L, mean=M)

The offset-exponential mean is not free: the paper specifies a 5% probability of the
divergence being older than the soft maximum, so for age = offset + Exp(mean),

    mean = (upper - offset) / -ln(0.05)

which is what compute_offset_exponential_mean.py prints. This script derives it the same
way rather than hard-coding the numbers.

Usage:
    python swap_calibration_prior.py IN.lphy --scheme uniform    -o OUT.lphy
    python swap_calibration_prior.py IN.lphy --scheme suggested  -o OUT.lphy
    python swap_calibration_prior.py IN.lphy --scheme both       # writes both, named
                                                                 # after the input folder

Note this preserves the clades of the input script. The checked-in uniformPrior and
suggestedPrior scripts calibrate a slightly different set (Colobinae and Cebidae in place
of Cercopithecidae and Callitrichidae_Cebidae), so output will not match them byte for
byte; it is the same model with the input's clade choices.
"""

import argparse
import math
import os
import re

TAIL_PROBABILITY = 0.05  # P(age > upper), per de Vries and Beck

# Clades the paper gives an offset-exponential shape to, by the name used for the
# calibration variable (cEuarchontoglires -> "Euarchontoglires"). Any clade not listed
# keeps its uniform bounds in the suggested scheme.
OFFSET_EXPONENTIAL_CLADES = {
    "Euarchontoglires",
    "Euarchonta",
    "Primates",
    "Cercopithecidae",
    "Hominoidea",
    "Hominidae",
}

CALIBRATION = re.compile(
    r"^(?P<indent>\s*)c(?P<clade>\w+)\s*=\s*calibration\("
    r"\s*taxa\s*=\s*(?P<taxa>\w+)\s*,"
    r"\s*upper\s*=\s*(?P<upper>[\d.eE+-]+)\s*,"
    r"\s*lower\s*=\s*(?P<lower>[\d.eE+-]+)\s*\)\s*;\s*$",
    re.MULTILINE,
)

CONDITIONED = re.compile(
    r"^(?P<indent>\s*)(?P<name>\w+)\s*~\s*ConditionedMRCAPrior\(calibrations\s*=\s*"
    r"(?P<list>\[[^\]]*\])\s*\)\s*;\s*$",
    re.MULTILINE,
)


def offset_exponential_mean(offset, upper, tail_probability=TAIL_PROBABILITY):
    """The Exp mean putting `tail_probability` of the mass above `upper`."""
    return (upper - offset) / -math.log(tail_probability)


def _fmt(x):
    """Print a float the way the LPhy scripts do: no trailing zeros, but keep one dp."""
    s = f"{x:.6f}".rstrip("0")
    return s + "0" if s.endswith(".") else s


def swap(source, scheme):
    """Return `source` with its calibration block rewritten into `scheme`."""
    if not CALIBRATION.search(source):
        raise SystemExit("no `cX = calibration(...)` lines found; is this a "
                         "calibrationPrior script?")

    def replace(m):
        indent, clade = m.group("indent"), m.group("clade")
        taxa = m.group("taxa")
        upper, lower = float(m.group("upper")), float(m.group("lower"))

        if scheme == "suggested" and clade in OFFSET_EXPONENTIAL_CLADES:
            mean = offset_exponential_mean(lower, upper)
            return (f"{indent}c{clade} ~ OffsetExponentialMRCA(taxa={taxa}, "
                    f"offset={_fmt(lower)}, mean={_fmt(mean)});")
        return (f"{indent}c{clade} ~ UniformMRCA(taxa={taxa}, "
                f"upper={_fmt(upper)}, lower={_fmt(lower)});")

    out = CALIBRATION.sub(replace, source)

    # the joint prior becomes a plain array of independent calibrations
    out, n = CONDITIONED.subn(
        lambda m: f"{m.group('indent')}{m.group('name')} = "
                  f"toArray(calibrations={m.group('list')});",
        out,
    )
    if n != 1:
        raise SystemExit(f"expected exactly one ConditionedMRCAPrior line, found {n}")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", help="calibrationPrior .lphy script")
    ap.add_argument("--scheme", choices=("uniform", "suggested", "both"), default="both")
    ap.add_argument("-o", "--output", help="output file (single scheme only); "
                                           "default writes next to the input")
    args = ap.parse_args()

    source = open(args.input).read()
    schemes = ("uniform", "suggested") if args.scheme == "both" else (args.scheme,)
    if args.output and len(schemes) > 1:
        raise SystemExit("-o takes a single --scheme")

    for scheme in schemes:
        out = args.output or os.path.join(os.path.dirname(args.input) or ".",
                                          scheme + "Prior.lphy")
        with open(out, "w") as fh:
            fh.write(swap(source, scheme))
        print("wrote", out)


if __name__ == "__main__":
    main()
