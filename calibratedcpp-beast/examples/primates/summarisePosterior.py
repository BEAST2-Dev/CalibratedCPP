#!/usr/bin/env python3
"""
Crown-Primates (TaxonSet3) age summary from BEAST 2 logs: mean age and the fraction
of post-burnin trees older than the K-Pg boundary (66 Ma), per run.

Only nogapN files are used (plain nogap is intentionally excluded).

Usage:
    Edit the RUNS dict below and run:
        python3 bayes_factor_age66.py
"""

import csv
import math
import os

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(BASE, "data")
OUT_CSV = os.path.join(BASE, "taxonset3_age_summary.csv")

RUNS = {
    # CPP, suggested priors (offset-exponential + uniform), conditionOnCalibrations=true (default)
    "CPP_suggested, with data":               os.path.join(DATA, "allCalibrations_suggestedPriors_nogapN.txt"),
    "CPP_suggested, prior only":              os.path.join(DATA, "sample-from-prior_allCalibrations_suggestedPriors.txt"),
    # CPP, suggested priors, conditionOnCalibrations=false
    "CPP_suggested_condCalFalse, with data":  os.path.join(DATA, "allCalibrations_suggestedPriors_nogapN_condCalFalse.txt"),
    "CPP_suggested_condCalFalse, prior only": os.path.join(DATA, "sample-from-prior_allCalibrations_suggestedPriors_condCalFalse.txt"),
    # CPP, uniform priors (all calibrations Uniform), conditionOnCalibrations=true (default)
    "CPP_uniform, with data":                 os.path.join(DATA, "allCalibrations_nogapN.txt"),
    "CPP_uniform, prior only":                os.path.join(DATA, "sample-from-prior_allCalibrations_uniform.txt"),
    # CPP, uniform priors, conditionOnCalibrations=false
    "CPP_uniform_condCalFalse, with data":    os.path.join(DATA, "allCalibrations_nogapN_condCalFalse.txt"),
    "CPP_uniform_condCalFalse, prior only":   os.path.join(DATA, "sample-from-prior_allCalibrations_uniform_condCalFalse.txt"),
    # BD (BirthDeathGernhard08Model), suggested priors -- no conditionOnCalibrations equivalent
    "BD_suggested, with data":                os.path.join(DATA, "allCalibrations_birthDeath_suggestedPriors_nogapN.txt"),
    "BD_suggested, prior only":               os.path.join(DATA, "sample-from-prior_allCalibrations_birthDeath_suggestedPriors.txt"),
    # BD, uniform priors (all calibrations Uniform) -- no conditionOnCalibrations equivalent
    "BD_uniform, with data":                  os.path.join(DATA, "allCalibrations_birthDeath_nogapN.txt"),
    "BD_uniform, prior only":                 os.path.join(DATA, "sample-from-prior_allCalibrations_birthDeath_uniform.txt"),
}

# Each pair is (H0 run name, H1 run name), both "with data" runs. BF = odds_older(H0) /
# odds_older(H1) -- a direct odds ratio between the two runs' own older/younger split,
# with no prior-normalization step (unlike the earlier within-model younger/older BF,
# both sides here already have data).
BF_PAIRS = [
    ("CPP_suggested, with data", "CPP_suggested_condCalFalse, with data"),
    ("CPP_suggested, with data", "BD_suggested, with data"),
    ("CPP_suggested, with data", "CPP_uniform, with data"),
    ("CPP_suggested_condCalFalse, with data", "CPP_uniform_condCalFalse, with data"),
    ("BD_suggested, with data", "BD_uniform, with data"),
]

# Supervisor's cross-model BF: M1 = the model that concludes younger, M2 = the model that
# concludes older. BF = [P(M1|D)/P(M2|D)] x [P(M2)/P(M1)], using RAW "count over 66"
# proportions (not odds) for all four terms -- P(Mi|D) from each model's own with-data run,
# P(Mi) from each model's own prior-only run. Each tuple is (label_m1, label_m2), where the
# labels index into RUNS via "<label>, with data" / "<label>, prior only".
MODEL_PAIRS = [
    ("BD_suggested", "CPP_suggested"),
    ("BD_uniform", "CPP_uniform"),
]

COLUMN = "mrca.age(TaxonSet3)"
BOUNDARY = 66.0
BURNIN = 0.10


def read_column(path, column, burnin):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = header.index(column)
        vals = []
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= idx:
                continue
            vals.append(float(fields[idx]))
    cut = int(len(vals) * burnin)
    return vals[cut:]


def hpd_interval(vals, mass=0.95):
    """Narrowest interval containing `mass` fraction of the (sorted) samples."""
    s = np.sort(np.array(vals))
    n = len(s)
    interval_idx = int(np.floor(mass * n))
    if interval_idx >= n:
        return s[0], s[-1]
    widths = s[interval_idx:] - s[: n - interval_idx]
    best = np.argmin(widths)
    return s[best], s[best + interval_idx]


def interpret(bf, name_h0, name_h1):
    """Kass & Raftery (1995) scale, stated for whichever hypothesis wins."""
    b = bf if bf >= 1 else 1.0 / bf
    favored = name_h0 if bf >= 1 else name_h1
    twoln = 2 * math.log(b)
    if twoln < 2:
        strength = "not worth more than a bare mention"
    elif twoln < 6:
        strength = "positive"
    elif twoln < 10:
        strength = "strong"
    else:
        strength = "very strong"
    return twoln, favored, strength


def odds(ratio):
    return ratio / (1 - ratio)


def main():
    print(f"TaxonSet3 (Primates) age summary | boundary={BOUNDARY} Ma | burn-in={BURNIN:.0%}\n")

    # Pass 1: read every run's own stats (mean, HPD, ratio_older). No BF yet -- BF needs a
    # with-data run paired against its own matching prior-only run, computed in pass 2.
    stats = {}
    for name, path in RUNS.items():
        if not os.path.exists(path):
            print(f"{name}: MISSING ({path})")
            continue
        vals = read_column(path, COLUMN, BURNIN)
        n = len(vals)
        mean_age = sum(vals) / n
        hpd_lo, hpd_hi = hpd_interval(vals)
        n_older = sum(1 for v in vals if v >= BOUNDARY)
        ratio_older = n_older / n
        stats[name] = {"n": n, "mean_age": mean_age, "hpd_lo": hpd_lo, "hpd_hi": hpd_hi,
                        "n_older": n_older, "ratio_older": ratio_older}

    # Pass 2: BF = odds_older(with data) / odds_older(matching prior only) -- the data/prior
    # comparison that's actually derivable from Bayes' rule (see the earlier derivation).
    # Only "with data" rows get a BF; "prior only" rows are just the reference point used
    # inside that division, so they don't get one of their own.
    rows = []
    results = {}
    for name, s in stats.items():
        results[name] = s
        bf, support = None, ""
        if name.endswith(", with data"):
            prior_name = name.replace(", with data", ", prior only")
            s_prior = stats.get(prior_name)
            if s_prior is None:
                support = "n/a (no prior-only run)"
            elif s["n_older"] in (0, s["n"]) or s_prior["n_older"] in (0, s_prior["n"]):
                support = "n/a (empty tail)"
            else:
                bf = odds(s["ratio_older"]) / odds(s_prior["ratio_older"])
                support = "over 66" if bf > 1 else "under 66"

        print(f"{name}:")
        print(f"  n = {s['n']}")
        print(f"  mean age = {s['mean_age']:.3f} Ma   95% HPD = [{s['hpd_lo']:.3f}, {s['hpd_hi']:.3f}]")
        bf_str = "" if bf is None else f"   BF (data/prior) = {bf:.4g}   supports: {support}"
        print(f"  older/total = {s['n_older']}/{s['n']} = {s['ratio_older']:.4f}{bf_str}\n")

        rows.append([name, s["n"], s["mean_age"], s["hpd_lo"], s["hpd_hi"], bf, support, s["n_older"], s["ratio_older"]])

    bf_rows = []
    for name_h0, name_h1 in BF_PAIRS:
        r0, r1 = results.get(name_h0), results.get(name_h1)
        print(f"=== BF: H0 = {name_h0}  vs  H1 = {name_h1} ===")
        if r0 is None or r1 is None:
            missing = [n for n, r in [(name_h0, r0), (name_h1, r1)] if r is None]
            print(f"  Missing: {', '.join(missing)} -- skipping.\n")
            continue
        if r0["n_older"] in (0, r0["n"]) or r1["n_older"] in (0, r1["n"]):
            print("  WARNING: an empty tail (0 or all trees older) -- odds/BF unreliable by counting.\n")
            continue

        odds_h0 = r0["ratio_older"] / (1 - r0["ratio_older"])
        odds_h1 = r1["ratio_older"] / (1 - r1["ratio_older"])
        bf = odds_h0 / odds_h1
        twoln, winner, strength = interpret(bf, name_h0, name_h1)

        print(f"  odds_older(H0) = {odds_h0:.4g}")
        print(f"  odds_older(H1) = {odds_h1:.4g}")
        print(f"  BF (H0/H1) = {bf:.4g}   2 ln BF = {twoln:.3f}")
        print(f"  Winner: {winner}, strength: {strength}\n")

        bf_rows.append([name_h0, name_h1, odds_h0, odds_h1, bf, twoln, winner, strength])

    model_rows = []
    for m1, m2 in MODEL_PAIRS:
        r1_data, r1_prior = results.get(f"{m1}, with data"), results.get(f"{m1}, prior only")
        r2_data, r2_prior = results.get(f"{m2}, with data"), results.get(f"{m2}, prior only")
        print(f"=== model BF: M1 = {m1} (younger)  vs  M2 = {m2} (older) ===")
        needed = {f"{m1}, with data": r1_data, f"{m1}, prior only": r1_prior,
                  f"{m2}, with data": r2_data, f"{m2}, prior only": r2_prior}
        missing = [n for n, r in needed.items() if r is None]
        if missing:
            print(f"  Missing: {', '.join(missing)} -- skipping.\n")
            continue

        p_m1_d, p_m1 = r1_data["ratio_older"], r1_prior["ratio_older"]
        p_m2_d, p_m2 = r2_data["ratio_older"], r2_prior["ratio_older"]
        bf = (p_m1_d / p_m2_d) * (p_m2 / p_m1)
        twoln, winner, strength = interpret(bf, m1, m2)

        print(f"  P(M1|D) [{m1}, with data, over 66]   = {p_m1_d:.4g}")
        print(f"  P(M2|D) [{m2}, with data, over 66]   = {p_m2_d:.4g}")
        print(f"  P(M1)   [{m1}, prior only, over 66]  = {p_m1:.4g}")
        print(f"  P(M2)   [{m2}, prior only, over 66]  = {p_m2:.4g}")
        print(f"  BF = (P(M1|D)/P(M2|D)) x (P(M2)/P(M1)) = {bf:.4g}   2 ln BF = {twoln:.3f}")
        print(f"  Winner: {winner}, strength: {strength}\n")

        model_rows.append([m1, m2, p_m1_d, p_m2_d, p_m1, p_m2, bf, twoln, winner, strength])

    with open(OUT_CSV, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["run", "n", "mean_age", "hpd95_lower", "hpd95_upper", "BF_data_over_prior", "support", "n_older", "ratio_older"])
        w.writerows(rows)
        w.writerow([])
        w.writerow(["H0", "H1", "odds_h0", "odds_h1", "BF", "twoln_bf", "winner", "strength"])
        w.writerows(bf_rows)
        w.writerow([])
        w.writerow(["M1", "M2", "P(M1|D)", "P(M2|D)", "P(M1)", "P(M2)", "BF", "twoln_bf", "winner", "strength"])
        w.writerows(model_rows)
    print(f"Wrote {OUT_CSV}")


if __name__ == "__main__":
    main()
