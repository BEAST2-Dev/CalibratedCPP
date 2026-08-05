# Calibration Prior Validation

Validates the `CalibrationPrior` BEAST likelihood by comparing MCMC samples from
the BEAST XML against ages drawn directly from the equivalent LPhy `ConditionedMRCAPrior`
simulator.  Agreement between the two confirms that the likelihood and the simulator
are targeting the same distribution.

The validation uses 10 calibrated clades on 11 taxa.  Clade numbering matches
`calibrationprior_simulation.r` (clade 1 = root):

| Clade | Taxa         | Lower | Upper | Node type            |
|------:|:-------------|------:|------:|:---------------------|
|     1 | leaf 1–11    | 10.0  | 10.5  | lognormal (root)     |
|     2 | leaf 1–8     |  9.4  |  9.6  | truncated lognormal  |
|     3 | leaf 9–11    |  4.8  |  5.0  | truncated lognormal  |
|     4 | leaf 9–10    |  4.0  |  5.0  | Beta (overlaps 3)    |
|     5 | leaf 7–8     |  8.0  |  9.5  | Beta (overlaps 2)    |
|     6 | leaf 1–6     |  6.8  |  8.0  | truncated lognormal  |
|     7 | leaf 5–6     |  1.8  |  2.5  | truncated lognormal  |
|     8 | leaf 1–4     |  2.0  |  3.2  | truncated lognormal  |
|     9 | leaf 1–3     |  1.7  |  2.5  | Beta (overlaps 8)    |
|    10 | leaf 1–2     |  1.6  |  2.5  | Beta (overlaps 9)    |

## Files

| File | Description |
|:-----|:------------|
| `test_calibration_prior.xml` | BEAST 2.8 XML — samples tree topology under the `CalibrationPrior` likelihood |
| `test_calibration_prior.lphy` | LPhy script — defines the same 10 calibrations via `ConditionedMRCAPrior` |
| `calibration_prior_comparison.py` | Python script — overlaid density plots comparing BEAST MCMC vs LPhy samples |
| `test_calibration_prior.log` | BEAST MCMC log (MRCA ages, produced by running the XML) |
| `lphy_wsim.tsv` | LPhy simulated ages (produced by running the Java sampler) |
| `mcmc_vs_lphy_comparison.pdf` | Output comparison plots |

## Usage

### 1 - Generate BEAST XML with LPhyBeast

```bash
mvn -pl calibratedcpp-lphybeast exec:exec \
  -Dlphybeast.args="convert -o test_calibration_prior.xml \
    ../calibratedcpp-beast/validation/calibrationprior/test_calibration_prior.lphy"
```

### 2 — Run BEAST MCMC

From the project root:

```bash
mvn -pl calibratedcpp-beast exec:exec \
  -Dbeast.args="-overwrite -prefix validation/calibrationprior/ validation/calibrationprior/test_calibration_prior.xml"
```

This produces `test_calibration_prior.log` with columns `mrca.clade1` … `mrca.clade10`.

### 3 — Generate LPhy simulated samples

From the `calibratedcpp-lphy` module directory:

```bash
cd calibratedcpp-lphy
mvn test-compile exec:java \
  -Dexec.mainClass=CalibrationPriorValidationSampler \
  -Dexec.classpathScope=test
```

This runs `ConditionedMRCAPrior.sample()` 50 000 times and writes
`../calibratedcpp-beast/validation/calibrationprior/lphy_wsim.tsv`.
Per-clade coverage is printed to stdout.

### 4 — Plot comparison

```bash
cd calibratedcpp-beast/validation/calibrationprior
python calibration_prior_comparison.py
```

Requires `numpy`, `scipy`, `pandas`, `matplotlib`.
Produces `mcmc_vs_lphy_comparison.pdf` with overlaid density plots for each clade,
dashed lines at the calibration bounds, and a console table of coverage and t-test results.

## Interpreting results

Each plot shows the marginal distribution of a clade's MRCA age under:
- **BEAST MCMC** (orange) — posterior sampled by MCMC with the `CalibrationPrior` likelihood and no data
- **LPhy sim** (blue) — direct draws from `ConditionedMRCAPrior`

The two distributions should overlap closely.  The coverage printed to stdout
(fraction of samples falling within the calibration bounds) should be near 0.90
for all clades, since the default confidence level is 90%.
