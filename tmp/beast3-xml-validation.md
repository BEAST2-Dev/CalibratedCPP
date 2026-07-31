# BEAST3 Example XML Validation

Validated all `*_b3.xml` files under `calibratedcpp-beast/src/test/resources/calibratedcpp/examples/`, excluding the `legacy/` folder (BEAST2 XMLs), using `BeastMain -validate` (parses and builds the model without running MCMC). Re-validated as of the latest edit: **40 pass, 10 fail**, unchanged from the prior run.

## Known issues: 10 currently-failing files, 3 independent root causes

**1. Missing class `orc.consoperators.UcldScalerOperator` (7 files)** — `b3test/cats_halfAlignment_b3.xml`, `lphyBeastXmls/cats_b3.xml`, `lphyBeastXmls/cats_halfAlignment_b3.xml`, `lphyBeastXmls/cats_halfAlignment_uncon_b3.xml`, `lphyBeastXmls/cats_uncon_b3.xml`, `lphyBeastXmls/cats_varSites_b3.xml`, `lphyBeastXmls/cats_varSites_uncon_b3.xml`. The single largest remaining blocker. Investigated by checking the sibling source checkout at `~/WorkSpace/ORC` (the `io.github.jordandouglas:beast-orc` package source):
- `git log --all -S"UcldScalerOperator"` in that repo shows the class was **deleted outright** in commit `41f6b51` ("preparation for v1.3.0 release for beast 2.8") — `src/orc/consoperators/UcldScalerOperator.java` (459 lines) was removed as part of migrating the whole package to Maven/beast 2.8, alongside several other files (old `build.xml`, some example XMLs). The diff shows no renamed replacement class added in the same commit — the current `orc` package (versions 1.3.0-SNAPSHOT/1.3.1, both checked) simply has no equivalent operator.
- The old class implemented a specialized joint scale-move on the UCLD relaxed-clock stdev (`ucldStdev`) together with branch rates, preserving proposal probability under the lognormal rate distribution — a purpose-built operator, not a generic one, so there's no drop-in equivalent already elsewhere in `beast-base` or the `orc` package.
- **Not yet fixed** — pending a decision on how to proceed (e.g. drop the wrapped operator and just keep the generic `ScaleOperator`/`RealRandomWalkOperator` siblings already present in the same `AdaptableOperatorSampler` block, which are structurally valid stand-ins, vs. some other resolution).

**2. Could not find object associated with idref `branchRates.prior` (2 files)** — `lphyBeastXmls/cats_sp_b3.xml`, `lphyBeastXmls/cats_sp_uncon_b3.xml`. A dangling reference, unrelated to the GTR/UcldScaler issues. **Not yet investigated.**

**3. `MCMC` has no input named `value` (1 file)** — `cats/xml/cats-sample_from_prior_consistent_prior_not_conditioned_b3.xml`. Likely a stray/malformed attribute directly under the `<run>` element. **Not yet investigated.**

All three are independent, pre-existing issues, unrelated to each other and to everything already fixed this session.

Update 1: the 26 files that failed with a removed `calibration.CalibrationClade` class were fixed (the model's `calibrations` input now takes `TaxonSet` directly, so the `<calibrations spec="calibration.CalibrationClade"><taxa spec="TaxonSet">...` wrapper was flattened to `<calibrations spec="TaxonSet">...`). 9 of those 26 now passed; the other 17 failed for unrelated, pre-existing reasons uncovered by fixing the CalibrationClade issue.

Update 2: two more targeted fixes applied, both of which resolved the *reported* error but uncovered a further, different pre-existing error in the same files:
- The 8 `lphyBeastXmls/*` files: the `pi` parameter's `spec` was changed from `RealVectorParam` to `SimplexParam` (a `RealVectorParam` subtype implementing `Simplex`), fixing the `frequencies` type mismatch. These now fail one step later — the `GTR` substitution model's `rates="@rates"` input no longer exists; GTR now expects six separate rate inputs (`rateAC`, `rateAG`, `rateAT`, `rateCG`, `rateCT`, `rateGT`) instead of one combined vector.
- `b3test/cats_halfAlignment_b3.xml`: added `io.github.beast2-dev:beast-labs:2.1.0-beta2` as a dependency of `calibratedcpp-beast` (built locally from `~/WorkSpace/BEASTLabs`, targets the same `beast-base` version this project uses), which resolves the missing `beastlabs.evolution.tree.RNNIMetric`. The file now fails one step later on a different missing class, `orc.consoperators.UcldScalerOperator`, which isn't in any ORC-related jar currently on the classpath.

Update 3: removed the obsolete `<log spec="calibrationprior.logger.MRCALogger">` entry from the 8 `cats/xml/*_consistent_prior*` files. That class was deleted from the codebase (commit `2af7ee4`, "Updated primates xmls") without a replacement, and the regenerated `primates` XMLs in that same commit simply dropped MRCA-age logging — this fix follows that precedent rather than inventing a new logger. **7 of the 8 now pass.** The 8th (`cats-sample_from_prior_consistent_prior_not_conditioned_b3.xml`) hits a different, pre-existing, unrelated error: `MCMC` has no input named `value`.

Update 4: `b3test/b3_b3.xml` used the wrong `spec` for its birth-death distribution — `spec="beast.base.spec.evolution.speciation.CalibratedBirthDeathModel"` points at a *different*, built-in BEAST-base class that happens to share the same simple name but has no `conditionOnRoot` input. This project has its own `calibratedcpp.CalibratedBirthDeathModel` (extends `CalibratedCoalescentPointProcess`, which is where `conditionOnRoot` actually lives), and the element's other attributes (`diversificationRate`, `turnover`, `rho`, nested `<calibrations spec="TaxonSet">`) all match that class's inputs exactly — confirming it was a namespace mix-up. Changed `spec` to `calibratedcpp.CalibratedBirthDeathModel`; the `conditionOnRoot` error is gone, but the file now hits the same pre-existing `GTR` `rates` issue as the `lphyBeastXmls/*` files (no net change in pass count).

Update 5: for the 8 `lphyBeastXmls/*` files, the `rates` vector was split into 5 independent `RealScalarParam(PositiveReal)` state nodes (`rateAC`,`rateAG`,`rateAT`,`rateCG`,`rateCT`; `rateGT` intentionally omitted so it defaults to `1.0`, mirroring the already-working `cats/xml/*_consistent_prior*` files), each given an independent `LogNormal(0,1)` prior in place of the old `Dirichlet` prior on the vector, wired into `GTR` via its 5 named rate attributes, with a `ScaleOperator` and log entries per rate replacing the old `DeltaExchangeOperator`/vector log entries.

For `b3test/b3_b3.xml` specifically, a different, narrower fix was used instead at the user's request: rather than splitting `rates`, the `substModel`'s `spec` was changed from `beast.base.spec.evolution.substitutionmodel.GTR` to `substmodels.nucleotide.GTR` (from the `io.github.beast2-dev:substmodels` package already on the classpath). That class's `rates` input is inherited from `GeneralSubstitutionModel` as `Input<RealVector<NonNegativeReal>>` — since Java generics are erased at runtime, the original 6-dim `SimplexParam` `rates` stateNode, its `Dirichlet` prior, and its `DeltaExchangeOperator` all still satisfy that input unchanged, so nothing else in the file needed to move. (Note: `b3test/b3_b3.xml` had reverted to its unmodified git state after an interrupted turn, losing the earlier `CalibrationClade`→`TaxonSet` and `CalibratedBirthDeathModel`-class fixes — both were reapplied before this GTR fix.) **`b3test/b3_b3.xml` now passes fully** (`-validate` prints `Done!`).

The `lphyBeastXmls/*` files still don't pass: splitting `rates` unblocked `GTR`, but 6 of the 8 then hit the same pre-existing `orc.consoperators.UcldScalerOperator` missing-class issue as `b3test/cats_halfAlignment_b3.xml`, and the other 2 hit a different pre-existing error (`Could not find object associated with idref branchRates.prior`) that was previously masked by the GTR error.

Update 6: at the user's request, all 8 `lphyBeastXmls/*` files were switched from `beast.base.spec.evolution.substitutionmodel.GTR` (split-scalar rates) to `substmodels.nucleotide.GTR` (single vector rates), matching the change made to `b3test/b3_b3.xml`. Since `substmodels.nucleotide.GTR` wants a single `rates` vector, not separate `rateAC`/`rateAG`/etc. inputs, the Update-5 scalar split was reversed for these 8 files: the original vector `rates` parameter, `Dirichlet` prior, and `DeltaExchangeOperator` were restored **exactly as they existed in git HEAD** (recovered via `git show`, not fabricated — including the true `rateGT` value that had been dropped during the scalar split), and only the `substModel spec` was changed. 16 other files (`cats/xml/*_consistent_prior*`, `*_uniform_mrca*`) also use the split-scalar style with the built-in `beast.base` GTR and currently pass — those were explicitly left untouched, per user instruction, and still use `beast.base.spec.evolution.substitutionmodel.GTR`. No change in pass/fail count: all 8 files were already blocked by the unrelated `orc` issues below, and swapping the GTR class doesn't touch those.

**Overall result: 40 pass, 10 fail.**

## lphyBeastXmls/ (8 files)

All 8 now use `substmodels.nucleotide.GTR` with the original vector-style `rates` (see Update 6). Each still hits one of two other pre-existing, unrelated errors that were previously masked by the GTR issue.

**Result: 0 pass, 8 fail.**

| XML | Status | Reason |
|---|---|---|
| lphyBeastXmls/cats_b3.xml | FAIL | Missing class `orc.consoperators.UcldScalerOperator` (GTR rates issue fixed) |
| lphyBeastXmls/cats_halfAlignment_b3.xml | FAIL | Missing class `orc.consoperators.UcldScalerOperator` (GTR rates issue fixed) |
| lphyBeastXmls/cats_halfAlignment_uncon_b3.xml | FAIL | Missing class `orc.consoperators.UcldScalerOperator` (GTR rates issue fixed) |
| lphyBeastXmls/cats_sp_b3.xml | FAIL | Could not find object associated with idref `branchRates.prior` (GTR rates issue fixed) |
| lphyBeastXmls/cats_sp_uncon_b3.xml | FAIL | Could not find object associated with idref `branchRates.prior` (GTR rates issue fixed) |
| lphyBeastXmls/cats_uncon_b3.xml | FAIL | Missing class `orc.consoperators.UcldScalerOperator` (GTR rates issue fixed) |
| lphyBeastXmls/cats_varSites_b3.xml | FAIL | Missing class `orc.consoperators.UcldScalerOperator` (GTR rates issue fixed) |
| lphyBeastXmls/cats_varSites_uncon_b3.xml | FAIL | Missing class `orc.consoperators.UcldScalerOperator` (GTR rates issue fixed) |

## Everything else (42 files)

**Result: 40 pass, 2 fail.**

### Failed (2)

| XML | Reason |
|---|---|
| b3test/cats_halfAlignment_b3.xml | Missing class `orc.consoperators.UcldScalerOperator` (beastlabs dependency fixed; this is the next missing class) |
| cats/xml/cats-sample_from_prior_consistent_prior_not_conditioned_b3.xml | `MCMC` has no input named `value` (MRCALogger issue fixed; this is a separate, unrelated error) |

### Passed (40)

- agedependent/cats_halfAlignment_ADB_b3.xml *(fixed)*
- b3test/b3_b3.xml *(fixed)*
- cats/xml/cats-full_consistent_prior_b3.xml *(fixed)*
- cats/xml/cats-full_consistent_prior_not_conditioned_b3.xml *(fixed)*
- cats/xml/cats-full_uniform_mrca_b3.xml *(fixed)*
- cats/xml/cats-full_uniform_mrca_not_conditioned_b3.xml *(fixed)*
- cats/xml/cats-half_consistent_prior_b3.xml *(fixed)*
- cats/xml/cats-half_consistent_prior_not_conditioned_b3.xml *(fixed)*
- cats/xml/cats-half_uniform_mrca_b3.xml *(fixed)*
- cats/xml/cats-half_uniform_mrca_not_conditioned_b3.xml *(fixed)*
- cats/xml/cats-sample_from_prior_consistent_prior_b3.xml *(fixed)*
- cats/xml/cats-sample_from_prior_uniform_mrca_b3.xml *(fixed)*
- cats/xml/cats-sample_from_prior_uniform_mrca_not_conditioned_b3.xml *(fixed)*
- cats/xml/cats-var_consistent_prior_b3.xml *(fixed)*
- cats/xml/cats-var_consistent_prior_not_conditioned_b3.xml *(fixed)*
- cats/xml/cats-var_uniform_mrca_b3.xml *(fixed)*
- cats/xml/cats-var_uniform_mrca_not_conditioned_b3.xml *(fixed)*
- primates/factorXml/allCalibrations_birthDeath_nogapN_b3.xml
- primates/factorXml/allCalibrations_birthDeath_suggestedPriors_nogap_b3.xml
- primates/factorXml/allCalibrations_birthDeath_suggestedPriors_nogapN_b3.xml
- primates/factorXml/allCalibrations_nogapN_b3.xml
- primates/factorXml/allCalibrations_nogapN_condCalFalse_b3.xml
- primates/factorXml/allCalibrations_suggestedPriors_nogap_b3.xml
- primates/factorXml/allCalibrations_suggestedPriors_nogapN_b3.xml
- primates/factorXml/allCalibrations_suggestedPriors_nogapN_condCalFalse_b3.xml
- primates/factorXml/sample-from-prior_allCalibrations_birthDeath_suggestedPriors_b3.xml
- primates/factorXml/sample-from-prior_allCalibrations_birthDeath_uniform_b3.xml
- primates/factorXml/sample-from-prior_allCalibrations_suggestedPriors_b3.xml
- primates/factorXml/sample-from-prior_allCalibrations_suggestedPriors_condCalFalse_b3.xml
- primates/factorXml/sample-from-prior_allCalibrations_uniform_b3.xml
- primates/primates_consistent_prior_b3.xml
- primates/primates_consistent_prior_not-conditioned_b3.xml
- primates/primates_noroot_node5_b3.xml
- primates/primates_remove8_b3.xml
- primates/primates_uniform_prior_b3.xml
- primates/primates_uniform_prior_not-conditioned_b3.xml
- primates/sample-from-prior_primates_consistent_prior_b3.xml
- primates/sample-from-prior_primates_consistent_prior_not-conditioned_b3.xml
- primates/sample-from-prior_primates_uniform_prior_b3.xml
- primates/sample-from-prior_primates_uniform_prior_not-conditioned_b3.xml
