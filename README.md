# The Calibrated Coalescent Point Process
[![Docker CI](https://github.com/BEAST2-Dev/CalibratedCPP/actions/workflows/maven-tests.yml/badge.svg)](https://github.com/BEAST2-Dev/CalibratedCPP/actions/workflows/maven-tests.yml)

This is an implementation of the Calibrated Coalescent Point Process (Calibrated CPP) in BEAST3.

The CPP is a model of ultrametric trees where node ages are i.i.d. random variables. This captures a general class of birth-death processes with time-dependent birth and death rates (as in BDSKY) and age dependent death rates.
- Individuals die at some rate $\mu(x,t)$ that depends on the age of the individual $x$ and the time $t$.
- Individuals give birth to new individuals at a rate $\lambda(t)$ depending only on time.

[Lambert and Stadler (2012)](https://doi.org/10.1016/j.tpb.2013.10.002) show that these models give a uniform distribution over ranked labelled (or oriented) tree topologies, and that the node ages are i.i.d. random variables. The node age density where $\lambda=\lambda(t)$ and $\mu=\mu(x,t)$ and each individual is sampled with probability $\rho$ is,

$$q(t) = \frac{\rho\lambda (\lambda-\mu)e^{-(\lambda-\mu)t}}{\rho\lambda+(\lambda(1-\rho)-\mu)e^{-(\lambda-\mu)t}}$$

the cumulative distribution function is,

$$Q(t) = \frac{\rho\lambda(1-e^{-(\lambda-\mu)t})}{\rho\lambda+(\lambda(1-\rho)-\mu)e^{-(\lambda-\mu)t}}.$$

The Calibrated CPP is a calibrated tree prior using the CPP. Calibrated tree priors are used for molecular clock dating by conditioning on the existence and ages of the most recent common ancestors of monophyletic clades.

## For developer

Build the package:
```bash
mvn clean package -DskipTests
```

Run an example XML with BEAST:
```bash
mvn -pl calibratedcpp-beast exec:exec -Dbeast.args="-validate src/test/resources/calibratedcpp/examples/b3test/b3_b3.xml"
```
Start BEAUti:
```bash
cd calibratedcpp-beast
mvn exec:exec -Dbeast.module=beast.fx -Dbeast.main=beastfx.app.beauti.Beauti
```

### Release

This repo does not yet have an automated `ci-publish.yml` / `release.sh` (unlike
the [gold-standard BEAST 3 package setup](https://github.com/CompEvol/beast3/blob/master/scripts/package-release-setup.md)),
so releases are done manually:

1. Tag and push the release:

```bash
git tag v0.0.1
git push origin v0.0.1
```

2. Build the BEAST3 package ZIP:

```bash
mvn clean package
```

This produces `calibratedcpp-beast/target/calibratedcpp-beast-<version>.zip`
(built from `calibratedcpp-beast/src/assembly/beast-package.xml`), containing
the package jar, `version.xml`, and examples.

Create a GitHub release for the tag and attach the ZIP.

Add/update the CBAN entry so the package is discoverable in BEAUti: submit
a PR to [CompEvol/CBAN/packages2.8.xml](https://github.com/CompEvol/CBAN/blob/master/packages2.8.xml)
pointing at the ZIP's download URL.

## Project structure
CalibratedCPP contains 5 subprojects:

### calibratedcpp-beast

The implementation has the following structure:
- `CalibratedCoalescentPointProcess` extends `SpeciesTreeDistribution` and takes a list of `calibrations`, `tree`, and the `origin` age OR `conditionOnRoot` as inputs.
- `CalibratedBirthDeathModel` extends `CalibratedCoalescentPointProcess` and implements `calculateLogNodeAgeDensity()` and `calculateLogNodeAgeCDF()` with node age density and CDF for the constant rate birth-death process.
- `CalibratedBirthDeathSkylineModel` extends `CalibratedCoalescentPointProcess` and implements `calculateLogNodeAgeDensity()` and `calculateLogNodeAgeCDF()` with node age density and CDF for the birth-death process with piecewise constant rates.
- `CalibratedAgeDependentExtinctionModel` extends `CalibratedCoalescentPointProcess` and implements `calculateLogNodeAgeDensity()` and `calculateLogNodeAgeCDF()` for the case when individuals have arbitrary lifetime distributions and give birth at a constant rate. The special case where individuals have Erlang distributed (Gamma distributed with integer shape parameter) lifetimes has a fast solution.

### calibratedcpp-lphy

The Lphy simulator is currently implemented within LinguaPhylo as a generative method called `CalibratedCPP`:
- The simulator takes birth rate, death rate, sampling probability for the birth-death process.
- Calibration information is taken from the output of `ConditionedMRCAPrior`, other leaf names are optional to pass in.
- Stem age can be passed in if there is no root calibration, but if root calibration and stemAge are both specified, the tree will still be root conditioned.

### calibratedcpp-lphy-studio

The Lphy studio is a visualiser that users can type in Lphy scripts directly or load local lphy scripts.
Calibrations can be viewed as calibration taxa names and the sampled ages.

### calibratedcpp-lphybeast

This subproject is converting Lphy simulators to XMLs for BEAST3 running:
- Construct BirthDeathSkylineModel with birth death parameters and conditions on the origin of the tree. The birth and death rates are SkylineParameter, which allow users to specify how times varies. Default false as BEAST3 model, users can manually change the XML to turn them on.
- For calibrations sourced from `ConditionedMRCAPrior`, the default is to build a single joint `CalibrationPrior` (preserving the nested/overlap calibration structure — LogNormal at each disjoint root, Beta on overlapping child/parent ratios, truncated LogNormal on nested non-overlapping children). Pass `-MRCAPrior` to `convert`/`run` to instead build one independent, per-clade `MRCAPrior(Uniform(lower,upper))` (no joint structure).
- `conditionOnCalibrations` on `CalibratedBirthDeathSkylineModel`/`CalibratedAgeDependentBirthDeathModel` defaults to `true` iff calibrations are provided. Pass `-conditionOnCalibrations true|false` to `convert`/`run` to override it, instead of hand-editing the output XML.

### calibratedcpp-lphybeast-launcher

This module is for launching lphybeast.
Run from the project root:
```angular2html
mvn -pl calibratedcpp-lphybeast-launcher exec:exec -Dlphybeast.args="convert ../calibratedcpp-lphy/examples/example.lphy"
mvn -pl calibratedcpp-lphybeast-launcher exec:exec -Dlphybeast.args="run -l 30000000 ../calibratedcpp-lphy/examples/example.lphy"
```

Two calibratedcpp-only flags are available on `convert`/`run`, in addition to lphybeast's own
options (paths inside `-Dlphybeast.args="..."` should be absolute, since the forked process's
working directory is `calibratedcpp-lphybeast-launcher/`, not the repo root):
```bash
# Default: joint CalibrationPrior for ConditionedMRCAPrior-sourced calibrations
mvn -pl calibratedcpp-lphybeast-launcher exec:exec -Dlphybeast.args="convert ../calibratedcpp-lphy/examples/example.lphy"

# -MRCAPrior: independent per-clade MRCAPrior(Uniform) instead of the joint CalibrationPrior
mvn -pl calibratedcpp-lphybeast-launcher exec:exec -Dlphybeast.args="convert -MRCAPrior ../calibratedcpp-lphy/examples/example.lphy"

# -conditionOnCalibrations: specify conditionOnCalibrations (true/false)
mvn -pl calibratedcpp-lphybeast-launcher exec:exec -Dlphybeast.args="convert -conditionOnCalibrations false ../calibratedcpp-lphy/examples/example.lphy"

# flags combine freely, and work with 'run' (convert + launch BEAST) the same way
mvn -pl calibratedcpp-lphybeast-launcher exec:exec -Dlphybeast.args="run -MRCAPrior -conditionOnCalibrations false ../calibratedcpp-lphy/examples/example.lphy"
```

## Well-Calibrated Study
Run LPhyBEAST using the configuration described in the [Developer Guide](#developer-guide).
Simulated parameter values will be written to the `.log` file, while the trees, alignment, and calibrations information will be saved to separate output files.
The simulated data are set as the starting values in the LPhyBEAST-generated XML, which can be passed directly to BEAST3 for Bayesian phylogenetic inference.
BEAST3 produces a posterior distribution that can be compared with the data generated by LPhyBEAST for well-calibrated study.

Run from the project root:
```bash
# Convert with 100 replicates
mvn -pl calibratedcpp-lphybeast-launcher exec:exec -Dlphybeast.args="convert -r 100 ../calibratedcpp-lphy/examples/wellCalibratedStudy.lphy"

# Run BEAST3 on each replicate
for i in {0..99}; do
    echo "Running wellCalibratedStudy_r${i}.xml"
    mvn -pl calibratedcpp-beast exec:exec -Dbeast.args="-overwrite examples/b3test/wellCalibratedStudy_r${i}.xml"
done
```

## License

CalibratedCPP is free software. It may be modified and distributed under the terms
of the GNU General Public License version 3 or, at your option, any later
version. A copy of this license should be found in the file COPYING located in
the root directory of this repository. If this file is absent for some reason,
it can also be retrieved from https://www.gnu.org/licenses.
