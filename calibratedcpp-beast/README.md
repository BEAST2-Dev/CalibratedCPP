# Calibrated Coalescent Point Process BEAST3 Implementation

Birth-death models are typically parameterised by a per-lineage birth-rate $(\lambda)$, a per-lineage death-rate $(\mu)$, and a probability with which each individual at the present time is sampled $(\rho)$. Other parameters of interest are:
- The reproductive number: $R=\lambda/\mu$.
- The diversification rate: $d=\lambda-\mu$.
- Turnover: $\tau=\mu/\lambda$.

We can specify any two of these parameters together except for turnover $(\tau)$ and the reproductive number $(R)$.

Both `CalibratedBirthDeathModel` and `CalibratedBirthDeathSkylineModel` extend `CalibratedCoalescentPointProcess` which itself extends `SpeciesTreeDistribution`. The  inputs of a `CalibratedCoalescentPointProcess` object are:
- `origin`: A `RealParameter` for the start time of the process. This is optional.
- `conditionOnRoot`: A `Boolean` which is true the process is conditioned on the root height. Either `origin` must be provided or `conditionOnRoot` should be set to true.
- `tree`: A phylogenetic tree.
- `calibrations`: A list of `TaxonSet` objects which we want to condition being monophyletic and the TMRCA of the clade. Note: providing the all the taxa in the tree as a `TaxonSet` is equivalent to setting `conditionOnRoot=true`.
- `calibrationForest`: A `CalibrationForest` supplying the same clades as an alternative to `calibrations` (provide exactly one of the two). It is typically a `CalibrationForestParser` built from an annotated Newick constraint tree. Because it *is* a `CalibrationForest`, a single parsed object can be shared by `idref` with a `CalibrationPrior` so both read identical clades. See the [documentation](https://beast2-dev.github.io/CalibratedCPP/) for the constraint-tree format.
- `conditionOnCalibrations`: A `Boolean` which is `true` if we want to condition on the `calibrations` provided and `false` otherwise. Default is `true`.

## `CalibratedBirthDeathModel`

The birth-death model takes as inputs any TWO of birth-rate $(\lambda)$, death-rate $(\mu)$, diversification rate $(d)$, reproductive number $(R)$, and turnover $(\tau)$ as `RealParameter` objects; and the sampling probability $(\rho)$ as a `RealParameter` in the interval $[0,1]$.

## `CalibratedBirthDeathSkylineModel`

The birth-death skyline model takes as inputs any TWO of the time dependent piecewise constant parameters: birth-rate $(\lambda(t))$, death-rate $(\mu(t))$, diversification rate $(d(t))$, effective reproductive number $(R_e(t))$, and turnover $(\tau(t))$ as `SkylineParameter` objects; and the sampling probability $(\rho)$ as a `RealParameter` in the interval $[0,1]$.

`SkylineParameter` objects take as inputs:
- `values`: A `RealParameter` object specifying the value of the parameter in each time interval ordered from root to tips.
- `changeTimes`: Optional, change times can be specified as `RealParameter`objects. If this is not provided, the rates are spread evenly over the age of the tree (origin/root height).
- `timesAreAges`: A `Boolean` which is `true` if the change times are given as ages before the present. The default is `false`.  
- `timesAreRelative`: A `Boolean` which is `true` if the change times are given as relative to the age of the tree (origin/root height). The default is `false`.

## `CalibratedAgeDependentExtinctionModel`

The age dependent extinction model describes a population where the lifetime of individuals are drawn from an arbitrary distribution and give birth to new individuals at a constant rate $\lambda$.

`CalibratedAgeDependentExtinctionModel` takes as inputs:
- `lifetimeDistribution`: A `ParametricDistribution` object describing the lifetimes of the individuals. The special case of an Erlang distribution (Gamma with integer shape parameter) is handeled separately from the general case which requires solving a Volterra integro-differential equation (VIDE). 
- `birthRate`: The rate at which individuals give birth to new individuals.
- `rho`: The probability with which each extant individual is sampled.



## Validation

In the [`validation`](../validation) directory we validate the `BirthDeathModel` against the Heled and Drummond implementation of the Calibrated Birth-Death process and benchmark the performance. We also benchmark the likelihood calculation time on 93 real calibration analyses available on [phylodata](https://pypi.org/project/phylodata/).
