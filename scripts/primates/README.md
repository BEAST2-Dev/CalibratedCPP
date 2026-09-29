# Primates analysis scripts

Everything behind the primates figures: does conditioning the tree prior on the
calibrations (`conditionOnCalibrations`) change the estimated node ages?

The `codon_deepCal` analysis (`--dataset=grid-codon-deepcal`) is the codon-partitioned
alignment with the `codon` calibration set plus de Vries and Beck node #6 (crown
Strepsirrhini, uniform 36.573-55.8 Ma), nine calibrations in all. Source:
`primates_codon_deepCal*.lphy`. It is run as 3 calibration schemes x 2 tree priors:

    primates_codon_deepCal[_<scheme>][_condFalse]

| scheme infix | calibration scheme |
|---|---|
| (none) | joint (Beta-LogNormal) calibration prior |
| `_suggested` | de Vries and Beck independent calibration priors |
| `_uniform` | uniform independent calibration priors |

`_condFalse` marks the regular (unconditioned) tree prior. Each analysis was run as
five replicates plus one sample-from-prior chain.

## Folders

- `data/` — the `codon_deepCal` runs: five replicates and the sample-from-prior chains
  (`*-fromPrior.log`/`.trees`).
- `data/combined/` — the LogCombiner output (`*-combined.log`/`.trees`, duplicate-seed
  replicates left out) and the `*_summary.tree` files written by TreeAnnotator.
- `figures/` — the figures and the RF table, named `<dataset>_*`.
- `clade_names.lphy` — clade definitions (including uncalibrated ones such as Primates)
  used to name clades in the figures.

The XMLs live in `calibratedcpp-beast/src/test/resources/calibratedcpp/examples/primates/xmls`
(`primatesGridNewCalibration/` for `codon_deepCal`).

## Figures

```bash
./make_figures.sh                    # DATASET=grid-codon-deepcal FMT=pdf by default
```

runs, for one dataset:

```bash
python3 plot_posterior_condCal.py --dataset=$d Primates   # <d>_Primates_age.pdf
python3 plot_clade_condCal.py --dataset=$d                # <d>_clade_ages.pdf
python3 plot_node_distributions.py --dataset=$d [--cond=false]   # <d>_node_distributions_cond.pdf
python3 compare_rf.py --dataset=$d [--unrooted]           # <d>_rf.csv
```

`plot_posterior_condCal.py` draws one clade's age (default `Primates`; any name from
`clade_names.lphy` or the XML's LPhy header works) under conditioning on vs off,
sample-from-prior above posterior, one column per calibration scheme.
`plot_clade_condCal.py` plots every clade's median age + 95% HPD, regular (x) vs
calibrated (y) tree prior. `plot_node_distributions.py` draws the posterior age of every
node under the three schemes, with the de Vries and Beck calibration priors beneath.
`compare_rf.py` prints the pairwise Robinson-Foulds distances among the six summary
trees and the clades on which any pair differs. All default to `grid-codon-deepcal`.

Clade ages are read from the `.trees` files rather than the traces, because BEAST only
logs `mrca.age()` for calibrated clades; they are cached as `.ages_*.npy` beside the tree
files, delete those to force a re-read.

`plot_posterior_condCal.py` is self-contained; the others share `primates_common.py`
(file naming, XML/TaxonSet parsing, summary-tree reader).

## Combining runs and summary trees

```bash
./run_logcombiner.sh [DIR]           # <stem>-rep{$REPS} -> DIR/combined/<stem>-combined (default data, REPS="1 2 3")
./run_treeannotator.sh [-f] [DIR]    # .trees -> <stem>_summary.tree (default data/combined)
```

`run_logcombiner.sh` and `run_treeannotator.sh` each point at a local install; override
with `$LOGCOMBINER` / `$TREEANNOTATOR`.
