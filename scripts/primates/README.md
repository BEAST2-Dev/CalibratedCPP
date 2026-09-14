# Primates analysis scripts

Everything behind the primates figures: does conditioning the tree prior on the
calibrations (`conditionOnCalibrations`) change the estimated node ages?

The grid is 3 alignments x 3 calibration schemes x 2 tree priors, run as

    primates_<alignment>[_<scheme>][_condFalse]

| alignment | `--dataset` | |
|---|---|---|
| `codon` | `grid-codon` | codon-partitioned, full alignment |
| `codon_noGapN` | `grid-codon-nogapN` | codon-partitioned, gap/N codons removed |
| `unpartition` | `grid-unpartition` | unpartitioned, gap/N columns removed |

| scheme infix | calibration scheme |
|---|---|
| (none) | joint (Beta-LogNormal) calibration prior |
| `_suggested` | de Vries and Beck independent calibration priors |
| `_uniform` | uniform independent calibration priors |

`_condFalse` marks the regular (unconditioned) tree prior. Each analysis was run as
three replicates plus one sample-from-prior chain.

## Folders

- `grid/data/` — the LogCombiner output (`*-combined.log`/`.trees`, three replicates
  pooled after a 10% burn-in each), the sample-from-prior chains (`*-fromPrior.txt`/
  `.trees`) and the `*_summary.tree` files written by TreeAnnotator. The per-replicate
  runs are not kept.
- `grid/figures/` — the SI figures.
- `iqtree/` — ML reference trees and their `.iqtree` reports, generated under the same
  models as the BEAST runs; `primates_400loci*.treefile` are the three grid alignments,
  the two codon ones run with `-p iqtree/<name>_codon.nex` so IQ-TREE partitions by
  codon position exactly as BEAST does.
- `refTrees/` — published reference topologies (de Vries and Beck, ASTRAL).

The XMLs live in
`calibratedcpp-beast/src/test/resources/calibratedcpp/examples/primates/xmls/primatesGrid`.

## Figures

```bash
for d in grid-codon grid-codon-nogapN grid-unpartition; do
    python3 plot_posterior_condCal.py --format=pdf --dataset=$d Primates   # <d>_Primates_age.pdf
    python3 plot_clade_condCal.py --format=pdf --dataset=$d                # <d>_clade_ages.pdf
done
python3 plot_clade_merged.py         # the three clade figures stacked -> grid_clade_ages.pdf

python3 compare_mrca_ages.py         # grid/data logs -> grid/mrca_age_comparison_grid.csv
python3 compute_grid_stats.py        # grid/data trees -> grid/topology_support.json
python3 plot_grid_bounds.py [--cond=false]   # SI: calibrated nodes vs fossil bounds
python3 plot_grid_topology.py                # SI: the two contested nodes, PP per analysis
```

`plot_posterior_condCal.py` draws one clade's age (default `Primates`; any name from
the XML's LPhy header works) under conditioning on vs off, sample-from-prior above
posterior, one column per calibration scheme. `plot_clade_condCal.py` does the same for
every calibrated clade at once. Both default to `grid-codon-nogapN`. Clade ages are read
from the `.trees` files rather than the traces, because BEAST only logs `mrca.age()`
for calibrated clades; they are cached in `grid/data/.ages_*.npy`, delete those to force
a re-read.

`plot_clade_merged.py` stacks the three per-alignment clade figures into
`grid_clade_ages.pdf`, each block under a subtitle naming its alignment, with one set
of axis limits shared across all nine panels so the rows can be compared directly.

`plot_posterior_condCal.py` is self-contained. `plot_clade_condCal.py`,
`plot_clade_merged.py` and `compare_mrca_ages.py` share `primates_common.py` (file
naming, XML/TaxonSet parsing, summary-tree reader). `grid_common.py` holds the grid
design, file locations and the Newick walker for `compute_grid_stats.py` and the two
`plot_grid_*.py` scripts.

`compare_mrca_ages.py` tabulates every `mrca.age()` column (mean + 95% HPD) across the
logs, labelling clades from each run's own XML. `compute_grid_stats.py` is pure Python
and reads each combined `.trees` once to get the posterior support for each resolution
of the two contested nodes, which is absent from the MCC tree. `plot_grid_bounds.py`
reads the CSV and the fossil bounds from the XML's LPhy header.

## Other checks

```bash
python3 check_ess.py [DIR]           # min ESS per combined run (default grid/data)
```

## Preparing inputs

```bash
./run_iqtree.sh                      # ML tree -> iqtree/
./run_logcombiner.sh DIR             # <stem>-rep{1,2,3} -> DIR/combined/<stem>-combined
./run_treeannotator.sh [-f] DIR      # .trees -> <stem>_summary.tree
python3 filter_gap_columns.py        # drop columns with any gap/N
python3 filter_codon_columns.py      # same, but whole codons, frame preserved
python3 sample_codon_alignment.py    # concatenate N random coding loci
python3 make_prior_variants.py       # derive the scheme / condFalse / fromPrior XMLs
python3 swap_calibration_prior.py    # rewrite a joint-prior LPhy script into an independent-prior scheme
python3 make_sample_from_prior.py    # older 8-variant generator for the pre-grid base XMLs
python3 compute_offset_exponential_mean.py   # de Vries and Beck Table 1 -> Exponential means
```

`run_iqtree.sh`, `run_logcombiner.sh` and `run_treeannotator.sh` each point at a local
install; override with `$IQTREE` / `$LOGCOMBINER` / `$TREEANNOTATOR`.
