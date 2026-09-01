# Primates analysis scripts

Everything behind the primates figures: does conditioning the tree prior on the
calibrations (`conditionOnCalibrations`) change the estimated node ages?

Three calibration schemes are compared, each run with conditioning on and off:

| file stem | scheme |
|---|---|
| `calibrationPrior` | joint (Beta-LogNormal) calibration prior |
| `suggestedPrior` | de Vries and Beck independent calibration priors |
| `uniformPrior` | uniform independent calibration priors |

A `-condFalse` suffix marks the regular (unconditioned) tree prior, and
`-fromPrior` marks a sample-from-prior run.

## Folders

- `data/` — BEAST output for the nogapN alignment: `.txt` traces, `.trees` posteriors,
  and the `*_summary.tree` files written by TreeAnnotator. `data/superseded/` holds
  chains replaced by later runs.
- `codons/` — the same, for the codon-partitioned analyses.
- `iqtree/` — ML reference trees (`primates.treefile` for the full alignment,
  `primates_nogapN.treefile` for the filtered one) and their `.iqtree` reports. These
trees are generated under the same models as BEAST3 runs.
- `refTrees/` — published reference topologies (de Vries and Beck, ASTRAL).

The XMLs themselves live in
`calibratedcpp-beast/src/test/resources/calibratedcpp/examples/primates/xmls`.

## Figures

```bash
python3 plot_posterior_condCal.py            # nogapN_condCal_Primates.png
python3 plot_posterior_condCal.py Colobinae  # any named clade
python3 plot_clade_condCal.py                # nogapN_cladeCondCal.png
```

Every script defaults to the `nogapN` dataset in `data/`; no arguments needed.

`plot_posterior_condCal.py` is self-contained. `plot_clade_condCal.py`,
`compare_rf.py` and `compare_mrca_ages.py` share `primates_common.py`, which holds the
file naming, the XML/TaxonSet parsing and the summary-tree reader. The two plotting
scripts and `compare_rf.py` take `--dataset=NAME` (`nogapN`, `codon`, `codon-nogapN`)
to switch alignments; the codon sets read `codons/`. `compare_mrca_ages.py` instead
takes folder arguments (`python3 compare_mrca_ages.py codons`), since it labels each
log from whichever XML matches its filename.

Note the `codon`/`codon-nogapN` runs currently in `codons/` were built on alignments
that are not in reading frame, so their codon partitions carry no positional signal
(r1 ~ r2 ~ r3). The in-frame XMLs are in `xmls/newRun/`.

Clade ages are read from the `.trees` files rather than the traces, because
BEAST only logs `mrca.age()` for calibrated clades. They are cached in
`data/.ages_*.npy`; delete those to force a re-read.

## Comparisons

```bash
python3 compare_rf.py [--unrooted]   # RF distances, summary trees vs references
python3 compare_mrca_ages.py [DIR]   # node ages (mean + 95% HPD); default data/
```

## Preparing inputs

```bash
./run_iqtree.sh                      # ML tree -> iqtree/
./run_treeannotator.sh [-f] [DIR]    # .trees -> <stem>_summary.tree (default: data)
python3 filter_gap_columns.py        # drop columns with any gap/N
python3 filter_codon_columns.py      # same, but whole codons, frame preserved
python3 sample_codon_alignment.py    # concatenate N random coding loci
python3 make_sample_from_prior.py    # derive the 8 XML variants from the base XMLs
python3 compute_offset_exponential_mean.py   # de Vries and Beck Table 1 -> Exponential means
```

`run_iqtree.sh` and `run_treeannotator.sh` each point at a local install; override with
`$IQTREE` / `$TREEANNOTATOR`.
