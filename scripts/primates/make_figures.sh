#!/usr/bin/env bash
# Regenerate the figures and the RF table for one dataset into figures/.
set -euo pipefail
cd "$(dirname "$0")"
DATASET=${DATASET:-grid-codon-deepcal}
FMT=${FMT:-pdf}
python3 plot_posterior_condCal.py --format=$FMT --dataset=$DATASET Primates   # <dataset>_Primates_age
python3 plot_clade_condCal.py --format=$FMT --dataset=$DATASET                # <dataset>_clade_ages
python3 plot_node_distributions.py --format=$FMT --dataset=$DATASET           # <dataset>_node_distributions_cond
python3 compare_rf.py --dataset=$DATASET                                      # <dataset>_rf.csv
