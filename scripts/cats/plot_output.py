import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from pathlib import Path
import numpy as np
from scipy.stats import gaussian_kde
from scipy.spatial.distance import jensenshannon

# --- 1. SETUP ---
SCRIPT_DIR = Path(__file__).parent
os.chdir(SCRIPT_DIR)
os.makedirs('data', exist_ok=True)

# --- 2. CONFIGURATION ---
# TaxonSet1/2/3 are the three calibrated clades; the uniform-MRCA runs log them via
# their MRCAPriors, the calibration-prior runs via CalibrationPrior.log(), which emits
# the same mrca.age(<taxonset id>) column names.
target_columns = [
    'mrca.age(TaxonSet1)',
    'mrca.age(TaxonSet2)',
    'mrca.age(TaxonSet3)'
]

# Everything this script reads and writes lives in data/: the BEAST runs put their
# $(filebase).log there, and the figures go back alongside them.
DATA_DIR = 'data'
OUT_DIR = DATA_DIR

PRIOR_TYPES = [
    ('Calibration prior', 'consistent_prior'),
    ('Uniform MRCA prior', 'uniform_mrca'),
]

ANALYSES = [
    ('Prior', 'cats-sample_from_prior'),
    ('Half Alignment', 'cats-half'),
    ('Full Alignment', 'cats-full'),
]

CONDITIONS = [
    ('Conditioned', ''),
    ('Not Conditioned', '_not_conditioned'),
]


def build_file_map():
    """(filename, analysis, condition, prior type) for every run we know how to name."""
    entries = []
    for prior_label, prior_tag in PRIOR_TYPES:
        for analysis, stem in ANALYSES:
            for condition, cond_tag in CONDITIONS:
                entries.append(("%s_%s%s_b3.log" % (stem, prior_tag, cond_tag),
                                analysis, condition, prior_label))
    return entries


def resolve(filename):
    candidate = os.path.join(DATA_DIR, filename)
    return candidate if os.path.exists(candidate) else None


file_map = build_file_map()

burnin_fraction = 0.1

# --- 3. HELPER FUNCTION ---
def load_all_data(file_map):
    all_data = []

    for filename, analysis_type, condition_type, prior_type in file_map:
        path = resolve(filename)
        if path is None:
            print(f"File not found: {filename}")
            continue
        try:
            df = pd.read_csv(path, sep='\t', comment='#')
            cols = [c for c in target_columns if c in df.columns]
            if not cols:
                print(f"Warning: Targets not found in {path}")
                continue
            drop_n = int(len(df) * burnin_fraction)
            df = df.iloc[drop_n:].copy()

            melted = df[cols].melt(var_name='TaxonSet', value_name='Age')
            melted['Analysis'] = analysis_type
            melted['Condition'] = condition_type
            melted['Prior'] = prior_type

            melted['TaxonSet'] = melted['TaxonSet'].str.replace('mrca.age(', '', regex=False).str.replace(')', '', regex=False)
            all_data.append(melted)
            print(f"Loaded: {path}  [{prior_type} / {analysis_type} / {condition_type}]")
        except Exception as e:
            print(f"Error reading {path}: {e}")

    if all_data:
        return pd.concat(all_data, ignore_index=True)
    return pd.DataFrame()

# --- 4. LOAD DATA ---
df = load_all_data(file_map)

if df.empty:
    print("No data loaded.")
    exit()

# --- 5. PLOTTING THE HISTOGRAMS ---
plt.rcParams.update({
    'axes.labelsize': 18,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'legend.fontsize': 16,
    'pdf.fonttype': 42,      # Ensures fonts embed correctly in the PDF
    'ps.fonttype': 42
})

taxon_sets = [c.replace('mrca.age(', '').replace(')', '') for c in target_columns]
analysis_types = [a for a, _ in ANALYSES]
prior_types = [p for p, _ in PRIOR_TYPES]
colors = {'Conditioned': 'tab:blue', 'Not Conditioned': 'tab:orange'}

title_map = {
    'TaxonSet1': 'Clade 1',
    'TaxonSet2': 'Clade 2',
    'TaxonSet3': 'Clade 3'
}

file_tag = {
    'Calibration prior': 'calibrationPrior',
    'Uniform MRCA prior': 'uniformMRCA',
}

# One 3x3 figure per prior type: rows are the analyses, columns the calibrated clades.
for prior in prior_types:
    subset_prior = df[df['Prior'] == prior]
    if subset_prior.empty:
        print(f"No runs loaded for {prior}, skipping its figure")
        continue

    fig, axes = plt.subplots(3, 3, figsize=(18, 14), constrained_layout=True)

    for row_idx, analysis in enumerate(analysis_types):
        subset_analysis = subset_prior[subset_prior['Analysis'] == analysis]

        for col_idx, taxon in enumerate(taxon_sets):
            ax = axes[row_idx, col_idx]
            data_to_plot = subset_analysis[subset_analysis['TaxonSet'] == taxon]

            if not data_to_plot.empty:
                sns.histplot(
                    data=data_to_plot,
                    x='Age',
                    hue='Condition',
                    fill=True,
                    palette=colors,
                    alpha=0.4,
                    ax=ax,
                    common_norm=False,
                    element="step",
                    stat="density",
                    linewidth=0.5
                )

                display_name = title_map.get(taxon, taxon)
                ax.set_title(f"{display_name} ({analysis})", fontsize=20, fontweight='bold')
                ax.set_xlabel("Age" if row_idx == 2 else "")

                if (row_idx == 0 and col_idx == 2):
                    sns.move_legend(ax, "upper right", title=None)
                else:
                    if ax.get_legend(): ax.get_legend().remove()
            else:
                ax.set_visible(False)

    fig.suptitle(prior, fontsize=24, fontweight='bold')
    out_pdf = os.path.join(OUT_DIR, f"cats_comparison_plot_{file_tag[prior]}.pdf")
    plt.savefig(out_pdf, format='pdf', bbox_inches='tight')
    plt.close(fig)
    print(f"Comparison plot saved as '{out_pdf}'")

# --- 6. JENSEN-SHANNON DIVERGENCE COMPUTATION ---
jsd_records = []

for prior in prior_types:
    subset_prior = df[df['Prior'] == prior]

    for analysis in analysis_types:
        subset_analysis = subset_prior[subset_prior['Analysis'] == analysis]

        for taxon in taxon_sets:
            display_name = title_map.get(taxon, taxon)
            data_taxon = subset_analysis[subset_analysis['TaxonSet'] == taxon]

            cond_data = data_taxon[data_taxon['Condition'] == 'Conditioned']['Age'].values
            not_cond_data = data_taxon[data_taxon['Condition'] == 'Not Conditioned']['Age'].values

            if len(cond_data) > 1 and len(not_cond_data) > 1:
                min_val = min(cond_data.min(), not_cond_data.min())
                max_val = max(cond_data.max(), not_cond_data.max())
                grid = np.linspace(min_val, max_val, 1000)

                kde_cond = gaussian_kde(cond_data)(grid)
                kde_not_cond = gaussian_kde(not_cond_data)(grid)

                p = kde_cond / np.sum(kde_cond)
                q = kde_not_cond / np.sum(kde_not_cond)

                js_dist = jensenshannon(p, q, base=2.0)
                js_div = js_dist ** 2

                jsd_records.append({
                    'Prior': prior,
                    'Analysis': analysis,
                    'Clade': display_name,
                    'JSD': js_div
                })

df_jsd = pd.DataFrame(jsd_records)

if df_jsd.empty:
    print("No conditioned/not-conditioned pairs available, skipping the JSD barplot")
    raise SystemExit

# --- 7. PLOT THE JSD BARPLOT ---
# One panel per prior type, sharing a y axis so the two are directly comparable.
present_priors = [p for p in prior_types if (df_jsd['Prior'] == p).any()]
fig2, axes2 = plt.subplots(1, len(present_priors), figsize=(12 * len(present_priors), 8),
                           sharey=True, squeeze=False)

# Define distinct colors for the Analysis Types
analysis_colors = {
    'Prior': 'tab:gray',
    'Half Alignment': 'tab:blue',
    'Full Alignment': 'tab:green'
}

for ax2, prior in zip(axes2[0], present_priors):
    sns.barplot(
        data=df_jsd[df_jsd['Prior'] == prior],
        x='Clade',          # <-- Groups the bars by Clade along the X-axis
        y='JSD',
        hue='Analysis',     # <-- Colors the bars by Prior, Half, Full
        hue_order=analysis_types,   # Keeps the colors in chronological order
        order=['Clade 1', 'Clade 2', 'Clade 3'],  # Keeps Clades in order
        palette=analysis_colors,
        ax=ax2,
        edgecolor='black'
    )

    ax2.set_title(prior, fontsize=20, fontweight='bold', pad=15)
    ax2.set_ylabel('JS Divergence', fontsize=18)
    ax2.set_xlabel('Clade', fontsize=18)

    if ax2 is axes2[0][-1]:
        sns.move_legend(ax2, "upper right", title="Analysis Type", title_fontsize=16, fontsize=14)
    elif ax2.get_legend():
        ax2.get_legend().remove()

fig2.suptitle('Jensen-Shannon Divergence: Conditioned vs. Not Conditioned',
              fontsize=22, fontweight='bold')
jsd_pdf = os.path.join(OUT_DIR, "cats_jsd_barplot.pdf")
plt.savefig(jsd_pdf, format='pdf', bbox_inches='tight')

print(f"JSD Barplot saved as '{jsd_pdf}'")
