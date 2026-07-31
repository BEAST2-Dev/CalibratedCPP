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

# --- 2. CONFIGURATION ---
target_columns = [
    'mrca.age(TaxonSet1)', 
    'mrca.age(TaxonSet2)', 
    'mrca.age(TaxonSet3)'
]

file_map = [
    # PRIOR FILES
    ('cats-sample_from_prior.log', 'Prior', 'Conditioned'),
    ('cats-sample_from_prior_not_conditioned.log', 'Prior', 'Not Conditioned'),
    
    # HALF FILES
    ('cats-half.log', 'Half Alignment', 'Conditioned'),
    ('cats-half_not_conditioned.log', 'Half Alignment', 'Not Conditioned'),
    
    # FULL FILES 
    ('cats-full.log', 'Full Alignment', 'Conditioned'),
    ('cats-full_not_conditioned.log', 'Full Alignment', 'Not Conditioned')
]

burnin_fraction = 0.1

# --- 3. HELPER FUNCTION ---
def load_all_data(file_map):
    all_data = []
    
    for filename, analysis_type, condition_type in file_map:
        if os.path.exists(filename):
            try:
                df = pd.read_csv(filename, sep='\t', comment='#')
                cols = [c for c in target_columns if c in df.columns]
                if cols:
                    drop_n = int(len(df) * burnin_fraction)
                    df = df.iloc[drop_n:].copy()
                    
                    melted = df[cols].melt(var_name='TaxonSet', value_name='Age')
                    melted['Analysis'] = analysis_type
                    melted['Condition'] = condition_type
                    
                    melted['TaxonSet'] = melted['TaxonSet'].str.replace('mrca.age(', '', regex=False).str.replace(')', '', regex=False)
                    all_data.append(melted)
                    print(f"Loaded: {filename}")
                else:
                    print(f"Warning: Targets not found in {filename}")
            except Exception as e:
                print(f"Error reading {filename}: {e}")
        else:
            print(f"File not found: {filename}")

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

fig, axes = plt.subplots(3, 3, figsize=(18, 14), constrained_layout=True)

taxon_sets = [c.replace('mrca.age(', '').replace(')', '') for c in target_columns]
analysis_types = ['Prior', 'Half Alignment', 'Full Alignment']
colors = {'Conditioned': 'tab:blue', 'Not Conditioned': 'tab:orange'}

title_map = {
    'TaxonSet1': 'Clade 1',
    'TaxonSet2': 'Clade 2',
    'TaxonSet3': 'Clade 3'
}

for row_idx, analysis in enumerate(analysis_types):
    subset_analysis = df[df['Analysis'] == analysis]
    
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

plt.savefig("cats_comparison_plot.pdf", format='pdf', bbox_inches='tight')

# --- 6. JENSEN-SHANNON DIVERGENCE COMPUTATION ---
jsd_records = []

for analysis in analysis_types:
    subset_analysis = df[df['Analysis'] == analysis]
    
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
                'Analysis': analysis,
                'Clade': display_name,
                'JSD': js_div
            })

df_jsd = pd.DataFrame(jsd_records)

# --- 7. PLOT THE JSD BARPLOT ---
fig2, ax2 = plt.subplots(figsize=(12, 8))

# Define distinct colors for the Analysis Types
analysis_colors = {
    'Prior': 'tab:gray', 
    'Half Alignment': 'tab:blue', 
    'Full Alignment': 'tab:green'
}

sns.barplot(
    data=df_jsd,
    x='Clade',          # <-- Groups the bars by Clade along the X-axis
    y='JSD',
    hue='Analysis',     # <-- Colors the bars by Prior, Half, Full
    hue_order=analysis_types, # Keeps the colors in chronological order
    order=['Clade 1', 'Clade 2', 'Clade 3'], # Keeps Clades in order
    palette=analysis_colors,
    ax=ax2,
    edgecolor='black'
)

# Formatting the plot
ax2.set_title('Jensen-Shannon Divergence:\nConditioned vs. Not Conditioned', fontsize=20, fontweight='bold', pad=15)
ax2.set_ylabel('JS Divergence', fontsize=18)
ax2.set_xlabel('Clade', fontsize=18)

# Customize the legend to reflect the new Hue
sns.move_legend(ax2, "upper right", title="Analysis Type", title_fontsize=16, fontsize=14)

plt.savefig("cats_jsd_barplot.pdf", format='pdf', bbox_inches='tight')

print("JSD Barplot saved as 'cats_jsd_barplot.pdf'")