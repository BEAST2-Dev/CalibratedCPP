#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# 1. Setup Relative Paths
SCRIPT_DIR = Path(__file__).parent
BENCH_CSV = SCRIPT_DIR / "benchmark_results.csv"
COMP_CSV = SCRIPT_DIR / "validation_results.csv"
README_PATH = SCRIPT_DIR / "README.md"

# Define separate output paths as PDFs
PLOT_BENCH_PATH = SCRIPT_DIR / "benchmark_plot.pdf"
PLOT_COMP_PATH = SCRIPT_DIR / "validation_plot.pdf"

# Style configuration
COLORS = {"Calibrated CPP": "black", "Heled & Drummond": "red"}

def update_readme(df):
    if not README_PATH.exists():
        return

    content = README_PATH.read_text(encoding='utf-8')

    # 1. Clean and Prepare Data
    df_clean = df[['Benchmark', 'Param: nCalibrations', 'Score']].copy()
    df_clean['Benchmark'] = df_clean['Benchmark'].apply(
        lambda x: "Calibrated CPP" if "measureCPP" in x else "Heled & Drummond"
    )

    # 2. Pivot to horizontal format
    pivot_df = df_clean.pivot(index='Param: nCalibrations', columns='Benchmark', values='Score')
    pivot_df = pivot_df.reset_index()
    
    # Ensure "Number of Calibrations" is an integer to prevent "1.00"
    pivot_df['Param: nCalibrations'] = pivot_df['Param: nCalibrations'].astype(int)
    
    # Rename columns for the final header
    pivot_df.columns = [
        "Number of Calibrations", 
        "Time (μs) Calibrated CPP", 
        "Time (μs) Heled & Drummond"
    ]

    # 3. Generate Markdown Table
    md_table = pivot_df.to_markdown(index=False, tablefmt="github", floatfmt=("", ".2f", ".2f"))

    # 4. Slice and Stitch
    start_anchor = "<!-- start table -->"
    end_anchor = "<!-- end table -->"
    start_idx = content.find(start_anchor)
    end_idx = content.find(end_anchor)

    if start_idx != -1 and end_idx != -1:
        prefix = content[:start_idx + len(start_anchor)]
        suffix = content[end_idx:]
        new_content = f"{prefix}\n\n{md_table}\n\n{suffix}"
        README_PATH.write_text(new_content, encoding='utf-8')
        print("README.md table updated with integer calibrations.")

def main():
    if not BENCH_CSV.exists() or not COMP_CSV.exists():
        print("Error: Required CSV files are missing.")
        return

    # Load Data
    df_bench = pd.read_csv(BENCH_CSV)
    df_comp = pd.read_csv(COMP_CSV)

    # ==========================================
    # Plot 1: Benchmark (Time)
    # ==========================================
    fig1, ax1 = plt.subplots(figsize=(6, 5))

    for label, group in df_bench.groupby('Benchmark'):
        name = "Calibrated CPP" if "measureCPP" in label else "Heled & Drummond"
        ax1.plot(group['Param: nCalibrations'], group['Score'], 
                 label=name, color=COLORS[name], marker='o', linewidth=1.5)
    
    ax1.set_yscale('log')
    ax1.set_xlabel('Number of calibrations (k)')
    ax1.set_ylabel('Time (μs)')
    ax1.set_title('Likelihood Calculation Time')
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='best', frameon=True)

    fig1.tight_layout()
    # Save as PDF (dpi is ignored for vector graphics like PDFs, but safe to leave)
    fig1.savefig(PLOT_BENCH_PATH)
    plt.close(fig1)
    print(f"Benchmark plot saved to {PLOT_BENCH_PATH}")

    # ==========================================
    # Plot 2: Likelihood (Validation)
    # ==========================================
    fig2, ax2 = plt.subplots(figsize=(6, 5))

    ax2.plot(df_comp['BirthRate'], df_comp['HeledAndDrummond_LogLikelihood'], 
             color=COLORS["Heled & Drummond"], label="Heled & Drummond", linewidth=1.5)
    ax2.scatter(df_comp['BirthRate'], df_comp['CPP_LogLikelihood'], 
                color=COLORS["Calibrated CPP"], label="Calibrated CPP", s=20, zorder=3)
    
    ax2.set_xlabel('Birth-rate')
    ax2.set_ylabel('Log-likelihood')
    ax2.set_title('Likelihood Validation')
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='best', frameon=True)

    fig2.tight_layout()
    fig2.savefig(PLOT_COMP_PATH)
    plt.close(fig2)
    print(f"Validation plot saved to {PLOT_COMP_PATH}")

    # Update Readme
    update_readme(df_bench)

if __name__ == "__main__":
    main()