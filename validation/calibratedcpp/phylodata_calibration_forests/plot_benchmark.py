#!/usr/bin/env python3
"""
Plot benchmark results from CalibratedCPP constraint tree benchmarks
and replace the regression/table info in README.md accurately.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re
from pathlib import Path

# Use script directory for relative paths
SCRIPT_DIR = Path(__file__).parent
CSV_PATH = SCRIPT_DIR / "benchmark_results.csv"
PLOT_PATH = SCRIPT_DIR / "benchmark_analysis_plot.pdf"
README_PATH = SCRIPT_DIR / "README.md"


def update_readme(a, b, c, r2, top_df):
    if not README_PATH.exists():
        return

    content = README_PATH.read_text(encoding='utf-8')

    # 1. Update Regression section (Simple string splitting)
    # Target the formula line specifically
    formula_start = "time(μs) ="
    if formula_start in content:
        parts = content.split(formula_start)
        # split[0] is everything before formula, split[1] is everything after
        # We find the end of the line in the second part
        line_end = parts[1].find("\n")
        new_formula = f" {a:.3f} × complexity + {b:.3f} × taxa + {c:.1f}"
        content = parts[0] + formula_start + new_formula + parts[1][line_end:]

    # Target the R2 line
    r2_anchor = "**Regression model** (R² ="
    if r2_anchor in content:
        parts = content.split(r2_anchor)
        line_end = parts[1].find("):")
        content = parts[0] + r2_anchor + f" {r2:.2f}" + parts[1][line_end:]

    # 2. Build the new table string
    table_header = "| File | Taxa | Calibrations | Complexity | Time (μs) |"
    table_divider = "|------|------|--------------|------------|-----------|"
    rows = [f"| {str(row['file']).replace('.newick', '')} | {int(row['taxa']):,} | {int(row['calibrations'])} | {int(row['complexity']):,} | {row['p50_us']:,.0f} |" for _, row in top_df.iterrows()]
    new_table_body = f"{table_header}\n{table_divider}\n" + "\n".join(rows)

    # 3. Slice and Stitch the Table
    start_anchor = "<!-- start of table -->"
    end_anchor = "<!-- end of table -->"
    
    start_idx = content.find(start_anchor)
    end_idx = content.find(end_anchor)

    if start_idx != -1 and end_idx != -1:
        # We keep everything up to the start anchor + the anchor itself
        prefix = content[:start_idx + len(start_anchor)]
        # We keep everything from the end anchor to the end of file
        suffix = content[end_idx:]
        
        # Stitch it together
        new_content = f"{prefix}\n{new_table_body}\n{suffix}"
        
        README_PATH.write_text(new_content, encoding='utf-8')
        print(f"Successfully updated: {README_PATH}")
    else:
        print(f"Error: Could not find '{start_anchor}' or '{end_anchor}' in README.")


def main():
    df = pd.read_csv(CSV_PATH)
    df = df[df['error'].isna() | (df['error'] == '')]

    print(f"Successful benchmarks: {len(df)}")

    # Fit multivariate: time ~ complexity + taxa
    X = np.column_stack([df['complexity'].values, df['taxa'].values, np.ones(len(df))])
    y = df['p50_us'].values
    coeffs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    a, b, c = coeffs

    y_pred = a * df['complexity'].values + b * df['taxa'].values + c
    ss_res = np.sum((y - y_pred)**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    r2 = 1 - ss_res / ss_tot

    print(f"Regression: time = {a:.5f}*complexity + {b:.3f}*taxa + {c:.1f}")
    print(f"R² = {r2:.4f}")

    # Create figure
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    mask_nonzero = df['complexity'] > 0
    axes[0].scatter(df.loc[mask_nonzero, 'complexity'], df.loc[mask_nonzero, 'p50_us'],
                    alpha=0.6, s=40, label=f'complexity>0 (n={mask_nonzero.sum()})')
    axes[0].scatter(df.loc[~mask_nonzero, 'complexity'] + 0.5, df.loc[~mask_nonzero, 'p50_us'],
                    alpha=0.6, s=40, c='orange', label=f'complexity=0 (n={(~mask_nonzero).sum()})')
    axes[0].set_xscale('log')
    axes[0].set_yscale('log')
    axes[0].set_xlabel('Complexity Score')
    axes[0].set_ylabel('Median Time (μs)')
    axes[0].set_title('Time vs Complexity (Log-Log)')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    # Plot 2: Time vs Taxa
    axes[1].scatter(df['taxa'], df['p50_us'], alpha=0.6, s=40)
    axes[1].set_yscale('log')
    axes[1].set_xscale('log')
    axes[1].set_xlabel('Number of Taxa')
    axes[1].set_ylabel('Median Time (μs)')
    axes[1].set_title('Time vs Taxa Count')
    axes[1].grid(True, alpha=0.3)

    # Plot 3: Actual vs Predicted
    axes[2].scatter(y_pred, y, alpha=0.6, s=40)
    axes[2].plot([0, max(y)], [0, max(y)], 'r--', label='Perfect fit', linewidth=2)
    axes[2].set_xscale('log')
    axes[2].set_yscale('log')
    axes[2].set_xlabel('Predicted Time (μs)')
    axes[2].set_ylabel('Actual Time (μs)')
    axes[2].set_title(f'Multivariate Fit (R²={r2:.4f})')
    axes[2].legend()
    axes[2].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(PLOT_PATH, dpi=150)
    print(f"Saved: {PLOT_PATH}")

    # Prepare Top 5 and update README
    top_5 = df.sort_values(by='p50_us', ascending=False).head(5)
    update_readme(a, b, c, r2, top_5)


if __name__ == "__main__":
    main()