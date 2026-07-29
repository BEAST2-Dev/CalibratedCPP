#!/usr/bin/env python3
"""
BEAST2 Log File Analyzer
Compare multiple runs to check for convergence and consistency
"""

import pandas as pd
import numpy as np
import sys
from pathlib import Path
from typing import Dict, List, Tuple
import warnings
warnings.filterwarnings('ignore')

class BEAST2LogAnalyzer:
    def __init__(self, log_files: List[str]):
        """
        Initialize analyzer with log files
        
        Args:
            log_files: List of paths to .log or .txt files
        """
        self.log_files = log_files
        self.burnin = 0.1  # Fixed at 10%
        self.data = {}
        self.stats = {}
        self.load_logs()
    
    def load_logs(self):
        """Load and parse BEAST2 log files"""
        print(f"Loading {len(self.log_files)} log files...\n")

        for i, log_file in enumerate(self.log_files):
            try:
                df = pd.read_csv(
                    log_file,
                    sep='\t',
                    comment='#',
                    dtype={'Sample': int}
                )

                # Drop TaxonSet columns except mrca.age(TaxonSet*)
                drop = [c for c in df.columns
                        if 'TaxonSet' in c and not c.startswith('mrca.age(')]
                df = df.drop(columns=drop)

                # Apply burnin
                burnin_samples = int(len(df) * self.burnin)
                df = df.iloc[burnin_samples:].reset_index(drop=True)

                self.data[log_file] = df
                print(f"✓ Run {i+1}: {Path(log_file).name}")
                print(f"  Loaded {len(df)} samples (after {self.burnin*100:.0f}% burnin)")

            except Exception as e:
                print(f"✗ Error loading {log_file}: {e}")
    
    def calculate_statistics(self):
        """Calculate summary statistics for each run"""
        print("\n" + "="*80)
        print("SUMMARY STATISTICS BY RUN")
        print("="*80 + "\n")
        
        for log_file, df in self.data.items():
            run_name = Path(log_file).name
            self.stats[log_file] = {}
            
            print(f"\n{run_name}:")
            print("-" * 60)
            
            for col in df.columns:
                if col == 'Sample':
                    continue
                
                try:
                    mean = df[col].mean()
                    std = df[col].std()
                    hpd_low = df[col].quantile(0.025)
                    hpd_high = df[col].quantile(0.975)
                    ess = self.calculate_ess(df[col].values)
                    
                    self.stats[log_file][col] = {
                        'mean': mean,
                        'std': std,
                        'hpd_low': hpd_low,
                        'hpd_high': hpd_high,
                        'ess': ess
                    }
                    
                except:
                    pass
    
    def calculate_ess(self, samples: np.ndarray, max_lag: int = None) -> float:
        """
        Calculate Effective Sample Size using autocorrelation
        """
        if max_lag is None:
            max_lag = min(len(samples) // 2, 500)
        
        n = len(samples)
        mean = np.mean(samples)
        c0 = np.sum((samples - mean) ** 2) / n
        
        if c0 == 0:
            return n
        
        acf = 1.0
        for lag in range(1, max_lag):
            c_lag = np.sum((samples[:-lag] - mean) * (samples[lag:] - mean)) / n
            acf += 2 * (c_lag / c0)
            if acf <= 0:
                break
        
        ess = n / max(acf, 1.0)
        return max(ess, 1)
    
    def compare_runs(self, tolerance: float = 0.05) -> Dict:
        """
        Compare statistics across runs
        
        Args:
            tolerance: Relative tolerance for comparing means (default 5%)
        """
        print("\n" + "="*80)
        print("COMPARISON ACROSS RUNS")
        print("="*80 + "\n")
        
        if len(self.data) < 2:
            print("Need at least 2 runs to compare!")
            return {}
        
        results = {
            'identical': [],
            'similar': [],
            'different': []
        }
        
        # Get all parameters
        all_params = set()
        for stats in self.stats.values():
            all_params.update(stats.keys())
        
        all_params = sorted(list(all_params))
        
        comparison_table = []
        
        for param in all_params:
            means = []
            esses = []
            file_names = []
            
            for log_file, stats_dict in self.stats.items():
                if param in stats_dict:
                    means.append(stats_dict[param]['mean'])
                    esses.append(stats_dict[param]['ess'])
                    file_names.append(Path(log_file).stem)
            
            if len(means) > 1:
                mean_values = np.array(means)
                grand_mean = np.mean(mean_values)

                # Relative difference from grand mean of all replicates
                if grand_mean != 0:
                    rel_diff = np.max(np.abs(mean_values - grand_mean)) / np.abs(grand_mean)
                else:
                    rel_diff = np.max(np.abs(mean_values - grand_mean))

                # Status
                if rel_diff < tolerance * 0.1:
                    status = "✓ IDENTICAL"
                    results['identical'].append(param)
                elif rel_diff < tolerance:
                    status = "≈ SIMILAR"
                    results['similar'].append(param)
                else:
                    status = "✗ DIFFERENT"
                    results['different'].append(param)

                # ESS check
                ess_min = min(esses)
                ess_status = "✓" if ess_min > 200 else "⚠" if ess_min > 50 else "✗"

                comparison_table.append({
                    'Parameter': param,
                    'Status': status,
                    'GrandMean': f"{grand_mean:.4e}",
                    'Rel.Diff': f"{rel_diff*100:.2f}%",
                    'Min_ESS': f"{ess_min:.0f}",
                    'ESS_Status': ess_status
                })
        
        # Print comparison table
        comparison_df = pd.DataFrame(comparison_table)
        print(comparison_df.to_string(index=False))
        
        return results
    
    def print_detailed_comparison(self, tolerance: float = 0.05):
        """Print detailed statistics only for parameters with differences"""
        print("\n" + "="*80)
        print("DETAILED COMPARISON - PARAMETERS WITH DIFFERENCES")
        print("="*80 + "\n")

        all_params = set()
        for stats in self.stats.values():
            all_params.update(stats.keys())

        all_params = sorted(list(all_params))

        params_with_diff = []

        for param in all_params:
            means = []
            for log_file, stats_dict in self.stats.items():
                if param in stats_dict:
                    means.append(stats_dict[param]['mean'])

            if len(means) > 1:
                mean_values = np.array(means)
                grand_mean = np.mean(mean_values)
                if grand_mean != 0:
                    rel_diff = np.max(np.abs(mean_values - grand_mean)) / np.abs(grand_mean)
                else:
                    rel_diff = np.max(np.abs(mean_values - grand_mean))

                # Only show if there are differences
                if rel_diff >= tolerance * 0.1:  # Show if relative diff > 0.5%
                    params_with_diff.append((param, rel_diff))

        if not params_with_diff:
            print("✓ No significant differences found - all parameters converged well!\n")
            return

        for param, rel_diff in params_with_diff:
            print(f"\n{param}  (Rel. Diff: {rel_diff*100:.2f}%)")
            print("-" * 70)

            param_data = []
            for log_file, stats_dict in self.stats.items():
                if param in stats_dict:
                    s = stats_dict[param]
                    run_name = Path(log_file).stem
                    param_data.append({
                        'Run': run_name,
                        'Mean': f"{s['mean']:.6e}",
                        'Std': f"{s['std']:.6e}",
                        'HPD_Low': f"{s['hpd_low']:.6e}",
                        'HPD_High': f"{s['hpd_high']:.6e}",
                        'ESS': f"{s['ess']:.0f}"
                    })

            if param_data:
                param_df = pd.DataFrame(param_data)
                print(param_df.to_string(index=False))

    def generate_report(self, output_file: str = None):
        """Generate a summary report"""
        print("\n" + "="*80)
        print("CONVERGENCE ASSESSMENT")
        print("="*80 + "\n")

        self.calculate_statistics()
        comparison = self.compare_runs()

        print(f"\n\nSummary:")
        print(f"  Identical parameters: {len(comparison['identical'])}")
        print(f"  Similar parameters:   {len(comparison['similar'])}")
        print(f"  Different parameters: {len(comparison['different'])}")

        if comparison['different']:
            print(f"\n⚠️  WARNING: {len(comparison['different'])} parameters differ significantly!")
            print(f"   Consider increasing chain length or checking model specification.")
        else:
            print(f"\n✓ All parameters converged well across runs!")

        self.print_detailed_comparison()

        # Optional: save to file
        if output_file:
            with open(output_file, 'w') as f:
                f.write(f"BEAST2 Log File Comparison Report\n")
                f.write(f"=" * 80 + "\n\n")
                f.write(f"Files analyzed: {len(self.log_files)}\n")
                f.write(f"Burnin: {self.burnin*100:.0f}%\n\n")
                f.write(f"Identical parameters: {len(comparison['identical'])}\n")
                f.write(f"Similar parameters: {len(comparison['similar'])}\n")
                f.write(f"Different parameters: {len(comparison['different'])}\n")

            print(f"\n✓ Report saved to {output_file}")


def main():
    # Default to 'log/' subdirectory next to this script
    script_dir = Path(__file__).parent
    if len(sys.argv) > 1:
        folder = Path(sys.argv[1])
    else:
        folder = script_dir / 'log'

    if not folder.exists():
        print(f"Error: Folder '{folder}' does not exist!")
        print("\nUsage: python checkIdentical.py /path/to/log/files")
        sys.exit(1)

    # b3_r*.log (replicates only, not b3.log) + cats *.txt — the 6 comparison logs
    log_files = sorted(folder.glob('b3_r*.log')) + sorted(folder.glob('*.txt'))

    if not log_files:
        print(f"No matching log files found in '{folder}'!")
        print("\nUsage: python checkIdentical.py /path/to/log/files")
        sys.exit(1)

    log_files = [str(f) for f in log_files]

    print(f"Found {len(log_files)} log files in '{folder}':\n")
    for f in log_files:
        print(f"  - {Path(f).name}")

    # Run analysis
    analyzer = BEAST2LogAnalyzer(log_files)
    analyzer.generate_report()


if __name__ == '__main__':
    main()