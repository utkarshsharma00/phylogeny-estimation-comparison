#!/usr/bin/env python3
"""
Python script for summarizing tree comparison results
Usage: python summarize_results.py <results_directory>
"""

import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def parse_comparison_file(file_path):
    """Parse a comparison file and extract RF, FN, and FP values."""
    results = {}
    try:
        with open(file_path, 'r') as f:
            for line in f:
                if "RF distance:" in line:
                    value = line.split(":")[1].strip()
                    results['RF'] = float(value) if value != 'N/A' else np.nan
                elif "FN rate:" in line:
                    value = line.split(":")[1].strip()
                    results['FN'] = float(value) if value != 'N/A' else np.nan
                elif "FP rate:" in line:
                    value = line.split(":")[1].strip()
                    results['FP'] = float(value) if value != 'N/A' else np.nan
    except Exception as e:
        print(f"Error parsing {file_path}: {e}")
        results = {'RF': np.nan, 'FN': np.nan, 'FP': np.nan}
    return results

def main():
    if len(sys.argv) != 2:
        print("Usage: python summarize_results.py <results_directory>")
        sys.exit(1)
    
    results_dir = sys.argv[1]
    
    # Create a dataframe to store results
    columns = ['Model', 'Method', 'Replicate', 'RF', 'FN', 'FP']
    results_df = pd.DataFrame(columns=columns)
    
    # Methods to analyze
    methods = {
        'fasttree_gtr': 'FastTree (GTR)',
        'fasttree_jc': 'FastTree (JC)',
        'nj_jc': 'NJ (JC)',
        'nj_logdet': 'NJ (LogDet)',
        'nj_pdist': 'NJ (p-distance)'
    }
    
    # Process each model condition
    for model in ['1000M1', '1000M4']:
        model_dir = os.path.join(results_dir, model)
        if not os.path.exists(model_dir):
            print(f"Warning: Model directory {model_dir} does not exist")
            continue
            
        # Process each replicate
        for rep_dir in os.listdir(model_dir):
            rep_path = os.path.join(model_dir, rep_dir)
            if not os.path.isdir(rep_path):
                continue
                
            # Extract replicate number
            rep = rep_dir  # e.g., "R1"
            
            # Process each method
            for method_key, method_name in methods.items():
                comparison_file = os.path.join(rep_path, f'{method_key}_comparison.txt')
                
                if os.path.exists(comparison_file):
                    results = parse_comparison_file(comparison_file)
                    
                    # Add to dataframe
                    new_row = {
                        'Model': model,
                        'Method': method_name,
                        'Replicate': rep,
                        'RF': results.get('RF', np.nan),
                        'FN': results.get('FN', np.nan),
                        'FP': results.get('FP', np.nan)
                    }
                    results_df = pd.concat([results_df, pd.DataFrame([new_row])], ignore_index=True)
    
    # Save raw results to CSV
    raw_file = os.path.join(results_dir, 'raw_results.csv')
    results_df.to_csv(raw_file, index=False)
    print(f"Raw results saved to {raw_file}")
    
    # Calculate average metrics for each method and model
    summary_df = results_df.groupby(['Model', 'Method']).agg({
        'FN': ['mean', 'std'],
        'FP': ['mean', 'std']
    }).reset_index()
    
    # Save summary to CSV
    summary_file = os.path.join(results_dir, 'summary_results.csv')
    summary_df.to_csv(summary_file)
    print(f"Summary results saved to {summary_file}")
    
    # Create plots
    create_plots(results_df, results_dir)

def create_plots(df, results_dir):
    """Create plots to visualize results."""
    # Skip plotting if no valid data
    if df.empty or df['FN'].isna().all() or df['FP'].isna().all():
        print("Warning: Not enough valid data for plotting")
        return
        
    # Average FN rates by method and model
    plt.figure(figsize=(12, 6))
    
    # Group by model and method, then calculate mean FN
    fn_means = df.groupby(['Model', 'Method'])['FN'].mean().reset_index()
    
    # Pivot to create a table suitable for grouped bar chart
    fn_pivot = fn_means.pivot(index='Method', columns='Model', values='FN')
    
    # Plot grouped bar chart
    ax = fn_pivot.plot(kind='bar', capsize=4)
    
    plt.title('Average False Negative Rates by Method and Model')
    plt.ylabel('False Negative Rate')
    plt.xlabel('Method')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(results_dir, 'fn_rates_by_method.png'))
    
    # Average FP rates by method and model
    plt.figure(figsize=(12, 6))
    
    # Group by model and method, then calculate mean FP
    fp_means = df.groupby(['Model', 'Method'])['FP'].mean().reset_index()
    
    # Pivot to create a table suitable for grouped bar chart
    fp_pivot = fp_means.pivot(index='Method', columns='Model', values='FP')
    
    # Plot grouped bar chart
    ax = fp_pivot.plot(kind='bar', capsize=4)
    
    plt.title('Average False Positive Rates by Method and Model')
    plt.ylabel('False Positive Rate')
    plt.xlabel('Method')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(results_dir, 'fp_rates_by_method.png'))

if __name__ == "__main__":
    main()