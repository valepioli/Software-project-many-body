"""
Created on Mon 20 Apr 17:35 2026
@author: Pioli Valeria

Analysis script: Loads simulation data, processes units, 
and generates visualizations.
"""

import numpy as np
import argparse
import os
from src.config import EF
from src.plotting import plot_bcs_bec_crossover, plot_physical_regimes

def main():
    # 1. CLI ARGUMENT PARSING
    parser = argparse.ArgumentParser(description="Analyze and plot BCS-BEC crossover data.")
    
    parser.add_argument("--data_dir", type=str, default="results", 
                        help="Folder where the simulation data is stored (default: results)")
    
    parser.add_argument("--data_file", type=str, default="crossover_data.txt", 
                        help="Name of the numerical data file (default: crossover_data.txt)")
    
    parser.add_argument("--output_dir", type=str, default="images", 
                        help="Folder where plots will be saved (default: images)")

    args = parser.parse_args()

    # 2. PATH CONSTRUCTION
    data_path = os.path.join(args.data_dir, args.data_file)
    
    # Check if data exists
    if not os.path.exists(data_path):
        print(f"Error: Data file '{data_path}' not found.")
        print("Please run main.py first or check your --data_dir and --data_file arguments.")
        sys.exit(1)

    # Ensure output directory exists
    os.makedirs(args.output_dir, exist_ok=True)

    # 3. DATA LOADING
    print(f"Loading data from: {data_path}...")
    try:
        # Expected format: [1/kFa, mu, Delta]
        data = np.loadtxt(data_path)
        interaction_range = data[:, 0]
        mu_vals = data[:, 1]
        delta_vals = data[:, 2]
    except Exception as e:
        print(f"Error reading file: {e}")
        sys.exit(1)

    # 4. DATA ANALYSIS & NORMALIZATION
    # Convert physical values to dimensionless units relative to Fermi Energy
    mu_normalized = mu_vals / EF
    delta_normalized = delta_vals / EF

    # 5. VISUALIZATION
    print(f"Generating plots and saving to '{args.output_dir}/'...")

    # Plot 1: Standard Crossover (Uses normalized values)
    plot_bcs_bec_crossover(interaction_range, mu_normalized, delta_normalized, 
                           save_path=args.output_dir)

    # Plot 2: Physical Regimes Infographic (Uses raw values for internal binding energy ratios)
    plot_physical_regimes(interaction_range, mu_vals, delta_vals, 
                          save_path=args.output_dir)

    print("Success: All visualizations updated.")

if __name__ == "__main__":
    main()
