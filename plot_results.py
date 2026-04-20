"""
Created on Mon 20 Apr 17:35 2026
@author: Pioli Valeria

Analysis script: Loads simulation data, processes units, 
and generates visualizations.
"""

import numpy as np
import os
from src.config import EF
from src.plotting import plot_bcs_bec_crossover, plot_physical_regimes

def main():
    data_path = "results/crossover_data.txt"
    images_dir = "images"
    
    if not os.path.exists(data_path):
        print(f"Error: {data_path} not found. Run main.py first.")
        return

    os.makedirs(images_dir, exist_ok=True)

    print(f"Loading data from {data_path}...")
    data = np.loadtxt(data_path)
    
    interaction_range = data[:, 0]
    mu_vals = data[:, 1]
    delta_vals = data[:, 2]

    # --- Analysis & Normalization ---
    # Convert physical values to dimensionless units relative to Fermi Energy
    mu_normalized = mu_vals / EF
    delta_normalized = delta_vals / EF

    print("Generating plots...")

    # Plot 1: Standard Crossover (Normalized)
    plot_bcs_bec_crossover(interaction_range, mu_normalized, delta_normalized, save_path=images_dir)

    # Plot 2: Physical Regimes Infographic (Uses raw values for internal ratios)
    plot_physical_regimes(interaction_range, mu_vals, delta_vals, save_path=images_dir)

    print(f"All plots saved in the '{images_dir}/' folder.")

if __name__ == "__main__":
    main()
