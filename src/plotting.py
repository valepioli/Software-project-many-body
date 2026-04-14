# plotting.py
"""
Created on Sun Apr 12 19:25 2026

@author: Pioli Valeria
"""

import matplotlib.pyplot as plt
import numpy as np

def plot_bcs_bec_crossover(interaction_range, mu_vals, delta_vals):
    """
    Generates  plot of mu and delta.

    Parameters:
    -----------
    interaction_range : array
        The values of 1/(kF * a) for the x-axis.
    mu_vals : array
        The solved chemical potential values.
    delta_vals : array
        The solved pairing gap values.
    """

    plt.figure(figsize=(10, 6))

    # Plot Mu and Delta
    plt.plot(interaction_range, mu_vals, 'o-', label=r'Chemical Potential $\mu / \epsilon_F$', markersize=4)
    plt.plot(interaction_range, delta_vals, 's-', label=r'Pairing Gap $\Delta / \epsilon_F$', markersize=4)

    # Reference lines for clarity
    plt.axhline(0, color='black', linestyle='--', linewidth=0.8, alpha=0.5)
    plt.axvline(0, color='gray', linestyle=':', label='Unitarity Limit')

    # Labeling
    plt.xlabel(r'Interaction Strength $1/(k_F a)$', fontsize=12)
    plt.ylabel(r'Energy [$\epsilon_F$ units]', fontsize=12)
    plt.title('BCS-BEC Crossover: Numerical Solution', fontsize=14)

    # Scientific formatting
    plt.grid(True, which='both', linestyle='--', alpha=0.6)
    plt.legend(fontsize=10)

    plt.tight_layout()
    # --- Save the plot ---
    if save_path:
        file_name = os.path.join(save_path, "crossover_plot.png")
        plt.savefig(file_name, dpi=300)
        print(f"Plot saved in {file_name}")

    plt.show()
