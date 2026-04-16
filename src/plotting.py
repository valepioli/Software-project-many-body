# plotting.py
"""
Created on Sun Apr 12 19:25 2026

@author: Pioli Valeria
"""

import matplotlib.pyplot as plt
import numpy as np
import os 

def plot_bcs_bec_crossover(interaction_range, mu_vals, delta_vals, save_path=None):
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

    fig, ax = plt.subplots(figsize=(10, 7))

    # 1. Plot Numerical Results
    ax.plot(interaction_range, mu_vals, 'b-o', label=r'Calculated $\mu / E_F$', 
            markersize=4, linewidth=1.5, zorder=3)
    ax.plot(interaction_range, delta_vals, 'r-s', label=r'Calculated $\Delta / E_F$', 
            markersize=4, linewidth=1.5, zorder=3)
    # 2. Benchmarks at Unitarity (x = 0)
    # We use short horizontal lines (hlines) as "target markers"
    bench_width = 0.15  # Width of the indicator bars
    
    # --- Mean Field Benchmarks (Theory) ---
    mf_mu, mf_delta = 0.59, 0.68
    ax.hlines([mf_mu, mf_delta], -bench_width, bench_width, colors='black', 
              linestyles='--', linewidth=1.2, alpha=0.8, zorder=4, label='Mean Field Theory')

    # --- Monte Carlo Benchmarks (Experimental/Refined) ---
    mc_mu, mc_delta = 0.37, 0.44
    ax.hlines([mc_mu, mc_delta], -bench_width, bench_width, colors='green', 
              linestyles='--', linewidth=1.2, alpha=0.8, zorder=4, label='Monte Carlo / Exp.')

    # 3. Annotations for the target bars
    # Using small text labels near the bars for immediate identification
    ax.text(bench_width+0.05, mf_mu, 'MF', color='black', va='center', fontsize=8, fontweight='bold')
    ax.text(bench_width+0.05, mf_delta, 'MF', color='black', va='center', fontsize=8, fontweight='bold')
    ax.text(bench_width+0.05, mc_mu, 'MC', color='green', va='center', fontsize=8, fontweight='bold')
    ax.text(bench_width+0.05, mc_delta, 'MC', color='green', va='center', fontsize=8, fontweight='bold')

    # 4. Reference Axes
    ax.axhline(0, color='black', linewidth=0.8, alpha=0.3)
    ax.axvline(0, color='gray', linestyle=':', linewidth=1.0, alpha=0.5)

    # 5. Styling and Labels
    ax.set_xlabel(r'Interaction Strength $1/(k_F a)$', fontsize=12)
    ax.set_ylabel(r'Energy [$\epsilon_F$ units]', fontsize=12)
    ax.set_title('BCS-BEC Crossover: Calculated vs Literature Target Values', fontsize=14)
    
    ax.grid(True, which='both', linestyle=':', alpha=0.3)
    ax.legend(loc='upper right', fontsize=9, frameon=True)
    
    # Text indications for regimes
    ax.annotate('BCS Regime', xy=(0.15, 0.75), xycoords='axes fraction', color='black', alpha=0.6)
    ax.annotate('BEC Regime', xy=(0.75, 0.75), xycoords='axes fraction', color='black', alpha=0.6)

    ax.set_ylim(-1.5, 2.0)
    ax.set_xlim(min(interaction_range), max(interaction_range))
    
    plt.tight_layout()

    # --- Save the plot ---
    if save_path:
        file_name = os.path.join(save_path, "crossover_plot.png")
        plt.savefig(file_name, dpi=300)
        print(f"Plot saved in {file_name}")

    plt.show()
