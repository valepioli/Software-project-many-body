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

    plt.close(fig)


def plot_physical_regimes(interaction_range, mu_vals, delta_vals, save_path=None):
    """
    Generates an infographic-style plot of the BCS-BEC crossover.
    X-axis: ratio of Chemical Potential to Gap (mu / Delta)
    Y-axis: normalized binding energy (epsilon_B / EF)
    """
    
    # 1. PHYSICAL CONSTANTS
    # In units where m=1, EF = 0.5 * kF^2. 
    # The vacuum binding energy epsilon_B = 1 / (m * a^2).
    # Thus, epsilon_B / EF = (1/a^2) / (kF^2 / 2) = 2 * (1 / kFa)^2.
    # Note: Binding energy in vacuum only exists for 1/kFa > 0 (BEC side).
    eb_over_ef = np.where(interaction_range > 0, 2.0 * (interaction_range)**2, 0)
    
    # 2. X-AXIS RATIO
    # The image plots the trajectory in the (mu/Delta) space.
    ratio_mu_delta = mu_vals / delta_vals

    plt.figure(figsize=(10, 7))

    # 3. CROSSOVER GRADIENT (The green 'glow' area)
    # We create a smooth green gradient centered at mu/Delta = 0
    for i in range(60):
        # Gaussian alpha for the glow effect
        alpha_val = 0.4 * np.exp(- (i-30)**2 / 150)
        plt.axvspan(-1.5 + i*0.05, -1.5 + (i+1)*0.05, color='green', alpha=alpha_val, zorder=1)

    # 4. MAIN CURVE
    # Black line showing the relation between binding energy and the mu/Delta ratio
    plt.plot(ratio_mu_delta, eb_over_ef, color='black', linewidth=2.5, zorder=5)

    # 5. REGIME LABELS (Rounded boxes)
    bbox_style = dict(boxstyle="round,pad=0.6", fc="white", ec="black", lw=1)
    # Positions match the three main regimes: BEC, Crossover, and BCS
    plt.text(-3.8, 18.5, "BEC", ha="center", bbox=bbox_style, fontsize=12, fontweight='bold', color='navy')
    plt.text(0, 18.5, "BCS-BEC\ncrossover", ha="center", bbox=bbox_style, fontsize=10, fontweight='bold', color='red')
    plt.text(3.8, 18.5, "BCS", ha="center", bbox=bbox_style, fontsize=12, fontweight='bold', color='navy')

    # 6. SCHEMATIC PAIRS (Visual representation of atoms)
    
    # BEC Side: Tightly bound, localized pairs
    plt.scatter([-4.2, -4.0, -3.2, -3.4], [13.5, 14.2, 12.5, 8.5], c='blue', s=80, zorder=10, edgecolors='gray')
    plt.scatter([-4.1, -3.9, -3.1, -3.3], [13.2, 13.9, 12.2, 8.2], c='red', s=80, zorder=10, edgecolors='gray')
    
    # BCS Side: Large, overlapping Cooper pairs (represented with transparent lines)
    for i in range(5):
        rx, ry = np.random.uniform(2.5, 4.5), np.random.uniform(8, 15)
        # Large faint bonds representing delocalized pairs
        plt.plot([rx, rx+1.2], [ry, ry-0.6], color='red', alpha=0.15, lw=8, solid_capstyle='round')
        plt.scatter([rx, rx+1.2], [ry, ry-0.6], c=['red', 'blue'], s=70, alpha=0.7, zorder=10)

    # 7. AXIS FORMATTING
    plt.xlabel(r"$\mu / \Delta$", fontsize=14)
    plt.ylabel(r"$\epsilon_B / \epsilon_F$", fontsize=14)
    
    # Limits set to match the reference image exactly
    plt.xlim(-5, 5)
    plt.ylim(0, 20)
    
    plt.axhline(0, color='black', lw=1.2)
    plt.grid(alpha=0.15, linestyle=':')
    plt.tight_layout()

    # 8. SAVE OUTPUT
    if save_path:
        out_file = os.path.join(save_path, "regimes_infographic.png")
        plt.savefig(out_file, dpi=300)
        print(f"Infographic saved to {out_file}")

