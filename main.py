#main.py
"""
Created on Sat Apr 4 18:15 2026

@author: Pioli Valeria
"""

print("BEC–BCS project started")

import numpy as np
import os
from src.config import k_max, N, n_target, EF
from src.physics import create_k_grid
from src.solver import solve_bcs_system
from src.plotting import plot_bcs_bec_crossover
from src.plotting import plot_physical_regimes

#1. SETUP
# k_max and N are now taken from config.py
k, dk = create_k_grid(k_max=k_max, N=N)
interaction_range = np.linspace(-2.0, 2.0, 40)

# =============================================================================
# 2. CROSSOVER SWEEP (ROOT-FINDING LOOP)
# =============================================================================

results = []

# Define the range for the dimensionless interaction parameter 1/(kF * a)
# Usually from -2 (BCS) to +2 (BEC)
interaction_range = np.linspace(-3.0, 3.0, 40) 
# initial_guess: We start with a guess close to the BCS limit (mu ~ EF, Delta small).
# This guess will be updated dynamically to improve convergence.
current_guess = [EF, 0.1]

for x in interaction_range:
    """
    We iterate over the dimensionless interaction parameter x = 1/(kF * a).
    BCS regime: x < 0 | Unitarity: x = 0 | BEC regime: x > 0
    """

    # 1. Map dimensionless parameter 'x' to scattering length 'a'.
    # For x = 0 (unitary limit), we use a large number (1e6) to approximate
    # the physical infinity 1/a -> 0.
    a = 1.0 / x if abs(x) > 1e-12 else 1e12

    # 2. Call the solver to find the [mu, Delta] pair that satisfies the system.
    sol = solve_bcs_system(a, n_target, k, dk, initial_guess=current_guess)

    results.append(sol)

    # 3. Update the initial guess.
    # If the solver was successful (not NaN), we use the current solution as the
    # starting point for the next interaction step. This ensures numerical
    # stability as the chemical potential and gap evolve smoothly across the crossover.
    if not np.isnan(sol[0]):
        current_guess = sol

# 3. DATA PREPROCESSING
# Separate the results into two arrays for plotting
results_array = np.array(results)
# NORMALIZATION: Divide by EF to get dimensionless units (mu/EF and Delta/EF)
# This ensures the BCS limit starts at 1.0, which is physically correct.
mu_vals= results_array[:, 0] 
delta_vals = results_array[:, 1] 

mu_vals_normalized = results_array[:, 0] / EF
delta_vals_normalized = results_array[:, 1] / EF

# --- Save results to folder ---
output_folder = "results"
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Combine interaction range, mu and delta into one matrix for saving
data_to_save = np.column_stack((interaction_range, mu_vals_normalized, delta_vals_normalized))

# Save as a text file
header = "1/kFa, mu/EF, Delta/EF"
np.savetxt(f"{output_folder}/crossover_data.txt", data_to_save, header=header, fmt="%.6f", delimiter="\t")
print(f"Numerical data saved in {output_folder}/crossover_data.txt")

# =============================================================================
# 4. VISUALIZATION
# =============================================================================
# Pass the folder path to the plot function to save the image
plot_bcs_bec_crossover(interaction_range, mu_vals_normalized, delta_vals_normalized, save_path=output_folder)
plot_physical_regimes(interaction_range, mu_vals, delta_vals)
