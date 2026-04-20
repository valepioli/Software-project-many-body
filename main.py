#main.py
"""
Created on Sat Apr 4 18:15 2026

@author: Pioli Valeria

Main simulation script: Solves the BCS-BEC crossover equations and 
saves the numerical results to a text file.
"""

print("BEC–BCS project started")

import argparse
import numpy as np
import os
from src.config import k_max, N, n_target, EF, start_x, end_x, steps
from src.physics import create_k_grid
from src.solver import solve_bcs_system
from src.plotting import plot_bcs_bec_crossover
from src.plotting import plot_physical_regimes

def main():
    # 1. CLI ARGUMENT PARSING
    parser = argparse.ArgumentParser(description="Simulate the BCS-BEC Crossover in a 3D Fermi Gas.")
    
    parser.add_argument("--k_max", type=float, default=k_max, help="Momentum cutoff")
    parser.add_argument("--n_points", type=int, default=N, help="Number of points in the k-grid")
    parser.add_argument("--start_x", type=float, default=start_x, help="Starting 1/(kF*a) (BCS)")
    parser.add_argument("--end_x", type=float, default=end_x, help="Ending 1/(kF*a) (BEC)")
    parser.add_argument("--steps", type=int, default=steps, help="Number of interaction steps")
    parser.add_argument("--output", type=str, default="results", help="Folder to save results")

    args = parser.parse_args()

    print(f"BEC–BCS project started: k_max={args.k_max}, N={args.n_points}, start_x={args.start_x}, end_x={args.end_x}, Steps={args.steps}")

    #1. SETUP
    os.makedirs(args.output, exist_ok=True)
    #create momentum space k grid
    k, dk = create_k_grid(k_max=args.k_max, N=args.n_points)
    # Define the range for the dimensionless interaction parameter 1/(kF * a)
    interaction_range = np.linspace(args.start_x, args.end_x, args.steps)

    # =============================================================================
    # 2. CROSSOVER SWEEP (ROOT-FINDING LOOP)
    # =============================================================================

    results = []

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

    # --- Save results to folder ---
    output_folder = "results"
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Combine interaction range, mu and delta into one matrix for saving
    data_to_save = np.column_stack((interaction_range, mu_vals, delta_vals))
   
    # Save as a text file
    header = "1/kFa, mu/EF, Delta/EF"
    np.savetxt(f"{output_folder}/crossover_data.txt", data_to_save, header=header, fmt="%.6f", delimiter="\t")
    print(f"Numerical data saved in {output_folder}/crossover_data.txt")


if __name__ == "__main__":
    main()
