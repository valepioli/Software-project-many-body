"""
Created on Tue Apr 24 14:30 2026
@author: Pioli Valeria

Convergence Study Script: 
Verifies how the solution (mu, Delta) behaves when changing numerical parameters.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root
from src.config import EF, n_target
from src.physics import create_k_grid
# Updated import to match your solver.py
from src.solver import solve_bcs_system, bcs_objective_functions

# --- SETUP ---
# Directory where all output plots and reports will be saved
output_dir = "convergence_study"

# Create the directory if it does not already exist
os.makedirs(output_dir, exist_ok=True)

# Physical setup: test the system at unitarity (1/k_F a = 0)
# A very large scattering length 'a' approximates the unitary limit
interaction_a = 1e12  

# Label used for plots and logging
x_label = "Unitarity (1/kFa = 0)"

print(f"--- Starting Convergence Study at {x_label} ---")


# =============================================================================
# TEST 1: Grid Density (N)
# =============================================================================
# Goal: Check how sensitive the solution is to the number of momentum grid points
# Increasing N should lead to convergence of μ and Δ
print("Running Test 1: Grid Density (N)...")

# Different grid resolutions to test (log-spaced manually)
n_list = [500, 1000, 2000, 5000, 10000, 20000, 40000]

# Keep momentum cutoff fixed while varying resolution
k_max_fixed = 100.0

# Store solutions for each N
results_n = []

for n in n_list:
    # Create momentum grid with N points up to k_max
    k, dk = create_k_grid(k_max=k_max_fixed, N=n)

    # Solve the BCS system for given grid
    # initial_guess = [chemical potential μ, gap Δ]
    sol = solve_bcs_system(
        interaction_a, n_target, k, dk,
        initial_guess=[EF, 0.5*EF]
    )

    # Store solution (μ, Δ)
    results_n.append(sol)

# Convert results to NumPy array for easier slicing
results_n = np.array(results_n)

# Normalize results by Fermi energy EF (dimensionless quantities)
mu_n = results_n[:, 0] / EF
delta_n = results_n[:, 1] / EF

# Plot convergence with respect to N
plt.figure(figsize=(8, 5))

# Chemical potential
plt.plot(n_list, mu_n, 'b-o', label=r'$\mu / E_F$')

# Gap parameter
plt.plot(n_list, delta_n, 'r-s', label=r'$\Delta / E_F$')

# Log scale helps visualize convergence across orders of magnitude
plt.xscale('log')

plt.xlabel('Number of points (N)')
plt.ylabel('Dimensionless Energy')

plt.title(f'Test 1: Convergence vs Grid Density ($k_{{max}}={k_max_fixed}$)')

# Grid styling for readability
plt.grid(True, which="both", ls="-", alpha=0.2)

plt.legend()

# Save plot to output directory
plt.savefig(os.path.join(output_dir, "convergence_N.png"))

plt.close()


# =============================================================================
# TEST 2: Momentum Cutoff (k_max)
# =============================================================================
# Goal: Check sensitivity to the upper momentum cutoff
# k_max must be large enough to capture all relevant physics
print("Running Test 2: Momentum Cutoff (k_max)...")

# Different cutoff values to test
k_max_list = [10, 50, 100, 200, 500, 1000]

# Keep grid density fixed while varying cutoff
n_fixed = 20000 

results_k = []

for km in k_max_list:
    # Create grid with fixed resolution but varying cutoff
    k, dk = create_k_grid(k_max=km, N=n_fixed)

    # Solve BCS equations
    sol = solve_bcs_system(
        interaction_a, n_target, k, dk,
        initial_guess=[EF, 0.5*EF]
    )

    results_k.append(sol)

# Convert to array
results_k = np.array(results_k)

# Normalize results
mu_k = results_k[:, 0] / EF
delta_k = results_k[:, 1] / EF

# Plot convergence with respect to k_max
plt.figure(figsize=(8, 5))

plt.plot(k_max_list, mu_k, 'b-o', label=r'$\mu / E_F$')
plt.plot(k_max_list, delta_k, 'r-s', label=r'$\Delta / E_F$')

plt.xlabel('Momentum Cutoff ($k_{max}$)')
plt.ylabel('Dimensionless Energy')

plt.title(f'Test 2: Convergence vs Momentum Cutoff ($N={n_fixed}$)')

plt.grid(True, alpha=0.3)
plt.legend()

plt.savefig(os.path.join(output_dir, "convergence_kmax.png"))

plt.close()


# =============================================================================
# TEST 3: Solver Residuals
# =============================================================================
# Goal: Verify that the numerical solver truly satisfies the equations
# by inspecting the residuals of the gap and number equations
print("Running Test 3: Solver Residuals...")

# Use a reasonably converged grid
k, dk = create_k_grid(k_max=100, N=10000)

# Slightly different initial guess to test solver robustness
initial_guess = [EF, 0.6*EF]

# Solve using root-finding directly on the objective functions
sol_root = root(
    bcs_objective_functions,
    initial_guess,
    args=(interaction_a, n_target, k, dk)
)

# Extract residuals of the two equations
res_gap, res_num = sol_root.fun

print(f"--- Solver Status ---")
print(f"Success: {sol_root.success}")

# Residuals should be close to zero if solution is accurate
print(f"Residual Gap Equation: {res_gap:.2e}")
print(f"Residual Number Equation: {res_num:.2e}")


# =============================================================================
# SAVE REPORT
# =============================================================================
# Write a summary of results to a text file
report_path = os.path.join(output_dir, "convergence_report.txt")

with open(report_path, "w") as f:
    f.write("CONVERGENCE STUDY REPORT\n")
    f.write("========================\n\n")

    f.write(f"Test at Unitarity (1/kFa = 0)\n")

    # Final solution (dimensionless form)
    f.write(
        f"Final Solution: mu/EF = {sol_root.x[0]/EF:.4f}, "
        f"Delta/EF = {sol_root.x[1]/EF:.4f}\n"
    )

    # Residuals as a measure of numerical accuracy
    f.write(
        f"Solver Residuals: Gap={res_gap:.2e}, "
        f"Number={res_num:.2e}\n"
    )

print(f"--- Study Complete. Results in '{output_dir}/' ---")