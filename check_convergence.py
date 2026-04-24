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
output_dir = "convergence_study"
os.makedirs(output_dir, exist_ok=True)

# Test at Unitarity (1/kFa = 0)
interaction_a = 1e12  
x_label = "Unitarity (1/kFa = 0)"

print(f"--- Starting Convergence Study at {x_label} ---")

# =============================================================================
# TEST 1: Grid Density (N)
# =============================================================================
print("Running Test 1: Grid Density (N)...")
n_list = [500, 1000, 2000, 5000, 10000, 20000, 40000]
k_max_fixed = 100.0
results_n = []

for n in n_list:
    k, dk = create_k_grid(k_max=k_max_fixed, N=n)
    sol = solve_bcs_system(interaction_a, n_target, k, dk, initial_guess=[EF, 0.5*EF])
    results_n.append(sol)

results_n = np.array(results_n)
mu_n = results_n[:, 0] / EF
delta_n = results_n[:, 1] / EF

plt.figure(figsize=(8, 5))
plt.plot(n_list, mu_n, 'b-o', label=r'$\mu / E_F$')
plt.plot(n_list, delta_n, 'r-s', label=r'$\Delta / E_F$')
plt.xscale('log')
plt.xlabel('Number of points (N)')
plt.ylabel('Dimensionless Energy')
plt.title(f'Test 1: Convergence vs Grid Density ($k_{{max}}={k_max_fixed}$)')
plt.grid(True, which="both", ls="-", alpha=0.2)
plt.legend()
plt.savefig(os.path.join(output_dir, "convergence_N.png"))
plt.close()

# =============================================================================
# TEST 2: Momentum Cutoff (k_max)
# =============================================================================
print("Running Test 2: Momentum Cutoff (k_max)...")
k_max_list = [10, 50, 100, 200, 500, 1000]
n_fixed = 20000 
results_k = []

for km in k_max_list:
    k, dk = create_k_grid(k_max=km, N=n_fixed)
    sol = solve_bcs_system(interaction_a, n_target, k, dk, initial_guess=[EF, 0.5*EF])
    results_k.append(sol)

results_k = np.array(results_k)
mu_k = results_k[:, 0] / EF
delta_k = results_k[:, 1] / EF

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
print("Running Test 3: Solver Residuals...")
k, dk = create_k_grid(k_max=100, N=10000)
initial_guess = [EF, 0.6*EF]

# Using bcs_objective_functions directly to check the residuals
sol_root = root(bcs_objective_functions, initial_guess, args=(interaction_a, n_target, k, dk))

res_gap, res_num = sol_root.fun
print(f"--- Solver Status ---")
print(f"Success: {sol_root.success}")
print(f"Residual Gap Equation: {res_gap:.2e}")
print(f"Residual Number Equation: {res_num:.2e}")

report_path = os.path.join(output_dir, "convergence_report.txt")
with open(report_path, "w") as f:
    f.write("CONVERGENCE STUDY REPORT\n")
    f.write("========================\n\n")
    f.write(f"Test at Unitarity (1/kFa = 0)\n")
    f.write(f"Final Solution: mu/EF = {sol_root.x[0]/EF:.4f}, Delta/EF = {sol_root.x[1]/EF:.4f}\n")
    f.write(f"Solver Residuals: Gap={res_gap:.2e}, Number={res_num:.2e}\n")

print(f"--- Study Complete. Results in '{output_dir}/' ---")
