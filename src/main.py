#main.py
"""
Created on Sat Apr 4 18:15 2026

@author: Pioli Valeria
"""

print("BEC–BCS project started")

from src.physics import create_k_grid
from src.utils import load_parameters

#1 Load the parameters
params = load_parameters("parameters.txt")

# This will pass N and k_max automatically
k, dk = create_k_grid(k_max=params['k_max'], N=int(params['N']))

n_target = params['n_target']

# =============================================================================
# 2. CROSSOVER SWEEP (ROOT-FINDING LOOP)
# =============================================================================

results = []

# initial_guess: We start with a guess close to the BCS limit (mu ~ EF, Delta small).
# This guess will be updated dynamically to improve convergence.
current_guess = [1.0, 0.5]

for x in interaction_range:
    """
    We iterate over the dimensionless interaction parameter x = 1/(kF * a).
    BCS regime: x < 0 | Unitarity: x = 0 | BEC regime: x > 0
    """

    # 1. Map dimensionless parameter 'x' to scattering length 'a'.
    # For x = 0 (unitary limit), we use a large number (1e6) to approximate
    # the physical infinity 1/a -> 0.
    a = 1.0 / x if abs(x) > 1e-10 else 1e6

    # 2. Call the solver to find the [mu, Delta] pair that satisfies the system.
    sol = solve_bcs_system(a, params['n_target'], k, dk, initial_guess=current_guess)

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
mu_vals = results_array[:, 0]
delta_vals = results_array[:, 1]

# 4. VISUALIZATION
plot_bcs_bec_crossover(interaction_range, mu_vals, delta_vals)
