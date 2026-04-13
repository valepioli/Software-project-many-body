# solver.py
"""
Created on Sat Apr 11 17:46 2026

@author: Pioli Valeria
"""

import numpy as np
from scipy.optimize import root
from src.physics import gap_integral, number_integral

def bcs_objective_functions(vars, a, n_target, k, dk):
    """
    Calculates the residuals for the Gap and Number equations.

    Parameters:
    -----------
    vars : list or tuple
        The unknowns [mu, Delta].
    a : float
        Scattering length.
    n_target : float
        The desired particle density.
    k, dk : array, float
        The momentum grid and its spacing.

    Returns:
    --------
    residuals : list
        The difference between calculated integrals and target values.
    """
    mu, Delta = vars
    # Use absolute value of Delta to prevent the solver from exploring 
    # unphysical negative gap values which cause square root errors.
    Delta_abs = np.abs(Delta)

    # Residual of the Gap Equation:
    # Theory (m=1, hbar=1): GapIntegral + 1/(4 * pi * a) = 0
    # The gap_integral function already includes the 1/(2*pi^2) factor.
    inv_a = 1.0 / a if a != 0 else 1e10 # Safety for 1/a
    res_gap = gap_integral(k, mu, Delta_abs, dk) + (inv_a / (4.0 * np.pi))

    # Residual of the Number Equation: Integral - n_target = 0
    # n_target should be kF^3 / (6 * pi^2) for a single spin species.
    res_num = number_integral(k, mu, Delta_abs, dk) - n_target
    return [res_gap, res_num]

def solve_bcs_system(a, n_target, k, dk, initial_guess=[1.0, 0.5]):
    """
    Finds the values of mu and Delta that satisfy the BCS-BEC crossover equations.

    Parameters:
    -----------
    a : float
        Scattering length.
    n_target : float
        Desired density.
    k, dk : array, float
        Grid parameters.
    initial_guess : list
        Starting point for the numerical solver.

    Returns:
    --------
    mu, Delta : floats
        The solved physical parameters. Returns [NaN, NaN] if solver fails.
    """
    sol = root(bcs_objective_functions, initial_guess, args=(a, n_target, k, dk), tol=1e-9)

    if sol.success:
        mu_sol, delta_sol = sol.x
        # Return the absolute value of Delta as the gap is physically positive
        return [mu_sol, np.abs(delta_sol)]
    else:
        # If solver fails, return NaNs so the main loop can handle the failure gracefully
        return [np.nan, np.nan]
