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
    
    # Residual of the Gap Equation: Integral - 1/(4 * pi * a)
    target_gap = 1.0 / (4.0 * np.pi * a)
    res_gap = gap_integral(k, mu, Delta, dk) - target_gap
    
    # Residual of the Number Equation: Integral - n_target
    res_num = number_integral(k, mu, Delta, dk) - n_target
    
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
    sol = root(bcs_objective_functions, initial_guess, args=(a, n_target, k, dk))
    
    if sol.success:
        return sol.x # [mu_solution, Delta_solution]
    else:
        # Return NaN to allow the loop to continue without crashing
        return [np.nan, np.nan]
