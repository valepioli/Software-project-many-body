import pytest
import numpy as np
from src.solver import bcs_objective_functions, solve_bcs_system
from src.physics import create_k_grid

def test_objective_functions_return_list():
    """
    WHAT: Verify that bcs_objective_functions returns a list of two values.
    WHY: The SciPy 'root' solver requires the objective function to return 
         a list or array of the same length as the input variables.
    """
    k, dk = create_k_grid()
    res = bcs_objective_functions([1.0, 0.5], a=1.0, n_target=1.0, k=k, dk=dk)

    assert len(res) == 2
    assert isinstance(res, list)

def test_solver_consistency_at_unitarity():
    """
    WHAT: Check if the solver can find a finite solution at Unitarity (1/a = 0).
    WHY: Unitarity is the most numerically challenging point of the crossover. 
         If the solver works here, the implementation is robust.
    """
    # At Unitarity, a is infinite (represented by a very large number)
    a_inf = 1e6
    n_target = 1.0
    k, dk = create_k_grid()

    mu, delta = solve_bcs_system(a_inf, n_target, k, dk)

    assert np.isfinite(mu), "Chemical potential should be a finite number."
    assert np.isfinite(delta), "Gap should be a finite number."
    assert delta > 0, "The pairing gap must be positive."

def test_solver_unsolvable_case():
    """
    WHAT: Ensure the solver handles failures by returning NaNs instead of crashing.
    WHY: To satisfy the requirement that the code should handle errors gracefully, 
         especially during automated loops over many parameters.
    """
    # Using an impossibly high density for a small grid
    k, dk = create_k_grid(k_max=1.0, N=10)
    mu, delta = solve_bcs_system(a=1.0, n_target=100.0, k=k, dk=dk)

    assert np.isnan(mu)
    assert np.isnan(delta)
