import pytest
import numpy as np
from src.physics import (
    create_k_grid, 
    compute_quasiparticle_energies, 
    gap_integral, 
    number_integral
)

# --- FIXTURES ---
# Fixtures allow us to reuse the same setup (the k-grid) for multiple tests
@pytest.fixture
def default_grid():
    """
    Provides a standard momentum grid to ensure test consistency.
    """
    return create_k_grid(k_max=10.0, N=1000)

# --- MOMENTUM GRID TESTS ---

def test_grid_size_matches_input():
    """
    WHAT: Check if create_k_grid returns an array of length N.
    WHY: To ensure the numerical resolution requested by the user is respected.
    """
    N_test = 500
    k, _ = create_k_grid(N=N_test)
    assert len(k) == N_test

def test_grid_avoids_zero():
    """
    WHAT: Ensure the first k-point is strictly positive.
    WHY: The gap equation has a 1/k^2 term; k=0 would cause a division by zero.
    """
    k, _ = create_k_grid()
    assert k[0] > 0

# --- ENERGY FUNCTION TESTS ---

@pytest.mark.parametrize("mu, Delta", [(1.0, 0.5), (0.0, 0.2), (-1.0, 1.0)])
def test_quasiparticle_energy_consistency(default_grid, mu, Delta):
    """
    WHAT: Verify E_k = sqrt(xi_k^2 + Delta^2) using np.allclose.
    WHY: To ensure the Bogoliubov dispersion relation is mathematically correct.
    """
    k, _ = default_grid
    _, xik, Ek = compute_quasiparticle_energies(k, mu, Delta)
    expected_Ek = np.sqrt(xik**2 + Delta**2)
    assert np.allclose(Ek, expected_Ek)

def test_quasiparticle_energy_is_above_gap(default_grid):
    """
    WHAT: Ensure E_k is always >= Delta.
    WHY: Physically, the gap Delta is the minimum energy for an excitation.
    """
    k, _ = default_grid
    Delta = 0.5
    _, _, Ek = compute_quasiparticle_energies(k, mu=1.0, Delta=Delta)
    assert np.all(Ek >= (Delta - 1e-9))

# --- INTEGRAL TESTS (PHYSICAL CONSISTENCY) ---

def test_number_density_is_always_positive(default_grid):
    """
    WHAT: Check if number_integral is positive.
    WHY: Density is a physical magnitude that cannot be negative.
    """
    k, dk = default_grid
    n = number_integral(k, mu=-1.0, Delta=0.5, dk=dk)
    assert n > 0

def test_density_increases_with_mu(default_grid):
    """
    WHAT: Verify that increasing mu increases density.
    WHY: Systems must have positive compressibility to be thermodynamically stable.
    """
    k, dk = default_grid
    n_low = number_integral(k, mu=0.5, Delta=0.5, dk=dk)
    n_high = number_integral(k, mu=1.0, Delta=0.5, dk=dk)
    assert n_high > n_low

def test_gap_integral_decreases_with_delta(default_grid):
    """
    WHAT: Verify the gap integral decreases as Delta increases.
    WHY: This monotonicity is required for the solver to converge to a unique solution.
    """
    k, dk = default_grid
    g_small = gap_integral(k, 1.0, 0.1, dk)
    g_large = gap_integral(k, 1.0, 1.0, dk)
    assert g_large < g_small

def test_integrals_are_finite(default_grid):
    """
    WHAT: Ensure integrals do not return NaN or Inf.
    WHY: The solver (fsolve) will crash if it encounters non-finite values.
    """
    k, dk = default_grid
    assert np.isfinite(gap_integral(k, 1.0, 0.5, dk))
    assert np.isfinite(number_integral(k, 1.0, 0.5, dk))
