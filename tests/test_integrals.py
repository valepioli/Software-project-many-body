import pytest
import numpy as np
from src.physics import create_k_grid, gap_integral, number_integral

# --- FIXTURES ---
# Fixtures allow us to reuse the same setup (the k-grid) for multiple tests
@pytest.fixture
def standard_grid():
    """
    Provides a standard grid for testing.
    """
    k, dk = create_k_grid()
    return k, dk

# --- UNIT TESTS ---

def test_grid_properties(standard_grid):
    """
    Check if the grid is generated correctly:
    - Should not contain zero (to avoid division by zero in the gap equation).
    - Should be monotonically increasing.
    """
    k, dk = standard_grid
    assert np.all(k > 0), "Grid points must be strictly positive."
    assert np.all(np.diff(k) > 0), "Grid points must be monotonically increasing."
    assert dk > 0, "Differential dk must be positive."

def test_integrals_finite_output(standard_grid):
    """
    Verify that the integrals return finite numerical values (no NaN or Inf).
    """
    k, dk = standard_grid
    mu, Delta = 1.0, 0.5
    
    gap_val = gap_integral(k, mu, Delta, dk)
    num_val = number_integral(k, mu, Delta, dk)
    
    assert np.isfinite(gap_val), "Gap integral returned NaN or Inf."
    assert np.isfinite(num_val), "Number integral returned NaN or Inf."

# --- PHYSICAL CONSISTENCY TESTS ---

def test_number_density_positivity(standard_grid):
    """
    Physical constraint: The number density (number_integral) must always 
    be positive, regardless of the chemical potential sign.
    """
    k, dk = standard_grid
    # Test for both BCS (positive mu) and BEC (negative mu) regimes
    for mu in [-2.0, 0.0, 2.0]:
        num_val = number_integral(k, mu, Delta=0.5, dk=dk)
        assert num_val > 0, f"Density must be positive even for mu={mu}."

def test_density_monotonicity_with_mu(standard_grid):
    """
    Physical constraint: Increasing the chemical potential (mu) 
    must increase the number density (n).
    """
    k, dk = standard_grid
    Delta = 0.5
    n_low = number_integral(k, mu=0.5, Delta=Delta, dk=dk)
    n_high = number_integral(k, mu=1.5, Delta=Delta, dk=dk)
    
    assert n_high > n_low, "Physical error: Density must increase with chemical potential."

@pytest.mark.parametrize("Delta", [0.1, 0.5, 1.0])
def test_gap_integral_trends(standard_grid, Delta):
    """
    The gap integral (regularized) should be a well-behaved float.
    Using parametrize allows testing multiple values of Delta automatically.
    """
    k, dk = standard_grid
    mu = 1.0
    val = gap_integral(k, mu, Delta, dk)
    assert isinstance(val, (float, np.float64))

# --- CONVERGENCE TESTS ---

def test_numerical_convergence():
    """
    Verify that doubling the grid density doesn't change the integral 
    significantly (Grid Convergence Study).
    """
    # Create a coarse grid and a fine grid
    # (Note: This assumes your create_k_grid can take N as an argument)
    k_coarse, dk_coarse = create_k_grid(N=1000) 
    k_fine, dk_fine = create_k_grid(N=2000)
    
    mu, Delta = 1.0, 0.5
    val_coarse = number_integral(k_coarse, mu, Delta, dk_coarse)
    val_fine = number_integral(k_fine, mu, Delta, dk_fine)
    
    relative_error = abs(val_fine - val_coarse) / val_fine
    assert relative_error < 1e-3, "Integral did not converge within 0.1% tolerance."
