# physics.py
"""
Created on Sat Apr 4 18:15 2026

@author: Pioli Valeria
"""
import numpy as np

# ========================================
# Momentum grid
# ========================================

def create_k_grid(k_max, N):
    """
    Creates a 1D momentum grid for numerical integration.

    Parameters:
    -----------
    k_max : float
        Maximum momentum (cutoff) in units of k_F.
    N : int
        Number of points in the grid.

    Returns:
    --------
    k : numpy.ndarray
        Array of momentum points from small value to k_max.
    dk : float
        Momentum spacing.
    """
    # avoid k=0 to prevent division by zero in gap integral
    k = np.linspace(1e-5, k_max, N)
    dk = k[1] - k[0]
    return k, dk

# ========================================
# Energy functions
# =======================================
def compute_quasiparticle_energies(k, mu, Delta):
    """
    Calculate all energy components for the BCS-BEC crossover.

    Parameters:
    -----------
    k : ndarray
        Momentum grid.
    mu : float
        Chemical potential.
    Delta : float
        Superconducting gap.

    Returns:
    --------
    eps_k : ndarray
        Kinetic energy (k^2 / 2m, here m=1 or 1/2 depending on convention).
    xi_k : ndarray
        Shifted energy (eps_k - mu).
    E_k : ndarray
        Quasiparticle excitation energy sqrt(xi_k^2 + Delta^2).
    """
    eps_k = k**2

    xi_k = eps_k - mu
    E_k = np.sqrt(xi_k**2 + Delta**2)

    return eps_k, xi_k, E_k

# ========================================
# Integrals for gap and number equations
# ========================================

def gap_integral(k, mu, Delta, dk):
    """
    Computes the regularized gap integral for the BCS equation.

    Integral form:
        I_gap = ∫ k^2 (1/(2 E_k) - 1/(2 k^2)) dk

    Parameters:
    -----------
    k : array
        Momentum grid
    mu : float
        Chemical potential
    Delta : float
        Pairing gap
    dk : float
        Momentum spacing

    Returns:
    --------
    float
        Value of the gap integral
    """
    # Unpack only what is needed (E_k is the third value)
    _, _, Ek = compute_quasiparticle_energies(k, mu, Delta)
    integrand = 1/(2*Ek) - 1/(2*k**2)
    return np.sum(k**2 * integrand) * dk

def number_integral(k, mu, Delta, dk):
    """
    Computes the particle number integral for the number equation.

    Integral form:
        n = ∫ k^2 * (1 - xi_k / E_k) dk

    Parameters:
    -----------
    k : array
        Momentum grid
    mu : float
        Chemical potential
    Delta : float
        Pairing gap
    dk : float
        Momentum spacing

    Returns:
    --------
    float
        Value of the number integral
    """
    # Unpack xi_k and E_k
    _, xik, Ek = compute_quasiparticle_energies(k, mu, Delta)
    integrand = 1 - xik/Ek
    return np.sum(k**2 * integrand) * dk
