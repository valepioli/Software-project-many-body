# physics.py
"""
Created on Sat Apr 4 18:15 2026

@author: Pioli Valeria
"""
import numpy as np

# ========================================
# Momentum grid
# ========================================

def create_k_grid(k_max=10.0, N=1000):
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
    eps_k =0.5*k**2

    xi_k = eps_k - mu
    E_k = np.sqrt(xi_k**2 + np.abs(Delta)**2)

    return eps_k, xi_k, E_k

# ========================================
# Integrals for gap and number equations
# ========================================

def gap_integral(k, mu, Delta, dk):
    """
    Computes the regularized gap integral for the BCS equation.

    To avoid UV divergence in 3D, we subtract the vacuum contribution:
    Integral form: ∫ [k^2 / 2π^2] * (1/(2 E_k) - 1/(2 eps_k)) dk
    Since eps_k = k^2/2, the second term simplifies to 1/k^2.
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
    # Regularized integrand: removes the linear divergence as k -> infinity
    integrand = (1.0 / (2.0 * Ek)) - (1.0 / (k**2))

    # (1 / 2*pi^2) comes from the 3D volume element 4*pi*k^2 / (2*pi)^3
    return (1.0 / (2.0 * np.pi**2)) * np.sum(k**2 * integrand) * dk

def number_integral(k, mu, Delta, dk):
    """

    Computes the particle density integral (per spin species).

    n_species = ∫ [k^2 / 2π^2] * v_k^2 dk
    where v_k^2 = 0.5 * (1 - xi_k / E_k) is the hole occupancy.

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
    # Occupancy probability v_k^2
    v_k_sq = 0.5 * (1.0 - xik / Ek)

    # Total density for one spin species in 3D
    return (1.0 / (2.0 * np.pi**2)) * np.sum(k**2 * v_k_sq) * dk
