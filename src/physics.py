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



