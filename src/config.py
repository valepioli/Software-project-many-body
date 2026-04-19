# config.py
"""
Created on Sun Apr 12 20:07 2026

@author: Pioli Valeria
"""

import numpy as np

# =============================================================================
# PHYSICAL AND NUMERICAL PARAMETERS
# =============================================================================

# Fermi Momentum (usually set to 1.0 as the unit of inverse length)
kF = 1.0

# Cutoff in k-space (units of kF)
# 100.0 is very large; ensure your solver remains stable
k_max = 100.0

# Number of points in the momentum grid
# Increase N if you increase k_max to maintain resolution
N = 10000

# Fermi Energy (EF = kF^2 / 2m). For m=1, EF = 0.5
EF = 0.5 * kF**2

# Target density for a single spin species (n_tot / 2)
# Formula: n = kF^3 / (6 * pi^2)
n_target = kF**3 / (6.0 * np.pi**2)

#Minimum of interaction range 1/kFa
int_min=-3

#Maximum of interaction range 1/kFa
int_max=3

#Number of points for the simulation(1/kFa steps)
steps=40
