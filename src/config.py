# config.py
"""
Created on Sun Apr 12 20:07 2026

@author: Pioli Valeria
"""

import numpy as np

# =============================================================================
# PHYSICAL AND NUMERICAL PARAMETERS
# =============================================================================

# Fermi Momentum (set to 1.0 as the unit of inverse length)
kF = 1.0

# Cutoff in k-space (units of kF)
k_max = 100.0

# Number of points in the momentum grid
N = 10000

# Fermi Energy (EF = kF^2 / 2m). For m=1, EF = 0.5
EF = 0.5 * kF**2

# Target density for a single spin species (n_tot / 2)
n_target = kF**3 / (6.0 * np.pi**2)

#Minimum of interaction range 1/kFa
start_x=-3

#Maximum of interaction range 1/kFa
end_x=3

#Number of points for the simulation(1/kFa steps)
steps=40
