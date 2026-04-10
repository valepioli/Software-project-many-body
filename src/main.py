#main.py
"""
Created on Sat Apr 4 18:15 2026

@author: Pioli Valeria
"""

print("BEC–BCS project started")

from src.physics import create_k_grid
from src.utils import load_parameters

# Load the parameters
params = load_parameters("parameters.txt")

# This will pass N=2000 and k_max=8.0 automatically
k, dk = create_k_grid(k_max=params['k_max'], N=int(params['N']))

n_target = params['n_target']
