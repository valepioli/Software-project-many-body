# BCS-BEC Crossover Simulation in 3D Fermi Gases

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue)](https://www.python.org/)

## Overview
This project simulates the **BCS-BEC Crossover** in a 3D ultracold Fermi gas at zero temperature ($T=0$). Using a mean-field approach, the software solves the coupled self-consistent equations to track the evolution of the system from weakly bound Cooper pairs (BCS limit) to a Bose-Einstein Condensate (BEC) of tightly bound dimers.

The transition is controlled by tuning the interaction strength via the dimensionless parameter $1/(k_F a)$, where $a$ is the s-wave scattering length.

---

## Physical Model

The simulation solves the mean-field equations for a continuum Fermi gas.

### 1. Excitation Spectrum
The energy of the quasiparticle excitations is given by:
$$E_k = \sqrt{(\epsilon_k - \mu)^2 + \Delta^2}$$
where $\epsilon_k = \frac{\hbar^2 k^2}{2m}$ is the single-particle kinetic energy.

### 2. The Regularized Gap Equation
In 3D, the contact interaction leads to a UV divergence. We implement a regularized version of the gap equation to ensure convergence:
$$\frac{m}{4\pi \hbar^2 a} = \int \frac{d^3k}{(2\pi)^3} \left( \frac{1}{2\epsilon_k} - \frac{1}{2E_k} \right)$$

### 3. The Number Equation
The total particle density $n$ is kept constant by solving for the chemical potential $\mu$:
$$n = \int \frac{d^3k}{(2\pi)^3} \left( 1 - \frac{\epsilon_k - \mu}{E_k} \right)$$

---

## Features
*   **Self-Consistent Solver**: Implements `scipy.optimize.root` with an iterative continuation method for high stability.
*   **UV Convergence**: Uses a subtraction scheme for the gap equation, making results independent of the high-momentum cutoff.
*   **Automatic Export**: Numerical data (`.txt`) and high-resolution plots (`.png`) are automatically saved to the `results/` folder.
*   **Dimensionless Units**: All results are normalized to the Fermi Energy $E_F$ and Fermi Momentum $k_F$.

---

## Project Structure
```text
├── src/
│   ├── config.py       # Physical constants (kF, EF, n)
│   ├── physics.py      # Integrals and Energy Spectrum
│   ├── solver.py       # Root-finding algorithm
│   └── plotting.py     # Visualization logic
├── results/            # Output: plots and numerical data
├── tests/              # Unit tests
├── main.py             # Main execution script
└── requirements.txt    # Dependencies (numpy, scipy, matplotlib)

## Goals

* Implement a numerical solver for the coupled equations
* Explore the crossover physics
* Produce plots of μ and Δ

---

## How to run

```bash
pip install -r requirements.txt
python src/main.py
```

---

## Status

Project started – work in progress.
