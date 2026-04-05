# BEC–BCS Crossover Simulation

## Overview

This project simulates the BEC–BCS crossover in ultracold Fermi gases using a mean-field approach.

The goal is to study how the system evolves from a Bose–Einstein condensate (BEC) of tightly bound molecules 
to a BCS superfluid of weakly bound Cooper pairs by tuning the interaction strength.

---

## Physical Model
The problem is described by a Fermi Hubbard many-body hamiltonian:
add picture

In the mean field approach we define a parameter pairing gap Δ, which is a measure of the strength of pairing:
add picture

Using the mean field approach we find the energy of the excitations of the system:
add picture

And the two fundamental equations:
* gap equation (coherence in pairing):
add picture

* number equation(mantains number of particles):
add picture 

We solve these mean-field equations at zero temperature, looking for the pairing gap and chemical potential so that both equations are 
satisfied, solving the non linear coupled problem we compute:

* the pairing gap Δ
* the chemical potential μ

as a function of the interaction parameter:

1/(k_F a)

---

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
