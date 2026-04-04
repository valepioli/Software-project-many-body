# BEC–BCS Crossover Simulation

## Overview

This project simulates the BEC–BCS crossover in ultracold Fermi gases using a mean-field approach.

The goal is to study how the system evolves from:

* a Bose–Einstein condensate (BEC) of tightly bound molecules
  to
* a BCS superfluid of weakly bound Cooper pairs

by tuning the interaction strength.

---

## Physical Model

We solve the mean-field BCS equations at zero temperature to compute:

* the pairing gap Δ
* the chemical potential μ

as a function of the interaction parameter:

1/(k_F a)

---

## Goals

* Implement a numerical solver for the coupled equations
* Explore the crossover physics
* Produce clear plots of μ and Δ

---

## How to run

```bash
pip install -r requirements.txt
python src/main.py
```

---

## Status

Project started – work in progress.
