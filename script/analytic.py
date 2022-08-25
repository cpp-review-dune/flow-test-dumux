#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sympy as sp
import numpy as np
import scipy as sc
from scipy.constants import g
from numpy.typing import ArrayLike


def show_version():
    sympy_version = f"SymPy version: {sp.__version__}\n"
    numpy_version = f"NumPy version: {np.__version__}\n"
    scipy_version = f"SciPy version: {sc.__version__}\n"

    return sympy_version + numpy_version + scipy_version


phi: float = 4e-1  # Porosity
rho_w: float = 1e3  # Impermeability
mu_w: float = 1e-3  # Water viscosity
K_int: float = 5e-12  # Intrinsic permeability
gravity_force: ArrayLike = np.array([0, 0, -g])  # Gravitational acceleration
p_n: float = 1e5  # Air pressure
pcsweLow: float = 1e4
pcsweHigh: float = 1e3
m: float = -5e5
sThres: float = 1e-2
EntryPc: float = 0
MaxPc: float = 1e10
p_d: float = 1e3
Lambda: float = 2

x, z, globalPos1, time, pwBottom, pwTop = sp.symbols(
    "x z globalPos1 time pwBottom pwTop"
)
(globalPos1 * globalPos1 * globalPos1 * time / 100) / 64
p_w = pwBottom + 0.5 * (pwTop - pwBottom) * (
    1 + sp.tan(5 * globalPos1 - 15 + time / 10)
)
p_c = p_n - p_w
S_w = lambda p_c: sThres + (p_c - pcsweLow) / m

if __name__ == "__main__":
    print(show_version())
    print(f"Φ = {phi}\nρ_w = {rho_w}\nμ_w = {mu_w}")
    print(gravity_force)
