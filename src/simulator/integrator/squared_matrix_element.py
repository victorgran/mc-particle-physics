"""
Definition of the squared matrix element for the electron-positron
to quark-antiquark cross-section to leading order.
"""

import numpy as np
from simulator.constants import alpha_QED, kappa, QCD_colors
from .particles import Electron, LightQuarks, ZBoson


e = Electron()
Q = LightQuarks()
Z = ZBoson()


def chi_funcs(s: float | np.ndarray):
    # Both chi functions have the same denominator.
    den = (s - Z.mass2) ** 2 + (Z.width2 * Z.mass2)
    chi1 = kappa * s * (s - Z.mass2) / den
    chi2 = (kappa ** 2) * (s ** 2) / den

    return chi1, chi2


def squared_matrix_element(s: float | np.ndarray, cos_theta: float | np.ndarray, q: int | np.ndarray):
    chi1, chi2 = chi_funcs(s)

    curly_brackets1 = (e.charge * Q.charge(q)) ** 2 + \
                      (2 * e.charge * e.V * Q.charge(q) * Q.V(q) * chi1) + \
                      (e.squared_coupling * Q.squared_coupling(q) * chi2)
    curly_brackets2 = (4 * e.charge * Q.charge(q) * e.A * Q.A(q) * chi1) + \
                      (8 * e.A * e.V * Q.A(q) * Q.V(q) * chi2)

    m_sqr = ((1. + cos_theta ** 2) * curly_brackets1) + (cos_theta * curly_brackets2)
    m_sqr *= (4. * np.pi * alpha_QED) ** 2
    m_sqr *= QCD_colors

    return m_sqr
