import numpy as np

from simulator import Electron, LightQuarks, ZBoson
from simulator.constants import *


e = Electron()
Q = LightQuarks()
Z = ZBoson()


def analytical_sigma(integral1: float, integral2: float, integral3: float) -> float:
    """
    Compute the cross-section for the electron-positron to quark-antiquark
    process to leading-order from the analytical results of the integrals.

    :param integral1: Integral of f(s)/s.
    :param integral2: Integral of chi1(s) * f(s) / s.
    :param integral3: Integral of chi2(s) * f(s) / s.
    :return: Cross-section in picobarns.
    """

    sigma = 4 * np.pi * (alpha_QED ** 2) * f_conv

    # Multiplies integral of f(s)/s.
    photon_amp = np.sum((e.charge * Q.charges) ** 2) * integral1

    # Multiplies integral of chi1(s) * f(s) / s.
    interference = np.sum(2 * e.charge * Q.charges * e.V * Q.vector_couplings) * integral2

    # Multiplies integral of chi2(s) * f(s) / s.
    zboson_amp = np.sum(e.squared_coupling * Q.squared_couplings) * integral3

    sigma *= (photon_amp + interference + zboson_amp)

    return sigma


def exact_value_1a() -> float:
    """
    Compute the exact cross-section for the case in section 1, part a),
    consisting on f(s) being a Dirac delta centered at the squared mass
    of the Z boson.

    :return: Cross-section in picobarns.
    """
    integral1 = 1 / Z.mass2  # Integral of f(s) / s.
    integral2 = 0.  # Integral of chi1(s) * f(s) / s.
    integral3 = (kappa ** 2) / Z.width2  # Integral of chi2(s) * f(s) / s.
    sigma = analytical_sigma(integral1, integral2, integral3)

    return sigma


def exact_value_1c() -> float:
    """
    Compute the exact cross-section for the case in section 1, part c),
    consisting on f(s) being a uniform distribution between (Z.mass - 3 * Z.width) ** 2
    and (Z.mass + 3 * Z.width) ** 2.

    :return: Cross-section in picobarns.
    """
    s_min = (Z.mass - 3 * Z.width) ** 2
    s_max = (Z.mass + 3 * Z.width) ** 2
    fs = 1 / (s_max - s_min)  # Flat beam spectrum.

    integral1 = fs * np.log(s_max / s_min)  # Integral of f(s) / s.

    log_term = 0.5 * np.log(((Z.mass ** 4) + Z.mass2 * (Z.width2 - 2 * s_max) + (s_max ** 2)) /
                            ((Z.mass ** 4) + Z.mass2 * (Z.width2 - 2 * s_min) + (s_min ** 2))
                            )

    integral2 = fs * kappa * log_term  # Integral of chi1(s) * f(s) / s.

    tan_term = - (Z.mass / Z.width) * (np.arctan((Z.mass2 - s_max) / (Z.mass * Z.width)) -
                                       np.arctan((Z.mass2 - s_min) / (Z.mass * Z.width)))

    integral3 = fs * (kappa ** 2) * (log_term + tan_term)  # Integral of chi2(s) * f(s) / s.

    sigma = analytical_sigma(integral1, integral2, integral3)

    return sigma


if __name__ == '__main__':
    exact_sigma_1a = exact_value_1a()
    print(f"Exact cross-section for point 1a: {exact_sigma_1a:,f}")
    exact_sigma_1c = exact_value_1c()
    print(f"Exact cross-section for point 1c: {exact_sigma_1c:,f}")
