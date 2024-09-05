"""
Definitions of Dirac, Uniform and Breit-Wigner distributions.

Each of these is implemented as a child class from the parent “Distribution” class,
which requires defining the sampling and evaluation properties for each of them.
"""
from __future__ import annotations

import numpy as np


class Distribution:
    """
    Basic distribution class with a sampling method
    and an expression to evaluate values on it.
    """
    def __init__(self, sample_method, evaluate_distro):
        self._sample_method = sample_method
        self._evaluate_distro = evaluate_distro

    def sample(self, sample_size):
        return self._sample_method(sample_size)

    def evaluate_distro(self, variables):
        return self._evaluate_distro(variables)


class Dirac(Distribution):
    """
    One-dimensional Dirac distribution with center at “x0”.
    Sampling it returns simply “x0” and evaluating it returns 1.
    """
    def __init__(self, x0: int | float):
        self._x0 = x0  # Root of the Dirac delta.
        def sample_method(sample_size: int): return np.ones(sample_size) * self._x0
        def evaluate_distro(_): return 1
        super().__init__(sample_method, evaluate_distro)

    def __repr__(self):
        distro_name = type(self).__name__
        return f"{distro_name}(x0={self.x0})"

    @property
    def x0(self):
        return self._x0


class Uniform(Distribution):
    """
    Uniform distribution between lower and upper values.
    Its sampling method draws values in the interval [lower, upper)
    and evaluating the distribution yields always 1/(upper - lower).
    """
    def __init__(self, lower: int | float, upper: int | float, rng=None):
        self._lower = lower
        self._upper = upper
        self._rng = rng if rng is not None else np.random.default_rng()

        def sample_method(sample_size: int):
            return self._rng.uniform(self._lower, self._upper, sample_size)

        def evaluate_distro(_): return 1. / (self._upper - self._lower)

        super().__init__(sample_method=sample_method, evaluate_distro=evaluate_distro)

    def __repr__(self):
        distro_name = type(self).__name__
        return f"{distro_name}(lower={self.lower}, upper={self.upper}, rng={self.rng})"

    @property
    def lower(self):
        return self._lower

    @property
    def upper(self):
        return self._upper

    @property
    def rng(self):
        return self._rng

    @rng.setter
    def rng(self, rng):
        self._rng = rng


class BreitWigner(Distribution):
    """
    Breit-Wigner distribution for a particle of given mass and decay_width,
    for values of s between s_min and s_max. The sampling method makes use
    of the Breit-Wigner mapping, generating uniform values of an auxiliary
    variable 'rho', and then transforming back to 's'.
    """
    def __init__(self, s_min: int | float, s_max: int | float,
                 mass: int | float, decay_width: int | float, rng=None):
        self._s_min = s_min
        self._s_max = s_max
        self._mass = mass
        self._width = decay_width
        self._rho_min = self.getRho(s_min)
        self._rho_max = self.getRho(s_max)
        self._rng = rng if rng is not None else np.random.default_rng()

        def sample_method(sample_size):
            rho_values = self._rng.uniform(self._rho_min, self._rho_max, sample_size)
            s_sample = (self._mass * self._width * np.tan(rho_values)) + (self._mass ** 2)
            return s_sample

        def evaluate_distro(s_values: int | float | np.ndarray):
            denominator = (s_values - (self._mass ** 2)) ** 2 + (self._mass * self._width) ** 2
            g_s = 1. / denominator  # Breit-Wigner propagator.
            g_s *= (self._mass * self._width) / (self._rho_max - self._rho_min)  # Normalization factor.
            return g_s

        super().__init__(sample_method=sample_method, evaluate_distro=evaluate_distro)

    def __repr__(self):
        distro_name = type(self).__name__
        return f"{distro_name}(s_min={self.s_min}, s_max={self.s_max}, mass={self.mass}, " \
               f"decay_width={self.decay_width}, rng={self.rng})"

    @property
    def s_min(self):
        return self._s_min

    @property
    def s_max(self):
        return self._s_max

    @property
    def mass(self):
        return self._mass

    @property
    def decay_width(self):
        return self._width

    @property
    def rng(self):
        return self._rng

    @rng.setter
    def rng(self, rng):
        self._rng = rng

    def getRho(self, s):
        rho = np.arctan((s - self._mass ** 2) / (self._mass * self._width))
        return rho
