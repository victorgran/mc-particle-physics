"""
Definitions of particles as classes, containing their relevant attributes for the project.
"""
from __future__ import annotations

import numpy as np

from simulator.constants import sin2w


class Electron:
    def __init__(self):
        self._charge = -1
        self._axial_coupling = -0.5  # Weak isospin.
        self._vector_coupling = self._axial_coupling - 2 * self._charge * sin2w
        self._squared_coupling = (self._axial_coupling ** 2) + (self._vector_coupling ** 2)

    @property
    def charge(self):
        return self._charge

    @property
    def A(self):
        return self._axial_coupling

    @property
    def V(self):
        return self._vector_coupling

    @property
    def squared_coupling(self):
        return self._squared_coupling


class LightQuarks:
    def __init__(self):
        self._flavours = ["down", "up", "strange", "charm", "bottom"]
        self._charges = np.array([-1 / 3, 2 / 3, -1 / 3, 2 / 3, -1 / 3])
        self._axial_couplings = np.array([-0.5, 0.5, -0.5, 0.5, -0.5])
        self._vector_couplings = self._axial_couplings - 2 * self._charges * sin2w
        self._squared_couplings = (self._axial_couplings ** 2) + (self._vector_couplings ** 2)

    @property
    def flavours(self):
        return self._flavours

    @property
    def charges(self):
        return self._charges

    @property
    def axial_couplings(self):
        return self._axial_couplings

    @property
    def vector_couplings(self):
        return self._vector_couplings

    @property
    def squared_couplings(self):
        return self._squared_couplings

    # Access single quark properties by index or by array.
    def flavour(self, flavour: int | np.ndarray):
        return self._flavours[flavour]

    def charge(self, flavour: int | np.ndarray):
        return self._charges[flavour]

    def A(self, flavour: int | np.ndarray):
        return self._axial_couplings[flavour]

    def V(self, flavour: int | np.ndarray):
        return self._vector_couplings[flavour]

    def squared_coupling(self, flavour: int | np.ndarray):
        return self._squared_couplings[flavour]


class ZBoson:
    def __init__(self):
        self._mass = 91.2
        self._width = 2.5  # Decay width.
        self._mass2 = self._mass ** 2
        self._width2 = self._width ** 2

    @property
    def mass(self):
        return self._mass

    @property
    def mass2(self):
        return self._mass2

    @property
    def width(self):
        return self._width

    @property
    def width2(self):
        return self._width2
