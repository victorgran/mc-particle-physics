import math as m


class Vec4:

    def __init__(self, E=0., px=0., py=0., pz=0.):
        self.E = E
        self.px = px
        self.py = py
        self.pz = pz

    def __getitem__(self, i):
        if i == 0:
            return self.E
        if i == 1:
            return self.px
        if i == 2:
            return self.py
        if i == 3:
            return self.pz
        raise Exception('Vec4D')

    def __repr__(self):
        return '({0},{1},{2},{3})'.format(self.E, self.px, self.py, self.pz)

    def __str__(self):
        return '({0},{1},{2},{3})'.format(self.E, self.px, self.py, self.pz)

    def __add__(self, v):
        return Vec4(self.E + v.E, self.px + v.px, self.py + v.py, self.pz + v.pz)

    def __sub__(self, v):
        return Vec4(self.E - v.E, self.px - v.px, self.py - v.py, self.pz - v.pz)

    def __neg__(self):
        return Vec4(-self.E, -self.px, -self.py, -self.pz)

    def __mul__(self, v):
        if isinstance(v, Vec4):
            return self.E * v.E - self.px * v.px - self.py * v.py - self.pz * v.pz
        return Vec4(self.E * v, self.px * v, self.py * v, self.pz * v)

    def __rmul__(self, v):
        if isinstance(v, Vec4):
            return self.E * v.E - self.px * v.px - self.py * v.py - self.pz * v.pz
        return Vec4(self.E * v, self.px * v, self.py * v, self.pz * v)

    def __div__(self, v):
        return Vec4(self.E / v, self.px / v, self.py / v, self.pz / v)

    def invariant_mass_squared(self):
        return self * self

    def invariant_mass(self):
        return m.sqrt(self.invariant_mass_squared())

    def length_3d_squared(self):
        return self.px * self.px + self.py * self.py + self.pz * self.pz

    def length_3d(self):
        return m.sqrt(self.length_3d_squared())

    def transverse_momentum_squared(self):
        return self.px * self.px + self.py * self.py

    def transverse_momentum(self):
        return m.sqrt(self.transverse_momentum_squared())

    def phi(self):
        if self.px == 0 and self.py == 0:
            return 0.0
        else:
            return m.atan2(self.py, self.px)

    def cross_product(self, v):
        return Vec4(0.0,
                    self.py * v.pz - self.pz * v.py,
                    self.pz * v.px - self.px * v.pz,
                    self.px * v.py - self.py * v.px)

    def boost(self, v):
        rsq = self.invariant_mass()
        v0 = (self.E * v.E - self.px * v.px - self.py * v.py - self.pz * v.pz) / rsq
        c1 = (v.E + v0) / (rsq + self.E)
        return Vec4(v0,
                    v.px - c1 * self.px,
                    v.py - c1 * self.py,
                    v.pz - c1 * self.pz)

    def boost_back(self, v):
        rsq = self.invariant_mass()
        v0 = (self.E * v.E + self.px * v.px + self.py * v.py + self.pz * v.pz) / rsq
        c1 = (v.E + v0) / (rsq + self.E)
        return Vec4(v0,
                    v.px + c1 * self.px,
                    v.py + c1 * self.py,
                    v.pz + c1 * self.pz)
