import math as m

from .vector import Vec4


class Particle:
    """A simple particle class."""

    def __init__(self, particle_num, momentum, color=None):
        """Initializes a particle given its particle number, a momentum, and a
        2-component list giving its color and its anti-color, where 0 stand for
        "no color"."""
        self.pid = particle_num
        self.mom = momentum
        self.color = color if color is not None else [0, 0]

    def __repr__(self):
        return "{0} {1} {2}".format(self.pid, self.mom, self.color)

    def __str__(self):
        return "{0} {1} {2}".format(self.pid, self.mom, self.color)

    def set(self, particle_num, momentum, color=None):
        self.pid = particle_num
        self.mom = momentum
        self.color = color if color is not None else [0, 0]

    def is_color_connected(self, other):
        """Checks if this and some other particle are "color-connected",
        i.e. a color in one of the two particles must find a corresponding
        anti-color in the other of the two particles."""
        return (self.color[0] > 0 and self.color[0] == other.color[1]) or \
               (self.color[1] > 0 and self.color[1] == other.color[0])


def check_event(event):
    """Checks momentum and color conservation in an event (= list of Particle
    instances)."""
    p_sum = Vec4()
    csum = {}
    for p in event:
        p_sum += p.mom
        if p.color[0] > 0: 
            csum[p.color[0]] = csum.get(p.color[0], 0) + 1
            if csum[p.color[0]] == 0:
                del csum[p.color[0]]
        if p.color[1] > 0:
            csum[p.color[1]] = csum.get(p.color[1], 0) - 1
            if csum[p.color[1]] == 0:
                del csum[p.color[1]]
    return (m.fabs(p_sum.E) < 1.e-12 and m.fabs(p_sum.px) < 1.e-12 and
            m.fabs(p_sum.py) < 1.e-12 and m.fabs(p_sum.pz) < 1.e-12 and
            len(csum) == 0
            )
