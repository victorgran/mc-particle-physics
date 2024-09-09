from .integrator.distributions import BreitWigner, Dirac, Uniform
from .integrator.particles import ZBoson

Z = ZBoson()


def getDistros(part: str):
    """
    Convenience function to retrieve the distributions of 's' and 'f(s)'
    for the corresponding part of the problem 1.

    Parts "a" and "b" correspond to a Dirac distribution in both cases.
    Parts "c" and "d" correspond to a uniform distribution for both.
    Parts "f" and "g" give a uniform distribution for 'f(s)' and a Breit-Wigner
    distribution for 's'.

    :param part: Part of the problem (a, b, c, d, f or g).
    :return: Distributions of 's' and 'f(s)', respectively.
    """

    parts = ["a", "b", "c", "d", "f", "g"]
    if part not in parts:
        exit(f"Part not valid. Options are: {parts}")

    if part in ["a", "b"]:
        s_distro = Dirac(Z.mass2)
        beam_distro = s_distro
    else:
        s_min = (Z.mass - 3 * Z.width) ** 2
        s_max = (Z.mass + 3 * Z.width) ** 2
        beam_distro = Uniform(s_min, s_max)
        if part in ["c", "d"]:
            s_distro = beam_distro
        else:
            s_distro = BreitWigner(s_min, s_max, Z.mass, Z.width)

    return s_distro, beam_distro
