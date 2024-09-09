import numpy as np

from simulator import Vec4, Particle, MonteCarloIntegrator


def getOutQuarksInfo(integrator: MonteCarloIntegrator, s: float):
    """
    Get the flavours of the outgoing quarks, according to the Monte
    Carlo Particle Numbering Scheme, as well as their momenta.

    :param integrator:
    :param s: Squared centre-of-mass energy.
    :return: Quarks flavour ids and their momenta.
    """
    E = np.sqrt(s) / 2  # Half centre-of-mass energy.

    # Get flavour particle ids from 1 to 5.
    flavour_ids = integrator.quark_flavours_samples + 1

    # Form the momenta of the outgoing quarks.
    cos_values = integrator.cos_samples
    sin_values = np.sqrt(1. - (cos_values ** 2))
    phi_values = integrator.phi_samples
    q_momentum = np.stack((-E * np.cos(phi_values) * sin_values,
                           -E * np.sin(phi_values) * sin_values,
                           -E * cos_values),
                          axis=1)

    return flavour_ids, q_momentum


def formEvents(integrator: MonteCarloIntegrator):
    """
    Form a list of events consisting on four particles: an incoming
    electron-positron pair and an outgoing quark-antiquark pair.

    :param integrator: Monte Carlo integrator to retrieve
    :return:
    """

    # Get values of kinematic variables.
    s_values = integrator.s_samples
    cos_values = integrator.cos_samples
    sin_values = np.sqrt(1. - (cos_values ** 2))
    phi_values = integrator.phi_samples

    # Get flavour particle ids from 1 to 5.
    flavour_ids = integrator.quark_flavours_samples + 1

    events = []

    for (s, cos, sin, phi, flv) in zip(s_values, cos_values, sin_values, phi_values, flavour_ids):
        E = np.sqrt(s) / 2  # Half centre-of-mass energy.

        # Incoming electron-positron pair.
        electron = Particle(11, Vec4(E, 0., 0., +E), color=None)
        positron = Particle(-11, Vec4(E, 0., 0., -E), color=None)

        # Outgoing quark momentum.
        p = np.array([-E * np.cos(phi) * sin, -E * np.sin(phi) * sin, -E * cos])

        # Outgoing quark-antiquark pair
        quark = Particle(flv, Vec4(E, *p), color=[1, 0])
        antiquark = Particle(-flv, Vec4(E, *(-p)), color=[0, 1])

        event = [electron, positron, quark, antiquark]
        events.append(event)

    return events
