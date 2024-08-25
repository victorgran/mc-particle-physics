import numpy as np
from tqdm import tqdm

from simulator import AlphaS, MonteCarloIntegrator, Shower, ZBoson, alpha_QCD_MZ
from partB import formEvents


def main() -> None:
    # Set QCD running coupling and energy scale.
    Z = ZBoson()
    alphas = AlphaS(Z.mass2, alpha_QCD_MZ)
    s = Z.mass2

    integrator = MonteCarloIntegrator(sum_quark_method="random")
    integrator.sampleDeltaSigma(sample_size)

    events = formEvents(integrator)
    fs_particles = []

    for event in tqdm(events):
        Shower(alphas).run(event, t=s)  # Running the shower modifies 'event' in-place.
        out_particles = (len(event) - 2)  # Don't count the 2 incoming particles.
        fs_particles.append(out_particles)

    fs_particles = np.array(fs_particles)

    # Weights for each event.
    d_sigmas = integrator.d_sigmas
    sigma, _ = integrator.integrateCrossSection()
    weights = d_sigmas * 4 * np.pi / (sample_size * sigma)  # Normalized weights.

    # Compute average and error estimate.
    avg_fs_particles = np.average(fs_particles, weights=weights)
    err_fs_particles = np.sqrt(np.sum(weights * (fs_particles - avg_fs_particles) ** 2) / sample_size)
    print(f"Average number of final-state particles: {avg_fs_particles:.2f} +/- {err_fs_particles:.2f}")

    return


if __name__ == '__main__':
    # NOTE: In order to ensure reproducibility, a seeded random generator has
    # been declared as a global variable within the utils/shower.py script.
    sample_size = 1_000
    main()
