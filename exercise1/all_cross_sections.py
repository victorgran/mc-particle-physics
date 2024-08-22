"""
Compute the cross-sections for the 3 different cases discussed in the project sheet.

Each cross-section is computed with two methods for the sum
over quark flavours: explicit sum and random flavour method.

The estimated values for the cross-section are given with their corresponding
error estimate and the exact value is provided as a reference.
"""

from simulator import MonteCarloIntegrator, getDistros, true_sigma_1a, true_sigma_1c


def showResults(results, true_value: float, title: str):
    sigma1, mc_err1, sigma2, mc_err2, sample_size = results
    print(title)
    print(f"Estimates for {sample_size:,d} Monte Carlo points.")
    print(f"Explicit sum result: {sigma1:,.2f} +/- {mc_err1:,.2f}")
    print(f"Rand flavour result: {sigma2:,.2f} +/- {mc_err2:,.2f}")
    print(f"Exact cross-section: {true_value:,.2f}")
    print()
    return


def compareQuarkSumMethods(s_distribution, beam_distribution, sample_size: int):
    integrator1 = MonteCarloIntegrator(s_distro=s_distribution,
                                       beam_distro=beam_distribution,
                                       sum_quark_method="explicit")
    integrator1.sample_diff_cross_section(sample_size)
    sigma1, mc_err1 = integrator1.integrate_cross_section()

    integrator2 = MonteCarloIntegrator(s_distro=s_distribution,
                                       beam_distro=beam_distribution,
                                       sum_quark_method="random")
    integrator2.sample_diff_cross_section(sample_size)
    sigma2, mc_err2 = integrator2.integrate_cross_section()

    return [sigma1, mc_err1, sigma2, mc_err2, sample_size]


def main() -> None:
    sample_size = 100_000
    print("\nCross-section results\n")
    s_distro, beam_distro = getDistros(part="b")
    results_b = compareQuarkSumMethods(s_distro, beam_distro, sample_size)
    showResults(results_b, true_sigma_1a, title="Fixed beam spectrum.")

    s_distro, beam_distro = getDistros(part="d")
    results_d = compareQuarkSumMethods(s_distro, beam_distro, sample_size)
    showResults(results_d, true_sigma_1c, title="Flat beam spectrum.")

    s_distro, beam_distro = getDistros(part="g")
    results_g = compareQuarkSumMethods(s_distro, beam_distro, sample_size)
    showResults(results_g, true_sigma_1c,
                title="Flat beam spectrum and Breit-Wigner distribution for s.")
    return


if __name__ == '__main__':
    main()
