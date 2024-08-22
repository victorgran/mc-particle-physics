import numpy as np

from simulator import MonteCarloIntegrator, true_sigma_1a
# from plotting import plot_mc_error


def main() -> None:
    ###############################################################################
    # MC estimates are statistically equivalent for both flavour summing options. #
    ###############################################################################
    integrator1 = MonteCarloIntegrator(sum_quark_method="explicit")
    integrator1.sample_diff_cross_section(sample_size)
    sigma1, mc_err1 = integrator1.integrate_cross_section()

    integrator2 = MonteCarloIntegrator(sum_quark_method="random")
    integrator2.sample_diff_cross_section(sample_size)
    sigma2, mc_err2 = integrator2.integrate_cross_section()

    print(f"Estimates for {sample_size:,d} Monte Carlo points.")
    print(f"Explicit sum result: {sigma1:,.2f} +/- {mc_err1:,.2f}")
    print(f"Rand flavour result: {sigma2:,.2f} +/- {mc_err2:,.2f}")
    print(f"Exact cross-section: {true_sigma_1a:,.2f}")

    ###########################
    # MC error-estimate plot. #
    ###########################

    # Make logarithmically-equally-spaced integer values for N.
    N_min_exp = 1  # Minimum power of 10 for N.
    N_max_exp = 7
    N_points = 30
    N_values = np.logspace(N_min_exp, N_max_exp, num=N_points, dtype=int)

    # Compute the Monte Carlo error at each value of N.
    mc_errors = np.array([])
    integrator = MonteCarloIntegrator(sum_quark_method="random")

    for N in N_values:
        integrator.sample_diff_cross_section(N)
        mc_error = integrator.integrate_cross_section()[1]
        mc_errors = np.append(mc_errors, mc_error)

    # plot_mc_error(N_values, mc_errors,
    #               save_fig=save_figure,
    #               save_name="./figures/exercise1b_plot.png"
    #               )

    return


if __name__ == '__main__':
    sample_size = 100_000  # Number of MC points to compute the cross-section.
    save_figure = False  # Save plot of Monte Carlo error estimate vs. sample size.
    main()
