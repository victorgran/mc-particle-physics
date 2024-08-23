import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from simulator import MonteCarloIntegrator, true_sigma_1a, getDistros, plotMonteCarloErrors
from all_cross_sections import compareQuarkSumMethods, showResults


def getMCErrors(sample_sizes: np.ndarray | list[int]):
    """
    Compute the Monte Carlo error estimate on the total cross-section, for different sample sizes N.

    :return:
    """

    N_sizes = sorted(sample_sizes)  # Make sure it's in ascending order.

    mc_errors = np.array([])
    integrator = MonteCarloIntegrator(sum_quark_method="random")

    for N in N_sizes:
        integrator.sample_diff_cross_section(N)
        mc_error = integrator.integrate_cross_section()[1]
        mc_errors = np.append(mc_errors, mc_error)

    return N_sizes, mc_errors


def main() -> None:
    ###############################################################################
    # MC estimates are statistically equivalent for both flavour summing options. #
    ###############################################################################
    print("\nCross-section results\n")
    s_distro, beam_distro = getDistros(part="b")
    results_b = compareQuarkSumMethods(s_distro, beam_distro, sample_size)
    showResults(results_b, true_sigma_1a, title="Fixed beam spectrum.")

    if load_data:
        df = pd.read_csv(data_path)
        sample_sizes = list(df["sample_sizes"])
        mc_errors = list(df["mc_errors"])
    else:
        sample_sizes, mc_errors = getMCErrors(N_values)
        if save_data:
            df = pd.DataFrame({"sample_sizes": sample_sizes, "mc_errors": mc_errors})
            df.to_csv(data_path, index=False)

    ###########################
    # MC error-estimate plot. #
    ###########################
    # TODO: Plot sample sizes vs. mc_errors.
    dpi = 400
    _ = plotMonteCarloErrors(sample_sizes, [mc_errors])

    if save_figure:
        plt.savefig(figure_path, dpi=dpi, bbox_inches='tight')
        print(f"Figure saved to {figure_path}.")

    plt.show()

    return


if __name__ == '__main__':
    sample_size = 100_000  # Number of MC points to compute the cross-section.

    #########################################
    # MC errors for different sample sizes. #
    #########################################
    load_data = True  # If true, the plot is generated from existing data.
    save_data = True  # If true, generated data for the plot is saved.
    data_path = "data/ex1_part_b.csv"  # Path to either load or save data.
    # Get logarithmically-equally-spaced points for sample sizes. Log is base 10.
    N_values = np.logspace(1, 7, num=30, dtype=int)

    save_figure = True  # Save plot of Monte Carlo error estimate vs. sample size.
    figure_path = "figures/ex1_part_b.png"
    main()
