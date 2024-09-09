"""
Create a single plot of the Monte Carlo error estimate vs. sample size for three different scenarios:
    1. fixed beam energy,
    2. flat beam spectrum with uniform sampling for s, and
    3. flat beam spectrum with Breit-Wigner importance sampling for s.

The plot includes two lines representing the decrease of the MC error estimate with 1/sqrt(N).
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path

from simulator import MonteCarloIntegrator, getDistros, plotMonteCarloErrors, setFontSizes


def getMCErrors(sample_sizes: np.ndarray | list[int], integrators: list[MonteCarloIntegrator]):
    """
    Compute the Monte Carlo error estimate on the total cross-section, for different sample sizes N.

    :param sample_sizes: List or numpy array of integers for the sample sizes.
    :param integrators: List of Monte Carlo integrators.
    :return: Ordered sample sizes and a list with the corresponding Monte Carlo errors for each integrator.
    """

    N_sizes = sorted(sample_sizes)  # Make sure it's in ascending order.

    mc_errors = [np.zeros(len(N_sizes)) for _ in range(len(integrators))]

    for (k, N) in enumerate(N_sizes):
        for (idx, integrator) in enumerate(integrators):
            integrator: MonteCarloIntegrator
            integrator.sampleDeltaSigma(N)
            mc_error = integrator.integrateCrossSection()[1]
            mc_errors[idx][k] = mc_error

    return N_sizes, mc_errors


def main() -> None:
    cases = ["Fixed", "Uniform", "Breit-Wigner"]

    if load_data:
        df = pd.read_csv(data_path)
        N_sizes = list(df["sample_sizes"])
        mc_errors = [list(df[key]) for key in cases]
    else:
        distros_b = getDistros("b")
        distros_d = getDistros("d")
        distros_g = getDistros("g")

        integrator_b = MonteCarloIntegrator(s_distro=distros_b[0], beam_distro=distros_b[1], sum_quark_method="random")
        integrator_d = MonteCarloIntegrator(s_distro=distros_d[0], beam_distro=distros_d[1], sum_quark_method="random")
        integrator_g = MonteCarloIntegrator(s_distro=distros_g[0], beam_distro=distros_g[1], sum_quark_method="random")

        N_sizes, mc_errors = getMCErrors(N_values, [integrator_b, integrator_d, integrator_g])

        if save_data:
            mc_error_data = {cases[k]: mc_errors[k] for k in range(len(cases))}
            mc_error_data["sample_sizes"] = N_sizes
            df = pd.DataFrame(mc_error_data)
            df.to_csv(data_path, index=False)

    setFontSizes(factor=font_size_div_factor, all_equal=save_figure)
    plotMonteCarloErrors(N_sizes, mc_errors,
                         labels=["Fixed", "Uniform", "Breit-Wigner"],
                         make_lines=[False, True, True])
    if save_figure:
        plt.savefig(figure_path, dpi=dpi, bbox_inches='tight')
        print(f"Figure saved to {figure_path}.")

    plt.show()
    return


if __name__ == '__main__':
    repo_dir = str(Path(__file__).parent.parent)  # Path to the repository directory.
    data_dir = repo_dir + "/data/"  # Path to data/
    figs_dir = repo_dir + "/figures/"  # Path to figures/

    ##################
    # Data variables #
    ##################
    load_data = True  # If true, the plot is generated from existing data.
    save_data = False  # If true, generated data for the plot is saved. Ignored if data is loaded.
    data_name = "ex1_all_mc_errors.csv"  # Filename for Monte Carlo error estimates data.
    data_path = data_dir + data_name  # Path to either load or save data. By default, saved or loaded from data/

    # Get logarithmically-equally-spaced points for sample sizes. Log is base 10.
    N_values = np.logspace(1, 7, num=30, dtype=int)

    ####################
    # Figure variables #
    ####################
    save_figure = False  # Save plot of Monte Carlo error estimate vs. sample size.
    figure_name = "ex1_all_mc_errors.png"
    font_size_div_factor = 0.8  # Adjust font sizes for visibility in document.
    dpi = 400  # Dots per inch for saving the figure.
    figure_path = figs_dir + figure_name  # Path to save the figure. By default, saved to figures/

    main()
