import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

from all_cross_sections import compareQuarkSumMethods, showResults
from all_mc_errors import getMCErrors
from simulator import getDistros, true_sigma_1c, MonteCarloIntegrator, setFontSizes, ZBoson, Dirac
from simulator.plotting import plotHistogram, plotMonteCarloErrors


def plotErrorEstimate() -> None:
    if load_data:
        df = pd.read_csv(data_path)
        sample_sizes = list(df["sample_sizes"])
        mc_errors = [np.array(df["mc_errors"])]
    else:
        s_distro, beam_distro = getDistros("d")
        integrator = MonteCarloIntegrator(s_distro=s_distro, beam_distro=beam_distro, sum_quark_method="random")
        sample_sizes, mc_errors = getMCErrors(N_values, [integrator])
        if save_data:
            df = pd.DataFrame({"sample_sizes": sample_sizes, "mc_errors": mc_errors[0]})
            df.to_csv(data_path, index=False)

    setFontSizes(factor=font_size_div_factor, all_equal=save_figure)
    _ = plotMonteCarloErrors(sample_sizes, mc_errors, make_lines=[True], labels=["Fixed beam"])

    if save_figure:
        plt.savefig(figure_path, dpi=dpi, bbox_inches='tight')
        print(f"Figure saved to {figure_path}.")

    plt.show()

    return


def sHistogram() -> None:
    Z = ZBoson()
    s_min = (Z.mass - 3 * Z.width) ** 2
    s_max = (Z.mass + 3 * Z.width) ** 2

    flat_distros = getDistros("d")
    integrator = MonteCarloIntegrator(s_distro=flat_distros[0], beam_distro=flat_distros[1], sum_quark_method="random")
    integrator.sampleDeltaSigma(histo_N)
    s_values = integrator.s_samples
    d_sigmas = integrator.d_sigmas

    weights = d_sigmas * (4 * np.pi * num_bins) / histo_N  # Area of histo == cross-section.
    weights *= (s_max - s_min)  # Height of histo == cross-section at s.

    # Scan 's' with fixed integrator.
    s_range = np.linspace(s_min, s_max, s_points)
    sigma_vals = []
    for s in s_range:
        s_distribution = Dirac(s)
        integrator = MonteCarloIntegrator(s_distro=s_distribution, sum_quark_method="random")
        integrator.sampleDeltaSigma(histo_N)
        sigma, _ = integrator.integrateCrossSection()
        sigma_vals.append(sigma)

    sigma_vals = np.array(sigma_vals)

    setFontSizes(factor=font_size_div_factor2, all_equal=save_histo)
    plotHistogram(s_values, weights, s_range, sigma_vals, num_bins, s_min, s_max)

    if save_histo:
        plt.savefig(histo_path, dpi=dpi, bbox_inches='tight')
        print(f"Figure saved to {histo_path}.")

    plt.show()

    return


def main() -> None:

    # Confirm that your Monte Carlo estimate for sigma with the flat beam spectrum is compatible with 9880(50) pb.
    print("\nCross-section results\n")
    s_distro, beam_distro = getDistros(part="d")
    results_b = compareQuarkSumMethods(s_distro, beam_distro, sample_size)
    showResults(results_b, true_sigma_1c, title="Fixed beam spectrum.")

    plotErrorEstimate()
    sHistogram()

    return


if __name__ == '__main__':
    repo_dir = str(Path(__file__).parent.parent)  # Path to the repository directory.
    data_dir = repo_dir + "/data/"  # Path to data/
    figs_dir = repo_dir + "/figures/"  # Path to figures/

    # ------------------------------------------------------------------------------------------------------------
    sample_size = 100_000  # Number of MC points to compute the cross-section.

    # ------------------------------------------------------------------------------------------------------------
    # MC errors for different sample sizes
    load_data = True  # If true, the plot is generated from existing data.
    save_data = False  # If true, generated data for the plot is saved. Ignored if data is loaded.
    data_name = "ex1_part_d.csv"  # Filename for Monte Carlo error estimates data.
    data_path = data_dir + data_name  # Path to either load or save data.

    # Get logarithmically-equally-spaced points for sample sizes. Log is base 10.
    N_values = np.logspace(1, 7, num=30, dtype=int)

    save_figure = False  # Save plot of Monte Carlo error estimate vs. sample size.
    figure_name = "ex1_part_d.png"
    font_size_div_factor = 0.7  # Adjust font sizes for visibility in document.
    dpi = 400  # Dots per inch for saving the figure.
    figure_path = figs_dir + figure_name  # Path to save the figure. By default, saved to figures/

    # ------------------------------------------------------------------------------------------------------------
    # s histogram
    num_bins = 51
    histo_N = 200_000  # Sample size to create histogram.
    s_points = 250  # Number of points to scan, with the fixed beam energy integrator, between s_min and s_max.

    save_histo = False
    histo_name = "ex1_part_d_histo.png"
    font_size_div_factor2 = 0.7  # Adjust font sizes for visibility in document.
    histo_path = figs_dir + histo_name  # Path to save the figure. By default, saved to figures/

    main()
