import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
import pandas as pd

from all_cross_sections import compareQuarkSumMethods, showResults
from all_mc_errors import getMCErrors
from simulator import getDistros, true_sigma_1c, MonteCarloIntegrator, setFontSizes, \
    plotMonteCarloErrors, ZBoson, Dirac


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


def plotHistogram() -> None:
    Z = ZBoson()
    s_min = (Z.mass - 3 * Z.width) ** 2
    s_max = (Z.mass + 3 * Z.width) ** 2

    flat_distros = getDistros("d")
    integrator = MonteCarloIntegrator(s_distro=flat_distros[0], beam_distro=flat_distros[1], sum_quark_method="random")
    integrator.sample_diff_cross_section(histo_N)
    s_values = integrator.s_samples
    d_sigmas = integrator.d_sigmas

    weights = d_sigmas * (4 * np.pi * num_bins) / histo_N  # Area of histo == cross-section.
    weights *= (s_max - s_min)  # Height of histo == cross-section at s.
    counts, bins = np.histogram(s_values, bins=num_bins, range=(s_min, s_max), weights=weights)

    # Scan 's' with fixed integrator.
    s_range = np.linspace(s_min, s_max, s_points)
    sigma_vals = []
    for s in s_range:
        s_distribution = Dirac(s)
        integrator = MonteCarloIntegrator(s_distro=s_distribution, beam_distro=s_distribution,
                                          sum_quark_method="random")
        integrator.sample_diff_cross_section(histo_N)
        sigma, _ = integrator.integrate_cross_section()
        sigma_vals.append(sigma)

    sigma_vals = np.array(sigma_vals)

    setFontSizes(factor=font_size_div_factor2, all_equal=save_histo)
    fig, ax = plt.subplots(layout='constrained', figsize=(7, 5))

    ax.plot(s_range, sigma_vals, color="#ED2939", label=r"$f(s) = \delta(s - s_{0})$", linewidth=2.1, zorder=1)
    ax.stairs(counts, bins, color="#0072bb", label=r"Uniform $f(s)$", linewidth=1.75, zorder=2)

    # Vertical line at Z mass squared.
    y_max = max(sigma_vals)
    ax.vlines(91.2 ** 2, 0, y_max, linestyles='dashed', colors="#058743", label=r"$M_{Z}^{2}$")

    ax.set_xlabel(r"$s \; [\text{GeV}^{2}]$")
    ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))

    ax.set_ylabel(r"$\sigma$ [pb]")
    ax.yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
    plt.legend()

    if save_histo:
        plt.savefig(histo_plot_path, dpi=dpi, bbox_inches='tight')
        print(f"Figure saved to {histo_plot_path}.")

    plt.show()

    return


def main() -> None:

    # Confirm that your Monte Carlo estimate for sigma with the flat beam spectrum is compatible with 9880(50) pb.
    print("\nCross-section results\n")
    s_distro, beam_distro = getDistros(part="d")
    results_b = compareQuarkSumMethods(s_distro, beam_distro, sample_size)
    showResults(results_b, true_sigma_1c, title="Fixed beam spectrum.")

    plotErrorEstimate()
    plotHistogram()

    return


if __name__ == '__main__':
    sample_size = 100_000  # Number of MC points to compute the cross-section.

    #########################################
    # MC errors for different sample sizes. #
    #########################################
    load_data = True  # If true, the plot is generated from existing data.
    save_data = False  # If true, generated data for the plot is saved. Ignored if data is loaded.
    data_path = "data/ex1_part_d.csv"  # Path to either load or save data.
    # Get logarithmically-equally-spaced points for sample sizes. Log is base 10.
    N_values = np.logspace(1, 7, num=30, dtype=int)

    save_figure = False  # Save plot of Monte Carlo error estimate vs. sample size.
    font_size_div_factor = 0.7  # Adjust font sizes for visibility in document.
    dpi = 400
    figure_path = "figures/ex1_part_d.png"

    ################
    # s histogram. #
    ################
    num_bins = 51
    histo_N = 200_000  # Sample size to create histogram.
    s_points = 250  # Number of points to scan, with the fixed beam energy integrator, between s_min and s_max.

    save_histo = False
    font_size_div_factor2 = 0.7  # Adjust font sizes for visibility in document.
    histo_plot_path = "figures/ex1_part_d_histo.png"
    main()
