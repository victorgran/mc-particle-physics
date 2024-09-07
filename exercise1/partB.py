import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

from simulator import MonteCarloIntegrator, true_sigma_1a, getDistros
from simulator.plotting import setFontSizes, plotMonteCarloErrors
from all_cross_sections import compareQuarkSumMethods, showResults
from all_mc_errors import getMCErrors


def main() -> None:
    ###############################################################################
    # MC estimates are statistically equivalent for both flavour summing options. #
    ###############################################################################
    s_distro, beam_distro = getDistros(part="b")
    results_b = compareQuarkSumMethods(s_distro, beam_distro, sample_size)
    showResults(results_b, true_sigma_1a,
                title="Cross-section results for the fixed beam spectrum.")

    ###########################
    # MC error-estimate plot. #
    ###########################

    if load_data:
        df = pd.read_csv(data_path)
        sample_sizes = list(df["sample_sizes"])
        mc_errors = [np.array(df["mc_errors"])]
    else:
        s_distro, beam_distro = getDistros("b")
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
    data_name = "ex1_part_b.csv"  # Filename for Monte Carlo error estimates data.
    data_path = data_dir + data_name  # Path to either load or save data.

    # Get logarithmically-equally-spaced points for sample sizes. Log is base 10.
    N_values = np.logspace(1, 7, num=30, dtype=int)

    save_figure = False  # Save plot of Monte Carlo error estimate vs. sample size.
    figure_name = "ex1_part_b.png"
    font_size_div_factor = 0.7  # Adjust font sizes for visibility in document.
    dpi = 400  # Dots per inch for saving the figure.
    figure_path = figs_dir + figure_name  # Path to save the figure. By default, saved to figures/

    main()
