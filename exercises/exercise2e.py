import matplotlib.pyplot as plt
from numpy import pi
from pathlib import Path
from tqdm import tqdm

from simulator.constants import alpha_QCD_MZ
from simulator.integrator import MonteCarloIntegrator, ZBoson
from simulator.plotting import plotJetHistograms, setFontSizes
from simulator.utils import AlphaS, Analysis, Shower

from exercise2b import formEvents


def runAnalysis():
    # Set QCD running coupling and energy scale.
    Z = ZBoson()
    alphas = AlphaS(Z.mass2, alpha_QCD_MZ)
    s = Z.mass2

    integrator = MonteCarloIntegrator(sum_quark_method="random")
    integrator.sampleDeltaSigma(sample_size)

    events = formEvents(integrator)

    # Weights for each event.
    d_sigmas = integrator.d_sigmas
    weights = d_sigmas * 4 * pi

    shower = Shower(alphas)
    analysis = Analysis()

    for event, weight in zip(tqdm(events), weights):
        shower.run(event, t=s)  # Running the shower modifies 'event' in-place.
        analysis.analyze(event, weight)

    analysis.finalize(file_name=data_path)

    return


def main():
    if not load_data:
        runAnalysis()

    sherpa_yoda = data_dir + "sherpa.yoda"  # Path to sherpa.yoda
    filenames = [sherpa_yoda, data_path + ".yoda"]
    colors = ["red", "blue"]
    labels = ["Sherpa", f"Simulator {sample_size:,.0f} events"]

    setFontSizes(font_size_div_factor, all_equal=True)

    # from simulator.utils import plot_jet_histograms
    # plot_jet_histograms(filenames)
    # plt.show()

    plots = plotJetHistograms(filenames, single_plot=False, colors=colors, labels=labels)

    for k, (fig, ax) in enumerate(plots):
        if save_figs:
            fig.savefig(figs_dir + f"jet_histo_{k + 1}.png", dpi=dpi, bbox_inches='tight')

        plt.show(block=True)

    return


if __name__ == '__main__':
    # NOTE: In order to ensure reproducibility, a seeded random generator has
    # been declared as a global variable within the utils/shower.py script.

    repo_dir = str(Path(__file__).parent.parent)  # Path to the repository directory.
    data_dir = repo_dir + "/data/"  # Path to data/
    figs_dir = repo_dir + "/figures/"  # Path to figures/

    sample_size = 1_000_000  # Number of events.

    load_data = True  # If true, plot(s) are generated from existing data.
    data_name = f"ex2e-analysis_{sample_size:_.0f}"  # By default, saved to data/
    data_path = data_dir + data_name

    save_figs = False
    font_size_div_factor = 0.7  # Adjust font sizes for visibility in document.
    dpi = 400

    main()
