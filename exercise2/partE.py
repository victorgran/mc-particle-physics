import matplotlib.pyplot as plt
from numpy import pi
from pathlib import Path
from tqdm import tqdm

from simulator.constants import alpha_QCD_MZ
from simulator.integrator import MonteCarloIntegrator, ZBoson
from simulator.utils import AlphaS, Analysis, Shower, plot_jet_histograms

from partB import formEvents


def runAnalysis(data_path: str):
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
    data_dir = str(Path(__file__).parent.parent) + "/data/"  # Path to data/
    data_path = data_dir + file_name
    if load_data:
        pass
    else:
        runAnalysis(data_path)

    sherpa_yoda = data_dir + "sherpa.yoda"  # Path to sherpa.yoda
    plot_jet_histograms([sherpa_yoda, data_path + ".yoda"])
    plt.show()

    return


if __name__ == '__main__':
    # NOTE: In order to ensure reproducibility, a seeded random generator has
    # been declared as a global variable within the utils/shower.py script.
    sample_size = 100_000
    load_data = True
    file_name = f"ex2e-analysis_{sample_size:_.0f}"  # By default, data is saved to data/
    save_fig = False
    main()
