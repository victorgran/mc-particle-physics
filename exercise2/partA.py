import matplotlib.pyplot as plt
import numpy as np

from simulator import alpha_QCD_MZ, ZBoson, setFontSizes, AlphaS


def main() -> None:
    Z = ZBoson()
    alpha = AlphaS(Z.mass2, alpha_QCD_MZ)

    alpha_vals = [alpha(t) for t in t_values]

    setFontSizes(factor=font_size_f, all_equal=False)
    fig, ax = plt.subplots(layout='constrained', figsize=(7, 5))
    ax.plot(t_values, alpha_vals, linewidth=1.7)
    ax.set_xlabel(r"$t \; [\mathrm{GeV}^{2}]$")
    ax.set_xscale('log')
    ax.set_ylabel(r"$\alpha_{s}$")
    plt.grid(linestyle='--', linewidth=0.8)

    if save_figure:
        plt.savefig(figure_path, dpi=dpi, bbox_inches='tight')
        print(f"Figure saved to {figure_path}.")

    plt.title("QCD coupling at different scales")
    plt.show()

    return


if __name__ == '__main__':
    t_min_exp = 1  # Minimum base 10 exponent in GeV2.
    t_max_exp = 4
    t_points = 300
    t_values = np.logspace(t_min_exp, t_max_exp, num=t_points)

    save_figure = True
    figure_path = "figures/exercise2a.png"
    font_size_f = 0.6  # Adjust font sizes for visibility in document.
    dpi = 400
    main()
