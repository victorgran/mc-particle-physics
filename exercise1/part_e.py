import gvar
import matplotlib.pyplot as plt
import numpy as np
import vegas

from simulator import squared_matrix_element, ZBoson, plotGrid, setFontSizes
from simulator.constants import *


def integrand(x):
    s, cos_theta, phi = x
    d_sigma = np.sum([squared_matrix_element(s, cos_theta, q) for q in range(N_q)])
    d_sigma *= f_conv / (64 * (np.pi ** 2) * s * (s_max - s_min))

    return d_sigma


def main() -> None:
    s_lims = [s_min, s_max]
    cos_lims = [-1, 1]
    phi_lims = [0, 2 * np.pi]

    gvar.ranseed(seed=42)  # Set the seed of the vegas' integrator.
    integrator = vegas.Integrator([s_lims, cos_lims, phi_lims])
    sigma_eval = integrator(integrand, nitn=10, neval=1000)

    print(sigma_eval.summary())

    adaptive_map = integrator.map
    grids = np.asarray(adaptive_map.grid)
    n_inc = np.asarray(adaptive_map.ninc)

    ids = [[0, 1], [0, 2], [1, 2]]
    labels = [r"$s$", r"$\cos{\theta}$", r"$\phi$"]

    for k in range(3):
        setFontSizes(factor=fig_div_factor, all_equal=save_figures)
        ax = plotGrid(grids, n_inc, ids[k], xlabel=labels[ids[k][0]], ylabel=labels[ids[k][1]])

        if save_figures:
            plt.savefig(f"figures/integration_grid_{k}.png", dpi=dpi, bbox_inches='tight')
            print(f"Figure saved to figures/integration_grid_{k}.png.")

        ax.set_title("Integration Grid")
        plt.show()

    return


if __name__ == '__main__':
    Z = ZBoson()
    s_min = (Z.mass - 3 * Z.width) ** 2
    s_max = (Z.mass + 3 * Z.width) ** 2

    save_figures = False
    fig_div_factor = 0.7
    dpi = 400
    main()
