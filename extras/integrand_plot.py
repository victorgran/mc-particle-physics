"""
Plot the cross-section density in kinematic space, i.e., the differential cross-section as a
function of the kinematic variables 's' and 'cos_theta' (there is no dependence on 'phi').
"""

import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
import numpy as np

from simulator import ZBoson, N_q, f_conv, squared_matrix_element


def main() -> None:
    s_range = np.linspace(s_min, s_max, 200)
    cos_rng = np.linspace(-1, 1, 200)

    s, cos = np.meshgrid(s_range, cos_rng)
    d_sigma = [squared_matrix_element(s, cos, q) for q in range(N_q)]
    d_sigma = np.sum(d_sigma, axis=0)
    d_sigma = np.divide(d_sigma, s)
    d_sigma *= f_conv / (64 * (np.pi ** 2) * (s_max - s_min))

    # Plot the surface
    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(projection='3d')
    ax.view_init(elev=34, azim=-60, roll=0)
    ax.plot_surface(s, cos, d_sigma, rstride=10, cstride=10, alpha=0.35, edgecolor="#0072bb")
    ax.set_xlabel(r"$s \; [\text{GeV}^{2}]$", fontsize=13)
    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=15)
    ax.tick_params(axis='z', labelsize=15)
    ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
    ax.set_ylabel(r"$ \cos{\theta} $", fontsize=15)
    ax.set_zlabel(r"$ \mathrm{d}\sigma / \mathrm{dV}$", fontsize=15)
    fig.tight_layout()
    plt.savefig("figures/integrand_plot.png", dpi=400)
    plt.show()

    return


if __name__ == '__main__':
    Z = ZBoson()
    s_min = (Z.mass - 3 * Z.width) ** 2
    s_max = (Z.mass + 3 * Z.width) ** 2
    main()
