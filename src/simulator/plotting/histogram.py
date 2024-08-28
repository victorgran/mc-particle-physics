import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
from numpy import histogram


def plotHistogram(s_values, weights, s_range, sigma_vals, num_bins: int, s_min: float, s_max: float):
    fig, ax = plt.subplots(layout='constrained', figsize=(7, 5))
    counts, bins = histogram(s_values, bins=num_bins, range=(s_min, s_max), weights=weights)
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
    return
