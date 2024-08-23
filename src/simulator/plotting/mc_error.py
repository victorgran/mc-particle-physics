import matplotlib.pyplot as plt
import numpy as np

from .colors import *

scatter_colors = [french_blue, midnight_blue, purple]
scatter_shapes = ["^", "x", None]


# TODO: Implement setting font sizes according to how the image will be resized in a document.
def plotMonteCarloErrors(sample_sizes: list[int],
                         mc_errors: list[float] | list[np.ndarray],
                         labels: list[str | None] = None):

    fig, ax = plt.subplots(figsize=(7, 5))

    labels = labels if labels is not None else [None for _ in range(len(mc_errors))]

    for (series, label, color, marker) in zip(mc_errors, labels, scatter_colors, scatter_shapes):
        ax.scatter(sample_sizes, series, label=label, color=color, zorder=5, marker=marker)
        linear_params = np.polyfit(np.log10(sample_sizes), np.log10(series), deg=1)
        line_equation = np.poly1d(linear_params)
        ax.plot(sample_sizes, np.power(10, line_equation(np.log10(sample_sizes))),
                label=fr"$Slope: {linear_params[0]:.2f} $", color=red, linewidth=2.2, zorder=2)

    ax.set_xlabel(r"Sample size $N$")
    ax.set_ylabel("MC error estimate")
    ax.set_xscale('log')
    ax.set_yscale('log')

    plt.legend()
    plt.grid(linestyle='--', linewidth=0.8, zorder=0)
    plt.tight_layout()

    return ax
