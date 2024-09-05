from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from .colors import *

scatter_colors = [french_blue, midnight_blue, purple]
scatter_shapes = ["^", "x", None]
line_colors = [red, orange, red]


def plotMonteCarloErrors(sample_sizes: list[int],
                         mc_errors: list[list[float]] | list[np.ndarray],
                         labels: list[str | None] = None, make_lines: list[bool] = None):
    """
    Plot the Monte Carlo error estimates at different sample sizes. Optionally, create a linear regression
    for any given series of MC error estimates and plot the resulting line labeled by its slope.

    :param sample_sizes: Sample sizes in ascending order.
    :param mc_errors: List of series of Monte Carlo error estimates.
    :param labels: Labels for each series of Monte Carlo error estimates.
    :param make_lines: Boolean array to specify if a linear regression is desired at a given series.
    :return: Axes of the figure.
    """

    fig, ax = plt.subplots(figsize=(7, 5))

    labels = labels if labels is not None else [None for _ in range(len(mc_errors))]
    make_lines = [False for _ in range(len(mc_errors))] if make_lines is None else make_lines

    for (series, label, color, marker, make_line, line_color) in zip(mc_errors, labels, scatter_colors,
                                                                     scatter_shapes, make_lines, line_colors):

        ax.scatter(sample_sizes, series, label=label, color=color, zorder=5, marker=marker)
        
        if not make_line:
            continue

        # Make linear regression and plot resulting line with its slope as label.
        linear_params = np.polyfit(np.log10(sample_sizes), np.log10(series), deg=1)
        line_equation = np.poly1d(linear_params)
        ax.plot(sample_sizes, np.power(10, line_equation(np.log10(sample_sizes))),
                label=fr"$Slope: {linear_params[0]:.2f} $", color=line_color, linewidth=2.2, zorder=2)

    ax.set_xlabel(r"Sample size $N$")
    ax.set_ylabel("MC error estimate")
    ax.set_xscale('log')
    ax.set_yscale('log')

    plt.legend()
    plt.grid(linestyle='--', linewidth=0.8, zorder=0)
    plt.tight_layout()

    return ax
