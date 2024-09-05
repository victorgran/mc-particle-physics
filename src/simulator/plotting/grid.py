import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
import numpy as np


def plotGrid(grid, n_inc, axes, n_grid=33, shrink=False, **kwargs):

    fig, ax = plt.subplots(figsize=(7, 5), layout='constrained')

    dx, dy = axes
    if dx is not None:
        n_skip = int(n_inc[dx] // n_grid)
        if n_skip < 1:
            n_skip = 1
        start = n_skip // 2
        xrange = [grid[dx, 0], grid[dx, n_inc[dx]]]
        xgrid = grid[dx, start::n_skip]
    else:
        xrange = [0., 1.]
        xgrid = None
    if dy is not None:
        n_skip = int(n_inc[dy] // n_grid)
        if n_skip < 1:
            n_skip = 1
        start = n_skip // 2
        yrange = [grid[dy, 0], grid[dy, n_inc[dy]]]
        y_grid = grid[dy, start::n_skip]
    else:
        yrange = [0., 1.]
        y_grid = None
    if shrink:
        if xgrid is not None:
            xrange = [min(xgrid), max(xgrid)]
        if y_grid is not None:
            yrange = [min(y_grid), max(y_grid)]
    if xgrid is not None:
        for i in range(len(xgrid)):
            plt.plot([xgrid[i], xgrid[i]], yrange, zorder=0, linewidth=0.85, color="#0072bb")
    if y_grid is not None:
        for i in range(len(y_grid)):
            plt.plot(xrange, [y_grid[i], y_grid[i]], zorder=0, linewidth=0.85, color="#0072bb")

    ax.set_xlim(*xrange)
    ax.set_ylim(*yrange)

    if 'xlabel' in kwargs:
        ax.set_xlabel(kwargs.get("xlabel"))
    if 'ylabel' in kwargs:
        ax.set_ylabel(kwargs.get("ylabel"))

    if dx == 0:
        ax.set_xlabel(r"$s \; [\mathrm{GeV}^{2}]$")
        ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
        # Line Z mass squared.
        ax.vlines(91.2 ** 2, *yrange, linewidth=2, colors="#ED2939", label=r"$M_{Z}^{2}$", zorder=2)

    if dy == 2:
        ax.set_yticks(np.arange(0, 2 * np.pi + 0.01, np.pi / 4))
        labels = ['$0$', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$', r'$\pi$',
                  r'$5\pi/4$', r'$3\pi/2$', r'$7\pi/4$', r'$2\pi$']
        ax.set_yticklabels(labels)

    return ax
