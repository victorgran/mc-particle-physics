import matplotlib.pyplot as plt


def setFontSizes(factor: float = 1., all_equal: bool = False):

    SMALL_SIZE = 8 / factor
    MEDIUM_SIZE = 10 / factor
    BIGGER_SIZE = 12 / factor

    if all_equal:
        SMALL_SIZE = BIGGER_SIZE
        MEDIUM_SIZE = BIGGER_SIZE

    plt.rc('font', size=SMALL_SIZE)  # Default text sizes.
    plt.rc('axes', titlesize=SMALL_SIZE)  # Axes title.
    plt.rc('axes', labelsize=MEDIUM_SIZE)  # x and y axes labels.
    plt.rc('xtick', labelsize=SMALL_SIZE)  # x ticks labels.
    plt.rc('ytick', labelsize=SMALL_SIZE)  # y ticks labels.
    plt.rc('legend', fontsize=SMALL_SIZE)  # Legend.
    plt.rc('figure', titlesize=BIGGER_SIZE)  # Figure title.

    return
