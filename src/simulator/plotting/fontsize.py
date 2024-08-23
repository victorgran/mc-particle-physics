import matplotlib.pyplot as plt


SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 12

plt.rc('font', size=SMALL_SIZE)          # Default text sizes.
plt.rc('axes', titlesize=SMALL_SIZE)     # Axes title.
plt.rc('axes', labelsize=MEDIUM_SIZE)    # x and y axes labels.
plt.rc('xtick', labelsize=SMALL_SIZE)    # x ticks labels.
plt.rc('ytick', labelsize=SMALL_SIZE)    # y ticks labels.
plt.rc('legend', fontsize=SMALL_SIZE)    # Legend.
plt.rc('figure', titlesize=BIGGER_SIZE)  # Figure title.


def setFontSizes():
    return
