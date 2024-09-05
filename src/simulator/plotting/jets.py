import matplotlib.pyplot as plt
import numpy as np


def data_objects(filenames: list[str], yoda_type="HISTO1D"):
    names = set()
    objects = {}
    for filename in filenames:
        with open(filename) as file:
            objects[filename] = {}
            current = None
            for line in file:
                words = line.split()
                if len(words) == 0:
                    continue
                if len(words) > 3 and words[2] == "YODA_" + yoda_type:
                    objects[filename][words[3]] = {"scale_factor": 1.0, "bins": []}
                    current = objects[filename][words[3]]
                    names.add(words[3])
                    continue
                if len(words) > 1 and "END" == words[1]:
                    current = None
                    continue
                if current is None:
                    continue
                if "ScaledBy" in words[0]:
                    current["scale_factor"] = float(line.split("=")[-1])
                    continue
                try:
                    if yoda_type == "HISTO1D":
                        xlow = float(words[0])
                        xhigh = float(words[1])
                        sumw = float(words[2])
                        sumw2 = float(words[3])
                    else:
                        xlow = float(words[0]) - float(words[1])
                        xhigh = float(words[0]) + float(words[2])
                        sumw = float(words[3])
                        sumw2 = 0.0
                    current["bins"].append([xlow, xhigh, sumw, sumw2])
                except ValueError:
                    pass

    return sorted(names), objects


def _singlePlot(filenames: list[str], **kwargs):

    histonames, histos = data_objects(filenames, yoda_type="HISTO1D")

    fig, ax = plt.subplots(figsize=(7, 5))
    sums = np.zeros(len(histos[filenames[0]][histonames[0]]["bins"]))
    sums2 = sums.copy()

    for i in range(len(sums)):
        sums[i] += sum([histos[filenames[0]][n]["bins"][i][2] for n in histonames])
        sums2[i] += sum([histos[filenames[0]][n]["bins"][i][3] for n in histonames])

    for histoname in histonames:
        histo = np.array(histos[filenames[0]][histoname]["bins"])

        x = np.append(histo[:, 0], histo[:, 1][-1])
        # widths = histo[:, 1] - histo[:, 0]
        centers = 0.5 * (x[:-1] + x[1:])
        y = histo[:, 2] / sums  # / widths
        y_err = np.sqrt(histo[:, 3]) / sums2  # / widths

        num_jets = int(histoname[-2])
        ax.step(x, np.append(y, y[-1]), where="post", linewidth=1.0, label=r"{}-rate".format(num_jets))
        ax.bar(centers, y, yerr=y_err, **kwargs)

    ax.legend()
    fig.tight_layout()

    plot = (fig, ax)

    return plot


def getBinsData(histo_data: dict) -> dict:

    bins_array = np.array(histo_data["bins"])
    xmin, xmax, areas, errs_area = bins_array.T
    widths = xmax - xmin
    centers = 0.5 * (xmax + xmin)
    heights = areas / widths
    errors = np.sqrt(errs_area) / widths
    edges = np.append(xmin, xmax[-1])

    bins_data = {"xmin": xmin, "xmax": xmax, "areas": areas, "errs_area": errs_area,
                 "widths": widths, "centers": centers, "heights": heights, "errors": errors, "edges": edges}

    return bins_data


def _fillAxes(histo_data, axes, is_reference: bool = False, ref_bins=None, color=None, label=None):

    histo_axis = axes[0]
    ratio_axis = axes[1]

    bins_data = getBinsData(histo_data)
    errors = bins_data["errors"]
    heights = bins_data["heights"]
    centers = bins_data["centers"]
    xmin = bins_data["xmin"]

    histo_axis.stairs(heights, bins_data["edges"], color=color, label=label)
    histo_axis.errorbar(centers, heights, yerr=errors, fmt='none', color=color)

    if is_reference:
        ratio_axis.step(xmin, [0] * len(xmin), label=label, color=color)
        ratio_axis.fill_between(xmin, [-1] * len(xmin), [1] * len(xmin), color='yellow', alpha=0.5)
        return bins_data
    else:
        if ref_bins is None:
            exit("Need reference data to plot error ratio.")

        heights = (heights - ref_bins["heights"]) / (ref_bins["errors"])
        errors /= ref_bins["errors"]

        ratio_axis.step(xmin, heights, where="post")
        ratio_axis.errorbar(centers, heights, yerr=errors, fmt='none')

    return


def plotJetHistogram(histo_files_data, histoname: str, colors: list[str] = None, labels: list[str] = None):

    # Start a figure composed of two subplots. Left: histogram / Right: Deviation ratio.
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(14, 5))
    histo_axis, ratio_axis = axes

    if colors is None:
        colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
    if labels is None:
        labels = [f"File {k + 1}" for k in range(len(histo_files_data))]

    # Take the first file as the reference and start filling the axis.
    ref_data = histo_files_data[0]
    ref_bins = _fillAxes(ref_data, axes, is_reference=True, color=colors[0], label=labels[0])

    if len(histo_files_data) >= 1:
        for (histo_data, color, label) in zip(histo_files_data[1:], colors[1:], labels[1:]):
            _fillAxes(histo_data, axes, ref_bins=ref_bins, color=color, label=label)

    # Formatting of the axes.
    histo_axis.set_yscale('log')
    ratio_axis.set_ylim(-5, 5)
    ratio_axis.set_ylabel("standard deviation")
    njet = int(histoname[-2])
    histo_axis.set_title(
        r"Differential ${}\to {}$ jet resolution".format(njet, njet + 1))
    jets = "{}{}".format(njet, njet + 1)
    histo_axis.set_xlabel(r"$\log_{10}(y_{" + jets + r"})$")
    histo_axis.set_ylabel(r"$\mathrm{d}\sigma/\mathrm{d}\log_{10}(y_{" + jets + r"})$ [pb]")

    x_lim = None
    if "y_23" in histoname:
        x_lim = (-4, -0.4)
        histo_axis.set_ylim(4e2, 4e4)
    elif "y_34" in histoname:
        x_lim = (-4, -0.8)
        histo_axis.set_ylim(8e-1, 8e4)
    elif "y_45" in histoname:
        x_lim = (-4, -1.2)
        histo_axis.set_ylim(8e-1, 8e4)
    elif "y_56" in histoname:
        x_lim = (-4, -1.6)
        histo_axis.set_ylim(8e-2, 8e4)

    histo_axis.set_xlim(x_lim)
    ratio_axis.set_xlim(x_lim)
    histo_axis.legend()
    fig.tight_layout()

    plot = (fig, axes)

    return plot


def plotJetHistograms(filenames: list[str], single_plot: bool = False,
                      colors: list[str] = None, labels: list[str] = None):
    """
    Plots jet rate histograms from YODA files. Pass a list of filenames to
    read histograms from. The first one is considered to be the reference
    histogram, against which all other files will be compared against. For this
    purpose, not only a nominal plot is drawn, but also a deviation plot, which
    indicates how many sigma are given bin deviates from the reference.
    """

    print("Plotting (might take a few moments) ...")
    with np.errstate(divide='ignore', invalid='ignore'):
        # There might be unfilled bins, therefore ignore numpy divide-by-0 errors,
        # cf. https://stackoverflow.com/questions/29950557/ignore-divide-by-0-warning-in-numpy

        if single_plot:
            plot_attrs = {"alpha": 0.0, "linewidth": 1.0}
            plot = _singlePlot(filenames, **plot_attrs)
            return [plot]

        # Multiple histograms.
        plots = []
        histonames, histos = data_objects(filenames, yoda_type="HISTO1D")

        for histoname in reversed(histonames):
            histo_data = []
            for filename, file_histos in histos.items():
                try:
                    data = file_histos[histoname]
                except KeyError:
                    print(f"File {filename} doesn't contain data for histogram {histoname}. Skipping...")
                    continue

                histo_data.append(data)

            if histo_data is None:
                print(f"Histogram {histoname} is registered but no data was found.")
                continue

            plot = plotJetHistogram(histo_data, histoname, colors=colors, labels=labels)
            plots.append(plot)

    # Figures are generated in reverse sequence to have proper display, but the
    # list needs to be reversed so that, e.g., saving the first figure corresponds
    # to the 2 -> 3 jets.
    plots.reverse()

    return plots
