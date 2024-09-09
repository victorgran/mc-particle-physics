import matplotlib.pyplot as plt
import numpy as np

lw = 1.0
colors = ['red', 'blue', 'green']

single_plot = False


def data_objects(filenames, yoda_type="HISTO1D"):
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


def plot_jet_histograms(filenames):
    """Plots jet rate histograms from YODA files. Pass a list of filenames to
    read histograms from. The first one is considered to be the reference
    histogram, against which all other files will be compared against. For this
    purpose, not only a nominal plot is drawn, but also a deviation plot, which
    indicates how many sigma are given bin deviates from the reference.
    """

    print("Plotting (might take a few moments) ...")

    with np.errstate(divide='ignore', invalid='ignore'):
        # There might be unfilled bins, therefore ignore numpy divide-by-0 errors,
        # cf. https://stackoverflow.com/questions/29950557/ignore-divide-by-0-warning-in-numpy

        histonames, histos = data_objects(filenames, yoda_type="HISTO1D")

        if single_plot:
            fig, ax = plt.subplots()
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
                yerr = np.sqrt(histo[:, 3]) / sums2  # / widths

                njets = int(histoname[-2])
                ax.step(x, np.append(y, y[-1]), where="post", linewidth=1.0, label=r"{}-rate".format(njets))
                ax.bar(centers, y, yerr=yerr, alpha=0.0, linewidth=lw)
            ax.legend()

        else:
            fig, axes = plt.subplots(nrows=len(histonames), ncols=2, figsize=[12, 15])
            axes = np.array(axes)
            main_axes = axes[:, 0]
            ratio_axes = axes[:, 1]
            for histoname, (main_ax, ratio_ax) in zip(histonames, zip(main_axes, ratio_axes)):
                try:
                    raw_central_histo = histos[filenames[0]][histoname]["bins"]
                    central_histo = np.array(raw_central_histo)
                except KeyError:
                    raw_central_histo = None
                    central_histo = None

                for (file, file_histos), color in zip(histos.items(), colors):
                    try:
                        raw_histo = file_histos[histoname]["bins"]
                        histo = np.array(raw_histo)
                    except KeyError:
                        continue

                    x = np.append(histo[:, 0], histo[:, 1][-1])
                    widths = histo[:, 1] - histo[:, 0]
                    centers = 0.5 * (x[:-1] + x[1:])
                    y = histo[:, 2] / widths
                    yerr = np.sqrt(histo[:, 3]) / widths

                    main_ax.step(x, np.append(y, y[-1]), where="post", linewidth=1.0, label=file, color=color)
                    main_ax.bar(centers, y, yerr=yerr, alpha=0.0, linewidth=lw, ecolor=color)

                    if raw_histo is raw_central_histo:
                        ratio_ax.step(x, [0] * len(x), color=color)
                        ratio_ax.fill_between(x, [-1] * len(x), [1] * len(x), color='yellow', alpha=0.5)
                    elif central_histo is not None:
                        y = (y - central_histo[:, 2] / widths) / (np.sqrt(central_histo[:, 3]) / widths)
                        yerr /= np.sqrt(central_histo[:, 3]) / widths
                        ratio_ax.step(x, np.append(y, y[-1]), where="post", linewidth=lw, color=color)
                        error_kw = {"elinewidth": lw}
                        ratio_ax.bar(centers, y, yerr=yerr, alpha=0.0, ecolor=color, error_kw=error_kw)

                main_ax.set_yscale('log')
                njet = int(histoname[-2])
                main_ax.set_title(
                    r"Differential ${}\to {}$ jet resolution".format(njet, njet + 1))
                jets = "{}{}".format(njet, njet + 1)
                main_ax.set_xlabel(r"$\log_{10}(y_{" + jets + r"})$")
                main_ax.set_ylabel(r"$\mathrm{d}\sigma/\mathrm{d}\log_{10}(y_{" + jets + r"})$ [pb]")
                ratio_ax.set_ylim(-5, 5)
                ratio_ax.set_ylabel("standard deviation")

                xlim = None
                if "y_23" in histoname:
                    xlim = (-4, -0.4)
                    main_ax.set_ylim(4e2, 4e4)
                    main_ax.legend()
                elif "y_34" in histoname:
                    xlim = (-4, -0.8)
                    main_ax.set_ylim(8e-1, 8e4)
                elif "y_45" in histoname:
                    xlim = (-4, -1.2)
                    main_ax.set_ylim(8e-1, 8e4)
                elif "y_56" in histoname:
                    xlim = (-4, -1.6)
                    main_ax.set_ylim(8e-2, 8e4)
                main_ax.set_xlim(xlim)
                ratio_ax.set_xlim(xlim)

            # fig.tight_layout()


def plot_scatters(filenames):
    histonames, histos = data_objects(filenames, yoda_type="SCATTER2D")

    fig, ax = plt.subplots()

    sums = np.zeros(len(histos[filenames[0]][histonames[0]]["bins"]))
    for i in range(len(sums)):
        sums[i] += sum([histos[filenames[0]][n]["bins"][i][2] for n in histonames])

    for histoname in histonames:
        histo = np.array(histos[filenames[0]][histoname]["bins"])

        x = np.append(histo[:, 0], histo[:, 1][-1])
        # widths = histo[:, 1] - histo[:, 0]
        y = histo[:, 2] / sums

        njets = int(histoname[-1])
        ax.step(x, np.append(y, y[-1]), where="post", linewidth=1.0, label=r"{}-rate".format(njets))
    ax.legend()
    ax.set_ylabel("relative n-jet event rates")
    ax.set_xlabel(r"$\log_{10}(y_\mathrm{cut})$")
