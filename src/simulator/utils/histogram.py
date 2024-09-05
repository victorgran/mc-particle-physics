import sys

"""This file reimplements the relevant classes of the YODA histogramming
framework."""


class Bin1D:
    """A single bin of a 1D histogram."""

    def __init__(self, xmin, xmax):
        self.xmin = xmin
        self.xmax = xmax
        self.w = 0.
        self.w2 = 0.
        self.wx = 0.
        self.wx2 = 0.
        self.n = 0.

    def __repr__(self):
        return str(self)

    def __str__(self):
        return "{0:10.6e}\t{1:10.6e}\t{2:10.6e}\t{3:10.6e}\t{4:10.6e}\t{5:10.6e}\t{6}".format(
            self.xmin, self.xmax, self.w, self.w2, self.wx, self.wx2, int(self.n))

    def format(self, tag):
        return "{0}\t{0}\t{1:10.6e}\t{2:10.6e}\t{3:10.6e}\t{4:10.6e}\t{5}".format(
            tag, self.w, self.w2, self.wx, self.wx2, int(self.n))

    def fill(self, x, w):
        """Add weight w."""
        self.w += w
        self.w2 += w * w
        self.wx += w * x
        self.wx2 += w * w * x
        self.n += 1.

    def scale(self, factor):
        """Scale current weight by `factor`."""
        self.w *= factor
        self.w2 *= factor * factor
        self.wx *= factor
        self.wx2 *= factor * factor


class Point2D:
    """A single point in a 2D scatter plot."""

    def __init__(self, xmin, xmax):
        self.x = (xmax + xmin) / 2
        self.xerr = [(xmax - xmin) / 2] * 2
        self.y = 0
        self.yerr = [0] * 2

    def __repr__(self):
        return str(self)

    def __str__(self):
        return "{0:10.6e}\t{1:10.6e}\t{2:10.6e}\t{3:10.6e}\t{4:10.6e}\t{5:10.6e}".format(
            self.x, self.xerr[0], self.xerr[1], self.y, self.yerr[0], self.yerr[1])

    def scale(self, factor):
        """Scale y dimension by `factor`."""
        self.y *= factor
        self.yerr = [err * factor for err in self.yerr]


class Histo1D:
    """A 1D histogram."""

    def __init__(self, nbin, xmin, xmax, name="/MC/untitled"):
        self.name = name
        self.bins = []
        self.uflow = Bin1D(-sys.float_info.max, xmin)
        width = (xmax - xmin) / nbin
        for i in range(nbin):
            self.bins.append(Bin1D(xmin + i * width, xmin + (i + 1) * width))
        self.oflow = Bin1D(xmax, sys.float_info.max)
        self.total = Bin1D(-sys.float_info.max, sys.float_info.max)
        self.scaled_by = 1.

    def __repr__(self):
        return str(self)

    def __str__(self):
        s = "# BEGIN YODA_HISTO1D {0}\n".format(self.name)
        s += "Path={0}\n".format(self.name)
        s += "ScaledBy={0}\n".format(self.scaled_by)
        s += "Title=\nType=Histo1D\n"
        s += "# ID\tID\tsumw\tsumw2\tsumwx\tsumwx2\tnumEntries\n"
        s += self.total.format("Total") + "\n"
        s += self.uflow.format("Underflow") + "\n"
        s += self.oflow.format("Overflow") + "\n"
        s += "# xlow\txhigh\tsumw\tsumw2\tsumwx\tsumwx2\tnumEntries\n"
        for i in range(len(self.bins)):
            s += "{0}\n".format(self.bins[i])
        s += "# END YODA_HISTO1D\n"
        return s

    def fill(self, x, w):
        """Fill a single point of weight w at a given x coordinate."""
        L = 0
        r = len(self.bins) - 1
        c = (L + r) // 2
        a = self.bins[c].xmin
        while r - L > 1:
            if x < a:
                r = c
            else:
                L = c
            c = (L + r) // 2
            a = self.bins[c].xmin
        if x > self.bins[r].xmin:
            self.oflow.fill(x, w)
        elif x < self.bins[L].xmin:
            self.uflow.fill(x, w)
        else:
            self.bins[L].fill(x, w)
        self.total.fill(x, w)

    def scale(self, factor):
        """Scale histogram weights (i.e. bin heights) by `factor`."""
        self.total.scale(factor)
        self.uflow.scale(factor)
        self.oflow.scale(factor)
        for i in range(len(self.bins)):
            self.bins[i].scale(factor)
        self.scaled_by = factor

    def plot(self):
        """Plots the histogram using the matplotlib library."""
        import matplotlib.pyplot as plt
        lefts = [b.xmin for b in self.bins]
        lefts.append(self.bins[-1].xmax)
        heights = [b.w / (lefts[1] - lefts[0]) for b in self.bins]
        heights.append(heights[-1])
        plt.step(lefts, heights, where='post')


class Scatter2D:

    def __init__(self, npoints, xmin, xmax, name="/MC/untitled"):
        self.name = name
        self.points = []
        width = (xmax - xmin) / npoints
        for i in range(npoints):
            self.points.append(Point2D(xmin + i * width, xmin + (i + 1) * width))
        self.scaled_by = 1.

    def __repr__(self):
        return str(self)

    def __str__(self):
        s = "# BEGIN YODA_SCATTER2D {0}\n".format(self.name)
        s += "Path={0}\n".format(self.name)
        s += "Title=\nType=Histo1D\n"
        s += "# xval\txerr-\txerr-\tyval\tyerr-\tyerr+\n"
        for i in range(len(self.points)):
            s += "{0}\n".format(self.points[i])
        s += "# END YODA_SCATTER2D\n"
        return s

    def scale(self, factor):
        """Scales the y coordinates by `factor`."""
        for i in range(len(self.points)):
            self.points[i].scale(factor)
        self.scaled_by = factor
