import math as m
import numpy as np

# from utils.vector import Vec4
# from utils.particle import Particle
from .histogram import Histo1D, Scatter2D


class Analysis:
    """Analyzer of 2->(n jet) scattering events, histogramming differential jet
    rate cross-sections

      d\\sigma / d\\log_10(y_{n,n+1})

    for various n. The integrated n-jet rates are also calculated and stored as
    scatter data.
    """

    def __init__(self):

        self.num_events = 0.

        # Construct list of histograms for the differential jet rates up to
        # y_{n_max,n_max+1} and for the corresponding integrated jet rates.
        n_bins = 100
        self.left_edge = -4.3
        self.right_edge = -0.3
        n_max = 4
        self.y_n = [
            Histo1D(n_bins, self.left_edge, self.right_edge, '/LL_JetRates/log10_y_{0}{1}'.format(i + 2, i + 3))
            for i in range(n_max)]
        self.y_n_integrated = [
            Scatter2D(n_bins, self.left_edge, self.right_edge, '/LL_JetRates/integ_log10_y_{0}'.format(i + 2))
            for i in range(n_max + 1)]

    def analyze(self, event, weight):
        """Adds a single event (= list of Particle instances)
        with corresponding Monte-Carlo weight to the histograms."""

        self.num_events += 1.

        # Fill differential j -> (j+1) splitting scale distributions if there
        # have not been a sufficient number of to cluster, we add the event to
        # the underflow of the histogram.
        y_ij_list = self.cluster(event)
        for j in range(len(self.y_n)):
            log_y = self.left_edge - 1
            if len(y_ij_list) > j:
                log_y = m.log10(y_ij_list[-1 - j])
            self.y_n[j].fill(log_y, weight)

        # Fill integrated j-jet rates.
        previous_logy = 1e20
        for j in range(len(self.y_n_integrated) - 1):
            j_jet_rate = self.y_n_integrated[j]
            log_y = self.left_edge - 1
            if len(y_ij_list) > j:
                log_y = m.log10(y_ij_list[-1 - j])
            for p in j_jet_rate.points:
                if log_y < p.x < previous_logy:
                    p.y += weight
            previous_logy = log_y
        for p in self.y_n_integrated[-1].points:
            if p.x < previous_logy:
                p.y += weight

    def finalize(self, file_name):
        """Scales the histograms properly and writes them out as a YODA file
        with the given file_name."""

        # Divide out the number of events to get the correct cross-section.
        for h in self.y_n:
            h.scale(1. / self.num_events)
        for s in self.y_n_integrated:
            s.scale(1. / self.num_events)

        # Write the histograms to a YODA file.
        file = open(file_name + ".yoda", "w")
        file.write("\n\n".join([str(h) for h in self.y_n]))
        file.write("\n\n")
        file.write("\n\n".join([str(s) for s in self.y_n_integrated]))
        file.close()

    @staticmethod
    def y_ij(p_i, p_j, q2):
        """Calculates the k_T-algorithm distance measure between four momenta
        p_i and p_j, and the reference scale Q^2"""
        pipj = p_i.px * p_j.px + p_i.py * p_j.py + p_i.pz * p_j.pz
        cos_theta = min(max(pipj / m.sqrt(p_i.length_3d_squared() * p_j.length_3d_squared()), -1.0), 1.0)
        return 2.0 * min(p_i.E ** 2, p_j.E ** 2) * (1.0 - cos_theta) / q2

    def cluster(self, event):
        """Applies the k_T clustering algorithm to an event (= list of
        Particle instances). A y_cut is not used, i.e. the clustering continues
        until only two jets are left.

        Returns a list of splitting scales y_ij, ordered from smallest y_ij to
        largest y_ij.
        """
        # TODO: Implement the k_T clustering algorithm described in
        #  the lecture here, and return all found splitting scales y_ij as a
        #  list, ordered from smallest y_ij to largest y_ij.

        # NOTE: The y_ij distance measure is already implemented, see above.
        # As a reference scale Q^2, use the invariant mass of two incoming (or
        # two outgoing) particles, which in our case equals the squared Z mass.

        # Setting the reference scale Q^2 as the invariant mass of the incoming particles.
        p1, p2 = [event[k].mom for k in range(2)]
        q2 = (p1 + p2).invariant_mass_squared()

        # Get the momenta of the final-state particles after showering.
        out_particles = event[2:]
        num_particles = len(out_particles)
        out_momenta = [particle.mom for particle in out_particles]

        splitting_scales = []

        # TODO: Start by creating a naive (stupid) implementation,
        #  check the results, and then come back to do it right.

        for k in range(num_particles - 2):
            # Compute the distance measures.
            y_distances = [self.y_ij(out_momenta[i], out_momenta[j], q2)
                           for i in range(num_particles - k - 1) for j in range(i+1, num_particles - k)]
            ij_values = [(i, j) for i in range(num_particles - 1 - k) for j in range(i+1, num_particles - k)]

            # Find the minimum and append it.
            min_id = np.argmin(y_distances)
            splitting_scales.append(y_distances[min_id])

            # Sum the momenta of the two particles giving the min. distance.
            i, j = ij_values[min_id]
            new_momentum = out_momenta[i] + out_momenta[j]

            # Re-do the list of outgoing momenta.
            out_momenta = [mom for (k, mom) in enumerate(out_momenta) if k != i and k != j]
            out_momenta.append(new_momentum)

        splitting_scales.sort()  # Sort the list from smallest to largest.

        # NOTE: There is no fixed y_cut, since we want to know at which point an
        # n-jet event starts looking like a (n+1)-jet event for a range of
        # different n (compare the kT jet fraction plot in the lecture, which
        # is a sort of integrated version of what we'd like to plot). Instead,
        # keep clustering until only two jets are left.

        # N particles have N(N-1)/2 distances. Each splitting leads from N to N-1.
        # In order to reach distances = 1, we require N - 2 splittings.

        # TODO:
        #  1. Start by computing all distance measures. How to arrange them?
        #  2. Find the minimum, append it to the list.
        #  3. Remove all the distance measures related to i and j.
        #  4. Sum momenta of particles i and j into particle k.
        #  5. Compute the distances of the remaining particles to particle k.

        return splitting_scales
