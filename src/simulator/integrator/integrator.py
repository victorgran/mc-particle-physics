import numpy as np

from simulator.constants import f_conv, N_q
from .distributions import Dirac, Uniform
from .particles import ZBoson
from .squared_matrix_element import squared_matrix_element


Z = ZBoson()


class MonteCarloIntegrator:
    """
    Monte Carlo Integrator for integrating the cross-section of the
    electron-positron to quark-antiquark process to leading order.

    The integrator requires specifying distributions for 's', 'cos_theta', 'phi', and the beam 'f(s)',
    defaulting to distributions Dirac(Z_mass2), Uniform(-1, 1), Uniform(0, 2 pi) and Dirac(Z_mass2), respectively.

    'sum_quark_method' can be 'explicit', to sum over all light-quark flavours, or 'random' which generates random
    a random flavour for each Monte Carlo point, takes the average, and multiplies by the number of flavours (here 5).

    'seed' specifies the seed of the random number generator.
    """
    def __init__(self, s_distro=None, cos_theta_distro=None, phi_distro=None, beam_distro=None,
                 sum_quark_method="explicit", seed: int = 42):

        self._seed = seed
        self.rng = np.random.default_rng(seed=self._seed)

        ##############
        # s sampling #
        ##############
        if s_distro is None:
            self._s_distro = Dirac(Z.mass2)
        else:
            self._s_distro = None
            self._setVariableDistro("_s_distro", s_distro)

        self._s_values = np.array([])

        ######################
        # cos theta sampling #
        ######################
        if cos_theta_distro is None:
            self._cos_distro = Uniform(-1., 1., rng=self.rng)
        else:
            self._cos_distro = None
            self._setVariableDistro("_cos_distro", cos_theta_distro)

        self._cos_vals = np.array([])

        ################
        # phi sampling #
        ################
        if phi_distro is None:
            self._phi_distro = Uniform(0., 2 * np.pi, rng=self.rng)
        else:
            self._phi_distro = None
            self._setVariableDistro("_phi_distro", phi_distro)

        self._phi_vals = np.array([])

        #####################
        # beam distribution #
        #####################
        if beam_distro is None:
            self._beam_distro = Dirac(Z.mass2)
        else:
            self._beam_distro = None
            self._setVariableDistro("_beam_distro", beam_distro)

        ###########################
        # sum over quark flavours #
        ###########################
        self._sum_quark_method = sum_quark_method
        self._sum_over_quarks = self._setQuarkMethod()
        self._quark_flavours = np.array([], dtype=int)

        self._d_sigmas = np.array([])

    def __repr__(self):
        class_name = type(self).__name__
        return f"{class_name}(s_distro={self._s_distro}, cos_theta_distro={self._cos_distro}, " \
               f"phi_distro={self._phi_distro}, beam_distro={self._beam_distro}, " \
               f"sum_quark_method={self._sum_quark_method}, seed={self._seed})"

    @property
    def s_samples(self):
        return np.copy(self._s_values)

    @property
    def cos_samples(self):
        return np.copy(self._cos_vals)

    @property
    def phi_samples(self):
        return np.copy(self._phi_vals)

    @property
    def quark_flavours_samples(self):
        if self._sum_quark_method == "explicit":
            exit("No flavours are sampled for 'explicit' quark sum method.")
        return np.copy(self._quark_flavours)

    @property
    def d_sigmas(self):
        return np.copy(self._d_sigmas)

    def _setVariableDistro(self, inst_variable, distro):
        if not hasattr(distro, 'sample'):
            exit(f"Given distribution for {inst_variable[1:]} does not have 'sample' method.")
        elif not hasattr(distro, 'evaluate_distro'):
            exit(f"Given distribution for {inst_variable[1:]} does not have 'evaluate_distro' method.")

        setattr(self, inst_variable, distro)
        setattr(getattr(self, inst_variable), "rng", self.rng)
        return

    def _setQuarkMethod(self):
        if self._sum_quark_method == "explicit":
            def sum_over_quarks(s_values, cos_vals, _):
                Mq = [squared_matrix_element(s_values, cos_vals, q) for q in range(N_q)]
                q_sum = np.sum(Mq, axis=0)
                return q_sum

        elif self._sum_quark_method == "random":
            def sum_over_quarks(s_values, cos_vals, _):
                sample_size = cos_vals.size
                q_flavours = self.rng.integers(0, N_q, sample_size)
                q_sum = N_q * squared_matrix_element(s_values, cos_vals, q_flavours)
                self._quark_flavours = np.append(self._quark_flavours, q_flavours)
                return q_sum

        else:
            exit(f"Method to sum quark flavours: '{self._sum_quark_method}' not implemented.")

        return sum_over_quarks

    def _divideByDistros(self, integrand):
        s_distro_vals = self._s_distro.evaluate_distro(self._s_values)
        cos_dist_vals = self._cos_distro.evaluate_distro(self._cos_vals)
        phi_dist_vals = self._phi_distro.evaluate_distro(self._phi_vals)

        ratios = np.divide(integrand, s_distro_vals)
        ratios = np.divide(ratios, cos_dist_vals)
        ratios = np.divide(ratios, phi_dist_vals)

        return ratios

    def sampleDeltaSigma(self, N):
        if self._d_sigmas.size > N:  # There are enough samples already.
            return

        sample_size = N if self._d_sigmas.size == 0 else N - self._d_sigmas.size

        s_values = self._s_distro.sample(sample_size)
        cos_vals = self._cos_distro.sample(sample_size)
        phi_vals = self._phi_distro.sample(sample_size)
        f_s = self._beam_distro.evaluate_distro(s_values)  # Beam spectrum distribution.

        d_sigmas = self._sum_over_quarks(s_values, cos_vals, phi_vals)
        d_sigmas *= (f_conv * f_s) / (64 * (np.pi ** 2) * s_values)

        # Add to lists to record the values.
        self._s_values = np.append(self._s_values, s_values)
        self._cos_vals = np.append(self._cos_vals, cos_vals)
        self._phi_vals = np.append(self._phi_vals, phi_vals)
        self._d_sigmas = np.append(self._d_sigmas, d_sigmas)

        return

    def integrateCrossSection(self) -> tuple[float, float]:
        N = self._d_sigmas.size
        if N == 0:
            exit("No samples for the differential cross-section have been generated.")

        # Cross-section estimator.
        d_sigmas = self._divideByDistros(self._d_sigmas)
        sigma_avg = np.sum(d_sigmas) / N

        # Monte Carlo error estimate.
        sigma2_avg = np.sum(d_sigmas ** 2) / N
        mc_err = np.sqrt((sigma2_avg - (sigma_avg ** 2)) / N)

        return sigma_avg, mc_err
