from numpy import pi, log

# QCD group constants
NC = 3.
TR = 1./2.
CA = NC
CF = (NC*NC-1.)/(2.*NC)


class AlphaS:
    """
    Calculator for the strong coupling \\alpha_s, implementing the running at
    the 0th and the 1st (=default) perturbative order.
    The starting point of the evolution is assumed to be \\alpha_s(\\m_Z**2).

    Examples usage, printing the strong coupling value at a scale of
    100 GeV**2:

      alphas = AlphaS(91.8**2, 0.118)
      print(alphas(100))

    """

    def __init__(self, squared_m_Z, alphas_at_squared_m_Z, perturbative_order=1,
                 squared_m_b=4.75**2, squared_m_c=1.3**2):
        self.order = perturbative_order
        self.mc2 = squared_m_c
        self.mb2 = squared_m_b
        self.mz2 = squared_m_Z
        self.at_sqr_mz = alphas_at_squared_m_Z
        self.at_sqr_mb = self(self.mb2)
        self.at_sqr_mc = self(self.mc2)

    @staticmethod
    def beta0(n_light_flavours):
        return 11. / 6. * CA - 2. / 3. * TR * n_light_flavours

    @staticmethod
    def beta1(n_light_flavours):
        return 17. / 6. * CA * CA - (5. / 3. * CA + CF) * TR * n_light_flavours

    def as0(self, t):
        if t >= self.mb2:
            t_ref = self.mz2
            alpha_ref = self.at_sqr_mz
            b0 = self.beta0(5)/(2.*pi)
        elif t >= self.mc2:
            t_ref = self.mb2
            alpha_ref = self.at_sqr_mb
            b0 = self.beta0(4)/(2.*pi)
        else:
            t_ref = self.mc2
            alpha_ref = self.at_sqr_mc
            b0 = self.beta0(3)/(2.*pi)
        return 1./(1./alpha_ref+b0*log(t/t_ref))

    def as1(self, t):
        if t >= self.mb2:
            t_ref = self.mz2
            alpha_ref = self.at_sqr_mz
            b0 = self.beta0(5)/(2.*pi)
            b1 = self.beta1(5)/pow(2.*pi, 2)
        elif t >= self.mc2:
            t_ref = self.mb2
            alpha_ref = self.at_sqr_mb
            b0 = self.beta0(4)/(2.*pi)
            b1 = self.beta1(4)/pow(2.*pi, 2)
        else:
            t_ref = self.mc2
            alpha_ref = self.at_sqr_mc
            b0 = self.beta0(3)/(2.*pi)
            b1 = self.beta1(3)/pow(2.*pi, 2)
        w = 1. + b0 * alpha_ref * log(t/t_ref)
        return alpha_ref/w*(1.-b1/b0*alpha_ref*log(w)/w)

    def __call__(self, t):
        if self.order == 0:
            return self.as0(t)
        return self.as1(t)
