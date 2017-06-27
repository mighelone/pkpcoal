"""
Chemical Percolation Devolatilization model.

New implementation in python of the Chemical Percolation
Devolatilization(CPD) model for coal devolatilization.

The CPD class is derived by the CPD class in the `PKP` module, which
calls the Fortran CPD code externally.

(c) Michele Vascellari 2016 Michele.Vascellari@vtc.tu-freiberg.de

Example
-------

Initialize the :class:`pkp.cpd.CPD` class using the proximate and
ultimate analysis::

    >>> import pkp.cpd
    >>> ua = {'C': 69,
              'H': 5,
              'O': 24.7,
              'N': 0.8,
              'S': 0.5}

    >>> pa = {'FC': 45.1,
              'VM': 50.6,
              'Ash': 4.3,
              'Moist': 19.0}
    >>> m = pkp.cpd.CPD(proximate_analysis=pa,
                        ultimate_analysis=ua,
                        pressure=101325,
                        name='Pittsburg')

Set the operating conditions::

    >>> op_cond = [[0, 500],
                   [0.001, 1400],
                   [0.01, 1400]]
    >>> m.operating_conditions = op_cond

Set the output directory where store the results::

    >>> m.path = './Results'

Set the parameters for the CPD run::

    >>> m.set_parameters(dt=1e-5, increment=2, dt_max=1e-5,
                         basename='test')

and finally run the model::

    >>> results = m.run()

The results are stored in a ``pandas`` dataframe.
"""

from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals


import pkp
import pkp.coal
import pkp.empirical_model
import pkp.triangle
import numpy as np
from autologging import logged
from scipy.integrate import ode
# from .binomial import bpmfln
from scipy.optimize import brentq, newton
import pandas as pd
import os
import warnings

from ._exceptions import ImportError

# Import Numba
try:
    from ._nb_functions import sum_x_n_calc, x_n_calc, fp, pstar_f
    from ._nb_functions import binomial, invernorm
    _use_numba = True
except ImportError:
    from ._np_functions import sum_x_n_calc, x_n_calc, fp, pstar_f
    from ._np_functions import binomial, invernorm
    _use_numba = False

# CPD constants
from ._cpd_correlations import CPD_CORRELATION, x_gas, y_gas

# define the binomial function
# binomial = bpmfln
# binomial = binom.pmf

Rgas = 1.987  # cal/mol-K


@logged
class CPD(pkp.coal.Coal, pkp.empirical_model.Model):
    """
    Chemical Percolation Devolatilization model (CPD).

    Run the Chemical Percolation Model (CPD) for evaluating the
    devolatilization behaviour of coal.

    It is based on the CPD_NLG version of the code:
    http://www.et.byu.edu/~tom/cpd/cpdcodes.html

    It is possible to use the Genetti correlation for estimating NMR
    parameters, or they can directly defined if they are known.
    """

    nmr_parameters = ['mdel', 'mw', 'p0', 'sig', 'c0']
    num_parameters = ['dt', 'increment', 'dt_max']

    ab = 2.602e15
    eb = 55400
    ebsig = 1800
    ac = 0.9
    ec = 0
    ag = 3.0e15
    eg = 69000
    egsig = 8100
    Acr = 3.0e15
    Ecr = 65000
    arad = 18.4
    erad = 6000
    fstable = 0.03
    an = 5.5e7
    en = 90000
    ensig = 0
    nmax = 20

    skip = 1

    def __init__(self, ultimate_analysis=None, proximate_analysis=None,
                 pressure=101325, name='CPD coal'):
        """
        Initialize CPD model.

        Parameters
        ----------
        proximate_analysis: dict
            Proximate analysis dict i.e:
            `{'FC': 45.1, 'VM': 50.6, 'Ash': 4.3, 'Moist': 19.0}`
        ultimate_analysis: dict
            Ultimate analysis dictionary i.e:
            `{'C': 80, 'H': 8, 'O': 12, 'N': 0, 'S': 0}`
        pressure: float
            Pressure of pyrolysis process
        name: str, unicode
            Reference name of the modelled coal

        """
        super(CPD, self).__init__(proximate_analysis=proximate_analysis,
                                  ultimate_analysis=ultimate_analysis,
                                  pressure=pressure,
                                  name=name)

        # check if they are in %
        self.fcar = self.ultimate_analysis['C']
        self.fhyd = self.ultimate_analysis['H']
        self.fnit = self.ultimate_analysis['N']
        self.foxy = self.ultimate_analysis['O']
        self.vm_daf = self.proximate_analysis_daf['VM']

        # set parameters -> this can be changed using
        # self.set_parameters
        self.solver = None
        self._set_NMR_parameters()
        self._set_numerical_parameters()

        # initialize fractions
        self.f = [[1, 0, 0, 0, 0]]
        self.n_frag = 20
        self.meta_n = np.zeros(self.n_frag)
        self.f_frag_n = np.zeros(self.n_frag)
        self.t_old = 0

    def set_parameters(self, **kwargs):
        """Set parameters for CPD.

        If a parameter is not defined it is not changed.
        If None the default value is setted.

        Parameters
        ----------
        nmr_parameters: dict
            NMR parameters dictionary. Keys are:
            ['mdel', 'mw', 'p0', 'sig', 'c0']
            If None values are calculated using Genetti correlation
        dt: float
            Time step
        increment: int
            Number of time step saved
        dt_max: float
            Max. time step
        solver: str
            CPD solver path
        basename: str
            Basename for CPD output files

        """
        # TODO set a warning for coal outside the parameters
        if 'nmr_parameters' in kwargs:
            self._set_NMR_parameters(
                nmr_parameters=kwargs['nmr_parameters'])

        num_parameters = ('dt', 'increment', 'dt_max')
        if any(p in kwargs for p in num_parameters):
            self._set_numerical_parameters(**kwargs)

        if 'solver' in kwargs:
            self.solver = kwargs['solver']

        if 'basename' in kwargs:
            self.basename = kwargs['basename']

    def get_parameters(self):
        """Get the parameters of CPD model."""
        nmr = {p: getattr(self, p)
               for p in self.nmr_parameters}
        par = {p: getattr(self, p)
               for p in self.num_parameters +
               ['basename', 'solver']}
        par['nmr_parameters'] = nmr
        return par

    def _set_numerical_parameters(self, dt=None, increment=None,
                                  dt_max=None, **kwargs):
        """Set numerical parameters."""
        # TODO get rid of it. Should be useless now
        if dt is None:
            dt = 1e-5
        if increment is None:
            increment = 1
        if dt_max is None:
            dt_max = 1e-5
        self.dt = dt
        self.increment = increment
        self.dt_max = dt_max

    def _set_NMR_parameters_from_correlation(self, nmr_parameters=None):
        """
        Set NMR parameters from correlation.

        Use the Genetti correlation for evaluating NMR input parameters for
        CPD.

        Parameters
        ----------
        parameters: dict, default=None
            Manually set parameters using a dictionary
            {'mdel': 0, 'mw': 0, 'p0': 0, 'sig': 0}

        """
        if nmr_parameters:
            [setattr(self, key, nmr_parameters[key])
             for key in self.nmr_parameters]
        else:
            c = CPD_CORRELATION.copy()
            self.c0 = (min(0.36,
                           max(0.118 * self.fcar * 100 - 10.1, 0.0)) +
                       min(0.15,
                           max(0.014 * self.foxy * 100 - 0.175, 0.0)))

            Y = (c[1] + c[2] * (self.fcar * 100.0) +
                 c[3] * (self.fcar * 100)**2 +
                 c[4] * (self.fhyd * 100) +
                 c[5] * (self.fhyd * 100)**2 +
                 c[6] * (self.foxy * 100) +
                 c[7] * (self.foxy * 100)**2 +
                 c[8] * (self.vm_daf * 100) +
                 c[9] * (self.vm_daf * 100)**2)
            [setattr(self, key, Y[i])
             for i, key in enumerate(self.nmr_parameters[:4])]

    def run(self, time=None, light_gas=True, n_frag=20):
        """
        Run CPD.

        Parameters
        ----------
        time: float
            End of calculation, if it is not specified use the last
            time in operating_conditions.
        light_gas: bool
            Calculate light gas using Genetti method
        n_frag: int (default=20)
            Number of fragments considered for tar formation

        Returns
        -------
        results: pandas.Dataframe
            Dataframe containg the results of CPD as a function of the
            residence time.

        """
        t, y, f = self._bridge_evolution(n_frag=n_frag, time_end=time)
        if self.increment > 1:
            t = t[::self.increment]
            y = y[::self.increment]
            f = f[::self.increment]
        data = np.concatenate([t[:, np.newaxis], y, f], axis=1)
        df = pd.DataFrame(data, index=t,
                          columns=['t', 'l', 'delta', 'c', 'char',
                                   'light_gas', 'tar', 'meta', 'cross'])
        df['T'] = [self.T(ti) for ti in t]
        df['p'] = self.intact_bridges(y.T)
        df['f'] = 1 - df['p']
        df['g1'], df['g2'] = self.gas(y.T)
        df['volatiles'] = df['tar'] + df['light_gas']

        if light_gas:
            X_gas = df['delta'] * 0.5 + df['l']
            X_gas = 1 - X_gas / X_gas.iloc[0]
            # self.find_triangle(plot='show')
            self.find_triangle()
            if self.triangle:
                # light gases are evaluated only if the coal
                # is inside one of the defined points
                X_gases = self.calc_lightgases(X_gas)
                for x, sp in zip(X_gases.T, ('CO', 'CO2', 'H2O', 'CH4')):
                    df[sp] = x * df['light_gas']
                df['others'] = 1 - X_gases.sum(axis=1)
                # f_gases = pd.DataFrame(
                #    (df['light_gas'].values * X_gases.T).T,
                #    columns=gas_species, index=df.index)
                # df = pd.concat([df, f_gases], axis=1)
                # df['others'] = 1 - (df['CO'] + df['CO2'] +
                #                    df['H2O'] + df['CH4'])

        df.to_csv(self._out_csv)
        return df

    def _set_NMR_parameters(self, nmr_parameters=None):
        """
        Calculate parameters using Genetti correlation.

        The NMR parameters are defined using the empirical correlation
        derived by Genetti. Input are the daf ultimate and proximate analyses.

        Parameters
        ----------
        parameters: dict, default=None
            Manually set parameters using a dictionary
            {'mdel': 0, 'mw': 0, 'p0': 0, 'sig': 0}

        """
        self._set_NMR_parameters_from_correlation(nmr_parameters)
        # test
        # self.mdel /= (1 - self.c0)
        # check this correction -> allow to obtain the same results as
        # the original code
        self.mdel /= (1 - self.c0)
        self.mdel -= 7

        self.sigma = self.sig - 1
        # average mass of the fused ring site
        self.ma = self.mw - self.sig * self.mdel
        # mass of bridges
        self.mb = 2 * self.mdel
        self.rba = self.mb / self.ma
        self.gasmw = self.rba * self.ma * 0.5
        self.solver = None

    def rate(self, t, y):
        """
        Integrand function for the ODE solver.

        Parameters
        ----------
        t: float
            Time
        y: array
            Number of bridges calculated respect to overall number of
            bridges:
            (labile bridges, side chains, char bridges)

        Return:
        dydt: array

        """
        l, delta, c, T = y
        # l, delta, c = np.maximum(y, 0)

        # Calculate the temperature.
        # It is valid only for prescribed particle temperatures
        kb, rho, kg = self._rates(y)
        f = 1 / (1 + rho)
        tol = 1e-6
        if l > tol:
            dldt = -kb * l
            dcdt = kb * f * l
        else:
            dldt, dcdt = 0.0, 0.0
        # if l < tol and delta < tol:
        #    ddeldt = 0.0
        # else:
        ddeldt = 2 * rho * kb * f * l - kg * delta
        # self.__log.warning('t: %s y: %s - Rates: %s - %s %s',
        #                   t, y, [dldt, ddeldt, dcdt], 2 * rho * kb * f, kg)
        return [dldt, ddeldt, dcdt]

    @property
    def y0(self):
        """Init solution."""
        return [self.p0 - self.c0,
                2 * (1 - self.p0),
                self.c0]

    def _bridge_evolution(self, n_frag=20, time_end=None):
        """
        Calculate the evolution of the bridges.

        Calculate the evolution of the bridges, solving the ODE.
        For each time step, the cross-linking, percolation and flash
        distillation functions are called in order to calculate the
        mass fractions of gas, tar and solid.

        Parameters
        ----------
        n_frag: int (20)
            Number of fragments considered
        time_end: float, None
            End of the calculation. If not specified use the maximum
            time in the operating_conditions.

        Return
        ------
        t, y, f
            Time, bridges and mass fraction arrays

        """
        # TODO needs more test!
        # variables are [l, d, c]
        backend = 'dopri5'
        # backend = 'vode'
        t0 = self.operating_conditions[0, 0]
        solver = ode(self._dydt).set_integrator(backend, nsteps=1,
                                                first_step=self.dt,
                                                # min_step=1e-6,
                                                atol=1e-6,
                                                rtol=1e-4,
                                                max_step=self.dt_max)

        # verbosity=1)
        solver._integrator.iwork[2] = -1
        y0 = [self.p0 - self.c0,
              2 * (1 - self.p0),
              self.c0]
        solver.set_initial_value(y0, t0)
        solver._integrator.iwork[2] = -1
        if not time_end:
            time_end = self.operating_conditions[-1, 0]
        else:
            assert time_end < self.operating_conditions[-1, 0], (
                'Set a time_end < operating_conditions[-1, 0]')

        t = [t0]
        y = [y0]
        f_solid, f_gas, f_tar, f_meta, f_cross = 1, 0, 0, 0, 0
        f = [[f_solid,
              f_gas,
              f_tar,
              f_meta,
              f_cross]]
        meta_n = np.zeros(n_frag)    # init metaplast to zeros
        f_frag_n = np.zeros(n_frag)  # init fragments to zeros
        warnings.filterwarnings("ignore", category=UserWarning)

        while solver.t < time_end:
            solver.integrate(time_end, step=True)
            self.__log.debug(
                '\n\nStart new time step\ntime=%s y=%s\n', solver.t,
                solver.y)
            self.__log.info('t=%s - y=%s', solver.t, solver.y)
            dt = solver.t - t[-1]
            T = self.T(solver.t)
            t.append(solver.t)
            y.append(solver.y)

            if f_meta > 1e-4:
                # cross-linking
                rate_cross = self._crosslinking(f_meta, T, dt)
                fract = 1 - rate_cross / f_meta
                f_meta -= rate_cross
                f_cross += rate_cross
            else:
                rate_cross = 0
                fract = 1
            self.__log.debug(
                'Crosslinking rate: %s / %s', rate_cross, fract)
            percolation = self._percolation(solver.y, f_tar,
                                            n_frag=n_frag, in_tar=True)
            # gas formed in the last step
            df_gas = max(percolation['f_gas'] - f_gas, 0)
            self.__log.debug('gas=%s, df_gas=%s', f_gas, df_gas)
            mw_n = percolation['mw_frag_n']
            # fragments formed in the last step
            df_n = percolation['f_frag_n'] - f_frag_n
            self.__log.debug('df_n=%s', df_n)

            tar_n, meta_n = self._flash_distillation(
                df_gas=df_gas, df_n=df_n, meta_n=meta_n,
                mw_n=mw_n, fracr=fract, T=T)

            f_frag_n = percolation['f_frag_n']
            f_gas = percolation['f_gas']
            f_tar += tar_n.sum()
            f_meta = meta_n.sum()
            f_solid = 1 - f_tar - f_gas
            f.append([f_solid, f_gas, f_tar, f_meta, f_cross])
            self.__log.debug('F=%s', f[-1])

        warnings.resetwarnings()
        t = np.array(t)
        y = np.array(y)
        f = np.array(f)

        return t, y, f

    def _rates(self, y):
        """Calculate _rates for the given temperature."""
        l, delta, c, T = y
        g1, g2 = self.gas(y)
        g = g1 + g2
        RT = T * Rgas
        # calculate bridge decomposition reaction rate
        eb = self.eb + self.ebsig * invernorm(
            1 - l / (self.p0 - self.c0))
        # self.__log.debug('Eb %s Eb0 %s', eb, self.eb)
        kb = self.ab * np.exp(-eb / RT)
        kc = self.ac * np.exp(-self.ec / RT)
        eg = self.eg + invernorm(0.5 * g / (1 - self.c0)) * self.egsig
        kg = self.ag * np.exp(-eg / RT)
        return kb, kc, kg
        # return rates(T, y, self.p0, self.c0, self.ab, self.eb, self.ebsig,
        #              self.ac, self.ec, self.ag, self.eg, self.egsig)

    def _percolation(self, y, f_tar=0, n_frag=20, in_tar=True):
        """
        Percolation statistic calculation.

        Parameters
        ----------
        y: list
            List of bridges parameters (l, delta, c)
        f_tar: float
            Fraction of tar already released
        n_frag: int, default=20
            Number of fragments calculated
        in_tar: bool

        Returns
        -------
        percolation, dict:
            Percolation dictionary {'f_gas', 'f_solid', 'f_frag',
            'f_frag_n', 'm_frag_n', 'pstar'}

        """
        self.__log.debug('\n\nStart Percolation\n')

        l, delta, c = y[:-1]
        p = l + c
        g1, g2 = self.gas(y)
        g = g1 + g2

        if in_tar:
            # Eq. 36, 37
            # Phi->a
            # Omega->b
            delta_fac = delta / (1 - p) if p < 0.9999 else 1
            self.__log.debug('delta/(1-p)=%s', delta_fac)
            a = 1 + self.rba * (l / p + (self.sigma - 1)
                                * 0.25 * delta_fac)
            b = 0.5 * delta_fac - l / p
            # p_threasold is the maxiumum value of pstar_eq
            # pstar is search from 0 to p_threasold
            p_threasold = 1. / self.sigma
            self.__log.debug('p thresold %s', p_threasold)
            if p > 0.999:
                pstar = 1
            elif p > p_threasold:
                # def fp(x): return x * (1 - x)**(self.sigma - 1)
                fpp = fp(p, self.sigma)

                # def pstar_f(x): return fp(x) - fpp
                # pstar = brentq(pstar_f, 0, p_threasold)
                pstar = newton(pstar_f, p_threasold * 0.5,
                               args=(self.sigma, fpp))
                self.__log.debug('Calc pstar with newton %s', pstar)
            else:
                pstar = p
            self.__log.debug('p %s, pstar %s', p, pstar)
            sfac = self.sig / (self.sigma - 1)
            # Eq. (5) fraction of bridges in finite fragments
            Fp = (pstar / p) ** sfac
            self.__log.debug(
                'Fraction of bridges in finite fragments=%s', Fp)
            Kp = Fp * (1 - self.sig * 0.5 * pstar)
            # Eq. (39) mass fraction of finite fragments
            f_frag = 2 * (a * Fp + b * Kp) / \
                (2 + self.rba * (1 - self.c0) * self.sig)
            self.__log.debug(
                'Mass fraction of finite fragments=%s', f_frag)

        # mass of gas released at time t Eq (31)
        # gas is produced only considering the remaining fragments in
        # the metaplast
        # gas is corrected with fraction of tar already released in the
        # gas phase. f_tar is from the previous time step
        mgas = self.mb * g * self.sig * 0.25 * (1 - f_tar)
        mtot = self.ma + self.mb * self.sig * 0.5 * (1 - self.c0)
        f_gas = mgas / mtot
        self.__log.debug(
            'fraction of gas (corrected with tar released) %s (%s)',
            f_gas, f_tar)
        f_solid = 1 - f_tar - f_gas
        self.__log.debug(
            ('fraction of remaining solid (includes finite'
             ' and inf. fragments) %s'), f_solid)

        n = np.arange(1, n_frag + 1)  # number of clusters in a fragment
        # broken bridges per cluster of size n
        tau = n * (self.sigma - 1) + 2
        s = n - 1  # intact bridges per cluster of size n
        n_bridges = tau + s
        # Eq. 32 mass of a finite fragment of size n
        mw_frag_n = (n * self.ma + (n - 1) * self.mb * l / p +
                     tau * self.mb * delta_fac * 0.25)

        # Eqs (1-4)
        Qn = self.sig / n_bridges / n * binomial(s, n_bridges, p)
        # Eq. (33) total mass of fragments of size
        m_frag_n = mw_frag_n * Qn
        self.__log.debug(
            'mass weight of finite fragments %s', mw_frag_n)

        # Eq 35 total mass of finite fragments
        m_frag = a * self.ma * Fp + b * self.mb * Kp

        # Eq 38 fraction of finite fragments
        # TODO this has to be corrected similarly to the f_gas equation
        f_frag = m_frag / mtot
        f_frag_n = m_frag_n / mtot
        self.__log.debug(
            'mass fraction of finite fragments %s', f_frag_n)

        self.__log.debug(
            'Total fraction of fragments sum %s / Eq.35 %s',
            f_frag_n.sum(), f_frag)

        return {'f_gas': f_gas,
                'f_solid': f_solid,
                'f_frag': f_frag,
                'f_frag_n': f_frag_n,
                'mw_frag_n': mw_frag_n,
                'pstar': pstar}

    def _tar_distribution(self):
        pass

    def intact_bridges(self, y):
        """
        Return the fraction of intact bridges.

        .. math::
            p = l + c

        Parameters
        ----------
        y: array
            Array containing (l, delta, c)

        Return
        ------
        p: float
            Fraction of intact bridges

        """
        l = y[0]
        c = y[2]
        return l + c

    def broken_bridges(self, y):
        """
        Return the fraction of broken bridges.

        .. math::
            f = 1 - (l + c)

        Parameters
        ----------
        y: array
            Array containing (l, delta, c)

        Return
        ------
        f: float
            Fraction of broken bridges

        """
        return 1 - self.intact_bridges(y)

    def gas(self, y):
        r"""
        Return the fraction of gas produced.

        .. math::
            g_1 = 2 \cdot f - \delta\\
            g_2 = 2 \cdot (c-c_0)

        Parameters
        ----------
        y: array
            Array containing (l, delta, c)

        Return
        ------
        g1, g2: (float, float)
            Gas fraction produced from side chain and charred bridge,
            respectively

        """
        f = self.broken_bridges(y)
        delta, c = y[1], y[2]
        g1 = 2 * f - delta
        g2 = 2 * (c - self.c0)
        return g1, g2

    def _crosslinking(self, f_meta, T, dt):
        """.Calculate the rate of cross-linking.

        Parameters
        ----------
        f_meta: float
            Fraction of metaplast
        T: float
            Temperature of the particle
        dt: float
            Integral delta time

        Return
        ------
        ratecr: float

        """
        return self.Acr * np.exp(-self.Ecr / Rgas / T) * f_meta * dt

    def _flash_distillation(self, df_gas, df_n, meta_n, mw_n, fracr, T):
        """
        Flash distillation.

        Calculate the flash distillation of tar species from the
        metaplast.

        Parameters
        ----------
        df_gas: float
            Incremental fraction of gas produced in the last time step
        df_n: array
            Incremental fraction of fragments produced in the last time
            step
        meta_n: array
            Fraction of metaplast from the previous time step
        mw_n: array
            Mass weight of the fragments
        fracr: float
            Fraction of metaplast reattached during cross-linking. This
            fraction does not partecipate to evaporation.
        T: float:
            Temperature of the particle

        Return
        ------
        tar_n_new: array
            Fraction of tar produced from flash distillation. This
            fraction is releases in the gas phase.
        meta_n_new: array
            Fraction of metaplast remaining in the particle.

        """
        self.__log.debug('\n\nStart flash_distillation\n')

        a = 87058.0
        b = 299.0
        g = 0.5903
        # mole fraction of n-mers contained in the metaplast
        self.__log.debug('Increment of fragments %s', df_n)
        self.__log.debug('Previous metaplast %s', meta_n)
        self.__log.debug('Cross-linking correction %s', fracr)
        F_n = np.append((df_n + meta_n * fracr) / mw_n, df_gas /
                        self.gasmw)
        if np.allclose(F_n, 0):
            self.__log.debug('F_n = 0 return tar, meta = 0')
            return F_n[:-1], F_n[:-1]
        F_n[F_n < 0] = 0
        self.__log.debug('F_n (mole) %s', F_n)

        F = F_n.sum()
        mw = np.append(mw_n, self.gasmw)
        # self.__log.debug('MW %s', mw)
        p_vap = a * np.exp(-b * mw ** g / T)
        self.__log.debug('p_vap %s', p_vap)
        k_n = p_vap * 101325 / self.pressure
        # self.__log.debug('kn %s', k_n)
        z_n = F_n / F
        # self.__log.debug('zn %s', z_n)
        # Eq. 52
        k_n_1 = k_n - 1

        # def x_n_calc(x): return z_n / (1 + k_n_1 * x)
        # Eq. 54

        # def funct(x): return np.sum(x_n_calc(x) * k_n_1)
        # gradf = lambda x: -np.sum(z_n * k_n_1 / (1 + k_n_1 * x)**2)
        if sum_x_n_calc(0, z_n, k_n_1) * sum_x_n_calc(0.9999, z_n, k_n_1) > 0:
            self.__log.debug('No vapor')
            fract_v = 0
            V = 0
            L = F
            x_n = F_n / F
            y_n = np.zeros_like(x_n)
        else:
            # fract_v = brentq(funct, 0, 0.999)
            fract_v = brentq(sum_x_n_calc, 0, 0.9999, args=(z_n, k_n_1))
            # np.testing.assert_almost_equal(fract_v, fract_v_n)
            # fract_v = newton(funct, 0.5)
            self.__log.debug('V/F = %s', fract_v)
            V = fract_v * F  # moles of tar
            L = F - V
            # mole fraction of n-mers in the metaplast
            x_n = x_n_calc(fract_v, z_n, k_n_1)
            # mole fraction of n-mers released as tar
            y_n = k_n * x_n if V > 0 else np.zeros_like(x_n)

        meta_n_new = x_n[:-1] * L * mw_n
        tar_n_new = y_n[:-1] * V * mw_n
        # assert np.allclose(
        #    meta_n_new + tar_n_new, F_n[:-1] * mw_n), \
        #    'Sum of xn+yn should be equal to Fn'
        self.__log.debug('metaplast fraction %s', x_n)
        self.__log.debug('tar fraction %s', y_n)

        self.__log.debug('Mass metaplast %s', meta_n_new.sum())
        self.__log.debug('Mass tar %s', tar_n_new.sum())

        return tar_n_new, meta_n_new

    def find_triangle(self, plot=False):
        """Find triangle for Genetti light gas correlation."""
        def distance(p1, p2):
            """Calc distance between two points."""
            d = p1 - p2
            return np.inner(d, d)

        points = np.array([[0.017773400000000002, 0.67172399999999999],
                           [0.020365399999999999, 0.58109549999999999],
                           [0.065940100000000001, 0.65505270000000004],
                           [0.077346499999999999, 0.80886970000000002],
                           [0.089362300000000006, 0.7575807],
                           [0.1077369, 0.85064280000000003],
                           [0.13058030000000001, 0.76711629999999997],
                           [0.13575690000000001, 0.85231900000000005],
                           [0.18034790000000001, 0.84992210000000001],
                           [0.20934410000000001, 0.78908880000000003],
                           [0.2603201, 0.85729379999999999],
                           [0.068699999999999997, 0.86299999999999999]])

        triangle_vertices = np.array([[0, 1, 2],
                                      [0, 4, 2],
                                      [0, 3, 4],
                                      [0, 11, 3],
                                      [11, 3, 5],
                                      [3, 5, 4],
                                      [4, 5, 6],
                                      [5, 7, 6],
                                      [7, 8, 6],
                                      [6, 8, 9],
                                      [8, 9, 10]])

        triangles = [pkp.triangle.Triangle(
            *(points[ti] for ti in t)) for t in triangle_vertices]

        # search triangle
        self.triangle = None
        for t, t_v in zip(triangles, triangle_vertices):
            if t.is_inside(self.van_kravelen):
                self.triangle = t
                self.triangle_coals = t_v
                self.triangle_weights = t.weights(self.van_kravelen)
                self.__log.debug('Find triangle %s %s', t, t_v)
                break
        if self.triangle is None:
            # add here a plot
            self.__log.error('Defined coal is outside of the coal triangles '
                             'defined in the Genetti correlation for light '
                             'gases\n'
                             'Light gases will not be calculated!')
            plot = 'show'
            # stop_calculation = True
            # we should check for the closest point

            distances = np.array(
                [distance(p, self.van_kravelen) for p in points])
            ref_coal = distances.argmin()
            self.__log.error('Closest coal is %s', ref_coal)

            for triangle, vertices in zip(triangles, triangle_vertices):
                if ref_coal in vertices:
                    self.triangle = triangle
                    self.triangle_coals = vertices
                    # self.triangle_weights = t.weights(self.van_kravelen)
                    self.triangle_weights = np.array([1 if ref_coal == v
                                                      else 0
                                                      for v in vertices])
                    break

        print('triangle='.format(self.triangle))
        if plot:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots()
            ax.scatter(points[:, 0], points[:, 1], s=100, c='black')
            for n, point in zip(np.arange(len(points)), points):
                ax.annotate(n, xy=point, xytext=(-5, 5),
                            textcoords='offset points')
            for t in triangles:
                t.plot(ax, color='black', linewidth=1)
            ax.set_xlabel('O:C')
            ax.set_ylabel('H:C')
            ax.set_xlim(xmin=0)
            if plot == 'show':
                ax.scatter(self.van_kravelen[0], self.van_kravelen[
                           1], c='red', s=100, label=self.name)
            name = os.path.join(self.path, self.basename + '-van_kravelen.png')
            ax.set_title('Van Kravelen diagram')
            ax.legend(loc='lower right')
            fig.savefig(name)
            plt.close(fig)

        # if stop_calculation:
        #    raise ValueError('Triangle not found for evaluating light gases\n'
        #                     'Check plot')

    def calc_lightgases(self, y):
        r"""
        Calculate lightgas distribution using Genetti method.

        Parameters
        ----------
        y_gas: array
            Fraction of light gases released during time. It is calculated
            using
            :math::`Y_g = 1-(\delta/2 + \pound)(\delta/2 + \pound)_0`

        """
        y_refs = np.array(
            [[np.interp(y, x_gas[i], y_gas[j, i]) for j in range(4)]
             for i in self.triangle_coals])
        return np.dot(y_refs.T, self.triangle_weights)

    def postprocess_step(self, t, y):
        """Post process the results after the time step."""
        solid, gas, tar, meta, cross = self.f[-1]
        T = y[-1]  # temperature
        dt = t - self.t_old

        if meta > 1e-4:
            # cross linking
            rate_cross = self._crosslinking(meta, T, dt)
            fract = 1 - rate_cross / meta
            meta -= rate_cross
            cross += rate_cross
        else:
            rate_cross = 0
            fract = 1
        self.__log.debug(
            'Crosslinking rate: %s / %s', rate_cross, fract)
        # TODO percolation use self.n_frag
        percolation = self._percolation(y, tar, n_frag=self.n_frag,
                                        in_tar=True)

        # gas formed in the last step
        df_gas = max(percolation['f_gas'] - gas, 0)
        self.__log.debug('gas=%s, df_gas=%s', gas, df_gas)
        mw_n = percolation['mw_frag_n']
        # fragments formed in the last step
        df_n = percolation['f_frag_n'] - self.f_frag_n
        self.__log.debug('df_n=%s', df_n)

        tar_n, self.meta_n = self._flash_distillation(
            df_gas=df_gas, df_n=df_n, meta_n=self.meta_n,
            mw_n=mw_n, fracr=fract, T=T)

        # store results
        self.f_frag_n = percolation['f_frag_n']
        gas = percolation['f_gas']
        tar += tar_n.sum()
        meta = self.meta_n.sum()
        solid = 1 - tar - gas
        self.f.append([solid, gas, tar, meta, cross])
        self.__log.debug('F=%s', self.f[-1])
        self.t_old = t

    def postprocess(self, t, y):
        """Postprocess results."""
        # stack y with f
        data = np.hstack([t[:, np.newaxis], y, self.f])
        columns = ['t', 'l', 'delta', 'c', 'T',
                   'solid', 'gas', 'tar', 'meta', 'cross']

        # calc light gas
        X_gas = y[:, 1] * 0.5 + y[:, 0]  # 1/2 delta + l
        X_gas = 1 - X_gas / X_gas[0]  # check if X_gas is different from f gas
        self.find_triangle()
        if self.triangle:
            X_gases = self.calc_lightgases(X_gas)
            data = np.hstack([data, X_gases])
            columns += ['CO', 'CO2', 'H2O', 'CH4']
        data = pd.DataFrame(data=data, columns=columns)
        if self.triangle:
            data['others'] = 1 - data[
                ['CO', 'CO2', 'H2O', 'CH4']].sum(axis=1)
        return data

    def get_yield(self, t, y):
        """Get the volatile yield."""
        return self.f[-1][1] + self.f[-1][2]
