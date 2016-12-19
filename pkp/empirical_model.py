'''
This module defines the class for the empirical pyrolysis models.
It contains classes for the following models:

* Single First Order Reaction (SFOR) model
    :class:`pkp.empirical_model.SFOR`
* Competing 2-Step Model (C2SM)
    :class:`pkp.empirical_model.C2SM`
* Distributed Activation Energy Model (DAEM)
    :class:`pkp.empirical_model.DAEM`
* Biagini-Tognotti model
    :class:`pkp.empirical_model.Biagini`
'''
from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals

import numpy as np

import pkp.detailed_model
import warnings
from autologging import logged

from scipy.integrate import ode
from scipy.misc import factorial
import collections

Rgas = 8314.33
sqrt2 = np.sqrt(2)
sqrtpi = np.sqrt(np.pi)


# http://stackoverflow.com/questions/11351032/named-tuple-and-optional-keyword-arguments
def namedtuple_with_defaults(typename, field_names, default_values=(),
                             units=None):
    T = collections.namedtuple(typename, field_names)
    T.__new__.__defaults__ = (None,) * len(T._fields)
    if isinstance(default_values, collections.Mapping):
        prototype = T(**default_values)
    else:
        prototype = T(*default_values)
    T.__new__.__defaults__ = tuple(prototype)
    if units is None:
        T.units = ('-') * len(T._fields)
    else:
        T.units = units
    return T


@logged
class EmpiricalModel(pkp.reactor.Reactor):
    '''
    Parent class for model.
    `y` is generally considered as the volatile yield released in the
    gas phase.
    '''
    _Parameters = namedtuple_with_defaults(typename='EmpiricalModel',
                                           field_names=('foo', 'bar'),
                                           default_values=(1, 1))
    _len_parameters = len(_Parameters._fields)
    _mask = np.array([True] * _len_parameters)

    # initial volatile yield
    y0 = 0

    def __init__(self, parameters=None):
        self.parameters = parameters

    @property
    def mask(self):
        return self._mask

    @classmethod
    def parameters_names(cls):
        return cls._Parameters._fields

    @classmethod
    def parameters_default(cls):
        return cls._Parameters.__new__.__defaults__

    @classmethod
    def parameters_units(cls):
        return cls._Parameters.units

    @property
    def len_parameters(self):
        '''
        Lenght of the parameters array/dictionary.
        '''
        return len(self._Parameters._fields)

    def _get_parameters(self):
        return self._parameters

    def _set_parameters(self, values):
        if values is None:
            self._parameters = self._Parameters()
        elif isinstance(values, dict):
            self._parameters = self._Parameters(**values)
        else:
            self._parameters = self._Parameters(*values)

    parameters = property(_get_parameters, _set_parameters,
                          doc=(
                              '_Parameters of the empirical models.'
                              ' They can be given as list/numpy array '
                              'or dictionary'))

    def run(self, t=None):
        '''
        Solve model using a ODE solver

        _Parameters
        ----------
        t: np.array, list, default=None
            Time array. This is used to take results from the ODE
            solver. If None times are automatically taken from
            the solver
        '''
        backend = 'dopri5'
        solver = ode(self.rate)
        t0 = self.operating_conditions[0, 0]
        solver.set_initial_value(self.y0, t0)

        # define the arguments for running the ODE solver
        args = [solver]
        ode_args = {
            'first_step': 1e-6,
            'max_step': 1e-4,
            'verbosity': 1
        }
        if t is None:
            ode_args['nsteps'] = 1
            ode_run = self._run_nostop
        else:
            ode_run = self._run_t
            ode_args['nsteps'] = 5000
            args.append(t)

        solver.set_integrator(backend, **ode_args)

        t, y = ode_run(*args)

        return t, np.squeeze(y)

    def _run_nostop(self, solver):
        '''
        Run the ODE solver stopping at then internal time step of the
        solver.

        _Parameters
        ----------
        solver: scipy.integrate.ode

        Returns
        -------
        t, y: np.ndarray
            Time and yields arrays.
        '''
        solver._integrator.iwork[2] = -1
        warnings.filterwarnings("ignore", category=UserWarning)
        time_end = self.operating_conditions[-1, 0]

        t = []
        y = []
        while solver.t < time_end:
            solver.integrate(time_end, step=True)
            # print(solver.t, solver.y)
            t.append(solver.t)
            y.append(solver.y)

        return np.array(t), np.array(y)

    def _run_t(self, solver, t):
        '''
        Run the ODE solver stopping at the prescribed time steps.

        _Parameters
        ----------
        solver: scipy.integrate.ode

        Returns
        -------
        t, y: np.ndarray
            Time and yields arrays.
        '''
        # self.__log.info('Solver backend %s', solver)
        y = []
        t_calc = []
        for ti in t:
            solver.integrate(ti)
            # print(solver.t, solver.y)
            y.append(solver.y)
            t_calc.append(solver.t)
            # print(solver.t)

        # if not np.allclose(t, t_calc):
        if not (t == t_calc).all():
            raise RuntimeError('t and t_calc not the same!')

        return np.array(t_calc), np.array(y)

    def rate(self, t, y):
        return 0

    @classmethod
    def unscale_parameters(cls, norm_parameters, parameters_min,
                           parameters_max):
        '''
        Unscale normalized parameters.
        _Parameters defined in `mask` are unscaled using the log values
        of the minimum and maximum parameters.

        .. math::
            p =  P (log_{10}(p_{max}) - log_{10}(p_{min})) +
            log_{10}(p_{min})

        _Parameters
        ----------
        norm_parameters: iterable
            List of normalized between 0 and 1 parameters
        parameters_min: iterable
            Minimum parameters
        parameters_max: iterable
            Maximum parameters

        Return
        ------
        unsc_par: array
            Unscaled paramters
        '''
        parameters_min = np.array(parameters_min)
        parameters_max = np.array(parameters_max)
        norm_parameters = np.array(norm_parameters)

        mask = np.array(cls._mask)
        parameters_min[mask] = np.log10(parameters_min[mask])
        parameters_max[mask] = np.log10(parameters_max[mask])

        unsc_par = (parameters_min + norm_parameters *
                    (parameters_max - parameters_min))

        # calculate A = 10^log10(A)
        unsc_par[mask] = np.power(10, unsc_par[mask])

        return unsc_par

    @classmethod
    def scale_parameters(cls, parameters, parameters_min,
                         parameters_max):
        '''
        Scale/normalize parameters using minimum and maximum values.
        _Parameters defined in `mask` are scaled using log values for
        the parameters.

        .. math::
            P = (log_{10}(p) - log_{10}(p_{min}))/
            (log_{10}(p_{max}) - log_{10}(p_{min}))

        _Parameters
        ----------
        parameters: iterable
            List of normalized between 0 and 1 parameters
        parameters_min: iterable
            Minimum parameters
        parameters_max: iterable
            Maximum parameters

        Return
        ------
        sc_par: array
            Scaled paramters
        '''
        if isinstance(parameters, dict):
            parameters = [parameters[p] for p in
                          cls.parameters_names()]
        parameters = np.array(parameters)
        parameters_min = np.array(parameters_min)
        parameters_max = np.array(parameters_max)

        mask = np.array(cls._mask)
        parameters_min[mask] = np.log10(parameters_min[mask])
        parameters_max[mask] = np.log10(parameters_max[mask])
        parameters[mask] = np.log10(parameters[mask])

        sc_par = ((parameters - parameters_min) /
                  (parameters_max - parameters_min))
        return sc_par


@logged
class SFOR(EmpiricalModel):
    '''
    The Single First Order Reaction (SFOR) model is characterized by a
    unique reaction, which converts raw coal into char and volatiles
    with a constant ratio.

    .. math::

        Raw -> y_0 volatiles + (1-y_0) char

    The reaction rate is given by:

    .. math::

        r = k(T) (y_0-y)

    with:

    .. math::

        k(T) = A e^{-E/Rg T}

    Where :math:`A` is the pre-exponential factor and :math:`E` is the
    activation energy. The model is characterized by a constant
    volatile yield :math:`y_0`, which does not depend on the conditions
    (heating rate and maximum temperature) of the devolatilization
    process.
    '''
    _Parameters = namedtuple_with_defaults(
        typename='SFOR',
        field_names=('A', 'E', 'y0'),
        default_values=(1e5, 50e6, 0.6),
        units=('1/s', 'J/kmol', '-'))
    _mask = np.array([True, False, False])

    def rate(self, t, y):
        '''
        Reaction rate used in the ODE solver.

        _Parameters
        ----------
        t: float
            Time
        y: float
            Volatile yield :math:`y(t)`

        Returns
        -------
        rate: float
            :math:`dy/dt`
        '''
        k = (self.parameters.A /
             np.exp(self.parameters.E / Rgas / self.T(t)))
        # return k * (1 - y - self.parameters['y0'])
        return k * (self.parameters.y0 - y)


@logged
class SFORT(SFOR):
    '''
    SFOR model with temperature threasold
    '''
    _Parameters = namedtuple_with_defaults(
        typename='SFORT',
        field_names=('A', 'E', 'y0', 'T'),
        default_values=(1e5, 50e6, 0.6, 500),
        units=('1/s', 'J/kmol', '-', 'K'))

    _mask = np.array([True, False, False, False])

    def rate(self, t, y):
        T = self.T(t)
        if T >= self.parameters.T:
            return super(SFORT, self).rate(t, y)
        else:
            return 0


@logged
class C2SM(EmpiricalModel):
    '''
    The Competing 2 Step Model(C2SM) is characterized by a
    two competing reaction with different activation energies and
    stoichiometry. The two reactions produce char and volatiles with
    stoichiometric fraction: math: `y_1` and: math: `y_2`:

    .. math::

        Raw -> y_1 volatiles + (1 - y_1) char

        Raw -> y_2 volatiles + (1 - y_2) char

    with rates: math: `R_1` and: math: `R_2`.
    The reaction rates are given by:

    .. math::

        r_1 = dy_1 / dt = y_1 k_1 s

        r_2 = dy_2 / dt = y_1 k_2 s

    and the overall release of volatiles is:

    .. math::

        r = (y_1 k_1 + y_2 k_2) s

    where: math: `s` is the remaining fraction of raw coal.
    The consumption of raw is given by:

    .. math::

        ds / dt = -(k_1 + k_2) s

    The reaction constants are given by:

    .. math::
        k_1(T) = A_1 e ^ {-E_1 / Rg T}

        k_2(T) = A_2 e^{-E_2/Rg T}

    Generally the first reactin is characterized by low activation
    energy with low release of volatiles, while the second by high
    activation energy and volatiles.
    '''
    _Parameters = namedtuple_with_defaults(
        typename='C2SM',
        field_names=('A1', 'E1', 'y1', 'A2', 'E2', 'y2'),
        default_values=(49e3, 34e6, 0.41, 7.2e7, 95e6, 0.58),
        units=('1/s', 'J/kmol', '-', '1/s', 'J/kmol', '-'))
    _mask = np.array([True, False, False, True, False, False])

    y0 = [0, 1]  # volatile yield, raw solid

    def rate(self, t, y):
        '''
        Reaction rate used in ODE solver.
        `y` is an array containing the overall volatile yield :math:`y`,
        and :math:`s` the raw coal fraction.

        _Parameters
        ----------
        t: float
            Time
        y: float
            Volatile yield :math:`y(t)`

        Returns
        -------
        rate: float
            :math:`dy/dt`
        '''
        RT = Rgas * self.T(t)
        k1, k2 = self._k(RT)
        dsdt = - (k1 + k2) * y[1]
        dydt = (self.parameters.y1 * k1 +
                self.parameters.y2 * k2) * y[1]
        return np.array([dydt, dsdt])

    def _k(self, RT):
        return (self.parameters.A1 / np.exp(
            self.parameters.E1 / RT),
            self.parameters.A2 / np.exp(
            self.parameters.E2 / RT))


@logged
class DAEM(EmpiricalModel):
    '''
    Calculates the devolatilization reaction using the Distributed
    Activation Energy Model (DAEM), using Hermit-Gaussian quadrature
    [Donskoi2000]_
    '''
    _Parameters = namedtuple_with_defaults(
        typename='DAEM',
        field_names=('A0', 'E0', 'sigma', 'y0'),
        default_values=(1e6, 100e6, 12e6, 0.6),
        units=('1/s', 'J/kmol', 'J/kmol', '-'))
    _mask = np.array([True, False, False, False])
    y0 = [0, 0, 0, 0, 0]

    n_quad = 4
    mt = 0.72
    x = np.array([np.sqrt((3 - np.sqrt(6)) / 2),
                  -np.sqrt((3 - np.sqrt(6)) / 2),
                  np.sqrt((3 + np.sqrt(6)) / 2),
                  -np.sqrt((3 + np.sqrt(6)) / 2)])
    H = 8 * pow(x, 3) - 12 * x  # hermite polynomial
    w = pow(2, n_quad - 1) * factorial(n_quad) * \
        np.sqrt(np.pi) / (pow(n_quad, 2) * pow(H, 2))
    Wm = w * np.exp(pow(x, 2))

    def rate(self, t, y):
        '''y, k0.., kn'''
        # TODO add with parameters

        # self.__log.debug('Em %s', Em)
        dIdt = (self.parameters.A0 *
                np.exp(-self._Em / Rgas / self.T(t)))
        # self.__log.debug('dkdt %s', dkdt)
        coeff1 = self.Wm * self.mt / sqrtpi
        coeff2 = np.exp(-pow((self._Em - self.parameters.E0) /
                             self.parameters.sigma, 2) / 2)
        coeff3 = np.exp(-y[1:]) * dIdt
        # self.__log.debug('coeff: %s %s %s', coeff1, coeff2, coeff3)
        # dydt = (self.parameters['y0'] - y[0]) * \
        #    np.sum(coeff1 + coeff2 + coeff3)
        dydt = self.parameters.y0 * np.sum(coeff1 * coeff2 * coeff3)
        # self.__log.debug('dydt %s', dydt)
        return np.append(dydt, dIdt)

    def _calc_Em(self):
        return (self.parameters.E0 +
                self.x * sqrt2 * self.parameters.sigma * self.mt)

    def _set_parameters(self, parameters):
        super(DAEM, self)._set_parameters(parameters)
        self._Em = self._calc_Em()

    def _get_parameters(self):
        return super(DAEM, self)._get_parameters()

    parameters = property(_get_parameters, _set_parameters)


@logged
class Biagini(EmpiricalModel):
    '''
    Calculates the devolatilization reaction using the Biagini model
    [Biagini2014]_

    It based on the :class:`SFOR` model:

    .. math::

        r = k(T) [y_0(T)-y]

    where the volatile yields depend on the temperature:

    .. math::

        y_0(T) = 1 - exp(-DI T/T_{st})

    Where :math:`DI` is the devolatilization index and
    :math:`T_{st}=1223` the standard temperature for estimating
    volatiles in ASTM.

    '''
    _Parameters = namedtuple_with_defaults(
        typename='Biagini',
        field_names=('A', 'E', 'k'),
        default_values=(1e6, 100e6, 0.5),
        units=('1/s', 'J/kmol', '-'))
    _mask = np.array([True, False, False])
    y0 = 0
    Tst = 1223

    def rate(self, t, y):
        T = self.T(t)
        y0 = 1 - np.exp(-self.parameters.k * T / self.Tst)
        k = (self.parameters.A /
             np.exp(self.parameters.E / Rgas / self.T(t)))
        return k * (y0 - y)
