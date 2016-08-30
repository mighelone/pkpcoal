'''
Empirical Model module
=============

Define the empirical models:

* SFOR
* C2SM
* DAEM
* ...
'''
from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals

import numpy as np

import pkp.detailed_model
import warnings
from autologging import logged

from scipy.integrate import ode
from scipy.misc import factorial

Rgas = 8314.33
sqrt2 = np.sqrt(2)
sqrtpi = np.sqrt(np.pi)


@logged
class EmpiricalModel(pkp.reactor.Reactor):
    '''
    Parent class for model.
    y is generally considered as the volatile yield released in the gas
    phase. y
    '''
    parameters_names = ['foo', 'bar']
    parameters_default = [1, 1]
    mask = np.array([True] * len(parameters_default))
    # initial volatile yield
    y0 = 0

    def __init__(self, parameters=None):
        self.parameters = parameters

    @property
    def len_parameters(self):
        return len(self.parameters_names)

    def _get_parameters(self):
        return self._parameters

    def _set_parameters(self, values):
        if values is None:
            par_values = self.parameters_default
        elif isinstance(values, (list, np.ndarray)):
            par_values = values
        elif isinstance(values, dict):
            if all(key in self.parameters_names for key in values):
                self._parameters = values
                return
            else:
                raise ValueError(
                    'Parameters should contains the'
                    ' following keys {}'.format(self.parameters_names))

        if len(par_values) > self.len_parameters:
            raise ValueError(
                'Length of parameters list/array should be {}'.format(
                    self.len_parameters))
        else:
            self._parameters = {k: v for k, v in
                                zip(self.parameters_names,
                                    par_values)}

    parameters = property(_get_parameters, _set_parameters)

    def run(self, t=None):
        '''
        Solve model using a ODE solver

        Parameters
        ----------
        t: np.array, list, default=None
            Time array. This is used to take results from the ODE
            solver. If None times are automatically taken from
            the solver
        '''
        backend = 'dopri5'
        # backend = 'vode'
        # vode_settings = {'first_step': 1e-6,
        #                 'max_step': 1e-4}
        t0 = self.operating_conditions[0, 0]
        solver = ode(self.rate).set_integrator(backend, nsteps=1,
                                               first_step=1e-6,
                                               max_step=1e-4,
                                               verbosity=1)
        # solver = ode(self.rate)
        # solver.set_integrator(backend)

        solver.set_initial_value(self.y0, t0)

        if t is None:
            t, y = self._run_nostop(solver)
        else:
            t_calc, y = self._run_t(solver, t)

            if not np.allclose(t, t_calc):
                #raise RuntimeError('t and t_calc not the same!')
                pass
        return t, np.squeeze(y)

    def _run_nostop(self, solver):
        solver._integrator.iwork[2] = -1
        warnings.filterwarnings("ignore", category=UserWarning)
        time_end = self.operating_conditions[-1, 0]

        #t = [self.operating_conditions[0, 0]]
        #y = [self.y0]
        t = []
        y = []
        while solver.t < time_end:
            solver.integrate(time_end, step=True)
            # print(solver.t, solver.y)
            t.append(solver.t)
            y.append(solver.y)

        return np.array(t), np.array(y)

    def _run_t(self, solver, t):
        y = []
        t_calc = []
        for ti in t:
            solver.integrate(ti)
            # print(solver.t, solver.y)
            y.append(solver.y)
            t_calc.append(solver.t)
            # print(solver.t)

        return np.array(t_calc), np.array(y)

    def rate(self, t, y):
        return 0

    @classmethod
    def unscale_parameters(cls, norm_parameters, parameters_min,
                           parameters_max):
        '''
        Unscale normalized parameters.
        A1 and A2 are stored as logA1, logA2

        Return
        ------
        unsc_par: array
            Unscaled paramters
        '''
        parameters_min = np.array(parameters_min)
        parameters_max = np.array(parameters_max)
        norm_parameters = np.array(norm_parameters)

        mask = np.array(cls.mask)
        parameters_min[mask] = np.log10(parameters_min[mask])
        parameters_max[mask] = np.log10(parameters_max[mask])

        unsc_par = (parameters_min + norm_parameters *
                    (parameters_max - parameters_min))

        # calculate A = 10^log10(A)
        unsc_par[mask] = np.power(10, unsc_par[mask])

        return unsc_par


@logged
class SFOR(EmpiricalModel):
    '''
    Single First Order Reaction (SFOR) model.
    The reaction rate is given by: :math:`r = k(T) (y0-y)` with
    :math `k(T) = A \exp(-E/Rg T)`.

    The class inherits the main methods from
    :meth:`empirical_model.EmpiricalModel`
    '''
    parameters_names = ['A', 'E', 'y0']
    parameters_default = [1e5, 50e6, 0.6]

    def rate(self, t, y):
        '''
        SFOR reaction rate, 1/s
        '''
        k = (self.parameters['A'] /
             np.exp(self.parameters['E'] / Rgas / self.T(t)))
        # return k * (1 - y - self.parameters['y0'])
        return k * (self.parameters['y0'] - y)

    @staticmethod
    def unscale_parameters(norm_parameters, parameters_min,
                           parameters_max):
        '''
        Unscale normalized parameters.
        A is stored as logA

        Return
        ------
        unsc_par: array
            Unscaled paramters
        '''
        # todo check if the conversion to array is too expensive!!
        parameters_min = np.array(parameters_min)
        parameters_min[0] = np.log10(parameters_min[0])
        parameters_max = np.array(parameters_max)
        parameters_max[0] = np.log10(parameters_max[0])
        norm_parameters = np.array(norm_parameters)
        unsc_par = (parameters_min + norm_parameters *
                    (parameters_max - parameters_min))

        # calculate A = 10^log10(A)
        unsc_par[0] = np.power(10, unsc_par[0])

        return unsc_par


@logged
class C2SM(EmpiricalModel):
    '''
    Competing 2 Step Model for pyrolysis
    '''
    parameters_names = ['A1', 'E1', 'y1', 'A2', 'E2', 'y2']
    parameters_default = [49e3, 34e6, 0.41, 7.2e7, 95e6, 0.58]
    mask = np.array([True, False, False, True, False, False])

    y0 = [0, 1]  # volatile yield, raw solid

    def rate(self, t, y):
        RT = Rgas * self.T(t)
        k1, k2 = self._k(RT)
        dsdt = - (k1 + k2) * y[1]
        dydt = (self.parameters['y1'] * k1 +
                self.parameters['y2'] * k2) * y[1]
        return np.array([dydt, dsdt])

    def _k(self, RT):
        return (self.parameters['A1'] / np.exp(
            self.parameters['E1'] / RT),
            self.parameters['A2'] / np.exp(
            self.parameters['E2'] / RT))


@logged
class DAEM(EmpiricalModel):
    '''alculates the devolatilization reaction using the Distributed
    Activation Energy Model (DAEM), using Hermit-Gaussian quadrature
    http://www.sciencedirect.com/science/article/pii/S0010218000001152'''

    parameters_names = ['A0', 'E0', 'sigma', 'y0']
    parameters_default = [1e6, 100e6, 12e6, 0.6]
    mask = np.array([True, False, False, False])
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
        dIdt = (self.parameters['A0'] *
                np.exp(-self._Em / Rgas / self.T(t)))
        # self.__log.debug('dkdt %s', dkdt)
        coeff1 = self.Wm * self.mt / sqrtpi
        coeff2 = np.exp(-pow((self._Em - self.parameters['E0']) /
                             self.parameters['sigma'], 2) / 2)
        coeff3 = np.exp(-y[1:]) * dIdt
        # self.__log.debug('coeff: %s %s %s', coeff1, coeff2, coeff3)
        # dydt = (self.parameters['y0'] - y[0]) * \
        #    np.sum(coeff1 + coeff2 + coeff3)
        dydt = self.parameters['y0'] * np.sum(coeff1 * coeff2 * coeff3)
        # self.__log.debug('dydt %s', dydt)
        return np.append(dydt, dIdt)

    def _calc_Em(self):
        return (self.parameters['E0'] +
                self.x * sqrt2 * self.parameters['sigma'] * self.mt)

    def _set_parameters(self, parameters):
        super(DAEM, self)._set_parameters(parameters)
        self._Em = self._calc_Em()

    def _get_parameters(self):
        return super(DAEM, self)._get_parameters()

    parameters = property(_get_parameters, _set_parameters)
