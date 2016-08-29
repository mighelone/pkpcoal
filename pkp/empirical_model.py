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


@logged
class EmpiricalModel(pkp.reactor.Reactor):
    '''
    Parent class for model.
    y is generally considered as the volatile yield released in the gas
    phase. y
    '''
    parameters_names = ['foo', 'bar']
    parameters_default = [1, 1]
    # initial volatile yield
    y0 = 0

    def __init__(self, parameters=None):
        self.parameters = parameters

    @property
    def len_parameters(self):
        return len(self.parameters_names)

    @property
    def parameters(self):
        return self._parameters

    @parameters.setter
    def parameters(self, values):
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
        #backend = 'vode'
        t0 = self.operating_conditions[0, 0]
        solver = ode(self.rate).set_integrator(backend, nsteps=1,
                                               first_step=1e-6,
                                               max_step=1e-4,
                                               verbosity=1)
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

    def unscale_parameters(self, norm_parameters,
                           parameters_min, parameters_max):
        return 0


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
             np.exp(self.parameters['E'] / 8314.33 / self.T(t)))
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
