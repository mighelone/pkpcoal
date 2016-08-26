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

from scipy.integrate import ode


class EmpiricalModel(pkp.detailed_model.Reactor):
    '''
    Parent class for model
    '''
    parameters_names = ['foo', 'bar']
    parameters_default = [1, 1]
    y0 = 1

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


class SFOR(EmpiricalModel):
    '''
    Single First Order Reaction (SFOR) model
    '''
    parameters_names = ['A', 'E', 'y0']
    parameters_default = [1e5, 50e6, 0.6]

    def rate(self, t, y):
        '''
        SFOR reaction rate, 1/s
        '''
        k = (self.parameters['A'] /
             np.exp(self.parameters['E'] / 8314.33 / self.T(t)))
        return k * (1 - y - self.parameters['y0'])
