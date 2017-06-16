'''
Define the :class:`pkp.reactor.Reactor`, which is used for prescribing
the operating conditions.
'''

from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals

from .empirical_model_t import *
from scipy.integrate import ode
import warnings


import numpy as np
from autologging import logged

# try:
#    from .interpolate import interp
# except:
#    print('Can not load numba interpolate, use numpy')
#    from numpy import interp

from .interpolate import interp


@logged
class Reactor(object):
    '''
    Base class for running devolatilization simulations
    '''

    def __init__(self, model=None, parameters=None):
        """
        Parameters
        ----------
        kinetic_model: pkp.empirical_model
            Model object

        """
        # super(Reactor, self).__init__()
        self._model = eval(model)(parameters)
        self.operating_conditions = None

    @property
    def operating_conditions(self):
        '''
        Operating conditions for devolatilization. They are defined as
        list of operating points [[t0, T0], [t1, T1], ..., [tn, Tn]]
        Each operating point is defined by the time in second and
        temperature in K.
        '''
        return self._operating_conditions

    @operating_conditions.setter
    def operating_conditions(self, conditions):
        if conditions is None:
            self.T = None
            self._operating_conditions = None
            return
        if not isinstance(conditions, (np.ndarray, list)):
            raise TypeError('Define conditions as list or numpy array')
        elif isinstance(conditions, list):
            conditions = np.array(conditions)
        if not conditions.ndim == 2:
            raise ValueError('Define conditions as array Nx2')
        if not conditions.shape[-1] == 2:
            raise ValueError('Define conditions as array Nx2')
        self._operating_conditions = conditions

        self._dTdt_array = (np.diff(self._operating_conditions[:, 1]) /
                            np.diff(self._operating_conditions[:, 0]))

        # def interp_tT(t):
        #     '''Interpolate time with temperature'''
        #     return interp(t, conditions[:, 0], conditions[:, 1])
        # self.T = interp_tT

    @property
    def y0(self):
        return self._model.y0 + [self.operating_conditions[0, 1]]

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
        # backend = 'dopri5'

        solver = ode(self.rate, jac=None)

        t0 = self.operating_conditions[0, 0]
        solver.set_initial_value(self.y0, t0)

        # define the arguments for running the ODE solver
        args = [solver]
        ode_args = {
            'first_step': 1e-5,
            'max_step': 1e-3,
        }
        if t is None:
            backend = 'dopri5'
            ode_args['nsteps'] = 1
            ode_args['verbosity'] = 2
            ode_run = self._run_nostop
        else:
            # backend = 'vode'
            backend = 'dopri5'
            ode_run = self._run_t
            ode_args['nsteps'] = 100000
            # ode_args['min_step'] = 1e-13
            args.append(t)

        solver.set_integrator(backend, **ode_args)
        warnings.filterwarnings("ignore", category=UserWarning)
        t, y = ode_run(*args)
        warnings.resetwarnings()

        return t, np.squeeze(y)

    def rate(self, t, y):
        return self._model.rate(t, y) + [self._dTdt(t, y)]

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

        t = [0]
        y = [np.array(self.y0)]
        while solver.t < time_end:
            solver.integrate(time_end, step=True)
            # print(solver.t, solver.y, self.rate(
            #    solver.t, solver.y), self.parameters.y0 - solver.y[0])
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
        if not np.allclose(t, t_calc):
            raise RuntimeError('t and t_calc not the same!')

        return np.array(t_calc), np.array(y)

    def _dTdt(self, t, y):
        t_array = self.operating_conditions[:, 0]
        if t < t_array[0] or t >= t_array[-1]:
            return 0.0
        index = next(
            (idx for idx, val in np.ndenumerate(t_array) if val > t))[0]

        return self._dTdt_array[index - 1]
