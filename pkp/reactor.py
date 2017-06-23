
"""
Reactor module.

Define the :class:`pkp.reactor.Reactor`, which is used for prescribing
the operating conditions.
"""

from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals

import numpy as np
from autologging import logged

from .empirical_model import SFOR, SFORT, C2SM, DAEM
from .polimi import Polimi
from .cpd import CPD
from scipy.integrate import ode
import warnings

from .interpolate import interp


@logged
class Reactor(object):
    """Base class for running devolatilization simulations."""

    _ode_parameters = {'first_step': 1e-5,
                       'max_step': 1e-3}

    def __init__(self, model=None, **kwargs):
        """
        Init reactor object.

        Parameters
        ----------
        kinetic_model: pkp.empirical_model
            Model object

        """
        # super(Reactor, self).__init__()
        self._model = eval(model)(**kwargs)
        self.operating_conditions = None

    @property
    def operating_conditions(self):
        """
        Define operating conditions of the reactor.

        Operating conditions for devolatilization. They are defined as
        list of operating points [[t0, T0], [t1, T1], ..., [tn, Tn]]
        Each operating point is defined by the time in second and
        temperature in K.
        """
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
        """Get initial solution vector for the reactor."""
        return np.append(self._model.y0, self.operating_conditions[0, 1])

    def run(self, t=None):
        """
        Run reactor for a given time.

        _Parameters
        ----------
        t: np.array, list, default=None
            Time array. This is used to take results from the ODE
            solver. If None times are automatically taken from
            the solver
        """
        # backend = 'dopri5'

        solver = ode(self.rate, jac=None)

        t0 = self.operating_conditions[0, 0]
        solver.set_initial_value(self.y0, t0)

        # define the arguments for running the ODE solver
        args = [solver]
        ode_args = dict(self._ode_parameters)
        # ode_args = {
        #    'first_step': 1e-5,
        #    'max_step': 1e-3,
        # }
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

        # return t, np.squeeze(y)
        return self.model.postprocess(t, y)

    def rate(self, t, y):
        """Rate for the ode integral."""
        return np.concatenate([self._model.rate(t, y), [self._dTdt(t, y)]])

    def _run_nostop(self, solver):
        """
        Run ODE solver with dense output.

        _Parameters
        ----------
        solver: scipy.integrate.ode

        Returns
        -------
        t, y: np.ndarray
            Time and yields arrays.

        """
        solver._integrator.iwork[2] = -1
        warnings.filterwarnings("ignore", category=UserWarning)
        time_end = self.operating_conditions[-1, 0]

        t = [0]
        y = [np.array(self.y0)]
        while solver.t < time_end:
            solver.integrate(time_end, step=True)
            # print(solver.t, solver.y, self.rate(
            #    solver.t, solver.y), self.parameters.y0 - solver.y[0])
            self.model.postprocess_step(solver.t, solver.y)
            t.append(solver.t)
            y.append(solver.y)

        return np.array(t), np.array(y)

    def _run_t(self, solver, t):
        """
        Run the ODE solver stopping at the prescribed time steps.

        _Parameters
        ----------
        solver: scipy.integrate.ode

        Returns
        -------
        t, y: np.ndarray
            Time and yields arrays.

        """
        # self.__log.info('Solver backend %s', solver)
        y = []
        t_calc = []
        for ti in t:
            solver.integrate(ti)
            y.append(solver.y)
            t_calc.append(solver.t)
            self.model.postprocess_step(solver.t, solver.y)
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

    @property
    def model(self):
        """Devolatilization Model used inside the reactor."""
        return self._model

    @property
    def reactor_parameters(self):
        """Reactor parameters dictionary."""
        return self._ode_parameters

    @property
    def model_parameters(self):
        """Model paramerers dictionary"""
        return self.model.parameters_dict

    def set_parameters(self, **kwargs):
        """
        Set the parameters.

        Keep the old values constant.

        Example
        -------
        ::
            # reactor.set_parameters(A=1e5)

        The command set the parameters A in the Model
        """
        model_parameters = self._model.parameters_dict
        # kwargs_model = {key: value for key, value in kwargs.items()
        #                if key in model_parameters}
        for key, value in kwargs.items():
            if key in model_parameters:
                model_parameters[key] = value
            elif key in self._ode_parameters:
                self._ode_parameters[key] = value

        self._model.set_parameters(**model_parameters)


@logged
class DTR(Reactor):
    """
    Class for Drop Tube Reactor.

    The gas temperature is assigned by a list/array of operating points (t, T),
    while particle temperature is evaluated from the heat transfer equation.
    """

    def __init__(self, model, **kwargs):
        """Init DTR."""
        super(DTR, self).__init__(model=model, **kwargs)
        self.density = 800
        self.dp = 100e-6
        self.daf = 0.9
        self.Vp = np.pi * self.dp ** 3 / 8
        self.Ap = np.pi * self.dp ** 2
        self.mp0 = self.density * self.Vp
        self.mash = self.mp0 * (1 - self.daf)
        self.mdaf = self.mp0 * self.daf
        self.cp = 1600
        self.h = 2000
        self.h_pyro = 0

    def calc_mass(self, y):
        """Calc mass of the particle for the given volatile yield y."""
        return self.mash + (1 - y) * self.mdaf

    @property
    def operating_conditions(self):
        """
        Operating conditions for devolatilization.

        They are defined as
        list of operating points [[t0, T0], [t1, T1], ..., [tn, Tn]]
        Each operating point is defined by the time in second and
        temperature in K.
        """
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

        def interp_tT(t):
            """Interpolate time with temperature."""
            return interp(t, conditions[:, 0], conditions[:, 1])
        self.Tg = interp_tT

    def _dTdt(self, t, yt):
        """Temperature time derivative."""
        T = yt[-1]
        # TODO y is not the correct volatile yield
        y = self.model.get_yield(t, yt)
        mp = self.calc_mass(y)
        Tg = self.Tg(t)
        h = 2 * 0.026 / self.dp
        qconv = h * self.Ap * (Tg - T)
        rate = self._model.rate(t, yt)[0]
        qpyro = -rate * self.mdaf * self.h_pyro
        qrad = 0
        qtot = qconv + qpyro + qrad
        return qtot / mp / self.cp
