"""
Empirical model.

This module defines the class for the empirical pyrolysis models.
It contains classes for the following models:

* Single First Order Reaction (SFOR) model
    :class:`pkp.empirical_model.SFOR`
* Competing 2-Step Model (C2SM)
    :class:`pkp.empirical_model.C2SM`
* Distributed Activation Energy Model (DAEM)
    :class:`pkp.empirical_model.DAEM`
* Biagini-Tognotti model
    :class:`pkp.empirical_model.BT`
"""
from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals

import numpy as np
import abc
from autologging import logged

from scipy.misc import factorial
import collections

Rgas = 8314.33
sqrt2 = np.sqrt(2)
sqrtpi = np.sqrt(np.pi)


def namedtuple_with_defaults(typename, field_names, default_values=(),
                             units=None):
    """
    Create a namedtuple with default values.

    see
    http://stackoverflow.com/questions/11351032/named-tuple-and-optional-keyword-arguments

    Parameters
    ----------
    typename: str
        Name of the namedtuple
    field_names: list
        List of names of the fields
    default_values: list
        Default values
    units: list
        Units of the parameters

    Returns
    -------
    namedtuple

    """
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
class Model(metaclass=abc.ABCMeta):
    """Abstract class for Model."""

    __metaclass__ = abc.ABCMeta
    # initial volatile yield
    y0 = [0]
    jacob = None

    @abc.abstractmethod
    def rate(self, t, y):
        """Rate for ODE."""
        return

    @abc.abstractmethod
    def set_parameters(self, *args, **kwargs):
        """Set parameters."""
        return

    @abc.abstractmethod
    def postprocess(self, t, y):
        """Post process data at the end of the run."""
        return

    @abc.abstractmethod
    def postprocess_step(self, t, y):
        """Post process data at the end of the time step."""
        return

    @abc.abstractmethod
    def get_yield(self, t, y):
        """Return the actual volatilization yield."""
        return

    @abc.abstractproperty
    def parameters_dict(self):
        """Return dictionary of the actual parameters."""
        return


@logged
class EmpiricalModel(Model):
    """
    Abstract class for empirical models.

    The derived class has to provide the `rate` method.
    """

    _Parameters = namedtuple_with_defaults(typename='EmpiricalModel',
                                           field_names=('foo', 'bar'),
                                           default_values=(1, 1))
    _len_parameters = len(_Parameters._fields)
    _mask = np.array([True] * _len_parameters)

    def __init__(self, *args, **kwargs):
        """Init empirical model."""
        self.set_parameters(*args, **kwargs)

    @property
    def mask(self):
        """Return the mask for the logaritmic property."""
        return self._mask

    def __str__(self):
        """Return description of the model."""
        return self.parameters.__str__()

    def __repr__(self):
        """Return description of the model."""
        return self.__str__()

    @classmethod
    def parameters_names(cls):
        """Return a list with the parameters names."""
        return cls._Parameters._fields

    @classmethod
    def parameters_default(cls):
        """Return a list with the default values of the parameters."""
        return cls._Parameters.__new__.__defaults__

    @classmethod
    def parameters_units(cls):
        """Return a list with the units of the parameters."""
        return cls._Parameters.units

    @property
    def len_parameters(self):
        """Lenght of the parameters array/dictionary."""
        return len(self._Parameters._fields)

    @property
    def parameters(self):
        """List of the actual parameters."""
        return self._parameters

    def set_parameters(self, *args, **kwargs):
        """
        Set the parameters of the model.

        Example
        -------

        Assuming that the model has parameters `x` and `y`::

            >>> model.set_parameters()

        Set the parameters to their default values.::

            >> model.set_parameters(1)

        Set the parameter `x` to 1. x is the first in the list of parameters.::

            >> model.set_parameters(1, 2)

        Set the parameter `x` to 1 and `y` to 2. It follows the order of the
        parameter list::

            >> model.set_parameters(x=1, y=2)

        Again set x to 1 and y to 2::

            >> model.set_parameters(y=2)

        Set y to 2 and x to its default values.
        """
        if len(args) > 0:
            if hasattr(args[0], '__iter__'):
                self._parameters = self._Parameters(*args[0])
            elif args[0] is None:
                self._parameters = self._Parameters()
            else:
                self._parameters = self._Parameters(*args)
        else:
            self._parameters = self._Parameters(**kwargs)

    @property
    def parameters_dict(self):
        """Return dictionary of the actual parameters."""
        return dict(zip(self.parameters_names(), self.parameters_list))

    @property
    def parameters_list(self):
        """Return a list of the actual parameters."""
        return [getattr(self.parameters, p) for p in self.parameters_names()]

    @classmethod
    def unscale_parameters(cls, norm_parameters, parameters_min,
                           parameters_max):
        """
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
        """
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
        """
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

        """
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

    def postprocess(self, t, y):
        """Post process results after ODE."""
        # TODO not sure if it is the right way to return data
        return t, y[:, [0, -1]]

    def postprocess_step(self, t, y):
        """Post process data at the end of the time step."""
        pass

    def get_yield(self, t, y):
        """Return volatile yield."""
        return y[0]


@logged
class SFOR(EmpiricalModel):
    """
    Single First Order Reaction (SFOR) model.

    The SFOR model is characterized by a unique reaction, which converts raw
    coal into char and volatiles with a constant ratio.

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

    """

    _Parameters = namedtuple_with_defaults(
        typename='SFOR',
        field_names=('A', 'E', 'y0'),
        default_values=(1e5, 50e6, 0.6),
        units=('1/s', 'J/kmol', '-'))
    _mask = np.array([True, False, False])
    y0 = [0]

    def rate(self, t, y):
        """
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

        """
        k = self._calc_k(y[1])
        # return k * (1 - y - self.parameters['y0'])
        dy = self.parameters.y0 - y[0]
        return [k * dy if dy > 1e-6 else 0]

    def _calc_k(self, T):
        return (self.parameters.A /
                np.exp(self.parameters.E / Rgas / T))

    # def jacob(self, t, y):
    #    return -self._calc_k(t)
    jacob = None


@logged
class SFORT(SFOR):
    """SFOR model with temperature threasold."""

    _Parameters = namedtuple_with_defaults(
        typename='SFORT',
        field_names=('A', 'E', 'y0', 'T'),
        default_values=(1e5, 50e6, 0.6, 500),
        units=('1/s', 'J/kmol', '-', 'K'))

    _mask = np.array([True, False, False, False])

    def rate(self, t, y):
        """Return reaction rate dy/dt."""
        if y[1] >= self.parameters.T:
            return super(SFORT, self).rate(t, y)
        else:
            return 0


@logged
class C2SM(EmpiricalModel):
    """
    Competing 2 Step Model (C2SM).

    The C2SM is characterized by a two competing reaction with different
    activation energies and stoichiometry. The two reactions produce char and
    volatiles with stoichiometric fraction: math: `y_1` and: math: `y_2`:

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

    """

    _Parameters = namedtuple_with_defaults(
        typename='C2SM',
        field_names=('A1', 'E1', 'y1', 'A2', 'E2', 'y2'),
        default_values=(49e3, 34e6, 0.41, 7.2e7, 95e6, 0.58),
        units=('1/s', 'J/kmol', '-', '1/s', 'J/kmol', '-'))
    _mask = np.array([True, False, False, True, False, False])

    y0 = [0, 1]  # volatile yield, raw solid

    def rate(self, t, y):
        """
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

        """
        k1, k2 = self._k(y[-1])
        if y[1] > 1e-6:
            dydt = [(self.parameters.y1 * k1 + self.parameters.y2 * k2) * y[1],
                    - (k1 + k2) * y[1]]
        else:
            dydt = [0, 0]
        return dydt

    # def jacob(self, t, y):
    #     k1, k2 = self._k(t)
    #     return np.array([[0, (self.parameters.y1 * k1 +
    #                           self.parameters.y2 * k2)],
    #                      [0, -(k1 + k2)]])

    def _k(self, T):
        """Calculate the reaction constants."""
        RT = Rgas * T
        return (self.parameters.A1 / np.exp(
            self.parameters.E1 / RT),
            self.parameters.A2 / np.exp(
            self.parameters.E2 / RT))


@logged
class DAEM(EmpiricalModel):
    """
    Distributed Activation Energy Model.

    Calculates the devolatilization reaction using the Distributed
    Activation Energy Model (DAEM), using Hermit-Gaussian quadrature
    [Donskoi2000]_

    """

    _Parameters = namedtuple_with_defaults(
        typename='DAEM',
        field_names=('A0', 'E0', 'sigma', 'y0'),
        default_values=(1e5, 50e6, 12e6, 0.6),
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

    def rate(self, t, yt):
        """
        Calculate rates.

        The array `yt` contains the overall volatile yields y, and the
        integral quadrature terms k0.., kn

        """
        # TODO add with parameters
        T = yt[-1]
        y = yt[:-1]
        # self.__log.debug('Em %s', Em)
        dIdt = (self.parameters.A0 *
                np.exp(-self._Em / Rgas / T))
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
        """Calculate activation energies for the quadrature points."""
        return (self.parameters.E0 +
                self.x * sqrt2 * self.parameters.sigma * self.mt)

    def set_parameters(self, *args, **kwargs):
        """Set parameters and recalculate energies."""
        super(DAEM, self).set_parameters(*args, **kwargs)
        self._Em = self._calc_Em()

    # parameters = property(_get_parameters, _set_parameters)


@logged
class BT(SFOR):
    """
    Biagini Tognotti (BT) model.

    Calculates the devolatilization reaction using the Biagini Tognotti (BT)
    model [Biagini2014]_.

    It based on the :class:`SFOR` model:

    .. math::

        r = k(T) [y_0(T)-y]

    where the volatile yields depend on the temperature:

    .. math::

        y_0(T) = 1 - exp(-DI T/T_{st})

    Where :math:`DI` is the devolatilization index and
    :math:`T_{st}=1223` the standard temperature for estimating
    volatiles in ASTM.

    """

    _Parameters = namedtuple_with_defaults(
        typename='Biagini',
        field_names=('A', 'E', 'k'),
        default_values=(1e5, 50e6, 0.5),
        units=('1/s', 'J/kmol', '-'))
    _mask = np.array([True, False, False])
    y0 = [0]
    Tst = 1223

    def rate(self, t, y):
        """Reaction rate."""
        T = y[-1]
        y0 = self._calc_y0(T)
        dy = (y0 - y[0])
        k = self._calc_k(T)
        return [k * dy]

    def _calc_y0(self, T):
        """Calculate y0."""
        return 1 - np.exp(-self.parameters.k * T / self.Tst)
