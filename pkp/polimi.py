"""
Polimi model module.

Define the class to work with Polimi model and for triangulation.
"""
from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals
import sys

import pkp.coal
import pkp.empirical_model
import numpy as np
import pandas as pd
import os
from autologging import logged

from pkp.coal import M_elements

from .triangle import Triangle
from ._exceptions import ImportError


try:
    import cantera
except ImportError:
    raise ImportError('Cantera not installed! Module cannot be used',
                      __name__)


def set_reference_coal(name, atoms):
    """
    Set the reference coals based on the atomic composition.

    Parameters
    ----------
    name: str
        Coal name
    atoms: dict
        Atom coefficient of the coal raw molecule

    Returns
    -------
    pkp.coal.Coal

    """
    atoms['N'] = 0
    atoms['S'] = 0
    ua = {el: (val * M_elements[el])
          for el, val in atoms.items()}
    return pkp.coal.Coal(
        name=name,
        ultimate_analysis=ua,
        proximate_analysis={'FC': 50,
                            'VM': 50,
                            'Ash': 0,
                            'Moist': 0})


# set the reference coals
coal1 = set_reference_coal('COAL1', atoms={'C': 12, 'H': 11, 'O': 0})
coal2 = set_reference_coal('COAL2', atoms={'C': 14, 'H': 10, 'O': 1})
coal3 = set_reference_coal('COAL3', atoms={'C': 12, 'H': 12, 'O': 5})
char = set_reference_coal('CHAR', atoms={'C': 1, 'H': 0, 'O': 0})


# Exceptions
class MechanismError(Exception):
    """Raise an exception if the mechanism for Polimi model is wrong."""

    pass


class CompositionError(Exception):
    """Raise an exception of composition is wrong."""

    pass


@logged
class TriangleCoal(Triangle):
    """Triangle class based on coal Van Kravelen diagram."""

    headers = ['O:C', 'H:C']

    def __init__(self, coal0, coal1, coal2):
        """
        Init Coal Triangle.

        Parameters
        ----------
        coal0, coal1, coal2: Coal, Polimi
            Coal vertices of the triangles, based on VK diagram

        """
        self.coal0 = coal0
        self.coal1 = coal1
        self.coal2 = coal2
        super(TriangleCoal, self).__init__(
            x0=coal0.van_kravelen,
            x1=coal1.van_kravelen,
            x2=coal2.van_kravelen)

    @staticmethod
    def _coal_to_x(coal):
        if isinstance(coal, (Polimi, pkp.coal.Coal)):
            return coal.van_kravelen
        else:
            return coal

    def is_inside(self, coal):
        """Check if the coal is inside the Triangles."""
        return super(
            TriangleCoal, self).is_inside(
                self._coal_to_x(coal))

    def weights(self, coal):
        """Return the weights of the given coal in the triangle."""
        return super(TriangleCoal, self).weights(
            self._coal_to_x(coal))

    def _coeff(self, coal):
        """Return the weights of the coal in the triangle."""
        return super(TriangleCoal, self)._coeff(
            self._coal_to_x(coal))

    def itercoals(self):
        """Iterate over coals returning coal vertices."""
        for c in [self.coal0, self.coal1, self.coal2]:
            yield c


# set the reference triangles
triangle_012 = TriangleCoal(char,
                            coal1,
                            coal2)

triangle_023 = TriangleCoal(char,
                            coal2,
                            coal3)

triangle_123 = TriangleCoal(coal1,
                            coal2,
                            coal3)


@logged
class Polimi(pkp.coal.Coal, pkp.empirical_model.Model):
    """
    Multi-Step Kinetic Devolatilizion Model (Polimi).

    Based on Sommariva (2010).
    """

    tar = ['VTAR1', 'VTAR2', 'VTAR3']
    light_gas = ['CO', 'CO2', 'H2O', 'H2', 'CH4', 'CH2', 'CH2', 'CH3O',
                 'BTX2']
    raw = ['COAL1', 'COAL2', 'COAL3']
    metaplast = ['GCH2', 'TAR1', 'GBTX2', 'GCH4', 'GCOH2',
                 'GCO2S', 'GH2O', 'GCOL', 'TAR2', 'GCO2TS',
                 'GCOAL3', 'GCO2', 'TAR3', 'GCOLS']
    char = ['CHAR', 'CHARH', 'CHARG']

    # define here the modificable parameters
    _parameters = ['mechanism']

    def __init__(self, proximate_analysis=None, ultimate_analysis=None,
                 pressure=101325, name='Coal', **kwargs):
        """
        Init Polimi model from ultimate analysis.

        Parameters
        ----------
        proximate_analysis: dict
        ultimate_analysis: dict
        pressure: float
        name: str

        See also
        --------
        :meth:`reference_coal`

        """
        super(Polimi, self).__init__(
            proximate_analysis=proximate_analysis,
            ultimate_analysis=ultimate_analysis,
            pressure=pressure,
            name=name)
        # this information should be setted in set_parameters
        self.mechanism = None
        self.skip = 1
        self.backend = None
        self._define_triangle()
        self.set_parameters(**kwargs)
        # self.parameters_dict = {}

    def set_parameters(self, **kwargs):
        """
        Set parameters for Polimi model.

        Parameters
        ----------
        mechanism: string_types
            Polimi mechanism in Cantera format

        """
        for key, value in kwargs.items():
            if key in self._parameters:
                setattr(self, key, value)

    def get_parameters(self):
        """Get the parameters dictionary."""
        return {key: getattr(self, key) for key in self._parameters}

    @property
    def parameters_dict(self):
        return self.get_parameters()

    # @property
    # def mechanism(self):
    def _get_mechanism(self):
        return self._mechanism

    # @mechanism.setter
    # def mechanism(self, value=None):
    def _set_mechanism(self, value=None):
        """Set mechanism. Default is COAL.xml."""
        if isinstance(value, cantera.Solution):
            self._mechanism = value
        else:
            if value is None:
                value = os.path.join(os.path.dirname(pkp.bins.__file__),
                                     'COAL.xml')
            try:
                self._mechanism = cantera.Solution(value)
            except:
                raise MechanismError('Cannot read {}'.format(value))
        self._mechanism.TP = 300, self.pressure
        self._calc_light_gas_index()

    mechanism = property(_get_mechanism, _set_mechanism,
                         doc='Mechanism in cantera format for Polimi')

    def _define_triangle(self):
        """Define coal triangle."""
        if triangle_012.is_inside(self):
            self.triangle = triangle_012
            self.inside = '012'
        elif triangle_023.is_inside(self):
            self.triangle = triangle_023
            self.inside = '023'
        elif triangle_123.is_inside(self):
            self.triangle = triangle_123
            self.inside = '123'
        else:
            raise CompositionError('Composition outside triangles')
        self.__log.debug('Coal is inside triangle %s', self.inside)
        self.triangle_weights = self.triangle.weights(self)
        self.composition = {c.name: self.triangle_weights[i]
                            for i, c in enumerate(
            self.triangle.itercoals())}
        self.mechanism.TPY = None, None, self.composition
        self.y0 = self.mechanism.Y

    def rate(self, t, y):
        """Volatilization rate."""
        self.mechanism.TPY = y[-1], self.pressure, y[:-1]
        return (self.mechanism.net_production_rates *
                self.mechanism.molecular_weights / self.mechanism.density)

    @classmethod
    def reference_coal(cls, ref_coal, proximate_analysis=None,
                       pressure=101325):
        """
        Init Polimi using reference coals.

        Define polimi model directly defining one of the reference
        coals.

        Parameters
        ----------
        ref_coal: str
            Name of the reference coals: COAL1, COAL2, COAL3, CHAR

        """
        try:
            c = getattr(sys.modules[__name__], ref_coal.lower())
        except AttributeError('Use one of the reference coal for Polimi'
                              ' model'):
            raise
        except:
            raise
        if proximate_analysis is None:
            proximate_analysis = c.proximate_analysis
        return cls(ultimate_analysis=c.ultimate_analysis,
                   proximate_analysis=proximate_analysis,
                   pressure=pressure,
                   name=c.name)

    def postprocess(self, t, y):
        """Postprocess results."""
        data = np.insert(y, 0, t, axis=1)[::self.skip]

        data = pd.DataFrame(data=data,
                            columns=['t'] + self.mechanism.species_names +
                            ['T'])
        for v in ('metaplast', 'char', 'raw', 'tar', 'light_gas'):
            data[v] = data[getattr(self, v)].sum(axis=1)
        data['solid'] = data[['metaplast', 'char', 'raw']].sum(axis=1)
        data['volatiles'] = data[['tar', 'light_gas']].sum(axis=1)
        return data

    def get_yield(self, t, y):
        """Get the volatile yield."""
        return y[self._light_gas_index].sum()

    def postprocess_step(self, t, y):
        """Post process at time step of the ODE."""
        pass

    def _calc_light_gas_index(self):
        """Get the list of the index of the light gas."""
        self._light_gas_index = [self._mechanism.species_index(sp)
                                 for sp in self.light_gas]
