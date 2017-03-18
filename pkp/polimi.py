from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals
from builtins import dict

import pkp.detailed_model
try:
    import cantera
except ModuleNotFoundError:
    raise ModuleNotFoundError('Polimi model needs Cantera to be used!')

import numpy as np
import warnings
import pandas as pd
import os
import sys
from autologging import logged

from pkp.detailed_model import M_elements
from scipy.integrate import ode

from .triangle import Triangle


def set_reference_coal(name, atoms):
    '''
    Set the reference coals based on the atomic composition

    Parameters
    ----------
    name: str
        Coal name
    atoms: dict
        Atom coefficient of the coal raw molecule

    Returns
    -------
    pkp.detailed_model.DetailedModel
    '''
    atoms['N'] = 0
    atoms['S'] = 0
    ua = {el: (val * M_elements[el])
          for el, val in atoms.items()}
    return pkp.detailed_model.DetailedModel(
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
    '''
    Raise an exception if the mechanism for Polimi model is wrong.
    '''
    pass


class CompositionError(Exception):
    '''
    Raise an exception of composition is wrong.
    '''
    pass


@logged
class TriangleCoal(Triangle):
    '''
    Triangle class based on coal Van Kravelen diagram
    '''

    headers = ['O:C', 'H:C']

    def __init__(self, coal0, coal1, coal2):
        '''
        Parameters
        ----------
        coal0, coal1, coal2: Coal, Polimi
            Coal vertices of the triangles, based on VK diagram
        '''
        self.coal0 = coal0
        self.coal1 = coal1
        self.coal2 = coal2
        super(TriangleCoal, self).__init__(
            x0=coal0.van_kravelen,
            x1=coal1.van_kravelen,
            x2=coal2.van_kravelen)

    @staticmethod
    def _coal_to_x(coal):
        if isinstance(coal, (Polimi,
                             pkp.detailed_model.DetailedModel)):
            return coal.van_kravelen
        else:
            return coal

    def is_inside(self, coal):
        '''
        Check if the coal is inside the Triangles.
        '''
        return super(
            TriangleCoal, self).is_inside(
                self._coal_to_x(coal))

    def weights(self, coal):
        return super(TriangleCoal, self).weights(
            self._coal_to_x(coal))

    def _coeff(self, coal):
        return super(TriangleCoal, self)._coeff(
            self._coal_to_x(coal))

    def itercoals(self):
        '''
        Iterate over coals returning coal vertices
        '''
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
class Polimi(pkp.detailed_model.DetailedModel):
    '''
    Polimi Multiple Step Kinetic Model for coal devolatilization
    Based on Sommariva (2010).
    '''
    tar = ['VTAR1', 'VTAR2', 'VTAR3']
    light_gas = ['CO', 'CO2', 'H2O', 'H2', 'CH4', 'CH2', 'CH2', 'CH3O',
                 'BTX2']
    raw = ['COAL1', 'COAL2', 'COAL3']
    metaplast = ['GCH2', 'TAR1', 'GBTX2', 'GCH4', 'GCOH2',
                 'GCO2S', 'GH2O', 'GCOL', 'TAR2', 'GCO2TS',
                 'GCOAL3', 'GCO2', 'TAR3', 'GCOLS']
    char = ['CHAR', 'CHARH', 'CHARG']

    def __init__(self, proximate_analysis, ultimate_analysis,
                 pressure=101325, name='Coal'):
        '''
        Parameters
        ----------
        proximate_analysis: dict
        ultimate_analysis: dict
        pressure: float
        name: str
        '''
        super(Polimi, self).__init__(
            proximate_analysis=proximate_analysis,
            ultimate_analysis=ultimate_analysis,
            pressure=pressure,
            name=name)
        self.mechanism = None
        self.skip = 1
        self.backend = None
        self._define_triangle()

    def set_parameters(self, **kwargs):
        '''
        Set parameters for Polimi model.

        Parameters
        ----------
        mechanism: string_types
            Polimi mechanism in Cantera format
        backend: string_types
            Backend for ODE solver.

        See also
        --------
        :meth:`backend`
        :meth:`mechanism`
        :meth:`skip`
        '''
        for key in ('mechanism', 'backend', 'skip'):
            if key in kwargs:
                setattr(self, key, kwargs[key])

    @property
    def backend(self):
        '''
        Backend for the ODE solver. Defauult is dopri. Possible values
        are 'dopri5', 'cvode', 'lsoda', 'dop853'.

        See also
        --------
        scipy.interp1d.ode
        '''
        return self._backend

    @backend.setter
    def backend(self, value=None):
        backend_keys = ['dopri5', 'cvode', 'lsoda', 'dop853']
        if value is None:
            self._backend = backend_keys[0]
        elif value in backend_keys:
            self._backend = value
        else:
            raise ValueError('Backend {} not allowed\n'
                             'Use {}'.format(value, backend_keys))

    # @property
    # def mechanism(self):
    def _get_mechanism(self):
        return self._mechanism

    # @mechanism.setter
    # def mechanism(self, value=None):
    def _set_mechanism(self, value=None):
        '''
        Set mechanism. Default is COAL.xml
        '''
        if value is None:
            value = os.path.join(os.path.dirname(pkp.bins.__file__),
                                 'COAL.xml')
        try:
            self._mechanism = cantera.Solution(value)
            self.mechanism.TP = 300, self.pressure
        except:
            raise MechanismError('Cannot read {}'.format(value))

    mechanism = property(_get_mechanism, _set_mechanism,
                         doc='Mechanism in cantera format for Polimi')

    def _define_triangle(self):
        '''
        Define in which triangle is the coal and calculate the
        composition based on the reference coals
        '''
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

    def run(self):
        '''
        Run Polimi model.

        Returns
        -------
        results: pandas.DataFrame
        '''
        mechanism = self.mechanism

        def dmidt(t, m):
            '''Calculate reaction rates'''
            mechanism.TPY = self.T(t), self.pressure, m
            return (mechanism.net_production_rates *
                    mechanism.molecular_weights / mechanism.density)

        backend = self.backend
        t0 = self.operating_conditions[0, 0]
        mechanism.Y = self.composition
        m0 = mechanism.Y

        solver = ode(dmidt).set_integrator(backend,
                                           nsteps=1)
        solver.set_initial_value(m0, t0)
        solver._integrator.iwork[2] = -1
        t = [t0]
        y = [m0]
        r = [dmidt(t0, m0)]
        warnings.filterwarnings("ignore", category=UserWarning)
        time_end = self.operating_conditions[-1, 0]

        while solver.t < time_end:
            # print(solver.t)
            # self.__log.debug('t=%s', solver.t)
            solver.integrate(time_end, step=True)
            t.append(solver.t)
            y.append(solver.y)
            r.append(dmidt(solver.t, solver.y))

        self.__log.debug('Set skip=%s', self.skip)
        t = np.array(t)[::self.skip]
        y = np.array(y)[::self.skip]
        data = pd.DataFrame(data=y,
                            columns=mechanism.species_names,
                            index=t)
        data.index.name = 'Time, s'
        data['T'] = self.T(t)
        for v in ('metaplast', 'char', 'raw', 'tar', 'light_gas'):
            data[v] = data[getattr(self, v)].sum(axis=1)

        data['solid'] = data[['metaplast', 'char', 'raw']].sum(axis=1)
        data['volatiles'] = data[['tar', 'light_gas']].sum(axis=1)

        self.__log.debug('Write %s', self._out_csv)
        data.to_csv(self._out_csv)

        return data

    @classmethod
    def reference_coal(cls, ref_coal, proximate_analysis=None,
                       pressure=101325):
        '''
        Define polimi model directly defining one of the reference
        coals.

        Parameters
        ----------
        ref_coal: str
            Name of the reference coals: COAL1, COAL2, COAL3, CHAR
        '''
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
