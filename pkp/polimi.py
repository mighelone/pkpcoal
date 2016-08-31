from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals

import pkp.detailed_model
import cantera
import numpy as np
import tabulate
import itertools
import warnings
import pandas as pd
import os
from autologging import logged

from pkp.detailed_model import M_elements
from scipy.integrate import ode


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
          for el, val in atoms.iteritems()}
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
    pass


class CompositionError(Exception):
    pass


class OutsideTriangleError(Exception):
    pass


@logged
class Triangle(object):
    '''
    Triangle class. Used for triangulation calculations
    '''
    headers = ['x', 'y']

    def __init__(self, x0=None, x1=None, x2=None):
        '''
        Parameters
        ----------
        x0, x1, x2, np.ndarray, list
            2D Array or list of the triangle vertices
        '''
        if x0 is None:
            x0 = np.array([0, 0])
        if x1 is None:
            x1 = np.array([1, 0])
        if x2 is None:
            x2 = np.array([0, 1])
        self.x0 = np.array(x0)
        self.x1 = np.array(x1)
        self.x2 = np.array(x2)

    def _coeff(self, x):
        '''
        Calculate coefficient of linear combination of x
        x-x0 = a*(x1-x0)+b*(x2-x0)
        '''
        v1 = self.x1 - self.x0
        v2 = self.x2 - self.x0
        v = x - self.x0
        matr = np.transpose(np.array([v1, v2]))
        return np.linalg.solve(matr, v)

    def is_inside(self, x):
        '''
        verify is point x is inside the triangle
        '''
        coeff = self._coeff(np.array(x))
        return all([
            coeff[0] >= 0,
            coeff[1] >= 0,
            coeff.sum() <= 1])

    def weights(self, x):
        '''
        Weights for the triangolation of vector x
        http://math.stackexchange.com/questions/1727200/compute-weight-of-a-point-on-a-3d-triangle

        Parameters
        ----------
        x: array, list
            Point for which weights are calculated

        Returns
        -------
        np.ndarray
            Weights array
        '''
        if not self.is_inside(x):
            raise OutsideTriangleError(
                'x={} is outside triangle\n{}'.format(x, self))
        w = np.cross(self.x0 - self.x1, self.x0 - self.x2)
        # note use the abs value for being sure that all areas are
        # negative
        return np.abs([np.cross(x - x0, x - x1) / w
                       for x0, x1 in itertools.combinations(
            self.__iter__(), 2)
        ])[::-1]

    def __iter__(self):
        for x in [self.x0, self.x1, self.x2]:
            yield x

    def __str__(self):
        s = tabulate.tabulate([x.tolist() for x in self.__iter__()],
                              headers=self.headers)
        return s

    def __repr__(self):
        s = super(Triangle, self).__repr__()
        s += '\n\n'
        s += self.__str__()
        return s


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
        mechanism: str, cantera.Solution
        '''
        super(Polimi, self).__init__(
            proximate_analysis=proximate_analysis,
            ultimate_analysis=ultimate_analysis,
            pressure=pressure,
            name=name)
        self.mechanism = None
        self.backend = None
        self._define_triangle()

    def set_parameters(self, **kwargs):
        '''
        '''
        for key in ('mechanism', 'backend'):
            if key in kwargs:
                setattr(self, key, kwargs[key])

    @property
    def backend(self):
        '''
        Return ODE backend
        '''
        return self._backend

    @backend.setter
    def backend(self, value=None):
        '''
        Set ODE solver backend.

        Parameters
        ----------
        value: str, default='dopri5'
            'dopri5', 'cvode', 'lsoda', 'dop853'
        Raise
        -----
        ValueError
            If backend does not exist
        '''
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

    mechanism = property(_get_mechanism, _set_mechanism)

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

        solver = ode(dmidt).set_integrator(backend, nsteps=1)
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

        t = np.array(t)
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
