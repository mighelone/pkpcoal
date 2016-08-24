from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals

import pkp.coalnew
import pkp.coalPolimi
import cantera
import numpy as np

from pkp.coalnew import M_elements


def set_reference_coal(name, atoms):
    atoms['N'] = 0
    atoms['S'] = 0
    ua = {el: (val * M_elements[el])
          for el, val in atoms.iteritems()}
    return pkp.coalnew.Coal(name=name, ultimate_analysis=ua,
                            proximate_analysis={'FC': 50,
                                                'VM': 50,
                                                'Ash': 0,
                                                'Moist': 0})

coal1 = set_reference_coal('Coal1', atoms={'C': 12, 'H': 11, 'O': 0})
coal2 = set_reference_coal('Coal2', atoms={'C': 14, 'H': 10, 'O': 1})
coal3 = set_reference_coal('Coal3', atoms={'C': 12, 'H': 12, 'O': 5})
char = set_reference_coal('Char', atoms={'C': 1, 'H': 0, 'O': 0})


class MechanismError(Exception):
    pass


class CompositionError(Exception):
    pass


class Triangle(object):

    def __init__(self, x0=None, x1=None, x2=None):
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
        calculate coefficient of linear combination of x
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


class TriangleCoal(Triangle):

    def __init__(self, coal0, coal1, coal2):
        self.coal0 = coal0
        self.coal1 = coal1
        self.coal2 = coal2
        super(TriangleCoal, self).__init__(
            x0=coal0.van_kravelen,
            x1=coal1.van_kravelen,
            x2=coal2.van_kravelen)

    def is_inside(self, coal):
        super(TriangleCoal, self).is_inside(coal.van_krevelen)


triangle_012 = TriangleCoal(x0=char,
                            x1=coal1,
                            x2=coal2)

triangle_023 = TriangleCoal(x0=char,
                            x1=coal2,
                            x2=coal3)

triangle_123 = TriangleCoal(x0=coal1,
                            x1=coal2,
                            x2=coal3)


class Polimi(pkp.coalnew.Coal):

    def __init__(self, proximate_analysis, ultimate_analysis,
                 pressure=101325, name='Coal', mechanism='COAL.xml'):
        super(Polimi, self).__init__(
            proximate_analysis=proximate_analysis,
            ultimate_analysis=ultimate_analysis,
            pressure=pressure,
            name=name)
        self.mechanism = mechanism
        self.set_triangle()

    @property
    def mechanism(self):
        return self._mechanism

    @mechanism.setter
    def mechanism(self, value):
        try:
            self._mechanism = cantera.Solution(value)
        except:
            raise MechanismError('Cannot read {}'.format(value))

    def set_triangle(self):
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

    def triangolate(self):
        pass
