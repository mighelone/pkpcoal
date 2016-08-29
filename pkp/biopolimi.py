from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals

import numpy as np
from autologging import logged

import pkp.polimi
from pkp.polimi import TriangleCoal
import cantera
import os


def set_reference_coal(name, c, h, o):
    return pkp.polimi.set_reference_coal(name,
                                         {'C': c,
                                          'H': h,
                                          'O': o})


def set_reference_biomass(name, comp):
    biomass.TPY = 300, 101325, comp
    return pkp.polimi.set_reference_coal(
        name,
        atoms={el: biomass.elemental_mass_fraction(el)
               for el in ('C', 'H', 'O')})


biomass_xml = os.path.join(os.path.dirname(pkp.bins.__file__),
                           'Biomass.xml')
biomass = cantera.Solution(biomass_xml)

bio1 = set_reference_coal('CELL', c=6, h=10, o=5)
bio2 = set_reference_coal('HCE', c=5, h=8, o=4)
bio3 = set_reference_coal(name='LIGC', c=15, h=14, o=4)
bio4 = set_reference_coal(name='LIGO', c=20, h=22, o=10)
bio5 = set_reference_coal(name='LIGH', c=22, h=28, o=9)


bioS1_comp = {'CELL': 0.6, 'HCE': 0.4}
bioS1 = set_reference_biomass(name='S1', comp=bioS1_comp)

bioS2_comp = {'LIGO': 0.8, 'LIGC': 0.2}
bioS2 = set_reference_biomass(name='S2', comp=bioS2_comp)

bioS3_comp = {'LIGH': 0.8, 'LIGC': 0.2}
bioS3 = set_reference_biomass(name='S3', comp=bioS3_comp)

bio_comp_species = list(set(bioS1_comp.keys() +
                            bioS2_comp.keys() +
                            bioS3_comp.keys()))

triangle_123 = TriangleCoal(bioS1, bioS2, bioS3)


@logged
class BioPolimi(pkp.polimi.Polimi):
    '''
    Biomass Polimi Multiple Step Kinetic Model
    '''
    light_gas = ['HAA',
                 'HMFU',
                 'LVG',
                 'XYLOSE',
                 'Glyoxal',
                 'Phenol',
                 'pCoumaryl',
                 'C11H12O4',
                 'C3H6O2',
                 'C3H4O2',
                 'ALD3',
                 'MECHO',
                 'C2H6OH',
                 'C2H4',
                 'CH3OH',
                 'CH4',
                 'CO2',
                 'CO',
                 'H2O',
                 'H2',
                 'GCO2',
                 'GCO',
                 'GCOH2',
                 'GH2',
                 'EtOH',
                 'HCOOH']
    raw = bio_comp_species
    char = ['Char']
    metaplast = ['CELLA',
                 'HCE',
                 'HCE1',
                 'HCE2',
                 'ACQUA',
                 'LIGOH',
                 'LIGCC',
                 'LIG',
                 'GCO2',
                 'GCO',
                 'GCOH2',
                 'GH2',
                 'GCH4',
                 'GCH3OH',
                 'GC2H4']
    tar = ['GH2']

    def _define_triangle(self):
        if triangle_123.is_inside(self):
            self.triangle = triangle_123
            w = self.triangle.weights(self)
            self.composition = {
                sp: (
                    np.dot(w,
                           [b.get(sp, 0)
                            for b in [bioS1_comp, bioS2_comp,
                                      bioS3_comp]]))
                for sp in bio_comp_species}
            self.weights = w
        else:
            raise pkp.polimi.CompositionError(
                'Biomass composition outside definition triangle')

    def _get_mechanism(self):
        return super(BioPolimi, self)._get_mechanism()

    def _set_mechanism(self, value=None):
        if value is None:
            value = biomass_xml
        super(BioPolimi, self)._set_mechanism(value)

    mechanism = property(_get_mechanism, _set_mechanism)
