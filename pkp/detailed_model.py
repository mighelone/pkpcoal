'''
Module contains the Detailed model base class. This class is used as
parent for detailed model classes, such :class:`pkp.cpd.CPD`,
:class:`pkp.polimi.Polimi` and :class:`pkp.biopolimi.BioPolimi`
'''
from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals
from six import string_types
from builtins import dict

import os
import numpy as np
import tabulate
import pkp
import pkp.reactor
import pkp.bins
import json
from autologging import logged
from distutils.dir_util import mkpath

from ._exceptions import PKPCompositionError, PKPConvertNumber

pa_keys = ['FC', 'VM', 'Ash', 'Moist']
pa_keys_daf = pa_keys[: 2]
ua_keys = ['C', 'H', 'O', 'N', 'S']

# calc M elements
M_elements = {'C': 12.010999999999999,
              'H': 1.0079400000000001,
              'N': 14.006740000000001,
              'O': 15.9994,
              'S': 32.065}

# heating value char
T_ref = 273

# entalpy of formation of main gas species
hf = {'CO2': -394427320.2748184,
      'CO': -111261940.94031826,
      'H2O': -242667933.54008716,
      'O2': -737319.0856519284,
      'SO2': -296840.0,
      'N2': -728930.9826824698,
      'char': -101.268}


with open(os.path.join(os.path.dirname(pkp.bins.__file__),
                       'el_fractions.json'), 'r') as f:
    el_fractions = json.load(f)

with open(os.path.join(os.path.dirname(pkp.bins.__file__),
                       'species.json'), 'r') as f:
    species = json.load(f)

lhv_char = (hf['char'] + hf['O2'] - hf['CO2']) / M_elements['C']

M_H2O = 2 * M_elements['H'] + M_elements['O']

# Latent Heat of Water in J/kg :
rH2O = 2263073


def normalize_dictionary(d):
    '''
    Normalize dictionary d and return a new one
    '''
    sum_d = sum(d.values())
    return {el: (val / sum_d) for el, val in d.items()}


@logged
class DetailedModel(pkp.reactor.Reactor):
    '''
    Detailed model class used as parent class for Devolatilization
    models
    '''

    def __init__(self, proximate_analysis=None, ultimate_analysis=None,
                 pressure=101325, hhv=None, name='Detailed model'):
        '''
        Parameters
        ----------
        proximate_analysis: dict
            Proximate analysis dict i.e:
            `{'FC': 45.1, 'VM': 50.6, 'Ash': 4.3, 'Moist': 19.0}`
        ultimate_analysis: dict
            Ultimate analysis dictionary i.e:
            `{'C': 80, 'H': 8, 'O': 12, 'N': 0, 'S': 0}`
        pressure: float
            Pressure of pyrolysis process in Pa
        name: str, unicode
            Reference name of the modelled coal
        '''
        super(DetailedModel, self).__init__()
        if ultimate_analysis is None:
            ultimate_analysis = {'C': 80, 'H': 8, 'O': 12, 'N': 0, 'S': 0}
        if proximate_analysis is None:
            proximate_analysis = {'FC': 45.1,
                                  'VM': 50.6, 'Ash': 4.3, 'Moist': 19.0}
        self.ultimate_analysis = ultimate_analysis
        self.proximate_analysis = proximate_analysis
        self.pressure = pressure
        self.__log.debug('Set pressure %s', self.pressure)
        self.name = name

        self.hhv = hhv

        self._operating_conditions = None
        self.T = None
        self.rho_dry = 1000.0

        self._basename = 'dummy'
        self._path = 'dummy'
        self.basename = None
        self.path = None

    @property
    def hhv(self):
        '''
        Higher heating value of coal as received in J/kg
        '''
        return self._hhv

    @hhv.setter
    def hhv(self, value):
        """Read HHV in MJ/kg"""
        if value:
            try:
                self._hhv = value * 1e6
            except TypeError:
                raise PKPConvertNumber(
                    'Define HHV as number or None: {}'.format(value))
            self._hhv_daf = self._hhv / self.daf
        else:
            self._hhv_daf = self.dulong()
            self._hhv = self.daf * self._hhv_daf

        self._lhv_daf = (
            self._hhv_daf - rH2O * self.ultimate_analysis['H'] *
            0.5 * M_H2O / M_elements['H'])
        self._lhv = (self._lhv_daf * self.daf -
                     self.proximate_analysis['Moist'] * rH2O)

    @property
    def lhv_char(self):
        '''
        Lower heating value of char in J/kg is calculated assuming char
        as graphite
        '''
        return lhv_char

    @property
    def hhv_daf(self):
        '''
        Higher heating value of coal on daf basis in J/kg
        '''
        return self._hhv_daf

    @property
    def lhv(self):
        '''
        Lower heating value of coal as received in J/kg
        '''
        return self._lhv

    @property
    def lhv_daf(self):
        '''
        Lower heating value of coal on daf basis in J/kg
        '''
        return self._lhv_daf

    def dulong(self):
        '''
        Calculate HHV_daf using the Dulong formula
        http://www.bti-europe.eu/downloads/CoalConversionFactsCalculations.pdf
        '''
        coeff = {
            'C': 33.3,
            'H': 144.2,
            'O': -18.025,
            'N': 0,
            'S': 9.3
        }
        return sum(c * self.ultimate_analysis[el]
                   for el, c in coeff.items()) * 1e6

    def postulate_species(self, y0, mw=200.0, include_nu=False):
        '''
        Calculate the volatile composition for the fitted empirical
        model assuming a unique *postulate* species. The composition is
        calculated assuming that char is composed only by carbon and
        the remaining carbon and other elements goes to the postulate
        species.

        Parameters
        ----------
        y0: float
            Final volatile yield
        mw: float
            Postulate species molecular weight
        include_nu: Bool
            Include stoichiometric coefficient in the output

        Returns
        -------
        dict:
            Dictionary containing the composition of the postulate
            species with molecular weight, enthalpy of formation and
            other information.

        See also
        --------
        :meth:`empirical_composition`
        '''
        assert 0 < y0 < 1, 'Define y0 between 0 and 1'
        ua_char = {el: (1 if el == 'C' else 0)
                   for el in self.ultimate_analysis}
        molecule = {el: ((val - (1 - y0) * ua_char[el]) * mw /
                         M_elements[el] / y0)
                    for el, val in self.ultimate_analysis.items()}
        molecule_name = ''.join('{}_{:4.3f} '.format(el, molecule[el])
                                for el in ['C', 'H', 'O', 'N', 'S'])

        lhv_vol = (self.lhv_daf - (1 - y0) * self.lhv_char) / y0

        nu = {
            'CO2': molecule['C'],
            'H2O': (molecule['H'] * 0.5),
            # reactant is negative
            'O2': -0.5 * (-molecule['O'] + 2 * molecule['C'] +
                          0.5 * molecule['H'] + 2 * molecule['S']),
            'SO2': molecule['S'],
            'N2': (molecule['N'] * 0.5)
        }

        nu_partial = nu.copy()
        nu_partial['CO'] = nu_partial.pop('CO2')
        nu_partial['O2'] = -0.5 * (-molecule['O'] + molecule['C'] +
                                   0.5 * molecule['H'] +
                                   2 * molecule['S'])

        full_oxidation = ('{} + {:4.3f} O2 -> '
                          '{:4.3f} CO2 + {:4.3f} H2O + '
                          '{:4.3f} N2 + {:4.3f} SO2'.format(
                              molecule_name,
                              -nu['O2'],
                              nu['CO2'],
                              nu['H2O'],
                              nu['N2'],
                              nu['SO2']))
        part_oxidation = ('{} + {:4.3f} O2 -> '
                          '{:4.3f} CO + {:4.3f} H2O + {:4.3f} N2'
                          ' + {:4.3f} SO2'.format(
                              molecule_name,
                              -nu_partial['O2'],
                              nu_partial['CO'],
                              nu_partial['H2O'],
                              nu_partial['N2'],
                              nu_partial['SO2']))

        hf_vol = np.sum(n * hf[el]
                        for el, n in nu.items()) + lhv_vol * mw

        postulate_dict = {
            'name': molecule_name,
            'formula': molecule,
            'molecular_weight': mw,
            'hf': hf_vol,
            'y0': y0,
            'reaction0': full_oxidation,
            'reaction1': part_oxidation
        }
        if include_nu:
            postulate_dict['nu'] = nu
            postulate_dict['nu_partial'] = nu_partial

        return postulate_dict

    @staticmethod
    def el_fraction(sp, el):
        '''Mass fraction of element el in species sp'''
        return el_fractions[sp][el]

    def empirical_composition(self, y0, tar, CO):
        '''
        Set the empirical composition of volatiles using the method by
        [Petersen2005]_
        This method is based on the element conservation balance and
        some assumption regarding tar and CO.

        The model assumes that C6H6 (tar), CO, CO2, N2, H2 

        Parameters
        ----------
        y0: float
            Final volatile yield
        tar: mass fraction of tar in volatiles
        CO: mass fraction of CO in volatiles

        Returns
        -------
        emp_dict: dict
            Empirical composition dictionary

        See also
        --------
        :meth:`postulate_species`
        '''

        def calc_remaining(comp):
            '''Remaining fraction of each elements'''
            return {el: (ua - tot_el_fraction(comp, el))
                    for el, ua in ultimate_analysis.items()}

        def tot_el_fraction(comp, element):
            '''Calc the total element fraction of the given element'''
            return np.sum([val * self.el_fraction(sp, element)
                           for sp, val in comp.items()])

        sum_ua = (sum(self.ultimate_analysis.values()) -
                  self.ultimate_analysis['S'])
        ultimate_analysis = {
            el: v / sum_ua
            for el, v in self.ultimate_analysis.items()
            if el != 'S'}
        self.__log.debug('Update ultimate_analysis %s',
                         ultimate_analysis)
        composition = {}
        composition['char'] = 1 - y0
        # assume tar as C6H6

        # assume N -> N2
        composition['N2'] = ultimate_analysis['N']
        # composition['CO'] = (CO * ultimate_analysis['O'] /
        #                     el_fraction('CO', 'O'))
        O_in_CO = self.el_fraction('CO', 'O')
        if ultimate_analysis['O'] > CO * O_in_CO:
            composition['CO'] = CO
            composition['CO2'] = ((ultimate_analysis['O'] - CO * O_in_CO) /
                                  self.el_fraction('CO2', 'O'))
        else:
            composition['CO'] = ultimate_analysis['O'] / O_in_CO
            composition['CO2'] = 0
        self.__log.debug('CO input %s', CO)
        self.__log.debug('CO set %s', composition['CO'])
        self.__log.debug('CO2 set %s', composition['CO2'])
        # composition['CO2'] = ((1 - CO) * ultimate_analysis['O'] /
        #                      el_fraction('CO2', 'O'))
        self.__log.debug('Vol composition %s', composition)

        remaining = calc_remaining(composition)
        self.__log.debug('Remaining element after CO/CO2/N2: %s',
                         remaining)

        C_in_tar = self.el_fraction('C6H6', 'C')
        self.__log.debug('C in TAR %s', C_in_tar)
        if C_in_tar * tar > remaining['C']:
            composition['C6H6'] = remaining['C'] / C_in_tar
            self.__log.debug('C in Tar > remaining C -> set tar: %s',
                             composition['C6H6'])
        else:
            composition['C6H6'] = tar
            self.__log.debug('Set TAR as %s', tar)

        # recalculate remaining
        remaining = calc_remaining(composition)
        self.__log.debug('Remaining element after tar: %s', remaining)

        c_to_h_mass = remaining['C'] / remaining['H']
        c_to_h_molar = c_to_h_mass * M_elements['H'] / M_elements['C']
        self.__log.debug('C/H molar: %s', c_to_h_molar)

        if 0 <= c_to_h_molar <= 0.5:
            # use C2H4 and H2
            composition['C2H4'] = (remaining['C'] /
                                   self.el_fraction('C2H4', 'C'))
            self.__log.debug('C2H4 %s', composition['C2H4'])
        elif 0.5 < c_to_h_molar < 1:
            # use C6H6
            composition['C6H6'] = (composition['C6H6'] +
                                   remaining['C'] /
                                   self.el_fraction('C6H6', 'C'))
            self.__log.debug('Update C6H6 %s', composition['C6H6'])

        self.__log.debug('Remaining element after C: %s', remaining)
        remaining = calc_remaining(composition)

        composition['H2'] = remaining['H']
        remaining = calc_remaining(composition)
        self.__log.debug('Remaining %s', remaining)
        self.__log.debug('Final Vol composition %s', composition)

        emp_dict = {'composition': composition,
                    'heat_pyro': self.heat_of_pyrolysis(
                        composition)}
        return emp_dict

    # def calc_element_fraction(self, element, species):
    #     '''
    #     Calculate element fraction of the given species.

    #     Parameters
    #     ----------
    #     element: str
    #         Name of element
    #     species: str
    #         Name of species

    #     Returns
    #     -------
    #     el_fraction: float
    #         fraction of element
    #     '''
    #     if species in ('Char', 'Solid'):
    #         if element == 'C':
    #             return 1.0
    #         else:
    #             return 0.0
    #     elif species == 'Other':
    #         return 0.0
    #     else:
    #         i = self.gas.element_index(element)
    #         return (self.gas.atomic_weights[i] *
    #                 self.gas.species(species).composition.get(element, 0) /
    #                 self.gas.molecular_weights[
    #                 self.gas.species_index(species)])

    def heat_of_volatiles(self, composition):
        '''
        Calculate the total heat of the volatile yields including char

        Parameters
        ----------
        composition: dict
            Volatile composition dictionary

        Returns
        -------
        hhv: float
            Heat of reaction of volatile gases mixture, MJ/kg
        '''
        return sum(val * self.heat_of_reaction_species(sp)
                   for sp, val in composition.items())

    def heat_of_pyrolysis(self, composition):
        '''
        Heat of pyrolysis. It is defined as the total heat released
        during pyrolysis per unit of volatiles.
        Positive heat of pyrolysis means that extra heat is removed to
        fullfill the energy balance of coal.

        Same convention used for latent heat of evaporation!!!

        .. math::
            \Delta H_{daf} = \Delta H_{vol} - \Delta H_{pyro}

        See also
        --------
        :meth:`heat_of_volatiles`
        :meth:`heat_of_reaction_species`

        '''
        heat_vol = self.heat_of_volatiles(composition)
        dh_pyro = self.lhv_daf - heat_vol
        return -dh_pyro / (1 - composition['char'])

    def heat_of_reaction_species(self, sp):
        '''
        Calculate the heat of reaction :math:`\Delta H` for the species
        sp:

        .. math::

            \sum_r \mu_r h_{f, r} = \sum_p \mu_p h_{f, p} + \Delta H

        where :math:`\mu` are the stoichiometric coefficient
        of the reactions, :math:`h_f` the enthalpy of formation
        and :math:`r` and :math:`p` are the reactants and products,
        respectively.

        Parameters
        ----------
        sp: str
            Species

        Returns
        -------
        Heat of reaction of the species, J/kg

        :meth:`heat_of_pyrolysis`
        :meth:`heat_of_reaction_species`
        '''
        sp_coeff = species[sp]
        if sp in ('N2'):
            return 0.0
        elif sp == 'char':
            return self.lhv_char
        hf = sp_coeff['hf']
        n_o2 = 0.5 * (
            2 * sp_coeff.get('C', 0) +
            0.5 * sp_coeff.get('H', 0) -
            sp_coeff.get('O', 0))
        n_co2 = sp_coeff.get('C', 0)
        # print 'n_co2', n_co2
        n_h2o = sp_coeff.get('H', 0) * 0.5
        # print 'n_h2o', n_h2o
        mw = sp_coeff['mw']

        # h_o2 = gas.species('O2').thermo.h(T_ref)
        # h_o2 = hf['O2']
        h_o2 = species['O2']['hf']
        # h_co2 = gas.species('CO2').thermo.h(T_ref)
        # h_co2 = hf['CO2']
        h_co2 = species['CO2']['hf']
        # h_h2o = gas.species('H2O').thermo.h(T_ref)
        # h_h2o = hf['H2O']
        h_h2o = species['H2O']['hf']

        heat_molar = (hf + n_o2 * h_o2 - n_co2 *
                      h_co2 - n_h2o * h_h2o)
        return heat_molar / mw

    @property
    def name(self):
        '''
        Reference name for the coal.
        '''
        return self._name

    @name.setter
    def name(self, value):
        if isinstance(value, string_types):
            self._name = value
        else:
            raise TypeError('Coal name should be a string')

    @property
    def ultimate_analysis(self):
        '''
        Ultimate analysis of coal in daf basis. It is defined as
        dictionary: {'C': 0.8, 'H':0.12, 'O': 0.05, 'N':0.2, 'S':0}.
        '''
        return self._ultimate_analysis

    @ultimate_analysis.setter
    def ultimate_analysis(self, ultimate_analysis):
        if not all((key in ultimate_analysis for key in ua_keys)):
            raise PKPCompositionError(
                'Ultimate analysis keys should be {}'.format(ua_keys))
        try:
            self._ultimate_analysis = normalize_dictionary(
                ultimate_analysis)
        except TypeError as e:
            raise PKPConvertNumber(
                'Error reading ultimate_analysis\n{}'.format(
                    ultimate_analysis))

    @property
    def rho_dry(self):
        '''
        Apparent density of dry coal in kg/m3.
        '''
        return self._rho_dry

    @rho_dry.setter
    def rho_dry(self, value):
        if isinstance(value, (float, int)):
            self._rho_dry = value
        else:
            raise TypeError('Define rho_dry as number')

    @property
    def proximate_analysis(self):
        '''
        Coal proximate analysis dictionary. It is defined as follows:
        `{'FC': 0.4, 'VM':0.4, 'Ash':0.05, 'Moist':0.15}`
        '''
        return self._proximate_analysis

    @proximate_analysis.setter
    def proximate_analysis(self, proximate_analysis):
        if not all((key in proximate_analysis for key in pa_keys)):
            raise PKPCompositionError(
                'Proximate analysis keys should be {}'.format(pa_keys))
        try:
            self._proximate_analysis = normalize_dictionary(
                proximate_analysis)
        except TypeError as e:
            raise PKPConvertNumber(
                'Error reading proximate_analysis\n{}'.format(
                    proximate_analysis))

        self._daf = sum((self._proximate_analysis[key]
                         for key in pa_keys_daf))
        self._proximate_analysis_daf = {
            key: (self.proximate_analysis[key] / self._daf)
            for key in pa_keys_daf}

    @property
    def proximate_analysis_daf(self):
        '''
        Coal proximate analysis on daf basis. It is defined as follows:
        `{'FC': 0.5, 'VM':0.5}`
        '''
        return self._proximate_analysis_daf

    @property
    def daf(self):
        '''
        Dry ash free fraction on coal as received. It correspond to the
        mass fractions of volatile and char.
        '''
        return self._daf

    @property
    def pressure(self):
        '''
        Operating pressure in Pa.
        '''
        return self._pressure

    @pressure.setter
    def pressure(self, value):
        self._pressure = value

    @property
    def path(self):
        '''
        Path where results are saved.
        '''
        return self._path

    @path.setter
    def path(self, value):
        if value is None:
            self._path = os.getcwd()
        else:
            self._path = os.path.abspath(value)
        self.__log.debug('Set path to %s', self._path)
        if not os.path.isdir(self._path):
            self.__log.debug('Create path %s', self._path)
            # os.mkdir(self._path)
            mkpath(self._path)
        # if you uodate the path update also the files
        self._set_basename(self._basename)

    @property
    def van_kravelen(self):
        '''
        Coordinate of van kravelen diagram, defined as (O/C, H/C) as
        received on molar basis.
        '''
        mol = {el: (self.ultimate_analysis[el] / M_elements[el])
               for el in ['C', 'H', 'O']}
        return np.array([mol['O'] / mol['C'], mol['H'] / mol['C']])

    # use explicit property instead of a decorator to redefine in
    # children classes
    def _get_basename(self):
        return self._basename

    def _set_basename(self, value):
        '''

        '''
        if value is None:
            value = (self.__class__.__name__ +
                     '_' + self.name.replace(' ', '_'))
        self._basename = value
        self.__log.debug('basename: %s', self.basename)
        self._out_csv = os.path.join(self.path, self.basename + '.csv')
        self.__log.debug('Out CSV %s', self._out_csv)

    basename = property(_get_basename, _set_basename,
                        doc='Basename for exporting results')

    def __str__(self):
        str = ('Coal: {}\n'.format(self.name))
        str += ''.join(['='] * (len(self.name) + 6))
        str += '\n\nUltimate Analysis\n'
        str += tabulate.tabulate(
            [[el, val]
             for el, val in self.ultimate_analysis.items()])
        str += '\n\nProximate Analysis\n'
        str += tabulate.tabulate(
            [[el, val]
             for el, val in self.proximate_analysis.items()])
        str += '\n'
        return str

    def __repr__(self):
        return self.__str__()

    # common interface for children classes
    def set_parameters(self, **kwargs):
        '''
        Set the calculation parameters
        '''
        pass

    def get_parameters(self):
        '''
        Return a dictionary parameters
        '''
        return {}

    def run(self, **kwargs):
        pass
