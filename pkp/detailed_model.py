'''
Module contains the Detailed model base class. This class is used as
parent for detailed model classes, such :class:`pkp.cpd.CPD`,
:class:`pkp.polimi.Polimi` and :class:`pkp.biopolimi.BioPolimi`
'''
from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals

import os
import numpy as np
import tabulate
import pkp
import pkp.reactor
from autologging import logged
from distutils.dir_util import mkpath

pa_keys = ['FC', 'VM', 'Ash', 'Moist']
pa_keys_daf = pa_keys[: 2]
ua_keys = ['C', 'H', 'O', 'N', 'S']
M_elements = {'C': 12.0, 'H': 1, 'O': 16.0, 'N': 28, 'S': 32}


def normalize_dictionary(d):
    '''
    Normalize dictionary d and return a new one
    '''
    sum_d = sum(d.itervalues())
    return {el: (val / sum_d) for el, val in d.iteritems()}


@logged
class DetailedModel(pkp.reactor.Reactor):
    '''
    Detailed model class used as parent class for Devolatilization
    models
    '''

    def __init__(self, proximate_analysis, ultimate_analysis,
                 pressure=101325, name='Detailed model'):
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
            Pressure of pyrolysis process
        name: str, unicode
            Reference name of the modelled coal
        '''
        super(DetailedModel, self).__init__()
        self.ultimate_analysis = ultimate_analysis
        self.proximate_analysis = proximate_analysis
        self.pressure = pressure
        self.name = name

        self._operating_conditions = None
        self.T = None
        self.rho_dry = 1000.0

        self._basename = 'dummy'
        self._path = 'dummy'
        self.basename = None
        self.path = None

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        if isinstance(value, (str, unicode)):
            self._name = value
        else:
            raise TypeError('Coal name should be a string')

    @property
    def ultimate_analysis(self):
        return self._ultimate_analysis

    @ultimate_analysis.setter
    def ultimate_analysis(self, ultimate_analysis):
        if not all((key in ultimate_analysis for key in ua_keys)):
            raise ValueError(
                'Ultimate analysis keys should be {}'.format(ua_keys))
        self._ultimate_analysis = normalize_dictionary(
            ultimate_analysis)

    @property
    def rho_dry(self):
        return self._rho_dry

    @rho_dry.setter
    def rho_dry(self, value):
        if isinstance(value, (float, int)):
            self._rho_dry = value
        else:
            raise TypeError('Define rho_dry as number')

    @property
    def proximate_analysis(self):
        return self._proximate_analysis

    @proximate_analysis.setter
    def proximate_analysis(self, proximate_analysis):
        if not all((key in proximate_analysis for key in pa_keys)):
            raise ValueError(
                'Proximate analysis keys should be {}'.format(pa_keys))
        self._proximate_analysis = normalize_dictionary(
            proximate_analysis)
        self._daf = sum((self._proximate_analysis[key]
                         for key in pa_keys_daf))
        self._proximate_analysis_daf = {
            key: (self.proximate_analysis[key] / self._daf)
            for key in pa_keys_daf}

    @property
    def proximate_analysis_daf(self):
        return self._proximate_analysis_daf

    @property
    def daf(self):
        return self._daf

    @property
    def pressure(self):
        return self._pressure

    @pressure.setter
    def pressure(self, value):
        self._pressure = value

    @property
    def path(self):
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
        Return coordinate of van kravelen diagram
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
        Define file base name for CPD results
        '''
        if value is None:
            value = (self.__class__.__name__ +
                     '_' + self.name.replace(' ', '_'))
        self._basename = value
        self.__log.debug('basename: %s', self.basename)
        self._out_csv = os.path.join(self.path, self.basename + '.csv')
        self.__log.debug('Out CSV %s', self._out_csv)

    basename = property(_get_basename, _set_basename)

    def __str__(self):
        str = ('Coal: {}\n'.format(self.name))
        str += ''.join(['='] * (len(self.name) + 6))
        str += '\n\nUltimate Analysis\n'
        str += tabulate.tabulate(
            [[el, val]
             for el, val in self.ultimate_analysis.iteritems()])
        str += '\n\nProximate Analysis\n'
        str += tabulate.tabulate(
            [[el, val]
             for el, val in self.proximate_analysis.iteritems()])
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
