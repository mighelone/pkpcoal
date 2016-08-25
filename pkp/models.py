'''
Models module
=============

Define the empirical models:

* SFOR
* C2SM
* DAEM
* ...
'''
from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals

import numpy as np


class Model(object):
    '''
    Parent class for model
    '''
    parameters_names = ['foo', 'bar']
    parameters_default = [1, 1]

    def __init__(self):
        self.parameters = None

    @property
    def len_parameters(self):
        return len(self.parameters_names)

    @property
    def parameters(self):
        return self._parameters

    @parameters.setter
    def parameters(self, values):
        if values is None:
            par_values = self.parameters_default
        elif isinstance(values, (list, np.ndarray)):
            par_values = values
        elif isinstance(values, dict):
            if all(key in self.parameters_names for key in values):
                self._parameters = values
                return
            else:
                raise ValueError(
                    'Parameters should contains the'
                    ' following keys {}'.format(self.parameters_names))

        if len(values) > self.len_parameters:
            raise ValueError(
                'Length of parameters list/array should be {}'.format(
                    self.len_parameters))
        else:
            self._parameters = {k: v for k, v in
                                zip(self.parameters_names,
                                    par_values)}

    @parameters
    def operating_conditions(self):


class SFOR(Model):
    '''
    Single First Order Reaction (SFOR) model
    '''
    parameters_names = ['A', 'E', 'Y0']
    parameters_default = [1e5, 50e6, 0.6]

    def __init__(self):
        self._parameters = {'A': 1e5, 'E': 50e6, 'y0': 0.6}
