'''
Empirical Model module
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

import pkp.detailed_model


class EmpiricalModel(pkp.detailed_model.Reactor):
    '''
    Parent class for model
    '''
    parameters_names = ['foo', 'bar']
    parameters_default = [1, 1]

    def __init__(self):
        self.parameters = None
        self.y0 = 1.0

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


class SFOR(EmpiricalModel):
    '''
    Single First Order Reaction (SFOR) model
    '''
    parameters_names = ['A', 'E', 'y0']
    parameters_default = [1e5, 50e6, 0.6]

    def __init__(self):
        self._parameters = {'A': 1e5, 'E': 50e6, 'y0': 0.6}

    def rate(self, t, y):
        '''
        SFOR reaction rate, 1/s
        '''
        k = (self.parameters['A'] /
             np.exp(self.parameters['E'] / 8314.33 / self.T(t)))
        return k * (y + self.parameters['y0'] - 1)
