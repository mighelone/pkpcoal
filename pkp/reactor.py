'''
Define the :class:`pkp.reactor.Reactor`, which is used for prescribing
the operating conditions.
'''

from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals

import numpy as np
from scipy.interpolate import interp1d
from autologging import logged


@logged
class Reactor(object):
    '''
    Base class for running devolatilization simulations
    '''

    def __init__(self):
        super(Reactor, self).__init__()
        self.operating_conditions = None

    @property
    def operating_conditions(self):
        return self._operating_conditions

    @operating_conditions.setter
    def operating_conditions(self, conditions):
        '''
        Define operating conditions for evaluating pyrolysis

        Parameters
        ----------
        conditions: np.ndarray, list
            [[t0, T0], ..., [tn, Tn]]
        '''
        if conditions is None:
            self.T = None
            self._operating_conditions = None
            return
        if not isinstance(conditions, (np.ndarray, list)):
            raise TypeError('Define conditions as list or numpy array')
        elif isinstance(conditions, list):
            conditions = np.array(conditions)
        if not conditions.ndim == 2:
            raise ValueError('Define conditions as array Nx2')
        if not conditions.shape[-1] == 2:
            raise ValueError('Define conditions as array Nx2')
        self._operating_conditions = conditions
        self.T = interp1d(conditions[:, 0], conditions[:, 1],
                          kind='linear')
