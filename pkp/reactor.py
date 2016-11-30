'''
Define the :class:`pkp.reactor.Reactor`, which is used for prescribing
the operating conditions.
'''

from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals

import numpy as np
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
        '''
        Operating conditions for devolatilization. They are defined as
        list of operating points [[t0, T0], [t1, T1], ..., [tn, Tn]]
        Each operating point is defined by the time in second and
        temperature in K.
        '''
        return self._operating_conditions

    @operating_conditions.setter
    def operating_conditions(self, conditions):
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

        def interp_tT(t):
            '''Interpolate time with temperature'''
            return np.interp(t, conditions[:, 0], conditions[:, 1])
        self.T = interp_tT
