'''
Minimization using scipy.
'''
from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals

import scipy.optimize
import pkp.evolution
import numpy as np

from autologging import logged


@logged
class Minimization(pkp.evolution.Evolution):

    def __init__(self):
        self._ntargets = 0
        self.ref_results = {}
        self._empirical_model = pkp.empirical_model.SFOR
        self._parameters_min = None
        self._parameters_max = None
        self._skip = 1

    def error(self, x):
        err = pkp.evolution.error(self, x)[0]
        self.__log.debug('x: %s - err: %s', x, err)
        return err

    def run(self, initial):
        '''
        Parameters
        ----------
        initial: dict, list
            Initial solution (non-scaled)
        '''
        self.__log.debug('Initial non-scaled: %s', initial)
        initial = self.empirical_model.scale_parameters(
            initial,
            self._parameters_min,
            self._parameters_max)
        self.__log.debug('Initial scaled: %s', initial)
        res = scipy.optimize.minimize(
            fun=self.error,
            x0=initial,
            args=(),
            method='BFGS',
            options={'disp': True}
            # maxiter=1000,
            # disp=True,
            # tol=1e-6)
        )
        best = self.unscale_parameters_final(res.x)
        self.results = res
        self.__log.debug('Best optimized parameters: %s', best)
        return {p: v for p, v in
                zip(self.empirical_model.parameters_names(), best)}
