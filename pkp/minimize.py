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

    def error(self, x):
        err = pkp.evolution.error(self, x)[0]
        print('Error=', err)
        return err

    def run(self, initial, **kwargs):
        if isinstance(initial, dict):
            initial = np.array(
                [initial[p]
                 for p in self.empirical_model.parameters_names])
        self.__log.debug('Initial non-scaled: %s', initial)
        # initial[0] = np.log10(initial[0])
        p_min = np.array(self._parameters_min)
        p_max = np.array(self._parameters_max)
        p_min[0] = np.log10(p_min[0])
        p_max[0] = np.log10(p_max[0])
        initial[0] = np.log10(initial[0])
        initial = ((initial - p_min) /
                   (p_max - p_min))
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
        self.__log.debug('Best optimized parameters: %s', best)
        # return {p: v for p, v in zip(
        #    self.empirical_model.parameters_names, best)}
        return res
