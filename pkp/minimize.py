'''
Minimization using scipy.
'''
from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals

import scipy.optimize
import pkp.evolution

from autologging import logged


@logged
class Minimization(pkp.evolution.Evolution):

    def __init__(self, initial):
        self.initial = initial

    def error(self, x):
        return pkp.evolution.error(self, x)[0]

    def run(self, **kwargs):
        res = scipy.optimize.minimize(
            fun=self.error,
            args=(),
            method='Nelder-Mead',
            tol=1e-6)
        best = self.unscale_parameters_final(res.x)
        self.__log.debug('Best optimized parameters: %s', best)
        return {p: v for p, v in zip(
            self.empirical_model.parameters_names, best)}
