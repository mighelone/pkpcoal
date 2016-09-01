from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals

import numpy as np
import pkp
import pkp.empirical_model
from pkp.empirical_model import Rgas


class C2SM2streams(pkp.empirical_model.C2SM):
    '''
    C2SM treating volatile yields as two stream :math:`y_1` and
    :math:`y_2`.
    '''
    y0 = [0, 0, 1]  # y1, y2, s

    def rate(self, t, y):
        RT = Rgas * self.T(t)
        k1, k2 = self._K(RT)

        dsdt = - (k1 + k2) * y[2]
        dy1dt = - self.parameters['y1'] * k1 * y[2]
        dy2dt = - self.parameters['y2'] * k2 * y[2]

        return np.array([dy1dt, dy2dt, dsdt])
