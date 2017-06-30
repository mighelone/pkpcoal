"""Triangle module."""
from autologging import logged
import numpy as np
import itertools
import tabulate


class OutsideTriangleError(Exception):
    """
    Raise an exception if the coal used is outside of the triangles
    defined.
    """

    pass


@logged
class Triangle(object):
    """
    Triangle class. Used for triangulation calculations
    """

    headers = ['x', 'y']

    def __init__(self, x0=None, x1=None, x2=None):
        """
        Parameters
        ----------
        x0, x1, x2, np.ndarray, list
            2D Array or list of the triangle vertices
        """
        if x0 is None:
            x0 = np.array([0, 0])
        if x1 is None:
            x1 = np.array([1, 0])
        if x2 is None:
            x2 = np.array([0, 1])
        self.x0 = np.array(x0)
        self.x1 = np.array(x1)
        self.x2 = np.array(x2)

    def _coeff(self, x):
        """
        Calculate coefficient of linear combination of x
        x-x0 = a*(x1-x0)+b*(x2-x0)
        """
        v1 = self.x1 - self.x0
        v2 = self.x2 - self.x0
        v = x - self.x0
        matr = np.transpose(np.array([v1, v2]))
        return np.linalg.solve(matr, v)

    def is_inside(self, x):
        """
        verify is point x is inside the triangle

        Returns
        -------
        bool
        """
        coeff = self._coeff(np.array(x))
        return all([
            coeff[0] >= 0,
            coeff[1] >= 0,
            coeff.sum() <= 1])

    def weights(self, x):
        """
        Weights for the triangolation of vector x
        http://math.stackexchange.com/questions/1727200/compute-weight-of-a-point-on-a-3d-triangle

        Parameters
        ----------
        x: array, list
            Point for which weights are calculated

        Returns
        -------
        np.ndarray
            Weights array
        """
        if not self.is_inside(x):
            raise OutsideTriangleError(
                'x={} is outside triangle\n{}'.format(x, self))
        w = np.cross(self.x0 - self.x1, self.x0 - self.x2)
        # note use the abs value for being sure that all areas are
        # negative
        return np.abs([np.cross(x - x0, x - x1) / w
                       for x0, x1 in itertools.combinations(
            self.__iter__(), 2)
        ])[::-1]

    def __iter__(self):
        for x in [self.x0, self.x1, self.x2]:
            yield x

    def __str__(self):
        s = tabulate.tabulate([x.tolist() for x in self.__iter__()],
                              headers=self.headers)
        return s

    def __repr__(self):
        s = super(Triangle, self).__repr__()
        s += '\n\n'
        s += self.__str__()
        return s

    def plot(self, ax, **kwargs):
        """Plot triangle"""
        xi, yi = ([self.x0[i], self.x1[i], self.x2[i], self.x0[i]]
                  for i in range(2))
        ax.plot(xi, yi, **kwargs)
