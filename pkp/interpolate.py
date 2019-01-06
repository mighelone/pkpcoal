"""Interpolation with numba."""
import sys

if sys.version_info >= (3, 6):
    NumbaError = ModuleNotFoundError
else:
    NumbaError = ImportError

try:
    from numba import jit

    # @jit("f8(f8, f8[:], f8[:])", nopython=True)
    @jit(nopython=True)
    def interp(xi, x, y):
        """
        Linear interpolation.

        Parameters
        ----------
        xi: float
            Value to interpolate
        x: array
            Grid of x points
        y: array
            Grid of y points

        Return
        ------
        floar: interpolated point

        """
        search = True
        x_min, x_max = x[0], x[-1]
        i_min, i_max = 0, len(x) - 1
        j = 0
        if xi > x_max:
            return y[i_max]
        elif xi < x_min:
            return y[i_min]
        while search:
            i = i_min + (i_max - i_min) // 2
            xs = x[i]
            if xi >= xs:
                i_min = i
            else:
                i_max = i

            if i_max - i_min == 1:
                search = False

            j += 1
            if j == 50:
                search = False

        return y[i_min] + (y[i_max] - y[i_min]) * (xi - x[i_min]) / (
            x[i_max] - x[i_min]
        )


except NumbaError:
    import numpy as np

    def interp(xi, x, y):
        """Interpolate with numpy."""
        return np.interp(xi, x, y, left=y[0], right=y[-1])
