import sys
if sys.version_info >= (3, 6):
    NumbaError = ModuleNotFoundError
else:
    NumbaError = ImportError

try:
    from numba import jit

    @jit("f8(f8, f8[:], f8[:])", nopython=True)
    def interp(xi, x, y):
        """
        Linear interpolation

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
        #print('i_min={}, i_max={}'.format(i_min, i_max))
        while search:
            # print('j={} i_min={}, i_max={}'.format(j, i_min, i_max))
            # i = int(i_min + ((i_max - i_min) / (x[i_max] - x[i_min])) * (xi - x[i_min]))
            i = i_min + (i_max - i_min) // 2
            xs = x[i]
            # print ('i={} xs={}'.format(i, xs))
            if xi >= xs:
                i_min = i
                #print('set i_min={} xs={} xi={}'.format(i_min, xs, xi))
            else:
                i_max = i
                #print('set i_max={} xs={} xi={}'.format(i_max, xs, xi))

            if i_max - i_min == 1:
                search = False

            j += 1
            if j == 50:
                search = False

            # print('')

        return y[i_min] + (y[i_max] - y[i_min]) * (xi - x[i_min]) / (x[i_max] - x[i_min])

except NumbaError:
    import numpy as np

    def interp(xi, x, y):
        return np.interp(xi, x, y, left=y[0], right=y[-1])
