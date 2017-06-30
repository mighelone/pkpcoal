"""
Implementing binomial pdf with gammaln.

see: Robert Kerns discussion in
http://groups.google.ca/group/comp.lang.python/browse_thread/thread/839009574397dc37

"""
import numpy as np
from scipy.misc import comb
from scipy import special
from ._exceptions import ImportError


def bpmf(k, n, p):
    """
    Binomial distribution using comb function in scipy.

    Parameters
    ----------

    """
    # this does not work for large n
    return comb(n, k) * (p**k) * ((1 - p)**(n - k))


try:
    import numba
    import math

    @numba.jit(nopython=True)
    def combinln(n, k):
        """Return combinln function with numba."""
        return math.lgamma(n + 1) - (math.lgamma(k + 1) +
                                     math.lgamma(n - k + 1))


    @numba.jit(nopython=True)
    def bpmfln(k, n, p):
        """Return bpmfln function with numba."""
        bnm = np.empty_like(n, dtype=np.float64)
        logp = math.log(p)
        one_logp = math.log(1 - p)
        for i in range(len(k)):
            bnm[i] = math.exp(combinln(n[i], k[i]) + k[i] *
                              logp + (n[i] - k[i]) * one_logp)
        return bnm

except ImportError:
    # proposed version using gammaln
    def combinln(n, k):
        """Return combinln function with numpy."""
        return (special.gammaln(n + 1) - (special.gammaln(k + 1) +
                                          special.gammaln(n - k + 1)))

    def bpmfln(k, n, p):
        """Return bpmfln function wuth numpy."""
        return np.exp(combinln(n, k) + k * np.log(p) + (n - k) * np.log(1 - p))


if __name__ == '__main__':

    n = 10
    p = 1.0e-5

    print("using gammaln")
    print(bpmfln(np.arange(11), 10, p))
    print("using comb")
    print(bpmf(np.arange(11), 10, p))

    # print(bpmfln(np.arange(11), 10, p) - bpmf(np.arange(11), 10, p))
    assert np.allclose(bpmfln(np.arange(11), 10, p),
                       bpmf(np.arange(11), 10, p))

    pmfnln = bpmfln(np.arange(5001), 5000, 0.99)
    print('n = 5000')
    print('nans', np.sum(np.isnan(pmfnln)))
    print('sum', np.sum(pmfnln))
    print('sum (repr)', repr(np.sum(pmfnln)))
    assert np.isclose(np.sum((pmfnln), np.sum(pmfnln)))
