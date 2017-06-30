"""
Here are defined the functions in numba
"""
import numba
import numpy as np
import math

from .interpolate import interp
from ._cpd_correlations import xx, yy

Rgas = 1.987  # cal/mol-K


@numba.jit(nopython=True)
def x_n_calc_i(x, z_n, k_n_1):
    """Calc xn_i from binomial distribution with numba."""
    return z_n / (1 + k_n_1 * x)


@numba.jit(nopython=True)
def x_n_calc(x, z_n, k_n_1):
    """Calc xn from binomial distribution with numba."""
    x_n = np.empty_like(z_n)
    for i in range(len(z_n)):
        x_n[i] = x_n_calc_i(x, z_n[i], k_n_1[i])
        # x_n[i] = z_n[i] / (1 + k_n_1[i] * x)
    return x_n


@numba.jit(nopython=True)
def sum_x_n_calc(x, z_n, k_n_1):
    """Calc sum of x_n."""
    _sum = 0
    for i in range(len(z_n)):
        # eq. 54
        x_n_calc = x_n_calc_i(x, z_n[i], k_n_1[i])
        _sum += k_n_1[i] * x_n_calc
    return _sum


@numba.jit(nopython=True)
def fp(x, sigma):
    """Fp for flash distillation."""
    return x * (1 - x)**(sigma - 1)


@numba.jit(nopython=True)
def pstar_f(x, sigma, fpp):
    """P star from percolation theory."""
    return fp(x, sigma) - fpp


@numba.jit(nopython=True)
def combinln(n, k):
    """
    Return the combined ln function.

    Parameters
    ----------
    n:
    k:

    """
    return math.lgamma(n + 1) - (math.lgamma(k + 1) +
                                 math.lgamma(n - k + 1))


@numba.jit(nopython=True)
def binomial(k, n, p):
    """
    Binomial function logarithmic.

    Parameters
    ----------
    k:
    n:
    p:

    """
    bnm = np.empty_like(n, dtype=np.float64)
    logp = math.log(p)
    one_logp = math.log(1 - p)
    for i in range(len(k)):
        bnm[i] = math.exp(combinln(n[i], k[i]) + k[i] *
                          logp + (n[i] - k[i]) * one_logp)
    return bnm


# FUNCTIONS

@numba.jit(nopython=True)
def invernorm(y):
    """
    Calculate the inverse normal distribution function.

    Calculate the inverse of the CDF of the normal distribution.
    It is a wrapper to scipy.stats.norm.ppf, which prevents to use
    values of the cumulative probability too small or too large.

    Parameters
    ----------
    y: float, array
        Cumulative probability

    Return
    ------
    float: inverse of the norm CDF

    """
    if y >= 0.5:
        return interp(y, yy, xx)
    else:
        return - interp(1.0 - y, yy, xx)
