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
    """
    """
    return z_n / (1 + k_n_1 * x)


@numba.jit(nopython=True)
def x_n_calc(x, z_n, k_n_1):
    """
    """
    x_n = np.empty_like(z_n)
    for i in range(len(z_n)):
        x_n[i] = x_n_calc_i(x, z_n[i], k_n_1[i])
        # x_n[i] = z_n[i] / (1 + k_n_1[i] * x)
    return x_n


@numba.jit(nopython=True)
def sum_x_n_calc(x, z_n, k_n_1):
    """
    """
    _sum = 0
    for i in range(len(z_n)):
        # eq. 54
        x_n_calc = x_n_calc_i(x, z_n[i], k_n_1[i])
        _sum += k_n_1[i] * x_n_calc
    return _sum


@numba.jit(nopython=True)
def fp(x, sigma):
    """
    """
    return x * (1 - x)**(sigma - 1)


@numba.jit(nopython=True)
def pstar_f(x, sigma, fpp):
    """
    """
    return fp(x, sigma) - fpp


@numba.jit(nopython=True)
def combinln(n, k):
    """

    Combined function

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
    Binomial function logarithmic

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
    '''
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
    '''
    # if y > 0.5:
    #    yp = y
    #    fac = 1
    # else:
    #    yp = 1 - y
    #    fac = -1
    #yp, fac = (y, 1.0) if y > 0.5 else (1 - y, -1.0)
    # return fac * np.interp(yp, yy, xx, right=3.4)
    # return fac * interp(yp, yy, xx)
    if y >= 0.5:
        return interp(y, yy, xx)
    else:
        return - interp(1.0 - y, yy, xx)


# @numba.jit(nopython=True)
# def rates(T, y, p0, c0, ab, eb, ebsig, ac, ec, ag, eg, egsig):
#     l, delta, c = y
#     f = 1.0 - (l + c)
#     g1 = 2 * f - delta
#     g2 = 2 * (c - c0)
#     g = g1 + g2
#     RT = Rgas * T

#     # bridges
#     eb = eb + ebsig * invernorm(
#         1 - l / (p0 - c0))
#     kb = ab * math.exp(-eb / RT)

#     # c
#     kc = ac * np.exp(-ec / RT)
#     eg = eg + invernorm(0.5 * g / (1 - c0)) * egsig
#     kg = ag * np.exp(-eg / RT)

#     return kb, kc, kg
