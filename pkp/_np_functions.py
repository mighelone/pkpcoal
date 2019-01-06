import numpy as np
from scipy import special

from .interpolate import interp
from ._cpd_correlations import xx, yy


def x_n_calc(x, z_n, k_n_1):
    return z_n / (1 + k_n_1 * x)


def sum_x_n_calc(x, z_n, k_n_1):
    """Eq 54"""
    return np.sum(z_n / (1 + k_n_1 * x) * k_n_1)


def fp(x, sigma):
    return x * (1 - x) ** (sigma - 1)


def pstar_f(x, sigma, fpp):
    return fp(x, sigma) - fpp


def combinln(n, k):
    """

    Combined function

    Parameters
    ----------
    n:
    k:
    """
    return special.gammaln(n + 1) - (
        special.gammaln(k + 1) + special.gammaln(n - k + 1)
    )


def binomial(k, n, p):
    """
    Binomial function logarithmic

    Parameters
    ----------
    k:
    n:
    p:
    """
    return np.exp(combinln(n, k) + k * np.log(p) + (n - k) * np.log(1 - p))


def invernorm(y):
    """
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
    # if y > 0.5:
    #    yp = y
    #    fac = 1
    # else:
    #    yp = 1 - y
    #    fac = -1
    yp, fac = (y, 1) if y > 0.5 else (1 - y, -1)
    # return fac * np.interp(yp, yy, xx, right=3.4)
    return fac * interp(yp, yy, xx)
