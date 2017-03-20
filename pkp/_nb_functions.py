"""
Here are defined the functions in numba
"""
import numba
import numpy as np


@numba.jit(nopython=True)
def x_n_calc_i(x, z_n, k_n_1):
    return z_n / (1 + k_n_1 * x)


@numba.jit(nopython=True)
def x_n_calc(x, z_n, k_n_1):
    x_n = np.empty_like(z_n)
    for i in range(len(z_n)):
        x_n[i] = x_n_calc_i(x, z_n[i], k_n_1[i])
        # x_n[i] = z_n[i] / (1 + k_n_1[i] * x)
    return x_n


@numba.jit(nopython=True)
def sum_x_n_calc(x, z_n, k_n_1):
    _sum = 0
    for i in range(len(z_n)):
        # eq. 54
        x_n_calc = x_n_calc_i(x, z_n[i], k_n_1[i])
        _sum += k_n_1[i] * x_n_calc
    return _sum


@numba.jit(nopython=True)
def fp(x, sigma):
    return x * (1 - x)**(sigma - 1)


@numba.jit(nopython=True)
def pstar_f(x, sigma, fpp):
    return fp(x, sigma) - fpp
