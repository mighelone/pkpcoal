import numpy as np


def x_n_calc(x, z_n, k_n_1):
    return z_n / (1 + k_n_1 * x)


def sum_x_n_calc(x, z_n, k_n_1):
    """Eq 54"""
    return np.sum(z_n / (1 + k_n_1 * x) * k_n_1)


def fp(x, sigma):
    return x * (1 - x)**(sigma - 1)


def pstar_f(x, sigma, fpp):
    return fp(x, sigma) - fpp
