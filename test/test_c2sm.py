from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals

import pytest
import pkp.empirical_model
from pkp.empirical_model import Rgas
import numpy as np

import matplotlib.pyplot as plt
try:
    plt.style.use('mystyle')
except:
    plt.style.use('ggplot')

par_default = pkp.empirical_model.C2SM.parameters_default()
op_conditions = [[0, 300], [0.1, 1300]]


@pytest.fixture(scope='module')
def c2sm():
    return pkp.empirical_model.C2SM()


def test_init(c2sm):
    for i, p in enumerate(c2sm.parameters_names()):
        assert getattr(c2sm.parameters, p) == par_default[i]


def test_operating_conditions(c2sm):
    c2sm.operating_conditions = op_conditions
    assert np.allclose(op_conditions, c2sm.operating_conditions)


def test_rate(c2sm):
    r = c2sm.rate(0, [0, 1])
    RT = Rgas * c2sm.T(0)
    k1 = par_default[0] * np.exp(-par_default[1] / RT)
    k2 = par_default[3] * np.exp(-par_default[4] / RT)
    k1m, k2m = c2sm._k(0)
    assert np.isclose(k1, k1m)
    assert np.isclose(k2, k2m)
    dsdt = -(k1 + k2)
    dydt = (par_default[2] * k1 + par_default[5] * k2)
    assert np.isclose(r[1], dsdt)
    assert np.isclose(r[0], dydt)


def test_run(c2sm):
    t, y = c2sm.run()
    assert y.ndim == 2
    assert y.shape[1] == 2
    assert y.shape[0] == len(t)

    fig, ax = plt.subplots()
    ax.plot(t, y[:, 1], label='Solid')
    ax.plot(t, y[:, 0], label='Volatile')
    ax.set_xlabel('t, s')
    ax.set_ylabel('Daf')
    ax.legend(loc='best')
    fig.savefig('test/test_c2sm.png')
    plt.close(fig)


def test_run_yields(c2sm):
    t, y = c2sm.run()

    oc_new = [[0, 300], [0.1, 1600]]
    c2sm.operating_conditions = oc_new
    t1, y1 = c2sm.run()

    fig, ax = plt.subplots()
    ax.plot(t, y[:, 0], label='Tmax={:4.3f}'.format(
        op_conditions[1][1]))
    ax.plot(t1, y1[:, 0], label='Tmax={:4.3f}'.format(
        oc_new[1][1]))
    ax.set_xlabel('t, s')
    ax.set_ylabel('Daf')
    ax.legend(loc='best')
    fig.savefig('test/test_c2sm_yields.png')
    plt.close(fig)


def test_unscale_parameters(c2sm):
    parameters_min = [1e3, 10e6, 0.1, 1e4, 50e6, 0.2]
    parameters_max = [1e6, 100e6, 0.5, 1e10, 200e6, 1]
    norm_parameters = [0.5] * 6
    parameters = c2sm.unscale_parameters(norm_parameters,
                                         parameters_min, parameters_max)

    for i in (0, 3):
        logA1min = np.log10(parameters_min[i])
        logA1max = np.log10(parameters_max[i])
        logA = logA1min + norm_parameters[i] * (logA1max - logA1min)
        assert np.isclose(logA, np.log10(parameters[i]))

    for i in (1, 2, 4, 5):
        A = norm_parameters[i] * (parameters_max[i] - parameters_min[i])
        A += parameters_min[i]
        assert np.isclose(A, parameters[i])
