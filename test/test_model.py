from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals

import pytest
import pkp.empirical_model
import numpy as np
import matplotlib.pyplot as plt
plt.style.use(['mystyle'])

operating_conditions = [[0, 500], [0.1, 1000], [0.2, 1500]]


@pytest.fixture
def sfor():
    return pkp.empirical_model.SFOR()


def test_init_model(sfor):
    assert sfor.parameters_default == [1e5, 50e6, 0.6]
    assert sfor.len_parameters == 3
    assert sfor.parameters_names == ['A', 'E', 'y0']


def test_set_parameters(sfor):
    # set parameters as list
    par_list = [1e6, 100e6, 0.5]
    sfor.parameters = par_list
    for key, value in sfor.parameters.iteritems():
        assert value == par_list[sfor.parameters_names.index(key)]

    # set as np.ndarray
    sfor.parameters = np.array(par_list)
    for key, value in sfor.parameters.iteritems():
        assert value == par_list[sfor.parameters_names.index(key)]

    # set as dict
    sfor.parameters = {k: v for k, v in zip(
        sfor.parameters_names, par_list)}
    for key, value in sfor.parameters.iteritems():
        assert value == par_list[sfor.parameters_names.index(key)]


def test_operating_conditions(sfor):
    sfor.operating_conditions = operating_conditions
    oc = operating_conditions
    assert np.alltrue(sfor.operating_conditions == np.array(oc))

    assert sfor.T(oc[0][0]) == oc[0][1]
    assert sfor.T(oc[1][0]) == oc[1][1]


def test_rate(sfor):
    par = sfor.parameters
    sfor.operating_conditions = operating_conditions

    t = 0
    T = sfor.T(t)
    rate = par['A'] * np.exp(-par['E'] / 8314.33 / T) * par['y0']

    assert np.isclose(sfor.rate(t=t, y=0), rate)


def test_run(sfor):
    parameters = {'A': 72e3,
                  'E': 55e6,
                  'y0': 0.5}
    sfor.parameters = parameters
    sfor.operating_conditions = [[0, 400],
                                 [0.005, 1400],
                                 [0.02, 1400]
                                 ]
    t, y = sfor.run()

    fig, ax = plt.subplots()
    ax.plot(t, y, label='Run Free t')

    t_max = sfor.operating_conditions[-1, 0]
    t0 = np.linspace(0, t_max, 100)

    _, y0 = sfor.run(t=t0)

    ax.plot(t0, y0, marker='x', linewidth=0, label='Run fix t')
    ax.set_xlabel('t, s')
    ax.set_ylabel('Daf solid mass fraction')
    ax.legend(loc='best')
    fig.savefig('test/test_run.pdf')
    plt.close(fig)


def test_unscale(sfor):
    par_min = np.array([1e3, 10e6, 0.4])
    par_max = np.array([1e8, 200e6, 0.6])
    scal_par = np.array([0.2, 0.5, 0.4])

    unsc_par = sfor.unscale_parameters(scal_par, par_min, par_max)
    unsc_par_c = (par_min + scal_par * (par_max - par_min))
    assert np.allclose(unsc_par[1:],
                       unsc_par_c[1:])
    assert np.isclose(np.log10(unsc_par[0]),
                      (np.log10(par_min) + scal_par * (
                          np.log10(par_max) - np.log10(par_min)))[0])
