from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals

import pytest
import pkp.empirical_model
import numpy as np

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

    assert np.isclose(sfor.rate(t=t, y=1), rate)
