import pytest

# from pkp.empirical_model_t import EmpiricalModel
from pkp.empirical_model import namedtuple_with_defaults, Rgas
from pkp.empirical_model import SFOR, C2SM
import numpy as np


def test_check_namedtuple():
    parameters = namedtuple_with_defaults('Test', field_names=('X', 'Y'),
                                          default_values=(0, 0))

    p = parameters()

    assert p.X == 0
    assert p.Y == 0
    assert len(p._fields) == 2
    assert p._fields == ('X', 'Y')


@pytest.fixture
def sfor():
    A, E, y0 = 1e7, 100e6, 0.5ebe
    return SFOR(A, E, y0)


@pytest.fixture
def c2sm():
    return C2SM()


def test_parameters(sfor):
    A, E, y0 = 1e5, 200e6, 0.6

    sfor.set_parameters(A)
    assert A == sfor.parameters.A
    assert sfor.parameters.E == sfor.parameters_default()[1]
    assert sfor.parameters.y0 == sfor.parameters_default()[2]

    sfor.set_parameters(A, E)
    assert A == sfor.parameters.A
    assert sfor.parameters.E == E
    assert sfor.parameters.y0 == sfor.parameters_default()[2]

    sfor.set_parameters(A, E, y0)
    assert A == sfor.parameters.A
    assert sfor.parameters.E == E
    assert sfor.parameters.y0 == y0

    sfor.set_parameters([A, E, y0])
    assert A == sfor.parameters.A
    assert sfor.parameters.E == E
    assert sfor.parameters.y0 == y0

    sfor.set_parameters(A=A, E=E, y0=y0)
    assert A == sfor.parameters.A
    assert sfor.parameters.E == E
    assert sfor.parameters.y0 == y0

    sfor.set_parameters(E=E)
    assert sfor.parameters.E == E
    assert sfor.parameters.y0 == sfor.parameters_default()[2]
    assert sfor.parameters.A == sfor.parameters_default()[0]


def test_parameters_list(sfor):
    par = sfor.parameters_list
    assert par[0] == sfor.parameters.A
    assert par[1] == sfor.parameters.E
    assert par[2] == sfor.parameters.y0
    par_dict = sfor.parameters_dict
    assert all(getattr(sfor.parameters, key) ==
               value for key, value in par_dict.items())


def test_sfor(sfor):
    A, E, y0 = sfor.parameters.A, sfor.parameters.E, sfor.parameters.y0
    model_1 = SFOR([A, E, y0])
    assert model_1.parameters.A == sfor.parameters.A


def test_sfor_rate(sfor):
    y, T = 0, 800
    rate = sfor.rate(0, [y, T])

    parameters = sfor.parameters

    k = sfor._calc_k(T)
    k_calc = parameters.A * np.exp(-parameters.E / Rgas / T)
    np.testing.assert_almost_equal(k, k_calc)

    rate = sfor.rate(0, [y, T])
    rate_calc = k_calc * (parameters.y0 - y)

    np.testing.assert_almost_equal(rate[0], rate_calc)
    # assert np.testing.assert_almost_equal(rate[0],
    #                                      (parameters.A *
    #                                       np.exp(-parameters.E / Rgas / T) *
    #                                       (parameters.y0 - y)))


def test_c2sm_rate(c2sm):
    y, s, T = 0.1, 0.6, 700
    parameters = c2sm.parameters
    k1, k2 = c2sm._k(T)

    k1c = parameters.A1 * np.exp(-parameters.E1 / Rgas / T)
    k2c = parameters.A2 * np.exp(-parameters.E2 / Rgas / T)

    np.testing.assert_almost_equal(k1, k1c)
    np.testing.assert_almost_equal(k2, k2c)

    rates = c2sm.rate(0, [y, s, T])

    dydt = (parameters.y1 * k1c + parameters.y2 * k2c) * s
    dsdt = - (k1c + k2c) * s

    np.testing.assert_almost_equal(dydt, rates[0])
    np.testing.assert_almost_equal(dsdt, rates[1])
