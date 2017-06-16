import pytest
import numpy as np

from pkp.reactor import Reactor
from pkp.empirical_model_t import SFOR, C2SM

operating_conditions = [[0, 400], [0.1, 400], [0.2, 1000]]


@pytest.fixture
def reactor():
    A, E, y0 = 1e7, 100e6, 0.5
    r = Reactor(model='SFOR', parameters=[A, E, y0])
    r.operating_conditions = operating_conditions
    return r


def test_reactor(reactor):
    assert isinstance(reactor._model, SFOR)


def test_operating_conditions(reactor):
    assert reactor.operating_conditions.shape == (3, 2)


def test_y0(reactor):
    assert reactor.y0 == [0, 400]


def test_rate(reactor):
    y, T = 0, 500
    rate = reactor.rate(0, [y, T])
    np.testing.assert_allclose(rate, reactor._model.rate(
        0, [y, T]) + [reactor._dTdt(0, [y, T])])


def test_dTdt(reactor):
    y = [0, 500]
    t = 0.05
    dTdt = reactor._dTdt(t, y)
    t1, T1 = operating_conditions[1]
    t0, T0 = operating_conditions[0]
    dTdt_c = (T1 - T0) / (t1 - t0)

    np.testing.assert_almost_equal(dTdt, dTdt_c)
    np.testing.assert_almost_equal(reactor._dTdt_array[0], dTdt_c)


def test_run(reactor):
    res = reactor.run()
    t, y = res
    np.testing.assert_almost_equal(t[0], reactor.operating_conditions[0, 0])
    np.testing.assert_almost_equal(t[-1], reactor.operating_conditions[-1, 0])
    np.testing.assert_almost_equal(
        y[-1, 1], reactor.operating_conditions[-1, 1], decimal=0)
    assert y[-1, 0] <= reactor._model.parameters.y0


def test_c2sm_run():
    r = Reactor('C2SM')

    r.operating_conditions = operating_conditions

    t, sol = r.run()

    y = sol[:, 0]
    s = sol[:, 1]
    T = sol[:, 2]

    assert np.alltrue(y + s <= 1)

    assert np.alltrue(np.diff(y) >= 0)
    assert np.alltrue(np.diff(s) <= 1)


def test_set_parameters(reactor):
    parameters_old = reactor.model_parameters

    A = 1e4

    reactor.set_parameters(A=A)
    assert reactor.model_parameters['A'] == A
    assert reactor.model_parameters['E'] == parameters_old['E']

    A, E, y0 = 1e5, 39e6, 0.9
    reactor.set_parameters(A=A, E=E, y0=y0)
    assert reactor.model_parameters['A'] == A
    assert reactor.model_parameters['E'] == E
    assert reactor.model_parameters['y0'] == y0

    dt = 1e-6
    reactor.set_parameters(A=A * 2, first_step=dt)
    assert reactor.model_parameters['A'] == A * 2
    assert reactor.model_parameters['E'] == E
    assert reactor.reactor_parameters['first_step'] == dt
