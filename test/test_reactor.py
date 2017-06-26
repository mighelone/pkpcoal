"""Test reactor module."""
import pytest
import numpy as np

from pkp.reactor import Reactor, DTR
from pkp.cpd import CPD
from pkp.empirical_model import SFOR

operating_conditions = [[0, 400], [0.1, 400], [0.2, 1000]]

sfor_parameters = {'A': 1e7, 'E': 50e6, 'y0': 0.5}
max_step = 1e-2


@pytest.fixture
def reactor():
    """Init the Reactor using the SFOR model."""
    r = Reactor(model='SFOR', max_step=max_step, **sfor_parameters)
    r.operating_conditions = operating_conditions
    return r


def test_reactor(reactor):
    """Test initialization of the reactor."""
    assert isinstance(reactor.model, SFOR)
    assert reactor.model_parameters == sfor_parameters
    assert reactor.reactor_parameters['max_step'] == max_step


def test_operating_conditions(reactor):
    """Test the operating conditions."""
    assert reactor.operating_conditions.shape == (3, 2)
    np.testing.assert_array_equal(reactor.operating_conditions,
                                  operating_conditions)


def test_y0(reactor):
    """Test intial solution of the ODE."""
    np.testing.assert_array_equal(reactor.y0, [0, operating_conditions[0][1]])


def test_rate(reactor):
    """Test the reactor rate used by the ODE solver."""
    y, T = 0, 500
    rate = reactor.rate(0, [y, T])
    dydt = reactor._model.rate(0, [y, T])
    np.testing.assert_allclose(rate, dydt + [reactor._dTdt(0, [y, T], dydt)])


def test_dTdt(reactor):
    """Test the temoerature gradient."""
    y = [0, 500]
    t = 0.05
    dydt = [0, 0]
    dTdt = reactor._dTdt(t, y, dydt)
    t1, T1 = operating_conditions[1]
    t0, T0 = operating_conditions[0]
    dTdt_c = (T1 - T0) / (t1 - t0)

    np.testing.assert_almost_equal(dTdt, dTdt_c)
    np.testing.assert_almost_equal(reactor._dTdt_array[0], dTdt_c)


def test_run(reactor):
    """Test run of the reactor."""
    res = reactor.run()
    t, y = res
    np.testing.assert_almost_equal(t[0], reactor.operating_conditions[0, 0])
    np.testing.assert_almost_equal(t[-1], reactor.operating_conditions[-1, 0])
    np.testing.assert_almost_equal(
        y[-1, 1], reactor.operating_conditions[-1, 1], decimal=0)
    assert y[-1, 0] <= reactor._model.parameters.y0


def test_c2sm_run():
    """Test reactor with C2SM."""
    r = Reactor('C2SM')

    r.operating_conditions = operating_conditions

    assert len(r.y0) == 3

    t, sol = r.run()

    y = sol[:, 0]
    # s = sol[:, 1]
    # T = sol[:, 1]
    parameters = r.model_parameters
    assert np.alltrue(y <= max(parameters['y1'], parameters['y2']))

    assert np.alltrue(np.diff(y) >= 0)
    # assert np.testing.assert_almost_equal(T[-1], operating_conditions[-1][1],
    #                                      decimal=0)


def test_set_parameters(reactor):
    """Test the parameters initialization."""
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


def test_reactor_cpd():
    """Test the reactor with the CPD model."""
    ua = {'C': 80, 'H': 8, 'O': 12, 'N': 0, 'S': 0}
    pressure = 5e5
    r = Reactor(model='CPD', ultimate_analysis=ua, pressure=pressure,
                name='test', max_step=1e-2)
    r.operating_conditions = operating_conditions

    assert isinstance(r.model, CPD)
    assert r.model.pressure == pressure
    assert r.model.ultimate_analysis['C'] == ua['C']/sum(ua.values())
    assert r.reactor_parameters['max_step'] == max_step

    sol = r.run()
    # check that the sum of solid, light_gas and tar is always ==1
    np.testing.assert_allclose(sol[['solid', 'gas', 'tar']].sum(axis=1),
                               1)


@pytest.fixture()
def dtr():
    """Init a Drop Tube Reactor."""
    r = DTR('SFOR')
    r.operating_conditions = operating_conditions
    return r


def test_gas_temperature(dtr):
    """Test the interpolation of the gas temperature."""
    assert all(dtr.Tg(point[0]) == point[1] for point in operating_conditions)


def test_calc_mass(dtr):
    """Test the evaluation of the particle mass."""
    assert dtr.mp0 == dtr.mash + dtr.mdaf
    y = 0.2
    assert dtr.calc_mass(0.2) == dtr.mash + (1-y) * dtr.mdaf


def test_run_dtr(dtr):
    """Test run of DTR."""
    t, y = dtr.run()
    assert all(y[:, 1] <= y[:, 2])
