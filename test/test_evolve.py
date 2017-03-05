import pytest
import numpy as np
import pkp
import pkp.evolution
import pkp.detailed_model

npop = 40
ngen = 30
t = np.linspace(0, 0.02, 100)
y = np.exp(-t * 1)
operating_conditions = [[0, 300],
                        [0.01, 1500],
                        [0.02, 1500]]

par_min = np.array([1e3, 50e6, 0.2])
par_max = np.array([1e10, 200e6, 0.8])

scaled_par = np.random.random(3)


@pytest.fixture
def ga():
    return pkp.evolution.Evolution(npop=npop, ngen=ngen)


def test_init(ga):
    assert ga._npop == npop
    assert ga._ngen == ngen

    assert ga.empirical_model == pkp.empirical_model.SFOR
    assert ga.n_targets == 0


def test_target(ga):

    ga.set_target(t=t, y=y, operating_conditions=[[0, 300], [0.1, 1500]])
    assert ga.n_targets == 1
    assert np.allclose(ga.ref_results['run0']['t'], t)
    assert np.allclose(ga.ref_results['run0']['y'], y)

    ga.set_target(t=t * 0.9, y=y * 1.1,
                  operating_conditions=[[0, 300], [0.1, 1500]])
    assert ga.n_targets == 2
    assert len(ga.ref_results['run1']['y']) == len(y)


def test_error(ga):
    ga.set_target(t=t, y=y, operating_conditions=[[0, 300], [0.1, 1500]])
    ga.operating_conditions = operating_conditions
    ga.parameters_range(par_min, par_max)

    assert np.allclose(par_min, ga._parameters_min)
    assert np.allclose(par_max, ga._parameters_max)

    unsc_par = ga.unscale_parameters(scaled_par)
    unsc_par_c = par_min + (par_max - par_min) * scaled_par
    assert np.allclose(unsc_par[1:], unsc_par_c[1:])
