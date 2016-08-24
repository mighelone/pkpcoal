from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals

import pkp.cpd
import pkp.coalnew
import numpy as np
import platform
import pytest

ua = {'C': 69,
      'H': 5,
      'O': 24.7,
      'N': 0.8,
      'S': 0.5}

pa = {'FC': 45.1,
      'VM': 50.6,
      'Ash': 4.3,
      'Moist': 19.0}

op_cond = [[0, 500],
           [0.001, 1400],
           [0.01, 1400]]


@pytest.fixture
def coal():
    return pkp.coalnew.Coal(proximate_analysis=pa, ultimate_analysis=ua)


@pytest.fixture
def cpd():
    return pkp.cpd.CPD(proximate_analysis=pa, ultimate_analysis=ua)


def test_normalize_dictionary():
    ua_norm = pkp.coalnew.normalize_dictionary(ua)
    assert sum(ua_norm.itervalues()) == 1
    assert ua['C'] / sum(ua.itervalues()) == ua_norm['C']


def test_coal_init(coal):
    assert ua['C'] / sum(ua.itervalues()) == coal.ultimate_analysis['C']
    assert coal.daf == (coal.proximate_analysis[
                        'VM'] + coal.proximate_analysis['FC'])
    assert coal.proximate_analysis_daf[
        'VM'] == coal.proximate_analysis['VM'] / coal.daf
    assert 'Ash' not in coal.proximate_analysis_daf

    ua_new = ua.copy()
    ua_new['C'] = 60
    coal.ultimate_analysis = ua_new
    assert pkp.coalnew.normalize_dictionary(
        ua_new) == coal.ultimate_analysis


def test_set_NMR(cpd):
    cpd._set_NMR_parameters()
    parameters_keys = ['mdel', 'mw', 'p0', 'sig']
    for key in parameters_keys:
        assert hasattr(cpd, key)

    assert np.isclose(cpd.fcar * 100, 69.0)
    assert np.isclose(cpd.VMdaf * 100, 52.87, atol=0.05)
    assert np.isclose(cpd.mdel, 46.5, atol=0.05)
    assert np.isclose(cpd.mw, 340.3, atol=0.05)
    assert np.isclose(cpd.p0, 0.620, atol=0.05)
    assert np.isclose(cpd.sig, 4.57, atol=0.05)
    assert np.isclose(cpd.c0, 0.15, atol=0.05)

    nmr = {'mdel': 45.5,
           'mw': 340,
           'p0': 0.6,
           'sig': 4.6,
           'c0': 0.1}
    cpd._set_NMR_parameters(nmr_parameters=nmr)
    for p in nmr:
        assert nmr[p] == getattr(cpd, p)


def test_set_parameters(cpd):
    par = cpd.get_parameters()
    print(par)

    par['basename'] = 'test'
    par['increment'] = 5

    cpd.set_parameters(**par)
    for key in par:
        if not key == 'nmr_parameters':
            assert getattr(cpd, key) == par[key]
        else:
            for p in par[key]:
                assert getattr(cpd, p) == par[key][p]


def test_set_numerical(cpd):
    par = cpd.get_parameters()
    par['increment'] = 5
    cpd._set_numerical_parameters(**par)
    assert cpd.increment == par['increment']


def test_operating_conditions(cpd):
    cpd.operating_conditions = op_cond
    assert cpd.operating_conditions.shape == (3, 2)


def test_run(cpd):
    cpd.path = './test'  # add path to set property
    cpd.operating_conditions = op_cond
    cpd.set_parameters(dt=1e-5, increment=2, dt_max=1e-5,
                       basename='test', nmr_parameters=None)
    res = cpd.run()
    res.index.name == 'Time(ms)'
