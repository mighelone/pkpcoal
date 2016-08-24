from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals

import pkp.cpd
import pkp.coalnew
import numpy as np
import platform

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

if platform.system() == 'Darwin':
    cpd_exe = './cpdnlg.x'
elif platform.system() == 'Linux':
    cpd_exe = './cpdnlg.x'
elif platform.system() == 'Windows':
    cpd_exe = './cpdnlg.exe'


def init_coal():
    return pkp.coalnew.Coal(proximate_analysis=pa, ultimate_analysis=ua)


def init_cpd():
    return pkp.cpd.CPD(proximate_analysis=pa, ultimate_analysis=ua)


def test_normalize_dictionary():
    ua_norm = pkp.coalnew.normalize_dictionary(ua)
    assert sum(ua_norm.itervalues()) == 1
    assert ua['C'] / sum(ua.itervalues()) == ua_norm['C']


def test_coal_init():
    coal = init_coal()
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


def test_calc_model_parameter():
    cpd = init_cpd()
    cpd.calc_model_parameters()
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


def test_operating_conditions():
    cpd = init_cpd()
    cpd.operating_conditions = op_cond
    assert cpd.operating_conditions.shape == (3, 2)


def test_write_input_file():
    cpd = init_cpd()
    cpd.calc_model_parameters()
    cpd.path = './test'
    cpd.set_numerical_parameters(dt=1e-5,
                                 increment=2,
                                 dt_max=1e-5)
    cpd.operating_conditions = op_cond
    cpd.write_input_files(fname='test')
    cpd.run(fname='test', solver=cpd_exe)
    cpd.read_results(fname='test')
    return cpd

results = test_write_input_file()
