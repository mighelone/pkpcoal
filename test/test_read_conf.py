import os
import sys
sys.path.insert(0, '../')
import pytest
import pkp
import ruamel_yaml as yaml
from pkp.detailed_model import M_elements, hf

import numpy as np

# from dictinput import settings

yml_file = 'input.yml'
with open(yml_file, 'r') as f:
    settings = yaml.load(f)


@pytest.fixture
def conf():
    return pkp.ReadConfiguration(settings)


def test_coal(conf):
    # for key in ('proximate_analysis', 'ultimate_analysis', 'rho_dry'):
    #    assert settings['Coal'][key] == getattr(conf, key)
    assert settings['Coal']['HHV'] == conf.HHV / 1e6

    op_cond = conf.operating_conditions


@pytest.fixture
def runner():
    return pkp.PKPRunner(yml_file)


def test_postulate(runner):
    # print(runner.ultimate_analysis)
    y0 = 0.6
    mw = 200
    post_dict = runner._postulate_species(y0, mw=mw)
    print(post_dict)
    mw = sum(post_dict['formula'][el] * mwi for el,
             mwi in M_elements.iteritems())
    # print(mw)
    assert np.isclose(post_dict['molecular_weight'], mw)

    formula = post_dict['formula']
    nu_CO2 = formula['C']
    nu_H2O = formula['H'] * 0.5
    nu_O2 = nu_CO2 + nu_H2O * 0.5 - formula['O'] * 0.5
    nu_SO2 = formula['S']

    lhv_vol = (runner.lhv_daf - (1 - y0) * runner.lhv_char) / y0

    # hfvol = nu_CO2 * hf['CO2'] + nu_H2O * \
    #    hf['H2O'] - nu_O2 * hf['O2'] + lhv_vol
    hf_vol = post_dict['hf']

    lhv_calc = (hf_vol + nu_O2 * hf['O2'] -
                nu_CO2 * hf['CO2'] + nu_SO2 * hf['SO2'] -
                nu_H2O * hf['H2O']) / mw
    assert np.isclose(lhv_calc, lhv_vol)


def test_dulong(runner):
    runner.ultimate_analysis = {
        'C': 0.8, 'H': 0.1, 'O': 0.1, 'S': 0, 'N': 0}
    hhv = (333 * 80 + 1442 * (10 - 10 / 8.)) * 1e3
    assert np.isclose(runner.dulong(), hhv)

    lhv_char = -(hf['CO2'] - hf['O2'] - hf['char']) / M_elements['C']
    assert np.isclose(lhv_char, runner.lhv_char)
