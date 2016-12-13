from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals
from builtins import dict

import os
import pytest
import pkp
import pkp.runner
import cantera

try:
    import ruamel_yaml as yaml
except:
    try:
        import ruamel.yaml as yaml
    except:
        import yaml
from pkp.detailed_model import M_elements, hf

import numpy as np

# from dictinput import settings

yml_file = 'input.yml'
with open(yml_file, 'r') as f:
    settings = yaml.load(f)


@pytest.fixture
def conf():
    return pkp.runner.ReadConfiguration(settings)


@pytest.fixture
def runner():
    return pkp.runner.PKPRunner(yml_file)


def test_postulate(runner):
    # print(runner.ultimate_analysis)
    y0 = 0.6
    mw = 200
    post_dict = runner.postulate_species(y0, mw=mw, include_nu=True)
    print(post_dict)
    mw = sum(post_dict['formula'][el] * mwi for el,
             mwi in M_elements.items())
    # print(mw)
    assert np.isclose(post_dict['molecular_weight'], mw)

    formula = post_dict['formula']
    nu_CO2 = formula['C']
    nu_H2O = formula['H'] * 0.5
    nu_SO2 = formula['S']
    nu_N2 = formula['N'] * 0.5
    nu_O2 = nu_CO2 + nu_H2O * 0.5 - formula['O'] * 0.5 + nu_SO2
    assert nu_CO2 == post_dict['nu']['CO2']
    assert nu_H2O == post_dict['nu']['H2O']
    assert np.isclose(nu_O2, -post_dict['nu']['O2'])
    assert nu_SO2 == post_dict['nu']['SO2']
    assert nu_N2 == post_dict['nu']['N2']

    lhv_vol = (runner.lhv_daf - (1 - y0) * runner.lhv_char) / y0

    # hfvol = nu_CO2 * hf['CO2'] + nu_H2O * \
    #    hf['H2O'] - nu_O2 * hf['O2'] + lhv_vol
    hf_vol = post_dict['hf']

    lhv_calc = (hf_vol + nu_O2 * hf['O2'] -
                nu_CO2 * hf['CO2'] - nu_SO2 * hf['SO2'] -
                nu_H2O * hf['H2O'] -
                nu_N2 * hf['N2']) / mw
    assert np.isclose(lhv_calc, lhv_vol)


def test_dulong(runner):
    runner.ultimate_analysis = {
        'C': 0.8, 'H': 0.1, 'O': 0.1, 'S': 0, 'N': 0}
    hhv = (333 * 80 + 1442 * (10 - 10 / 8.)) * 1e3
    assert np.isclose(runner.dulong(), hhv)

    lhv_char = -(hf['CO2'] - hf['O2'] - hf['char']) / M_elements['C']
    assert np.isclose(lhv_char, runner.lhv_char)


def test_empirical(runner):
    '''Test empirical composition'''
    comp_dict = runner.empirical_composition(y0=0.55, tar=0.3, CO=0.1)
    # pytest.set_trace()
    comp = comp_dict['composition']
    gas = cantera.Solution(
        os.path.join(os.path.dirname(pkp.bins.__file__),
                     '52.xml'))
    heat_vol = sum(val * runner.heat_of_reaction_species(sp, gas)
                   for sp, val in comp.items())

    species = ['CH4', 'O2', 'CO2', 'H2O']
    hf = [gas.species(sp).thermo.h(273) for sp in species]
    coeff = [1, 2, -1, -2]
    deltaH = np.dot(hf, coeff) / gas.molecular_weights[
        gas.species_index('CH4')]
    assert np.isclose(deltaH, runner.heat_of_reaction_species('CH4',
                                                              gas))
    print('LHV CH4', deltaH / 1e6)
    # check the heat of volatiles
    assert np.isclose(heat_vol, runner.heat_of_volatiles(comp, gas))

    # check heat of pyrolysis
    heat_of_pyrolysis = ((heat_vol - runner.lhv_daf) /
                         (1 - comp['char']))
    assert np.isclose(heat_of_pyrolysis,
                      runner.heat_of_pyrolysis(comp, gas))

    # check composition == element balance
    sum_ua = (sum(runner.ultimate_analysis.values()) -
              runner.ultimate_analysis['S'])
    print(sum_ua)
    print(runner.ultimate_analysis)
    for el in ('C', 'H', 'O', 'N'):
        sum_c = 0
        M_c = gas.atomic_weight(el)
        for sp, val in comp.items():
            if sp == 'char':
                if el == 'C':
                    sum_c += val
            else:
                n_c = gas.species(sp).composition.get(el, 0)
                mw = gas.molecular_weights[gas.species_index(sp)]
                sum_c += val * n_c * M_c / mw
        print('Sum ', el, sum_c)

        assert np.isclose(sum_c, runner.ultimate_analysis[el] / sum_ua)
