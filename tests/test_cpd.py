import numpy as np

import pkp.src.CPD as cpdsl
from test_main import linear_preProc
from test_main import cpd_linear
from test_main import cpd_run_single
from test_main import cpd_run_multi

from mock_data import *

import matplotlib.pyplot as plt

class TestCPDLauncher():

    def test_calcC0_zeroTest(self):
        # Test if c0=0.0 for no carbon and no Ox
        assert cpdsl.CPD.calcC0(0.0, 0.0) == 0.0
        # Test if for 100% carbon the maximum c0=0.36 is returned
        assert cpdsl.CPD.calcC0(1.0, 0.0) == 0.36
        # Test if for 100% oxygen the maximum c0=0.15 is returned
        assert cpdsl.CPD.calcC0(0.0, 1.0) == 0.15
        # Test if for 15% oxygen and 85% carbon a value of 0.035 is returned.
        # Value has been calculated with the online calculator tool
        # from http://www.et.byu.edu/~tom/cpd/correlation.html
        assert cpdsl.CPD.calcC0(0.85, 0.15) == 0.035

    def test_CPDResults_timeconversion(self, linear_preProc):
        results = linear_preProc
        assert np.array_equal(
                    results['time'],
                    cpd_linear['file'][0]/1000.0)
        assert np.array_equal(
                    results['time(ms)'],
                    cpd_linear['file'][0])

    def test_calcCoalParam(self):
        """ Test if for
                10% oxygen and
                70% carbon
                3%  hydrogen
                2%  nitrogen
                15% vm
            mdel == 43.3
            MW per cluster == 551.2
            Po == 0.65
            sigma+1 ==  5.91
            c0 == 0.0
            Value has been calculated with the online calculator tool
            from http://www.et.byu.edu/~tom/cpd/correlation.html """

        ret = cpdsl.CPD.CalcCoalParam(
            {'Carbon': 0.7, 'Hydrogen': 0.03,
             'Oxygen': 0.1, 'Nitrogen': 0.2},
            {'Volatile Matter': 0.15})

        expected = { 'c0'   : 0.0,
                     'mdel' : 43.27,
                     'mw'   : 551.16,
                     'sig'  : 5.91,
                     'p0'   : 0.65, }
        for key, value in expected.iteritems():
            assert abs(ret[key] - value) < 0.01

def test_cpd_results(cpd_run_single):
    for runNr, (opCond, res) in cpd_run_single.iteritems():
        MW_POST = 30
        print res.Qfactor(mock_pa)
        print res.VolatileCompositionMass(mock_pa,mock_ua)
        print res.VolatileCompositionMol(mock_pa,mock_ua,MW_POST)
        print res.ProductCompositionMol(mock_pa,mock_ua,MW_POST)
        print res.EnthalpyOfFormation(mock_pa,mock_ua,MW_POST,28220.0)

def test_cpd_results_multi(cpd_run_multi):
    res = cpd_run_multi
    heating_rates, qfactor = [], []

    for runNr,(operCond, result) in res.iteritems():
        heating_rates.append(2200/operCond[1][0])
        qfactor.append(result.Qfactor(mock_pa))
    
    fig, axs = plt.subplots()
    axs.scatter(heating_rates, qfactor, label="qFactor")

    plt.legend()
    axs.set_xlim([0,None])
    axs.set_ylim([1.4,1.6])
    axs.set_ylabel('qfactor [-]')
    axs.set_xlabel('heating rate [K/s]')
    fig.savefig('tests/qfactor_from_CPD.png')
