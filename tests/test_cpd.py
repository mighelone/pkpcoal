import numpy as np
import pytest

import pkp.src.CPD as cpdsl
from test_main import (linear_preProc, 
    cpd_linear, cpd_run_single, 
    cpd_run_multi, cpd_run_ifrf)
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
    from pkp.src.CoalThermoPhysics import Coal, R
    from pkp.src.CoalThermoPhysics import PostulateSubstance
    for runNr, (opCond, res) in cpd_run_single.iteritems():
        q = res.Qfactor()
        print "qFactor", q
        ps = PostulateSubstance(coal=res.coal,qFactor=q) 
        print ps.VolatileCompositionMass()
        print ps.VolatileCompositionMol()
        print ps.ProductCompositionMol()
        print ps.mechanism() 
        print "a6: ", ps.EnthalpyOfFormation() / R 

def test_cpd_results_ifrf(cpd_run_ifrf):
    from pkp.src.CoalThermoPhysics import Coal, R
    from pkp.src.CoalThermoPhysics import PostulateSubstance
    for runNr, (opCond, res) in cpd_run_ifrf.iteritems():
        q = res.Qfactor()
        print "qFactor", q
        ps = PostulateSubstance(coal=res.coal, qFactor=1.6) 
        print ps.VolatileCompositionMass()
        print ps.VolatileCompositionMol()
        print ps.ProductCompositionMol()
        print ps.mechanism() 
        print "ifrf, a6: ", ps.EnthalpyOfFormation() / R 

def test_cpd_results_ifrf(cpd_run_ifrf):
    from pkp.src.CoalThermoPhysics import Coal, R
    from pkp.src.CoalThermoPhysics import PostulateSubstance
    print "qFactor=1.45"
    for runNr, (opCond, res) in cpd_run_ifrf.iteritems():
        q = res.Qfactor()
        print "qFactor", q
        ps = PostulateSubstance(coal=res.coal, qFactor=1.45) 
        print ps.VolatileCompositionMass()
        print ps.VolatileCompositionMol()
        print ps.ProductCompositionMol()
        print ps.mechanism() 
        print "ifrf, a6: ", ps.EnthalpyOfFormation() / R 


def test_cpd_results_multi(cpd_run_multi):
    res = cpd_run_multi
    heating_rates, qfactor = [], []

    for runNr,(operCond, result) in res.iteritems():
        heating_rates.append(2200/operCond[1][0])
        qfactor.append(result.Qfactor())
    
    fig, axs = plt.subplots()
    axs.scatter(heating_rates, qfactor, label="qFactor")

    plt.legend()
    axs.set_xlim([0,None])
    axs.set_ylim([1.4,1.6])
    axs.set_ylabel('qfactor [-]')
    axs.set_xlabel('heating rate [K/s]')
    fig.savefig('tests/qfactor_from_CPD.png')
