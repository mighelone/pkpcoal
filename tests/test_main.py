import os

import pkp.src.CPD as cpdsl
import pkp.src.Models as mdls
import numpy as np

def test_calcC0_zeroTest():
    """ Test if c0=0.0 for no carbon and no Ox """
    assert cpdsl.CPD.calcC0(0.0, 0.0) == 0.0

def test_calcC0_maximumTestCarbon():
    """ Test if for 100% carbon the maximum c0=0.36 is returned """
    assert cpdsl.CPD.calcC0(1.0, 0.0) == 0.36

def test_calcC0_maximumTestOxygen():
    """ Test if for 100% oxygen the maximum c0=0.15 is returned """
    assert cpdsl.CPD.calcC0(0.0, 1.0) == 0.15

def test_calcC0_website_8515():
    """ Test if for 15%oxygen and 85% carbon a value of 0.035 is returned.
        Value has been calculated with the online calculator tool 
        from http://www.et.byu.edu/~tom/cpd/correlation.html """
    assert cpdsl.CPD.calcC0(0.85, 0.15) == 0.035

def test_calcCoalParam():
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

def test_calcMassCR():
    """ 

    """
    const_rate = mdls.constantRate({"k":1.0,"tstart":0.0,"finalYield":1.0})
    # Test if at t=0 no mass is released
    assert const_rate.calcMass(0.0,0.0) == 0.0
 


def test_full_main(tmpdir):
    from pkp.pkpcli import generate
    from pkp.pkpcli import fit
    fold = "/home/go/documents/code/pkp.git/inputs/"
    res = generate(folder=fold)
    fit = fit(folder=fold, results=res, selectPyrolModel="constantRate")
    print fit._tsv
    # mp.plotResults(res, fitsCR)
    # fitsArr = mp.startFittingProcedure(res, selectPyrolModel="arrheniusRate")
    # mp.plotResults(res, fitsArr)



