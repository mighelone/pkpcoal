import os

import pkp.src.CPD as cpdsl
import pkp.src.Models as mdls
import numpy as np

def test_calcC0_zeroTest():
    # Test if c0=0.0 for no carbon and no Ox
    assert cpdsl.CPD.calcC0(0.0, 0.0) == 0.0
    # Test if for 100% carbon the maximum c0=0.36 is returned
    assert cpdsl.CPD.calcC0(1.0, 0.0) == 0.36
    # Test if for 100% oxygen the maximum c0=0.15 is returned
    assert cpdsl.CPD.calcC0(0.0, 1.0) == 0.15
    # Test if for 15%oxygen and 85% carbon a value of 0.035 is returned.
    # Value has been calculated with the online calculator tool
    # from http://www.et.byu.edu/~tom/cpd/correlation.html
    assert cpdsl.CPD.calcC0(0.85, 0.15) == 0.035

def test_minimisation():
    from pkp.src.Fitter import LeastSquaresEstimator
    from pkp.src.Models import constantRate
    from pkp.src.CPD import CPDResult
    from numpy import array
    # # NOTE we mock to have more control
    # # and be independet from actual input files
    opt_params = {
         'GradBasedOpt':'fmin',
         'maxIter':1000,
         'scaleFactor': False, # What is this used for
         'Tolerance': False,   # Currently opt default
         'weightMass': 1.0,
         'weightRate': 0.0,
    }
    estimator = LeastSquaresEstimator(opt_params,finalYield=1.0)
    # # constant rate mock
    model_params = {'k':1.0, 'tstart':0.0, 'final_yield':1.0}
    constRate = constantRate(model_params)
    species = "MockSpecies"
    mock_res_dct = {
            "header": {"file": ["time(ms)",species,"temp"]},
            "data":   {"file": [ array([0.0,0.5,1.0]),
                                 array([0,0.5,1.0]),
                                 array([100.0,200.0,300.0])]}
            }
    results = CPDResult(dct=mock_res_dct)
    fit = estimator.estimate([results], constRate, species)
    print fit
    mock_res_dct['data']['file'][0] = array([0.0,0.25,0.5]) # double rate
    results = CPDResult(dct=mock_res_dct)
    constRate = constantRate(model_params)
    fitDoubleRate = estimator.estimate([results], constRate, species)
    print fitDoubleRate
    fitDoubleRateII = estimator.estimate([results], constRate, species)
    print fitDoubleRateII
    fitDoubleRateIII = estimator.estimate([results], constRate, species)
    print fitDoubleRateIII


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
    assert const_rate.calcMass(0.0,1.0) > 0.0
    # assert const_rate.calcMass(1.0,0.0) == 1.0

def test_full_main(tmpdir):
    pass
    # from pkp.pkpcli import generate
    # from pkp.pkpcli import fit
    # fold = "/home/go/documents/code/pkp.git/inputs/"
    # res = generate(folder=fold)
    # fit1run = fit(folder=fold, results=res,
    #     selectPyrolModel="constantRate")
    # fit1run.inputs['OperatingConditions']['runs'] = 1
    # fit2run = fit(folder=fold, results=res,
    #     selectPyrolModel="constantRate")
    # fit2run.inputs['OperatingConditions']['runs'] = 2
    # fit1 = fit1run.startFittingProcedure(res)
    # fit2 = fit2run.startFittingProcedure(res)
    # print fit1.res
    # print fit1._tsv
    # print fit2.res
    # print fit2._tsv

    # mp.plotResults(res, fitsCR)
    # fitsArr = mp.startFittingProcedure(res, selectPyrolModel="arrheniusRate")
    # mp.plotResults(res, fitsArr)



