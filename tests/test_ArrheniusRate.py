from test_main import *
import pkp.src.CPD as cpdsl
import pkp.src.Models as mdls
import numpy as np

import pytest
import matplotlib.pyplot as plt

import copy

@pytest.fixture
def arrhenius_fit(exp_preProc):
    """ initialize and estimate rates from mocked preproc data"""
    from pkp.src.Models import arrheniusRate
    print opt_params
    return arrheniusRate(opt_params, [exp_preProc], "Species")

def test_arrhenius_init(arrhenius_fit):
    """ test if instantiation of 
        constantRate model throws any errors """
    assert arrhenius_fit.name     == "ArrheniusRate"
    assert arrhenius_fit.species  == "Species"
    assert arrhenius_fit.calcMass == arrhenius_fit.calcMassArrhenius

def test_arrheniusRate_calcMass(arrhenius_fit):
    """ test if the constant rate calc mass function returns
        trivial results """
    fit        = arrhenius_fit
    zeroRate   = [0.0, 0.0, 0.0]
    result     = linear_yield
    zeroResult = result*0.0
    # test if for zero rate the yielded mass is always zero
    assert np.allclose(zeroResult,
                fit.calcMass(zeroRate, 0.0, time_array, temp_array))

def test_fit_arrheniusRate_model(arrhenius_fit):
    fit = arrhenius_fit
    optimizedParameter = fit.fit().x
    time = np.arange(0.0, 1.01, 0.01)
    print optimizedParameter
    fittedYield = fit.calcMass(optimizedParameter, 0.0, time, temp_array)
    
    fig, axs = plt.subplots()
    axs.plot(time, fittedYield, color='k', label="modeled yield")
    axs.plot(time, exp_yield, color='r', label="target yield", linestyle='--')

    plt.legend()
    axs.set_ylabel('yield [-]')
    axs.set_xlabel('time [s]')
    fig.savefig('tests/fit_arrhenius_rate_single.png')
