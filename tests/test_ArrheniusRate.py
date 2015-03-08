from test_main import *
import pkp.src.CPD as cpdsl
import pkp.src.Models as mdls
import numpy as np
import scipy as sp

import pytest
import matplotlib.pyplot as plt

import copy

@pytest.fixture
def arrhenius_fit(exp_preProc):
    """ initialize and estimate rates from mocked preproc data"""
    from pkp.src.Models import arrheniusRate
    print opt_params
    return arrheniusRate(opt_params, [exp_preProc], "Species")

@pytest.fixture
def arrhenius_fit_multi(exp_preProc):
    """ initialize and estimate rates from mocked preproc data"""
    from pkp.src.Models import arrheniusRate
    print opt_params
    return arrheniusRate(opt_params, [exp_preProc, exp_preProc], "Species")

def test_arrhenius_init(arrhenius_fit):
    """ test if instantiation of
        constantRate model throws any errors """
    assert arrhenius_fit.name     == "ArrheniusRate"
    assert arrhenius_fit.species  == "Species"
    assert arrhenius_fit.calcMass == arrhenius_fit.calcMassArrhenius

@pytest.mark.xfail
def test_arrheniusRate_calcMass(arrhenius_fit):
    """ test if the constant rate calc mass function returns
        trivial results """
    import scipy as sp
    fit        = arrhenius_fit
    zeroRate   = [0.0, 0.0, 0.0]
    result     = linear_yield
    zeroResult = result*0.0
    # test if for zero rate the yielded mass is always zero
    temp_array_interp = sp.interpolate.interp1d(temp_array, time_array)
    assert np.allclose(zeroResult,
                fit.calcMass(zeroRate, 0.0, time_array, temp_array_interp))

@pytest.mark.xfail
def test_fit_arrheniusRate_model(arrhenius_fit):
    fit = arrhenius_fit
    optimizedParameter = fit.fit().x
    time = np.arange(0.0, 1.01, 0.01)
    print "optparams_single", optimizedParameter
    temp_array_interp = sp.interpolate.interp1d(temp_array,
            time,
            bounds_error=False,
            fill_value=temp_array[-1])
    fittedYield = fit.calcMass(optimizedParameter, time, time, temp_array_interp)
    print fittedYield

    fig, axs = plt.subplots()
    axs.plot(time, fittedYield, color='k', label="modeled yield")
    axs.plot(time, exp_yield, color='r', label="target yield", linestyle='--')

    plt.legend()
    axs.set_ylabel('yield [-]')
    axs.set_xlabel('time [s]')
    fig.savefig('tests/fit_arrhenius_rate_single.png')

@pytest.mark.xfail
def test_fit_arrheniusRate_model_multi(arrhenius_fit_multi):
    fit = arrhenius_fit_multi
    optimizedParameter = fit.fit().x
    print "optparams_multi", optimizedParameter
    time = np.arange(0.0, 1.01, 0.01)
    temp_array_interp = sp.interpolate.interp1d(temp_array,
            time,
            bounds_error=False,
            fill_value=temp_array[-1])
    fittedYield = fit.calcMass(optimizedParameter, time, time, temp_array_interp)
    print fittedYield
    fig, axs = plt.subplots()
    axs.plot(time, fittedYield, color='k', label="modeled yield")
    axs.plot(time, exp_yield, color='r', label="target yield", linestyle='--')

    plt.legend()
    axs.set_ylabel('yield [-]')
    axs.set_xlabel('time [s]')
    fig.savefig('tests/fit_arrhenius_rate_multi.png')
