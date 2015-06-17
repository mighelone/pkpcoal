from test_main import *
import pkp.src.CPD as cpdsl
import pkp.src.Models as mdls
import numpy as np

import pytest
import matplotlib.pyplot as plt

import copy

@pytest.fixture
def linear_cr_fit(linear_preProc):
    """ initialize and estimate rates from mocked preproc data"""
    from pkp.src.Models import constantRate
    return constantRate(opt_params, [linear_preProc], "Species")

@pytest.fixture
def exp_cr_fit(exp_preProc):
    """ initialize and estimate rates from mocked preproc data"""
    from pkp.src.Models import constantRate
    return constantRate(opt_params, {'run0': exp_preProc}, "Species")

@pytest.fixture
def exp_cr_fit_bounded(exp_preProc):
    """ initialize and estimate rates from mocked preproc data"""
    from pkp.src.Models import constantRate
    return constantRate(opt_params_bounded, {'run0': exp_preProc}, "Species")

@pytest.fixture
def exp_cr_fit_bounded_multi(exp_preProc):
    """ initialize and estimate rates from mocked preproc data"""
    from pkp.src.Models import constantRate
    return constantRate(
        opt_params_bounded, {'run0':exp_preProc, 'run1':exp_preProc}, "Species")

def test_cr_init(linear_cr_fit):
    """ test if instantiation of 
        constantRate model throws any errors """
    fit =  linear_cr_fit
    assert fit.name     == "ConstantRate"
    assert fit.species  == "Species"
    assert fit.calcMass == fit.calcMassConstRate

def test_cr_calcMass(linear_cr_fit):
    """ test if the constant rate calc mass function returns
        trivial results """
    fit        = linear_cr_fit
    zeroRate   = [0.0, 0.0, 1.0]
    result     = linear_yield
    zeroResult = result*0.0
    zeroYield  = [1.0, 0.0, 0.0]
    unityYield = [1.0, 0.0, 1.0]
    # test if for zero rate the yielded mass is always zero
    assert np.array_equal(zeroResult,
                fit.calcMass(zeroRate, 0.0, time_array))
    assert np.array_equal(zeroResult,
                fit.calcMass(zeroYield, 0.0, time_array))
    #assert np.array_equal(unityYield[2]-np.exp(-result),
    #            fit.calcMass(unityYield, 0.0, time_array))

def test_cr_compute_model_error(exp_cr_fit):
    """ test if model error computation gives sensible results

    """        
    fit = exp_cr_fit
    unityYield = [1.0, 0.0, 1.0]
    # first test if model and preProc have the same yield
    # over time
    # assert np.array_equal(exp_yield,
    #     fit.calcMass(unityYield, 0.0, time_array))


    # Test the cumulative error functions
    from pkp.src.Models import Model
    def model_cumalative_error(arr):
        return Model.cumulative_error(arr, arr, arr, arr, 1.0, 1.0, 1.0)

    def model_cumalative_error_parallel_lines(delta):
        arr = np.zeros(10) 
        arr2 = arr + delta
        return Model.cumulative_error(arr, arr2, arr, arr2, 0.0, 1.0, 1.0)

    # two identical lines should result in zero error
    # no matter what the weights are
    assert model_cumalative_error(np.zeros(10)) == 0.0
    assert model_cumalative_error(np.ones(10)) == 0.0

    assert Model.cumulative_error(
            np.ones(10), np.zeros(10),
            np.ones(10), np.zeros(10),
            1.0, 1.0, 1.0) == 2.0

    assert Model.cumulative_error(
            np.zeros(10), np.zeros(10),
            np.ones(10), np.zeros(10),
            1.0, 1.0, 1.0) == 1.0

    assert Model.cumulative_error(
            np.ones(10), np.zeros(10),
            np.ones(10), np.zeros(10),
            0.0, 1.0, 1.0) == 1.0

    fig, axs = plt.subplots()

    # TODO Axis labels
    deltas = np.arange(0.0, 1.01, 0.01)
    errors = [model_cumalative_error_parallel_lines(delta)
                for delta in deltas]

    # assume linear dependency
    assert np.allclose(deltas, errors)

    shifted_deltas = deltas - 0.5
    shifted_errors = [model_cumalative_error_parallel_lines(delta)
        for delta in shifted_deltas]
    assert np.allclose(abs(shifted_deltas),shifted_errors)

    axs.plot(deltas, errors, color = 'k', label = "model error", linewidth = 2)
    axs.plot(deltas, shifted_errors, color = 'r',
        label = "model error shift", linewidth = 2)
    axs.set_ylabel('return values of error function')
    axs.set_xlabel('parallel lines delta')
    plt.legend()
    fig.savefig('tests/cummulative_error.png')

def test_fit_cr_model(exp_cr_fit, exp_cr_fit_bounded):
    fit = exp_cr_fit
    fitBounded = exp_cr_fit_bounded
    fit.fit()
    optimizedParameter = fit.parameter
    time = np.arange(0.0, 1.01, 0.01)
    print optimizedParameter
    fittedYield = fit.calcMass(optimizedParameter, 0.0, time)
    
    fig, axs = plt.subplots()
    axs.plot(time, fittedYield, color = 'k', label = "yield fit unbounded")
    axs.plot(time, exp_yield, color = 'r', label = "target yield")

    fitBounded.fit()
    optimizedParameter2 = fitBounded.parameter
    print optimizedParameter2
    fittedYieldBounded = fitBounded.calcMass(optimizedParameter2, 0.0, time)

    axs.plot(time, fittedYieldBounded,
        color = 'g', label = "yield fit", linestyle='--')

    plt.legend()
    fig.savefig('tests/fit_constant_rate_single.png')

@pytest.mark.xfail
def test_genetic_model(exp_cr_fit_bounded_multi):
    fit = exp_cr_fit_bounded_multi
    time = np.arange(0.0, 1.01, 0.01)
    fit.fit()
    optimizedParameter = fit.parameter
    print optimizedParameter
    fittedYield = fit.calcMass(optimizedParameter, 0.0, time)
    
    fig, axs = plt.subplots()
    axs.plot(time, fittedYield, color = 'k', label = "yield fit unbounded")
    axs.plot(time, exp_yield, color = 'r', label = "target yield")

    plt.legend()
    fig.savefig('tests/fit_constant_rate_multi.png')
