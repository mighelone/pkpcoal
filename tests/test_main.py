import os

import pkp.src.CPD as cpdsl
import pkp.src.Models as mdls
import numpy as np

import pytest
import matplotlib.pyplot as plt

import copy

# NOTE Mock data
opt_params = {
     'GradBasedOpt':'fmin',
     'maxIter':1000,
     'scaleFactor': False, # What is this used for
     'Tolerance': False,   # Currently opt default
     'weightMass': 1.0,
     'weightRate': 0.0,
     'constantRate': {
            'k':1.0,
            'tstart':0.0,
            'finalYield': 1.0},
}

opt_params_bounded = copy.deepcopy(opt_params)
opt_params_bounded['constantRate'].update(
    {'tstartBounds': (1e-6, 1e-5), 'finalYieldBounds': (1.0, 1.0)})

time_array   = np.arange(0.0, 1.001, 0.01)
time_array_ms= np.arange(0.0, 1.001, 0.01)/1000.0
linear_yield = np.arange(0.0, 1.001, 0.01)
exp_yield    = 1.0 - np.exp(-time_array)
temp_array   = time_array

cpd_header   = {"file": ["time(ms)", "Species", "temp"]}
cpd_linear   = {"file": [time_array_ms, linear_yield, temp_array]}
cpd_exp      = {"file": [time_array_ms, exp_yield, temp_array]}
#1.0-np.exp(-np.arange(0.0, 1.01,    0.01)),


# FIXTURES
@pytest.fixture
def linear_preProc():
    """ initialise a mocked cpd result with 
        linearily increasing  yield """
    from pkp.src.CPD import CPDResult
    res_dct = {"header": cpd_header, "data": cpd_linear}
    return CPDResult(dct=res_dct)

@pytest.fixture
def linear_cr_fit(linear_preProc):
    """ initialize and estimate rates from mocked preproc data"""
    from pkp.src.Models import constantRate
    return constantRate(opt_params, [linear_preProc], "Species")

@pytest.fixture
def exp_preProc():
    """ initialise a mocked cpd result with 
        exponentioal increasing  yield """
    from pkp.src.CPD import CPDResult
    res_dct = {"header": cpd_header, "data": cpd_exp}
    return CPDResult(dct=res_dct)

@pytest.fixture
def exp_cr_fit(exp_preProc):
    """ initialize and estimate rates from mocked preproc data"""
    from pkp.src.Models import constantRate
    return constantRate(opt_params, [exp_preProc], "Species")

@pytest.fixture
def exp_cr_fit_bounded(exp_preProc):
    """ initialize and estimate rates from mocked preproc data"""
    from pkp.src.Models import constantRate
    return constantRate(opt_params_bounded, [exp_preProc], "Species")

class TestFittingProcedures():
    """ Test class to test preformance of fitting 
        procedures
    """

    def test_cr_init(self, linear_cr_fit):
        """ test if instantiation of 
            constantRate model throws any errors """
        fit =  linear_cr_fit
        assert fit.name     == "ConstantRate"
        assert fit.species  == "Species"
        assert fit.calcMass == fit.calcMassConstRate

    def test_cr_calcMass(self, linear_cr_fit):
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
        assert np.array_equal(unityYield[2]-np.exp(-result),
                    fit.calcMass(unityYield, 0.0, time_array))

    def test_cr_compute_model_error(self, exp_cr_fit):
        """ test if model error computation gives sensible results

        """        
        fit = exp_cr_fit
        unityYield = [1.0, 0.0, 1.0]
        # first test if model and preProc have the same yield
        # over time
        assert np.array_equal(exp_yield,
            fit.calcMass(unityYield, 0.0, time_array))


        # Test the cumulative error functions
        from pkp.src.Models import Model
        def model_cumalative_error(arr):
            return Model.cumulative_error(arr, arr, arr, arr, 1.0, 1.0, 1.0)

        def model_cumalative_error_parallel(delta):
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
        errors = [model_cumalative_error_parallel(delta) for delta in deltas]

        # assume linear dependency
        assert np.allclose(deltas, errors)

        shifted_deltas = deltas - 0.5
        shifted_errors = [model_cumalative_error_parallel(delta)
            for delta in shifted_deltas]
        assert np.allclose(abs(shifted_deltas),shifted_errors)

        axs.plot(deltas, errors, color = 'k', label = "model error", linewidth = 2)
        axs.plot(deltas, shifted_errors, color = 'r', label = "model error shift", linewidth = 2)
        plt.legend()
        fig.savefig('tests/cummulative_error.png')

    def test_fit_cr_model(self, exp_cr_fit, exp_cr_fit_bounded):
        fit = exp_cr_fit
        fitBounded = exp_cr_fit_bounded
        optimizedParameter = fit.fit().x
        time = np.arange(0.0, 1.01, 0.01)
        print optimizedParameter
        fittedYield = fit.calcMass(optimizedParameter, 0.0, time)
        
        fig, axs = plt.subplots()
        axs.plot(time, fittedYield, color = 'k', label = "yield fit unbounded")
        axs.plot(time, exp_yield, color = 'r', label = "target yield")

        optimizedParameter2 = fitBounded.fit().x
        print optimizedParameter2
        fittedYieldBounded = fitBounded.calcMass(optimizedParameter2, 0.0, time)

        axs.plot(time, fittedYieldBounded,
            color = 'g', label = "yield fit", linestyle='--')

        plt.legend()
        fig.savefig('tests/fit.png')

    # def test_minimisation_linear():
    #     from pkp.src.CPD import CPDResult
    #     # # NOTE we mock to have more control
    #     # # and be independet from actual input files
    #         
    #     def constant_rate_checker(model_params):
    #         from copy import deepcopy
    #         mock_res_dct_local = deepcopy(mock_res_dct) # make a local copy
    #         fit  = init_and_est(model_params, mock_res_dct_local)
    #         fit2 = init_and_est(model_params, mock_res_dct_local)
    #         # Test if two runs with identical input 
    #         # produce identical output
    #         assert not (fit is fit2)
    #         assert np.array_equal(fit.parameter,fit2.parameter)
    #
    #         # lets double the rate
    #         mock_res_dct_double_rate = deepcopy(mock_res_dct_local)
    #         mock_res_dct_double_rate['data']['file'][0] = np.arange(0.0, 500.1, 5)
    #         fit_double = init_and_est(
    #                 model_params,
    #                 mock_res_dct_double_rate
    #             )
    #         assert not (fit is fit_double) # should be two different instances
    #         # but initial params should not be affected   
    #         assert fit_double.initialParameter == fit.initialParameter 
    #         # test if for doubled heating rate the rate constant is 
    #         # approximatley doubled and start_time is not affected 
    #         ratio = fit_double.k/fit.k
    #         # assert abs(fit.start_time - fit_double.start_time) < 0.1
    #         # assert abs(ratio - 2.0)/2.0 < 0.1
    #         return fit
    #
    #     def plot_and_save(fits, params, name):
    #         import matplotlib.pyplot as plt
    #         fig, axs = plt.subplots()
    #         colors = ['c', 'm', 'y', 'k','b','g','r']
    #         for i,fit in enumerate(fits):
    #             axs.scatter(
    #                 x=fit.time,
    #                 y=fit.mass,
    #                 marker='.',
    #                 label = str(params[i]),
    #                 color = colors[i],
    #              )
    #         results = CPDResult(dct=mock_res_dct)
    #         axs.plot(
    #              results['time'],
    #              results[species],
    #              color = 'k',
    #              label = "orig",
    #              linewidth = 2
    #         )
    #         plt.legend()
    #         fig.savefig('tests/' + name)
    #
    #     ks = [0.0, 1.0, 2.0, 10.0, 20.0, 50.0, 100.0]
    #     initial_params_k = [{'k':k, 'tstart':0.0, 'final_yield':1.0}
    #              for k in ks]
    #     
    #     tstarts = [-0.5, 0.0, 0.25, 0.5, 1.0]
    #     initial_params_t  = [{'k':10, 'tstart': tstart, 'final_yield':1.0}
    #              for tstart in tstarts]
    #
    #     fits_k = [constant_rate_checker(init) 
    #              for init in initial_params_k]
    #         
    #     fits_t = [constant_rate_checker(init) 
    #              for init in initial_params_t]
    #
    #     plot_and_save(fits_k, ks, "initial_k.png")
    #     plot_and_save(fits_t, tstarts, "initial_t.png")
    
   
# def test_modelError():
#     """ tests if computation of the deviation of
#         kinetic model to pre proc results works 
#         as expected    
#     """
#     pass
#
#
# def test_minimisation_exp():
#     model_params  = {'k': 0, 'tstart':0.0, 'final_yield':1.0}
#     fit = init_and_est(model_params, mock_res_dct_exp)
#     print "exponential fit"
#     print fit
 


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

    # Test if same CPD run gives same results
    # results = CPDResult(dct=mock_res_dct)
    # results2 = CPDResult(dct=mock_res_dct)
    # assert not (results is results2)
    # for species,field in results.iterspecies():
    #     assert results2[species] == field


