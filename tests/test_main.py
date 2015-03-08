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
     'arrheniusRate': {
            'preExp':1.0,
            'beta':0.0,
            'activationEnergy': 1.0},
}

mock_ua = {
       "Carbon": 69.20,
       "Hydrogen": 4.40,
       "Nitrogen": 1.42,
       "Oxygen": 9.98,
       "Sulphur": 0.58,
}

mock_pa = {
       "Fixed Carbon": 50.80,
       "Volatile Matter": 34.80,
       "Moisture": 5.80,
       "Ash": 8.60,
}

mock_cpd={
    'fit'    : 'constantRate',
    'active' : 'True',
    'deltaT' : '5e-06',
    'MW_TAR' : '130',
}

mock_opcond = {
    'pressure' : 1.0,
    'runs' : 1,
    'run1' : [ [ 0, 300], [ 0.034, 2000] ]
}

mock_input={'Coal': {
    'Proximate Analysis': mock_pa,
    'Ultimate Analysis': mock_ua
    },
    'CPD': mock_cpd,
    'OperatingConditions': mock_opcond
}

opt_params_bounded = copy.deepcopy(opt_params)
opt_params_bounded['constantRate'].update(
    {'tstartBounds': (1e-6, 1e-5), 'finalYieldBounds': (1.0, 1.0)})

time_array   = np.arange(0.0, 1.001, 0.01)
time_array_ms= np.arange(0.0, 1.001, 0.01)/1000.0
linear_yield = np.arange(0.0, 1.001, 0.01)
exp_yield    = 1.0 - np.exp(-time_array*5)
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
def exp_preProc():
    """ initialise a mocked cpd result with
        exponentioal increasing  yield """
    from pkp.src.CPD import CPDResult
    res_dct = {"header": cpd_header, "data": cpd_exp}
    return CPDResult(dct=res_dct)


@pytest.fixture
def cpd_run():
    from pkp.pkpcli import Generate
    return Generate(mock_input)

def test_cpd_results(cpd_run):
    res = cpd_run.executeSolver()
    print res[0].Qfactor(mock_pa)
    print res[0].VolatileCompositionMass(mock_pa,mock_ua)
    print res[0].VolatileCompositionMol(mock_pa,mock_ua,150)

class TestFittingProcedures():
    """ Test class to test preformance of fitting
        procedures
    """
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


