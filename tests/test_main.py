import os

import pkp.src.CPD as cpdsl
import pkp.src.Models as mdls
import pkp.src.CoalThermoPhysics as ctp
import numpy as np

import pytest
import matplotlib.pyplot as plt

import copy

from mock_data import *
# NOTE Mock data
#1.0-np.exp(-np.arange(0.0, 1.01,    0.01)),


# FIXTURES
@pytest.fixture
def linear_preProc():
    """ initialise a mocked cpd result with
        linearily increasing  yield """
    from pkp.src.CPD import CPDResult
    res_dct = {"header": cpd_header, "data": cpd_linear}
    return CPDResult(dct=res_dct, coal=ctp.Coal(mock_input['Coal']), runNr='run0')

@pytest.fixture
def exp_preProc():
    """ initialise a mocked cpd result with
        exponentioal increasing  yield """
    from pkp.src.CPD import CPDResult
    res_dct = {"header": cpd_header, "data": cpd_exp}
    return CPDResult(dct=res_dct, coal=ctp.Coal(mock_input['Coal']), runNr='run0')


@pytest.fixture(scope='session')
def cpd_run_single():
    from pkp.pkpcli import Generate
    return Generate(mock_input).executeSolver()

@pytest.fixture(scope='session')
def cpd_run_ifrf():
    from pkp.pkpcli import Generate
    return Generate(mock_input_ifrf).executeSolver()

@pytest.fixture(scope='session')
def cpd_run_multi():
    from pkp.pkpcli import Generate
    multi = mock_input
    multi['Coal']['OperatingConditions']['runs'] = 19
    return Generate(multi).executeSolver()

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


def almost(val,target):
    return abs((val-target)/target) < 0.05
