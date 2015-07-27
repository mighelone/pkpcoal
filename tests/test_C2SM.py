__author__ = 'vascella'

from test_main import *
import pkp.src.CPD as cpdsl
import pkp.src.Models as mdls
import numpy as np

import pytest
import matplotlib.pyplot as plt

import copy


@pytest.fixture
def c2sm_fit(exp_preProc):
    """ initialize and estimate rates from mocked preproc data"""
    from pkp.src.Models import C2SM
    print opt_params
    return C2SM(opt_params, {'run0': exp_preProc}, "Species")


def test_c2sm_init(c2sm_fit):
    """ test if instantiation of
        constantRate model throws any errors """
    assert c2sm_fit.name     == "C2SM"
    assert c2sm_fit.species  == "Species"
    #assert c2sm_fit.calcMass == 0

def test_c2sm_calcMass(c2sm_fit):
    '''
    test if the mass is calculated
    :param c2sm_fit:
    :return:
    '''
    import scipy as sp
    parameter = c2sm_fit.parameter
    temp_array_new = 400.+time_array*1000.
    temp_array_interp = sp.interpolate.interp1d(time_array, temp_array_new,
        bounds_error=False, fill_value=temp_array_new[-1])
    #print parameters[0]
    y = c2sm_fit.calcMassC2SM(parameter,init_mass=0,time=time_array,temp=temp_array_interp)
    assert y.size == time_array.size  # check that y size is equal to time
    assert all(y >= 0) and all(y <=1) # check that y is bounded between 0 and 1