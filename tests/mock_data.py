import copy 
import numpy as np

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
