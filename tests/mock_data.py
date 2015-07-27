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
            'kBounds': [1e-6, 0.99],
            'tstart':0.0,
            'tstartBounds': [0.0, 1.0],
            'finalYield': 1.0,
            'finalYieldBounds': [0.0, 1.0]},
     'arrheniusRate': {
            'preExp':1.0,
            'preExpBounds': [0.0, 1.0e10],
            'activationEnergy': 1.0,
            'activationEnergyBounds': [0.0, 1.0e6]},
     'C2SM': {
        'alpha1': 0.2,
        'alpha1Bounds': [0.0, 0.4],
        'A1': 1000.0,
        'A1Bounds': [0, 10000.],
        'E1': 5000.0,
        'E1Bounds': [1000.0, 10000.0],
        'alpha2': 0.8,
        'alpha2Bounds': [0.4, 1],
        'A2': 8000.0,
        'A2Bounds': [5000, 30000.],
        'E2': 9000.0,
        'E2Bounds': [1500.0, 20000.0]
     }
}

mock_ua = {"Carbon": 69.20, "Hydrogen": 4.40,
    "Nitrogen": 1.42, "Oxygen": 9.98, "Sulphur": 0.58}

mock_pa = {"Fixed Carbon": 50.80, "Volatile Matter": 34.80, 
           "Moisture": 5.80, "Ash": 8.60}

mock_cpd = {
    'fit'    : 'constantRate',
    'active' : 'True',
    'deltaT' : '1e-08',
    'MW_TAR' : '130',
}


T_init = [0.0, 300.0]
mock_opcond = {
    'run' + str(nr): [T_init, [0.000005*i, 2500]] for nr,i in enumerate(range(1,40,1))
}

mock_opcond['pressure'] = 1.0
mock_opcond['runs'] = 1

mock_input={'Coal': {
    'Proximate Analysis': mock_pa,
    'Ultimate Analysis': mock_ua,
    'hhv': 28220000,
    'MW_PS': 40,
    'OperatingConditions': mock_opcond,
    },
    'CPD': mock_cpd,
}

methane_as_coal = {'Coal': {
    'Ultimate Analysis': {"Carbon": 0.75 , "Hydrogen": 0.25,
        "Nitrogen": 0.0, "Oxygen": 0.0, "Sulphur": 0.0},
    'Proximate Analysis': {"Fixed Carbon": 0.0, "Volatile Matter": 100, 
        "Moisture": 0.0, "Ash": 0.0},
    'hhv': 50016000,
    'MW_PS': 16,
    'OperatingConditions': mock_opcond,
    },
}

co_as_coal = {'Coal': {
    'Ultimate Analysis': {"Carbon": 12.0/28.0 , "Hydrogen": 0.0,
        "Nitrogen": 0.0, "Oxygen": 16.0/28.0, "Sulphur": 0.0},
    'Proximate Analysis': {"Fixed Carbon": 0.0, "Volatile Matter": 100, 
        "Moisture": 0.0, "Ash": 0.0},
    'hhv': 10159968,
    'MW_PS': 28.0,
    'OperatingConditions': mock_opcond,
    },
}

propane_as_coal = {'Coal': {
    'Ultimate Analysis': {"Carbon": 36.0/44.0 , "Hydrogen": 1.0-36.0/44.0,
        "Nitrogen": 0.0, "Oxygen": 0.0, "Sulphur": 0.0},
    'Proximate Analysis': {"Fixed Carbon": 0.0, "Volatile Matter": 100, 
        "Moisture": 0.0, "Ash": 0.0},
    'hhv': 50016000,
    'MW_PS': 44.0,
    'OperatingConditions': mock_opcond,
    },
}

mock_ps  = {'Coal': {
    'Ultimate Analysis': {"Carbon": 0.5 , "Hydrogen": 0.5,
        "Nitrogen": 0.0, "Oxygen": 0.0, "Sulphur": 0.0},
    'Proximate Analysis': {"Fixed Carbon": 0.5, "Volatile Matter": 0.5, 
        "Moisture": 0.0, "Ash": 0.0},
    'hhv': 50016000,
    'MW_PS': 44.0,
    'OperatingConditions': mock_opcond,
    },
}

opt_params_bounded = copy.deepcopy(opt_params)
opt_params_bounded['constantRate'].update(
    {'tstartBounds': (1e-6, 1e-5), 'finalYieldBounds': (1.0, 1.0)})

time_array   = np.arange(0.01, 1.011, 0.01)
time_array_ms= np.arange(0.01, 1.011, 0.01)/1000.0
linear_yield = np.arange(0.0, 1.001, 0.01)
exp_yield    = 1.0 - np.exp(-time_array*5)
temp_array   = time_array

cpd_header   = {"file": ["time(ms)", "Species", "temp"]}
cpd_linear   = {"file": [time_array_ms, linear_yield, temp_array]}
cpd_exp      = {"file": [time_array_ms, exp_yield, temp_array]}
