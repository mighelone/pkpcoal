"""Correlation data for CPD."""

import numpy as np

CPD_CORRELATION = np.array(
    [
        [0.0, 0.0, 0.0, 0.0],
        [421.957, 1301.41, 0.489809, -52.1054],
        [-8.64692, 16.3879, -0.00981566, 1.63872],
        [0.0463894, -0.187493, 0.000133046, -0.0107548],
        [-8.47272, -454.773, 0.155483, -1.23688],
        [1.18173, 51.7109, -0.0243873, 0.0931937],
        [1.15366, -10.0720, 0.00705248, -0.165673],
        [-0.0434024, 0.0760827, 0.000219163, 0.00409556],
        [0.556772, 1.36022, -0.0110498, 0.00926097],
        [-0.00654575, -0.0313561, 0.000100939, -0.0000826717],
    ]
)

xx = np.array(
    [
        3.4,
        3.2,
        3.0,
        2.8,
        2.6,
        2.4,
        2.2,
        2.0,
        1.8,
        1.6,
        1.4,
        1.2,
        1.0,
        0.8,
        0.6,
        0.4,
        0.2,
        0.0,
    ]
)
yy = np.array(
    [
        0.9997,
        0.9993,
        0.9987,
        0.9974,
        0.9953,
        0.9918,
        0.9861,
        0.9772,
        0.9641,
        0.9452,
        0.9192,
        0.8849,
        0.8413,
        0.7881,
        0.7257,
        0.6554,
        0.5793,
        0.5,
    ]
)
xx = xx[::-1]
yy = yy[::-1]

n_ref_coals = 12
n_gas_species = 4
x_gas = np.array(
    [
        [0.0, 0.04, 0.11, 0.14, 0.21, 0.27, 0.34, 0.675, 0.9, 1.0],
        [0.0, 0.161, 0.442, 0.663, 0.777, 0.874, 0.921, 0.967, 1.0],
        [0.0, 0.022, 0.20, 0.430, 0.526, 0.64, 0.787, 0.875, 0.927, 0.955, 1.0],
        [0.0, 0.04, 0.12, 0.15, 0.23, 0.29, 0.36, 0.68, 0.9, 1.0],
        [0.0, 0.018, 0.058, 0.21, 0.417, 0.572, 0.696, 0.778, 0.821, 0.883, 0.932, 1.0],
        [
            0.0,
            0.052,
            0.144,
            0.291,
            0.498,
            0.639,
            0.746,
            0.859,
            0.925,
            0.949,
            0.966,
            1.0,
        ],
        [0.0, 0.063, 0.178, 0.33, 0.506, 0.612, 0.706, 0.813, 0.895, 0.94, 1.0],
        [0.0, 0.04, 0.12, 0.15, 0.23, 0.29, 0.36, 0.68, 0.9, 1.0],
        [0.0, 0.061, 0.146, 0.374, 0.535, 0.622, 0.714, 0.8, 0.883, 0.931, 0.964, 1.0],
        [
            0.0,
            0.034,
            0.087,
            0.179,
            0.316,
            0.472,
            0.585,
            0.694,
            0.777,
            0.872,
            0.935,
            1.0,
        ],
        [0.0, 0.04, 0.12, 0.16, 0.25, 0.31, 0.37, 0.68, 0.9, 1.0],
        [0.0, 0.02, 0.055, 0.17, 0.313, 0.434, 0.546, 0.716, 0.874, 0.935, 0.973, 1.0],
    ]
)


y_gas = np.array(
    [
        [
            [0.772, 0.772, 0.738, 0.455, 0.371, 0.304, 0.290, 0.273, 0.218, 0.218],
            [0.699, 0.632, 0.299, 0.269, 0.247, 0.249, 0.236, 0.225, 0.226],
            [0.0, 0.0, 0.35, 0.297, 0.301, 0.299, 0.284, 0.291, 0.306, 0.297, 0.283],
            [0.636, 0.636, 0.646, 0.550, 0.436, 0.320, 0.186, 0.199, 0.195, 0.195],
            [
                1.0,
                0.983,
                0.754,
                0.488,
                0.413,
                0.385,
                0.373,
                0.382,
                0.377,
                0.362,
                0.367,
                0.348,
            ],
            [
                0.665,
                0.636,
                0.604,
                0.508,
                0.435,
                0.409,
                0.383,
                0.362,
                0.351,
                0.343,
                0.342,
                0.339,
            ],
            [
                0.763,
                0.737,
                0.698,
                0.572,
                0.527,
                0.470,
                0.438,
                0.411,
                0.411,
                0.396,
                0.378,
            ],
            [0.748, 0.748, 0.637, 0.704, 0.490, 0.446, 0.348, 0.268, 0.266, 0.266],
            [
                0.0,
                0.0,
                0.385,
                0.461,
                0.396,
                0.369,
                0.344,
                0.323,
                0.292,
                0.277,
                0.266,
                0.257,
            ],
            [
                0.0,
                0.0,
                0.197,
                0.267,
                0.26,
                0.333,
                0.361,
                0.369,
                0.346,
                0.306,
                0.285,
                0.267,
            ],
            [0.521, 0.521, 0.55, 0.523, 0.511, 0.46, 0.414, 0.388, 0.313, 0.313],
            [
                0.0,
                0.0,
                0.291,
                0.335,
                0.264,
                0.271,
                0.261,
                0.211,
                0.171,
                0.160,
                0.153,
                0.149,
            ],
        ],
        [
            [0.0, 0.0, 0.0, 0.174, 0.174, 0.167, 0.129, 0.102, 0.071, 0.071],
            [0.259, 0.234, 0.113, 0.086, 0.097, 0.109, 0.116, 0.118, 0.122],
            [0.333, 0.327, 0.070, 0.052, 0.057, 0.06, 0.059, 0.062, 0.066, 0.08, 0.115],
            [0.194, 0.194, 0.152, 0.117, 0.116, 0.122, 0.081, 0.092, 0.065, 0.065],
            [
                0.0,
                0.0,
                0.0,
                0.122,
                0.103,
                0.086,
                0.083,
                0.082,
                0.085,
                0.086,
                0.093,
                0.128,
            ],
            [
                0.332,
                0.318,
                0.165,
                0.141,
                0.120,
                0.108,
                0.105,
                0.119,
                0.120,
                0.122,
                0.125,
                0.130,
            ],
            [0.229, 0.221, 0.125, 0.09, 0.07, 0.073, 0.083, 0.133, 0.132, 0.13, 0.147],
            [0.111, 0.111, 0.142, 0.175, 0.149, 0.155, 0.136, 0.122, 0.133, 0.133],
            [
                0.98,
                0.984,
                0.55,
                0.345,
                0.317,
                0.285,
                0.286,
                0.277,
                0.273,
                0.264,
                0.254,
                0.255,
            ],
            [
                0.993,
                0.989,
                0.786,
                0.572,
                0.519,
                0.416,
                0.375,
                0.345,
                0.335,
                0.32,
                0.303,
                0.299,
            ],
            [0.363, 0.363, 0.353, 0.325, 0.321, 0.35, 0.318, 0.251, 0.249, 0.249],
            [
                1.0,
                0.983,
                0.448,
                0.179,
                0.104,
                0.09,
                0.104,
                0.151,
                0.166,
                0.160,
                0.158,
                0.154,
            ],
        ],
        [
            [0.203, 0.203, 0.078, 0.160, 0.180, 0.219, 0.258, 0.294, 0.320, 0.320],
            [0.041, 0.037, 0.388, 0.389, 0.359, 0.332, 0.323, 0.307, 0.299],
            [
                0.667,
                0.655,
                0.42,
                0.454,
                0.444,
                0.419,
                0.382,
                0.353,
                0.331,
                0.321,
                0.306,
            ],
            [0.055, 0.055, 0.073, 0.088, 0.116, 0.124, 0.170, 0.15, 0.189, 0.189],
            [
                0.0,
                0.0,
                0.188,
                0.195,
                0.234,
                0.243,
                0.224,
                0.21,
                0.2,
                0.186,
                0.177,
                0.167,
            ],
            [
                0.0,
                0.0,
                0.11,
                0.155,
                0.176,
                0.172,
                0.185,
                0.173,
                0.163,
                0.159,
                0.156,
                0.151,
            ],
            [0.0, 0.0, 0.075, 0.136, 0.159, 0.178, 0.174, 0.157, 0.143, 0.141, 0.132],
            [0.02, 0.02, 0.026, 0.042, 0.045, 0.049, 0.064, 0.1, 0.128, 0.128],
            [
                0.0,
                0.0,
                0.0,
                0.029,
                0.048,
                0.067,
                0.069,
                0.072,
                0.069,
                0.066,
                0.063,
                0.061,
            ],
            [0.0, 0.0, 0.0, 0.0, 0.035, 0.05, 0.061, 0.058, 0.057, 0.053, 0.049, 0.046],
            [0.01, 0.01, 0.011, 0.016, 0.011, 0.021, 0.023, 0.035, 0.06, 0.06],
            [
                0.0,
                0.0,
                0.216,
                0.262,
                0.362,
                0.327,
                0.307,
                0.25,
                0.203,
                0.189,
                0.182,
                0.177,
            ],
        ],
        [
            [0.0, 0.0, 0.157, 0.121, 0.141, 0.112, 0.139, 0.085, 0.145, 0.145],
            [0.0, 0.0, 0.0, 0.057, 0.097, 0.109, 0.124, 0.15, 0.153],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.024, 0.078, 0.097, 0.099, 0.104, 0.099],
            [0.083, 0.083, 0.038, 0.066, 0.032, 0.168, 0.286, 0.324, 0.313, 0.313],
            [
                0.0,
                0.0,
                0.0,
                0.0,
                0.055,
                0.091,
                0.124,
                0.131,
                0.142,
                0.171,
                0.168,
                0.162,
            ],
            [
                0.0,
                0.0,
                0.0,
                0.028,
                0.093,
                0.129,
                0.142,
                0.162,
                0.181,
                0.191,
                0.193,
                0.195,
            ],
            [0.0, 0.0, 0.0, 0.075, 0.099, 0.122, 0.139, 0.133, 0.148, 0.167, 0.177],
            [0.101, 0.101, 0.173, 0.054, 0.219, 0.247, 0.335, 0.349, 0.28, 0.280],
            [
                0.0,
                0.0,
                0.055,
                0.115,
                0.151,
                0.168,
                0.172,
                0.2,
                0.236,
                0.264,
                0.287,
                0.298,
            ],
            [
                0.0,
                0.0,
                0.0,
                0.133,
                0.142,
                0.150,
                0.15,
                0.173,
                0.206,
                0.265,
                0.307,
                0.331,
            ],
            [0.096, 0.096, 0.066, 0.113, 0.123, 0.13, 0.2, 0.281, 0.334, 0.334],
            [
                0.0,
                0.0,
                0.0,
                0.084,
                0.078,
                0.115,
                0.130,
                0.191,
                0.262,
                0.294,
                0.311,
                0.322,
            ],
        ],
    ]
)

gas_species = ["H2O", "CO2", "CH4", "CO"]
