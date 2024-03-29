"""Correlation data for CPD."""

import numpy as np

CPD_CORRELATION = np.array([[0.0, 0.0, 0.0, 0.0],
                            [421.957, 1301.41, 0.489809, -52.1054],
                            [-8.64692, 16.3879, -0.00981566, 1.63872],
                            [0.0463894, -0.187493, 0.000133046,
                             -0.0107548],
                            [-8.47272, -454.773, 0.155483, -1.23688],
                            [1.18173, 51.7109, -0.0243873, 0.0931937],
                            [1.15366, -10.0720, 0.00705248, -0.165673],
                            [-0.0434024, 0.0760827, 0.000219163,
                             0.00409556],
                            [0.556772, 1.36022, -0.0110498, 0.00926097],
                            [-0.00654575, -0.0313561, 0.000100939,
                             -0.0000826717]], dtype=object)

xx = np.array([3.4, 3.2, 3., 2.8, 2.6, 2.4, 2.2, 2., 1.8,
               1.6, 1.4, 1.2, 1., .8, .6, .4, .2, 0.])
yy = np.array([.9997, .9993, .9987, .9974, .9953, .9918, .9861, .9772,
               .9641, .9452, .9192, .8849, .8413, .7881, .7257, .6554,
               .5793, .5])
xx = xx[::-1]
yy = yy[::-1]

n_ref_coals = 12
n_gas_species = 4
x_gas = np.array([[0., .04, .11, .14, .21, .27, .34, .675, .9, 1.],
                  [.0, .161, .442, .663, .777, .874, .921, .967, 1.],
                  [.0, .022, .20, .430, .526, .64, .787, .875, .927, .955, 1.],
                  [.0, .04, .12, .15, .23, .29, .36, .68, .9, 1.],
                  [.0, .018, .058, .21, .417, .572, .696, .778, .821, .883,
                   .932, 1.],
                  [.0, .052, .144, .291, .498, .639, .746, .859, .925, .949,
                   .966, 1.],
                  [.0, .063, .178, .33, .506, .612, .706, .813, .895, .94, 1.],
                  [.0, .04, .12, .15, .23, .29, .36, .68, .9, 1.],
                  [.0, .061, .146, .374, .535, .622, .714, .8, .883, .931,
                   .964, 1.],
                  [.0, .034, .087, .179, .316, .472, .585, .694, .777, .872,
                   .935, 1.],
                  [.0, .04, .12, .16, .25, .31, .37, .68, .9, 1.],
                  [.0, .02, .055, .17, .313, .434, .546, .716, .874, .935,
                   .973, 1.]], dtype=object)


y_gas = np.array([[[.772, .772, .738, .455, .371, .304, .290, .273, .218,
                    .218],
                   [.699, .632, .299, .269, .247, .249, .236, .225, .226],
                   [.0, .0, .35, .297, .301, .299, .284, .291, .306, .297,
                    .283],
                   [.636, .636, .646, .550, .436, .320, .186, .199, .195,
                    .195],
                   [1., .983, .754, .488, .413, .385, .373, .382, .377,
                       0.362, .367, .348],
                   [.665, .636, .604, .508, .435, .409, .383, .362, .351, .343,
                    .342, .339],
                   [.763, .737, .698, .572, .527, .470, .438, .411, .411, .396,
                    .378],
                   [.748, .748, .637, .704, .490, .446, .348, .268, .266,
                    .266],
                   [.0, .0, .385, .461, .396, .369, .344, .323, .292, .277,
                    .266, .257],
                   [.0, .0, .197, .267, .26, .333, .361, .369, .346, .306,
                    .285, .267],
                   [.521, .521, .55, .523, .511, .46, .414, .388, .313, 0.313],
                   [.0, .0, .291, .335, .264, .271, .261, .211, .171, .160,
                    .153, .149]],
                  [[.0, .0, .0, .174, .174, .167, .129, .102, .071, .071],
                   [.259, .234, .113, .086, .097, .109, .116, .118, .122],
                   [.333, .327, .070, .052, .057, .06, .059, .062, .066, .08,
                    0.115],
                   [.194, .194, .152, .117, .116, .122, .081, .092, .065,
                    .065],
                   [.0, .0, .0, .122, .103, .086, .083, .082, .085, .086, .093,
                    .128],
                   [.332, .318, .165, .141, .120, .108, .105, .119, .120, .122,
                    .125, .130],
                   [.229, .221, .125, .09, .07, .073, .083, .133, .132, .13,
                    .147],
                   [.111, .111, .142, .175, .149, .155, .136, .122, .133,
                    .133],
                   [.98, .984, .55, .345, .317, .285, .286, .277, .273, .264,
                    .254, .255],
                   [.993, .989, .786, .572, .519, .416, .375, .345, .335, .32,
                    .303, .299],
                   [.363, .363, .353, .325, .321, .35, .318, .251, .249, .249],
                   [1., .983, .448, .179, .104, .09, .104, .151, .166, .160,
                    .158, .154]],
                  [[.203, .203, .078, .160, .180, .219, .258, .294, .320,
                    .320],
                   [.041, .037, .388, .389, .359, .332, .323, .307, .299],
                   [.667, .655, .42, .454, .444, .419, .382, .353, .331, .321,
                    .306],
                   [.055, .055, .073, .088, .116, .124, .170, .15, .189, .189],
                   [.0, .0, .188, .195, .234, .243, .224, .21, .2, .186, .177,
                    .167],
                   [.0, .0, .11, .155, .176, .172, .185, .173, .163, .159,
                    .156, .151],
                   [.0, .0, .075, .136, .159, .178, .174, .157, .143, .141,
                    .132],
                   [.02, .02, .026, .042, .045, .049, .064, .1, .128, .128],
                   [.0, .0, .0, .029, .048, .067, .069, .072, .069, .066, .063,
                    .061],
                   [.0, .0, .0, .0, .035, .05, .061, 0.058, .057, .053, .049,
                    .046],
                   [.01, .01, .011, .016, .011, .021, .023, .035, .06, .06],
                   [.0, .0, .216, .262, .362, .327, .307, .25, .203, .189,
                    .182, .177]],
                  [[.0, .0, .157, .121, .141, .112, .139, .085, .145, .145],
                   [.0, .0, .0, .057, .097, .109, .124, .15, .153],
                   [.0, .0, .0, .0, .0, .024, .078, .097, .099, .104, .099],
                   [.083, .083, .038, .066, .032, .168, .286, .324, .313,
                    .313],
                   [.0, .0, .0, .0, .055, .091, .124, .131, .142, .171, .168,
                    .162],
                   [.0, .0, .0, .028, .093, .129, .142, .162, .181, .191, .193,
                    .195],
                   [.0, .0, .0, .075, .099, .122, .139, .133, .148, .167,
                    .177],
                   [.101, .101, .173, .054, .219, .247, .335, .349, .28, .280],
                   [.0, .0, .055, .115, .151, .168, .172, .2, .236, .264, .287,
                    .298],
                   [.0, .0, .0, .133, .142, .150, .15, .173, .206, .265, .307,
                    .331],
                   [.096, .096, .066, .113, .123, .13, .2, .281, .334, .334],
                   [.0, .0, .0, .084, .078, .115, .130, .191, .262, .294, .311,
                    .322]]], dtype=object)

gas_species = ['H2O', 'CO2', 'CH4', 'CO']
