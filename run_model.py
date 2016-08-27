import pkp.empirical_model
import numpy as np
par_min = np.array([1e3, 10e6, 0.4])
par_max = np.array([1e8, 200e6, 0.6])
scal_par = np.array([0.2, 0.5, 0.4])

unsc_par = pkp.empirical_model.SFOR.unscale_parameters(
    scal_par, par_min, par_max)
