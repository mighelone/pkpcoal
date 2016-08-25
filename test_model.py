from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals

import pytest
import pkp.models
import numpy as np


@pytest.fixture
def sfor():
    return pkp.models.SFOR()


def test_init_model(sfor):
    assert sfor.parameters_default == [1e5, 50e6, 0.6]
    assert sfor.len_parameters == 3
    assert sfor.parameters_names == ['A', 'E', 'Y0']

    # set parameters as list
    par_list = [1e6, 100e6, 0.5]
    sfor.parameters = par_list
    for key, value in sfor.parameters.iteritems():
        assert value == par_list[sfor.parameters_names.index(key)]

    # set as np.ndarray
    sfor.parameters = np.array(par_list)
    for key, value in sfor.parameters.iteritems():
        assert value == par_list[sfor.parameters_names.index(key)]

    # set as dict
    sfor.parameters = {k: v for k, v in zip(
        sfor.parameters_names, par_list)}
    for key, value in sfor.parameters.iteritems():
        assert value == par_list[sfor.parameters_names.index(key)]
