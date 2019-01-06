"""Test module for CPD."""
from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals
import pytest
import numpy as np

import pkp.cpd

ua = {"C": 74.12, "H": 4.96, "O": 13.18, "N": 1.45, "S": 0.0}

ua["S"] = 100 - sum(ua.values())

daf = 90.0
vm_daf = 43.37
pa = {
    "FC": (100 - vm_daf) / 100 * daf,
    "VM": vm_daf / 100 * daf,
    "Ash": 100 - daf,
    "Moist": 0,
}

pressure = 101325


@pytest.fixture
def cpd():
    """Init the CPD model."""
    return pkp.cpd.CPD(
        ultimate_analysis=ua, proximate_analysis=pa, pressure=101325, name="CPD coal"
    )


def test_cpd_init(cpd):
    """
    Test if the property are correctly initialized.

    Reference values are taken from
    https://www.et.byu.edu/~tom/cpd/correlation.html
    """
    mdel = 42.4
    mw = 383.0
    p0 = 0.501
    sigma = 5.20
    c0 = 0.00952

    # adjust mdel as in the fortran code
    # this is done in cpd._after_set_NMR
    # mdel /= (1-c0)
    # mdel -= 7

    np.testing.assert_almost_equal(cpd.ultimate_analysis["C"], 0.7412)
    np.testing.assert_almost_equal(cpd.ultimate_analysis["H"], 0.0496)

    assert cpd.pressure == pressure
    np.testing.assert_almost_equal(mdel, cpd.mdel, decimal=1)
    np.testing.assert_almost_equal(mw, cpd.mw, decimal=2)
    np.testing.assert_almost_equal(p0, cpd.p0, decimal=2)
    np.testing.assert_almost_equal(sigma, sigma, decimal=2)
    np.testing.assert_almost_equal(c0, cpd.c0, decimal=2)

    # check if the delta mw is correctly corrected
    mdel_corr = cpd.mdel / (1 - cpd.c0)
    mdel_corr -= 7
    ma = cpd.mw - cpd.sig * mdel_corr
    np.testing.assert_almost_equal(ma, cpd.ma)

    assert cpd.name == "CPD coal"


def test_set_parameters(cpd):
    """Test set parameters."""
    parameters = cpd.get_parameters()

    assert all(p in parameters for p in cpd.nmr_parameters + cpd.kin_parameters)

    mw = 360
    ab = parameters["ab"] * 2
    cpd.set_parameters(mw=mw, ab=ab)

    assert cpd.mw == mw
    assert cpd.ab == ab

    # check that afer set NMR was correctly executed
    mdel_corr = cpd.mdel / (1 - cpd.c0)
    mdel_corr -= 7
    ma = cpd.mw - cpd.sig * mdel_corr

    assert cpd.ma == ma
