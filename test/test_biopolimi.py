"""Test bioPolimi model."""
from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals

import pytest
import pkp.biopolimi
import numpy as np

ultimate_analysis = {"C": 1, "H": 8.33e-3, "O": 0.8, "N": 0, "S": 0}

proximate_analysis = {"FC": 0.1, "VM": 0.7, "Ash": 0.2, "Moist": 0.2}


@pytest.fixture
def bio():
    """Create a BioPolimi object."""
    return pkp.biopolimi.BioPolimi(
        proximate_analysis=proximate_analysis,
        ultimate_analysis=pkp.biopolimi.bioS1.ultimate_analysis,
        pressure=101325,
        name="biomass",
    )


def test_init(bio):
    """Test initialization."""
    print(bio.composition)
    assert "CELL" in bio.mechanism.species_names


def test_rate(bio):
    """Test rate."""
    T = 1000
    y = np.append(bio.y0, T)
    rate = bio.rate(0, y)

    assert len(rate) == len(y) - 1
