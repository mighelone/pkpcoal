from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals

import pytest
import pkp.biopolimi

ultimate_analysis = {'C': 1,
                     'H': 8.33e-3,
                     'O': 0.8,
                     'N': 0,
                     'S': 0}

proximate_analysis = {'FC': 0.1,
                      'VM': 0.7,
                      'Ash': 0.2,
                      'Moist': 0.2}


@pytest.fixture
def bio():
    return pkp.biopolimi.BioPolimi(
        proximate_analysis=proximate_analysis,
        ultimate_analysis=pkp.biopolimi.bioS1.ultimate_analysis,
        pressure=101325,
        name='biomass')


def test_init(bio):
    print(bio.composition)
    assert 'CELL' in bio.mechanism.species_names
