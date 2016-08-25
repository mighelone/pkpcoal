import pytest
import pkp
from dictinput import settings

yml_file = 'input.yml'


@pytest.fixture
def conf():
    return pkp.ReadConfiguration(settings)


def test_coal(conf):
    for key in ('proximate_analysis', 'ultimate_analysis', 'rho_dry'):
        assert settings['Coal'][key] == getattr(conf, key)
    assert settings['Coal']['HHV'] == conf.HHV / 1e6

    op_cond = conf.operating_conditions
