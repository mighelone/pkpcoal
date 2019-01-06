"""Test module for Polimi model."""

from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals

import pkp.polimi
import matplotlib.pyplot as plt
import itertools
import pytest
import numpy as np

ua0 = {"C": 1, "H": 0, "O": 0, "N": 0, "S": 0}

ua = {"C": 69, "H": 5, "O": 24.7, "N": 0.8, "S": 0.5}

pa = {"FC": 45.1, "VM": 50.6, "Ash": 4.3, "Moist": 19.0}

op_cond = [[0, 500], [0.001, 1400], [0.01, 1400]]


@pytest.fixture
def coal():
    """Init a Polimi instance."""
    return pkp.polimi.Polimi(
        ultimate_analysis=ua, proximate_analysis=pa, pressure=101325, name="Polimi test"
    )


@pytest.fixture
def triangle():
    """Define a triangle object."""
    return pkp.polimi.Triangle()


def test_init(coal):
    """Test initialization of Polimi model."""
    assert coal.name == "Polimi test"
    assert hasattr(coal, "van_kravelen")

    # check light_gas index
    mech = coal.mechanism
    light_gas = [mech.species_index(sp) for sp in coal.light_gas]
    assert light_gas == coal._light_gas_index


def test_triangle(triangle):
    """Test triangle."""
    point_inside = [0.2, 0.2]
    assert triangle.is_inside(point_inside)
    assert not triangle.is_inside([2, 2])

    for i, point in enumerate(triangle):
        assert triangle.weights(point)[i] == 1

    # assert (triangle.weights(point_inside) == 0).all()


def test_reference_coal():
    """Test reference coal."""
    fig, ax = plt.subplots()
    ref_coals = [pkp.polimi.coal1, pkp.polimi.coal2, pkp.polimi.coal3, pkp.polimi.char]

    triangle_lower = pkp.polimi.Triangle(
        x0=pkp.polimi.char.van_kravelen,
        x1=pkp.polimi.coal2.van_kravelen,
        x2=pkp.polimi.coal3.van_kravelen,
    )

    coal_test_true = [0.2, 0.6]
    coal_test_false = [0.3, 0.2]
    ax.plot(
        coal_test_true[0], coal_test_true[1], marker="<", markersize=10, color="red"
    )
    ax.plot(
        coal_test_false[0], coal_test_false[1], marker="<", markersize=10, color="blue"
    )
    for c in ref_coals:
        # ax.plot(c.van_kravelen[0], c.van_kravelen[1],
        #        marker='o', markersize=10, label=c.name,
        #        color='black')
        ax.annotate(c.name, xy=c.van_kravelen)
    for couple in itertools.combinations(ref_coals, 2):
        x = [c.van_kravelen[0] for c in couple]
        y = [c.van_kravelen[1] for c in couple]
        ax.plot(x, y, color="black", marker="o", markersize=10)
    ax.set_xlabel("O:C")
    ax.set_ylabel("H:C")
    ax.set_xlim([-0.05, 0.45])
    ax.set_ylim([-0.05, 1.05])
    # ax.legend(loc='best')

    fig.savefig("van_kravelen.png")
    plt.close(fig)

    assert triangle_lower.is_inside(coal_test_true)
    assert not triangle_lower.is_inside(coal_test_false)


def test_rate(coal):
    """Test rate calculation of Polimi."""
    T = 1000
    y = np.append(coal.y0, T)
    rate = coal.rate(0, y)

    mech = coal.mechanism
    mech.TPY = T, coal.pressure, coal.y0
    rate_c = mech.net_production_rates * mech.molecular_weights / mech.density

    np.testing.assert_allclose(rate, rate_c)

    assert coal.get_yield(0, coal.y0) == 0
