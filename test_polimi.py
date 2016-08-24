from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals

import pkp.polimi
from test_cpd import ua, pa, op_cond
import matplotlib.pyplot as plt
import itertools

try:
    plt.style.use(['mystyle'])
except:
    plt.style.use(['ggplot'])


def init_polimi():
    return pkp.polimi.Polimi(ultimate_analysis=ua,
                             proximate_analysis=pa,
                             pressure=101325,
                             name='Polimi test',
                             mechanism='COAL.xml')


def test_init():
    coal = init_polimi()


def test_triangle():
    triangle = pkp.polimi.Triangle()
    assert triangle.is_inside([0.2, 0.2])
    assert not triangle.is_inside([2, 2])


def test_reference_coal():
    fig, ax = plt.subplots()
    ref_coals = [pkp.polimi.coal1, pkp.polimi.coal2,
                 pkp.polimi.coal3, pkp.polimi.char]

    triangle_lower = pkp.polimi.Triangle(
        x0=pkp.polimi.char.van_kravelen,
        x1=pkp.polimi.coal2.van_kravelen,
        x2=pkp.polimi.coal3.van_kravelen)

    coal_test_true = [0.2, 0.6]
    coal_test_false = [0.3, 0.2]
    ax.plot(coal_test_true[0], coal_test_true[1], marker='<',
            markersize=10, color='red')
    ax.plot(coal_test_false[0], coal_test_false[1], marker='<',
            markersize=10, color='blue')
    for c in ref_coals:
        # ax.plot(c.van_kravelen[0], c.van_kravelen[1],
        #        marker='o', markersize=10, label=c.name,
        #        color='black')
        ax.annotate(c.name, xy=c.van_kravelen)
    for couple in itertools.combinations(ref_coals, 2):
        x = [c.van_kravelen[0] for c in couple]
        y = [c.van_kravelen[1] for c in couple]
        ax.plot(x, y, color='black', marker='o', markersize=10)
    ax.set_xlabel('O:C')
    ax.set_ylabel('H:C')
    ax.set_xlim([-0.05, 0.45])
    ax.set_ylim([-0.05, 1.05])
    # ax.legend(loc='best')

    fig.savefig('van_kravelen.png')
    plt.close(fig)

    assert triangle_lower.is_inside(coal_test_true)
    assert not triangle_lower.is_inside(coal_test_false)
