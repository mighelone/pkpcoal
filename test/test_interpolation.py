from pkp.interpolate import interp
import numpy as np


def test_ip():
    x = np.linspace(0, 1, 101)
    y = x ** 2

    xi = 0.3456

    np.testing.assert_almost_equal(interp(xi, x, y), np.interp(xi, x, y))
