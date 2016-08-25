from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals

from pkp.empirical_model import SFOR
import numpy as np
import matplotlib.pyplot as plt
plt.style.use(['mystyle'])

parameters = {'A': 72e3,
              'E': 55e6,
              'y0': 0.5}

sfor = SFOR(parameters=parameters)
print(sfor.parameters)
sfor.operating_conditions = [[0, 400],
                             [0.005, 1400],
                             [0.02, 1400]
                             ]
t, y = sfor.run()

fig, ax = plt.subplots()

ax.plot(t, y, marker='x')

t_max = sfor.operating_conditions[-1, 0]
t0 = np.linspace(0, t_max, 100)

_, y0 = sfor.run(t=t0)

ax.plot(t0, y0, marker='x')
plt.show()
