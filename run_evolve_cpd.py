from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals

from deap import algorithms
from deap import base
from deap import creator
from deap import tools

import array
import random
import numpy as np

import pkp
import pkp.empirical_model

import matplotlib.pyplot as plt
plt.style.use(['mystyle'])

runner = pkp.PKPRunner('input.yml')
runner.operating_conditions['runs'] = 1
runner.Polimi['active'] = False

res = runner.run(results_dir='./optdir')

res0 = res['CPD']['run0']
y_cpd = np.array(res0['fsolid'])
# t = np.array(res0['time(ms)']) * 1e-3
t = np.array(res0.index) * 1e-3
operating_conditions = runner.operating_conditions['run0']

# set a SFOR simulation
sfor_parameters = {'A': 72e3,
                   'E': 55e6,
                   'y0': 0.5}

m = pkp.empirical_model.SFOR(parameters=sfor_parameters)
m.operating_conditions = operating_conditions

_, y = m.run(t=t)

fig, ax = plt.subplots()
ax.plot(t, y_cpd, label='CPD')
ax.plot(t, y, label='SFOR')
ax.set_xlabel('t, s')
ax.set_ylabel('y, daf')
ax.legend(loc='best')

# plt.show()


def error(parameters):
    m = pkp.empirical_model.SFOR(
        parameters=parameters)
    m.operating_conditions = operating_conditions
    _, y = m.run(t=t)
    return np.mean(np.power((y - y_cpd), 2))

print('Error: {}'.format(error(sfor_parameters)))
