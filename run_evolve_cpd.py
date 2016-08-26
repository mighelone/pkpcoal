from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals

from deap import algorithms
from deap import base
from deap import creator
from deap import tools

import random
import numpy as np

import pkp
import pkp.empirical_model

import matplotlib.pyplot as plt
plt.style.use(['mystyle'])

# calc CPD solution
runner = pkp.PKPRunner('input.yml')
runner.operating_conditions['runs'] = 1
runner.Polimi['active'] = False
res = runner.run(results_dir='./optdir')
res0 = res['CPD']['run0']
y_cpd = np.array(res0['fsolid'])
# t = np.array(res0['time(ms)']) * 1e-3
t_cpd = np.array(res0.index) * 1e-3
operating_conditions = runner.operating_conditions['run0']

# set a SFOR simulation
# sfor_parameters = {'A': 72e3,
#                   'E': 55e6,
#                   'y0': 0.5}

# m = pkp.empirical_model.SFOR(parameters=sfor_parameters)
# m.operating_conditions = operating_conditions

# _, y = m.run(t=t)

# fig, ax = plt.subplots()
# ax.plot(t, y_cpd, label='CPD')
# ax.plot(t, y, label='SFOR')
# ax.set_xlabel('t, s')
# ax.set_ylabel('y, daf')
# ax.legend(loc='best')

parameters_min = np.array([2., 10, 0.3])  # logA, E MJ/kg, y0
parameters_max = np.array([5., 200, 0.7])
delta_parameters = parameters_max - parameters_min

parameters_0 = [3, 50, 0.5]


def scale_parameters(parameters):
    par_array = np.array(parameters)
    par_scaled = parameters_min + par_array * delta_parameters
    par_scaled[0] = np.power(10.0, par_scaled[0])
    par_scaled[1] = par_scaled[1] * 1e6

    return par_scaled


def error(parameters):
    # scale parameters and convert log
    par_scaled = scale_parameters(parameters)
    m = pkp.empirical_model.SFOR(
        parameters=par_scaled)
    m.operating_conditions = operating_conditions
    _, y = m.run(t=t_cpd)
    return np.mean(np.power((y - y_cpd), 2)),

#print('Error: {}'.format(error(parameters_0)))

# GA parameters
IND_SIZE = 3
NPOP = 30
CXPB, MUTPB, NGEN = 0.5, 0.2, 30

creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
creator.create("Individual", list, fitness=creator.FitnessMin)

toolbox = base.Toolbox()
# Attribute generator
#toolbox.register("attr_float", random.randrange, -100, 100)
toolbox.register("attr_float", random.random)
# Structure initializers
toolbox.register("individual", tools.initRepeat,
                 creator.Individual, toolbox.attr_float, n=IND_SIZE)
toolbox.register("population", tools.initRepeat, list,
                 toolbox.individual)

toolbox.register('mate', tools.cxTwoPoint)
toolbox.register('mutate', tools.mutGaussian, mu=0, sigma=1, indpb=0.2)
toolbox.register('select', tools.selTournament, tournsize=3)
#toolbox.register('evaluate', rosenbrock)
# toolbox.register('evaluate', bohachevsky)
toolbox.register('evaluate', error)

#id = toolbox.individual()
# print(id)

if __name__ == '__main__':
    pop = toolbox.population(n=NPOP)
    hof = tools.HallOfFame(1)
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", np.mean)
    stats.register("std", np.std)
    stats.register("min", np.min)
    stats.register("max", np.max)

    try:
        pop, log = algorithms.eaSimple(pop, toolbox, cxpb=CXPB,
                                       mutpb=MUTPB, ngen=NGEN,
                                       stats=stats, halloffame=hof,
                                       verbose=True)
    except KeyboardInterrupt:
        print('Stop evolution!')
