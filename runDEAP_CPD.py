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


def calc_CPD():
    '''Calc CPD'''
    runner = pkp.PKPRunner('input.yml')
    runner.operating_conditions['runs'] = 1
    runner.Polimi['active'] = False
    res = runner.run(results_dir='./optdir')
    return res['CPD']['run0']


def read_CPD():
    '''Read CPD results'''
    import pandas as pd
    return pd.read_csv('optdir/CPD_CPD-Run0.csv', index_col=0)


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
    return np.mean((y - y_cpd)**2),


calc = False

res0 = calc_CPD() if calc else read_CPD()


y_cpd = np.array(res0['fsolid'])
t_cpd = np.array(res0.index) * 1e-3

runner = pkp.PKPRunner('input.yml')
operating_conditions = runner.operating_conditions['run0']

parameters_min = np.array([3., 20, 0.4])  # logA, E MJ/kg, y0
parameters_max = np.array([14, 300, 0.8])
delta_parameters = parameters_max - parameters_min

parameters_0 = [4, 150, 0.5]


# GA parameters
IND_SIZE = 3
NPOP = 40
CXPB, MUTPB, NGEN = 0.7, 0.2, 30

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
        # pop, log = algorithms.eaSimple(pop, toolbox, cxpb=CXPB,
        #                               mutpb=MUTPB, ngen=NGEN,
        #                               stats=stats, halloffame=hof,
        #                               verbose=True)
        pop, log = algorithms.eaMuPlusLambda(pop, toolbox, mu=NPOP,
                                             lambda_=30,
                                             cxpb=CXPB,
                                             mutpb=MUTPB,
                                             ngen=NGEN,
                                             stats=stats,
                                             halloffame=hof,
                                             verbose=True)

    except KeyboardInterrupt:
        print('Stop evolution!')

    # plot the results
    fitnesses = np.array([p.fitness.values for p in pop])
    best = pop[fitnesses.argmin()]

    best_parameters = scale_parameters(best)

    print('Best population', best, best_parameters)

    m = pkp.empirical_model.SFOR(parameters=best_parameters)
    m.operating_conditions = operating_conditions
    _, y_sfor = m.run(t=t_cpd)

    fig, ax = plt.subplots()

    ax.plot(t_cpd, y_cpd, label='CPD')
    ax.plot(t_cpd, y_sfor, label='SFOR')

    ax.legend()

    plt.show()
