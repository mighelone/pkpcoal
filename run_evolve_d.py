from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals

from deap import base
from deap import creator
from deap import tools
from deap import algorithms

from deap.benchmarks import rosenbrock, bohachevsky, schaffer

import random
import numpy as np

IND_SIZE = 2
NPOP = 100

CXPB, MUTPB, NGEN = 0.2, 0.2, 1000

creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
creator.create("Individual", list, fitness=creator.FitnessMin)

toolbox = base.Toolbox()
# Attribute generator

toolbox.register("attr_float", random.randrange, -100, 100)
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
toolbox.register('evaluate', schaffer)


if __name__ == '__main__':
    pop = toolbox.population(n=300)
    hof = tools.HallOfFame(1)
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", np.mean)
    stats.register("std", np.std)
    stats.register("min", np.min)
    stats.register("max", np.max)

    pop, log = algorithms.eaSimple(pop, toolbox, cxpb=1.0, mutpb=0.3,
                                   ngen=40, stats=stats, halloffame=hof,
                                   verbose=True)
