from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals

from deap import algorithms
from deap import base
from deap import creator
from deap import tools

import array
import random
import numpy as np

creator.create("FitnessMax", base.Fitness, weights=(1.0,))
creator.create("Individual", array.array, typecode='b',
               fitness=creator.FitnessMax)

toolbox = base.Toolbox()

# Attribute generator
toolbox.register("attr_bool", random.randint, 0, 1)

# Structure initializers
toolbox.register("individual", tools.initRepeat,
                 creator.Individual, toolbox.attr_bool, 100)
toolbox.register("population", tools.initRepeat, list,
                 toolbox.individual)


def evalOneMax(individual):
    return sum(individual),


toolbox.register("evaluate", evalOneMax)
toolbox.register("mate", tools.cxTwoPoint)
toolbox.register("mutate", tools.mutFlipBit, indpb=0.05)
toolbox.register("select", tools.selTournament, tournsize=3)


def main():
    random.seed(64)

    pop = toolbox.population(n=300)
    hof = tools.HallOfFame(1)
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", np.mean)
    stats.register("std", np.std)
    stats.register("min", np.min)
    stats.register("max", np.max)

    pop, log = algorithms.eaSimple(pop, toolbox, cxpb=0.5, mutpb=0.2,
                                   ngen=40, stats=stats, halloffame=hof,
                                   verbose=True)

    return pop, log, hof

if __name__ == "__main__":
    main()
