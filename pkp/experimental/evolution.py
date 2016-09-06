from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals

from deap import base
from deap import creator
from deap import tools
from deap.benchmarks import binary
import numpy as np

import random

import pkp
import pkp.evolution


def checkBounds(min, max):
    def decorator(func):
        def wrappper(*args, **kargs):
            offspring = func(*args, **kargs)
            for child in offspring:
                for i in range(len(child)):
                    if child[i] > max:
                        child[i] = max
                    elif child[i] < min:
                        child[i] = min
            return offspring
        return wrappper
    return decorator


def error1(cls_, individual):
    '''
    Calculate the error for the given individual
    '''
    def error_run(y, y_t):
        return np.mean((y - y_t)**2, axis=0)
    err = 0
    parameters = cls_.unscale_parameters(individual)
    for run, results in cls_.ref_results.iteritems():
        m = cls_.empirical_model(parameters)
        m.operating_conditions = results['operating_conditions']
        _, y = m.run(results['t'])
        y = y[:, 0:2]
        err += error_run(y, results['y'])
        # del m
    return tuple(err.tolist())


class Evolve2Stream(pkp.evolution.Evolution):
    '''
    Optimize a empirical model with two streams y1 and y2
    '''

    def register(self):
        '''
        Register settings for the Evolution algorithm using DEAP
        Check if this can be done inside a function
        '''
        creator.create("FitnessMin", base.Fitness, weights=(-1.0, -1.0))
        creator.create("Individual", list, fitness=creator.FitnessMin)

        toolbox = base.Toolbox()
        # Attribute generator
        # toolbox.register("attr_float", random.randrange, -100, 100)

        # Structure initializers
        toolbox = self._individual(toolbox=toolbox)

        # toolbox.register('evaluate', self.error)

        self.toolbox = toolbox

    def _individual(self, toolbox):
        '''
        Set individual enconding. This function can be used for
        defining settings for specific evolution strategies
        '''
        # Random number generation
        # random.random generates float between 0 and 1
        toolbox.register("attr_float", random.random)
        # individual uses n chromosomes (as many model parameters)
        # and attr_float
        toolbox.register("individual", tools.initRepeat,
                         creator.Individual, toolbox.attr_float,
                         n=len(self.empirical_model.parameters_names))
        # define the fit function
        toolbox.register('evaluate', error1, self)
        # define the population as list of individuals
        toolbox.register("population", tools.initRepeat, list,
                         toolbox.individual)

        # define the mate algorithm
        # cxTwoPoint is good for integer/binary chromosomes
        # toolbox.register('mate', tools.cxTwoPoint)
        # blend crossover extending of 0.1 respect to the parameters
        # range. This can produce values out of the range 0-1.
        toolbox.register('mate', tools.cxBlend, alpha=0.1)
        # define the mutate algorithm
        toolbox.register('mutate', tools.mutGaussian, mu=0, sigma=1,
                         indpb=0.2)
        toolbox.decorate("mate", checkBounds(0, 1))
        toolbox.decorate("mutate", checkBounds(0, 1))
        # define the select algorithm
        # toolbox.register('select', tools.selTournament, tournsize=3)
        # pareto front selector
        toolbox.register("select", tools.selNSGA2)
        return toolbox

    def _set_stats(self):
        stats = tools.Statistics(lambda ind: ind.fitness.values)
        stats.register("avg", np.mean, axis=0)
        stats.register("std", np.std, axis=0)
        stats.register("min", np.min, axis=0)
        stats.register("max", np.max, axis=0)
        return stats
