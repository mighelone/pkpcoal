'''
Evolution module
================
Manage genetic evolution using DEAP
'''
from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals

import pkp.detailed_model
import pkp.empirical_model
import numpy as np
import random
from autologging import logged

from deap import base
from deap import creator
from deap import tools
from deap.benchmarks import binary
# from deap import algorithms
from pkp import algorithms


# import multiprocessing
from pathos.multiprocessing import ProcessPool
# from scoop import futures


def error(cls_, individual):
    '''
    Calculate the error for the given individual
    '''
    err = 0
    parameters = cls_.unscale_parameters(individual)
    for run, results in cls_.ref_results.iteritems():
        m = cls_.empirical_model(parameters)
        m.operating_conditions = results['operating_conditions']
        _, y = m.run(results['t'])
        if y.ndim == 2:
            # for multivariables case take only the first solution
            y = y[:, 0]
        err += cls_.error_run(y, results['y'])
        # del m
    return err,


# @binary.bin2float(0, 1, 16)

def error_binary(cls_, individual):
    @binary.bin2float(0, 1, 16)
    def f(individual, cls_):
        return error(cls_, individual)
    return f(individual, cls_)


#@binary.bin2float(0, 1, n_decoding)
# def error_binary(cls_, individual):


@logged
class Evolution(object):
    '''
    Evolution manager based on DEAP
    '''

    def __init__(self, npop=40, ngen=30, cxpb=0.6, mutpb=0.2, mu=None,
                 lambda_=None):
        '''
        Parameters
        ----------
        npop: int
            Size of population
        ngen: int
            Number of evolve
        cxpb: float
            Crossover probability (<1).
        mutpb: float
            Mutation probability (<1).
        mu: float
            The number of individuals to select for the next generation.
        lambda_: float
            The number of children to produce at each generation.
        '''

        # GA parameters
        self.__log.debug('Init Evolution')
        self._npop = npop
        self._ngen = ngen
        self._cxpb = cxpb
        self._mutpb = mutpb
        if mu is None:
            mu = npop
        if lambda_ is None:
            lambda_ = int(2 / 3 * npop)
        self._mu = mu
        self._lambda = lambda_
        self.__log.debug('Set npop=%s', self._npop)
        self.__log.debug('Set ngen=%s', self._ngen)
        self.__log.debug('Set cxpb=%s', self._cxpb)
        self.__log.debug('Set mutpb=%s', self._mutpb)

        self._ntargets = 0
        self.ref_results = {}

        self._empirical_model = pkp.empirical_model.SFOR
        self._parameters_min = None
        self._parameters_max = None

    def set_target(self, t, y, operating_conditions, every=1):
        '''
        Set the target conditions.
        This operation has to be done as many times as necessary

        Parameter
        ---------
        t: array
            Time vector
        y: array
            Yield vector
        '''
        if not len(t) == len(y):
            raise ValueError('Length of t and y should be the same')
        self.ref_results['run{}'.format(self.n_targets)] = {
            't': np.array(t)[::every],
            'y': np.array(y)[::every],
            'operating_conditions': operating_conditions
        }
        self._ntargets += 1
        self.__log.debug('Set target run(%s)', self._ntargets)

    @property
    def n_targets(self):
        return self._ntargets

    @property
    def empirical_model(self):
        return self._empirical_model

    @empirical_model.setter
    def empirical_model(self, model):
        '''
        Set the empirical model for the calibration
        '''
        # check attributes using the EmpiricalModel attributes
        if not issubclass(model, pkp.empirical_model.EmpiricalModel):
            raise TypeError('model has to be child of EmpiricalModel!')
        self.__log.debug('Set empirical_model to %s', model)
        self._empirical_model = model

    # def __call__(self, individual):
    #    '''
    #    Calculate the error for the given individual
    #    '''
    #    err = 0
    #    parameters = self.unscale_parameters(individual)
    #    for run, results in self.ref_results.iteritems():
    #        m = self.empirical_model(parameters)
    #        m.operating_conditions = results['operating_conditions']
    #        _, y = m.run(results['t'])
    #        err += self.error_run(y, results['y'])
    #        # del m
    #    return err,

    # def __call__(self, individual):
    #    return sum(individual),

    @staticmethod
    def error_run(y, y_t):
        return np.mean((y - y_t)**2)

    def evolve(self, n_p=1, verbose=True):
        '''
        Evolve the population using the MuPlusLambda evolutionary
        algorith.
        '''
        toolbox = self.toolbox

        # Process Pool of 4 workers
        if n_p > 1:
            # pool = multiprocessing.Pool(processes=n_p)
            # pool = Pool(2)
            pool = ProcessPool(nodes=n_p)
            toolbox.register("map", pool.map)
            # toolbox.register("map", futures.map)
            self.__log.warning('n_p not supported at the moment!\n'
                               'Only serial run supported!')

        pop = toolbox.population(n=self._npop)
        hof = tools.HallOfFame(1)
        stats = tools.Statistics(lambda ind: ind.fitness.values)
        stats.register("avg", np.mean)
        stats.register("std", np.std)
        stats.register("min", np.min)
        stats.register("max", np.max)

        # pop, log = algorithms.eaSimple(pop, toolbox, cxpb=CXPB,
        #                               mutpb=MUTPB, ngen=NGEN,
        #                               stats=stats, halloffame=hof,
        #                               verbose=True)
        pop, log = algorithms.eaMuPlusLambda(pop, toolbox,
                                             mu=self._mu,
                                             lambda_=self._lambda,
                                             cxpb=self._cxpb,
                                             mutpb=self._mutpb,
                                             ngen=self._ngen,
                                             stats=stats,
                                             halloffame=hof,
                                             verbose=verbose)

        if n_p > 1:
            pool.close()

        self.pop = pop
        self.log = log
        self.stats = stats

        fitnesses = np.array([p.fitness.values for p in pop])
        best = pop[fitnesses.argmin()]

        best_parameters = self.unscale_parameters(best)

        # print('Best population', best, best_parameters)
        return best_parameters

    def register(self):
        '''
        Register settings for the Evolution algorithm using DEAP
        Check if this can be done inside a function
        '''
        creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
        creator.create("Individual", list, fitness=creator.FitnessMin)

        toolbox = base.Toolbox()
        # Attribute generator
        # toolbox.register("attr_float", random.randrange, -100, 100)

        # Structure initializers
        toolbox = self._individual(toolbox=toolbox)
        toolbox.register("population", tools.initRepeat, list,
                         toolbox.individual)

        toolbox.register('mate', tools.cxTwoPoint)
        toolbox.register('mutate', tools.mutGaussian, mu=0, sigma=1,
                         indpb=0.2)
        toolbox.register('select', tools.selTournament, tournsize=3)
        # toolbox.register('evaluate', self.error)

        self.toolbox = toolbox

    def _individual(self, toolbox):
        '''Set individual enconding'''
        toolbox.register("attr_float", random.random)
        toolbox.register("individual", tools.initRepeat,
                         creator.Individual, toolbox.attr_float,
                         n=len(self.empirical_model.parameters_names))
        toolbox.register('evaluate', error, self)
        return toolbox

    def parameters_range(self, parameters_min, parameters_max):
        self.__log.debug('par min %s len %s', parameters_min,
                         len(parameters_min))
        self.__log.debug('par min %s len %s', parameters_max,
                         len(parameters_max))
        len_model = len(self.empirical_model.parameters_names)
        self.__log.debug('len_model %s is %s',
                         self.empirical_model,
                         len_model)
        if (len(parameters_min) != len_model or
                len(parameters_max) != len_model):
            raise ValueError(
                'Define parameters min and'
                ' max with length {}'.format(len_model))
        self._parameters_min = parameters_min
        self._parameters_max = parameters_max

    def unscale_parameters(self, norm_parameters):
        '''
        Unscale parameters for the given optimization.

        Note first define min and max parameters using
        `parameters_range`
        '''
        if (self._parameters_min is None or
                self._parameters_max is None):
            raise AssertionError(
                'Define at first the range of parameters')
        return self.empirical_model.unscale_parameters(
            norm_parameters, self._parameters_min,
            self._parameters_max)


class EvolutionBinary(Evolution):
    '''
    Evolution class using binary representation
    '''
    n_decoding = 16

    def _individual(self, toolbox):
        '''Set individual enconding using binary'''
        toolbox.register("attr_int", random.randint, 0, 1)
        toolbox.register("individual", tools.initRepeat,
                         creator.Individual, toolbox.attr_int,
                         n=self.n_decoding * len(
                             self.empirical_model.parameters_names))
        # toolbox.register('evaluate', error_binary, self)
        toolbox.register('evaluate', error_binary, self)
        return toolbox
