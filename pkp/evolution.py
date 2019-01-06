r"""
Genetic evolution module.

This module contains the class manager for the fitting of the detailed
model results using empirical models, as described by [Vascellari2013]_.

The class :class:`pkp.evolution.Evolution` allows to set different
results from detailed models, given in terms of time, yield arrays,
which corresponds to different operating conditions.

Different empirical models can be used for the calibration.
The calibration is based on the :math:`(\mu+\lambda)` genetic algorithm.

Example
-------
Define a new instance called `ga` of :class:`pkp.evolution.Evolution`::

    >>> import pkp.evolution
    >>> import pkp.empirical_model
    >>> ga = pkp.evolution.Evolution(npop=50, ngen=100, cxpb=0.2,
                                  mutpb=0.8, mu=50, lambda_=40)

The, set the reference conditions, assuming that `y` and `t` are arrays
evaluated from a detailed model::

    >>> ga.set_target(t=t, y=y,
                      operating_conditions=[[0, 400], [0.1, 1500]])

The `operating_conditions` is a list of time, temperature points.
Note that the maximum time `t` has to be smaller or equal to the maximum
time in the `operating_conditions`.
The same operation can be repeated as many time as wanted, defining
multiple operating conditions for the fit.

Then, set the empirical model to fit (i.e. the
:class:`pkp.empirical_model.SFOR` model)::

    >>> ga.empirical_model = pkp.empirical_model.SFOR

and the range of the parameters::

    >>> ga.parameters_range(parameters_min=[1e4, 50e6, 0.4],
                            parameters_max=[1e9, 200e6, 0.8])

Finally, register the parameters of the genetic algorithm using **DEAP**
toolbox and generate new `n` generations::

    >>> ga.register()
    >>> best = ga.run(n_p=4, verbose=True)
    gen	nevals	avg     	std     	min      	max
    0  	40    	0.249499	0.213491	0.0540713	0.980968
    1  	30    	0.118864	0.0777908	0.0105721	0.444615
    2  	30    	0.070285	0.0688494	0.00962298	0.432799

The best parameters are returned in `best`::

    >>> print(best)
    {u'A1': 18643.283889684601, u'A2': 574700893.25216532,
    u'y1': 0.37877698331076126, u'y2': 0.69828180155346786,
    u'E1': 48206052.660060763, u'E2': 121681065.13675486}

Information about the population and log of the evolution are
available::

    >>> ga.pop
    >>> ga.log

"""
from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals
from builtins import dict


import numpy as np
import random
import array
from autologging import logged

from deap import base
from deap import creator
from deap import tools
from deap.benchmarks import binary

# from deap import algorithms

from . import empirical_model
from . import algorithms
from . import reactor


# import multiprocessing
from pathos.multiprocessing import ProcessPool

# from scoop import futures

from ._exceptions import PKPModelError, PKPParametersError


def check_bounds(min, max):
    """
    Check if the individual exists in the given boundaries.

    Decorator which fix the limit of the function arguments to min and
    max values.

    Parameters
    ----------
    min: minimum value of the parameter
    max: minimum value of the parameter

    """

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


@logged
def error(cls_, individual):
    r"""
    Calculate the error function.

    Calculate the error for the given individual.
    The error is given by:

    .. math::
        E = \sum_j \sum_k (y_{ref, k}^j - y_{fit, k}^j)

    where the :math:`j` is the index of the run and :math:`k`
    the index of the single points for the given run.

    Parameters
    ----------
    cls_: Evolution
        Evolution instance is passed to the function for allowing
        multiprocessing.
    individual: iterable
        Iterable list containing the chromosome of the given individual.

    Return
    ------
    err: tuple
        Tuple containing the error for a multi-objective optimization.

    """
    err = 0
    parameters = cls_.unscale_parameters(individual)
    error._log.debug("Parameters:%s", parameters)
    for run, results in cls_.ref_results.items():
        # m = cls_.empirical_model(parameters)
        # m = pkp.reactor.Reactor(cls_.empirical_model,
        #                         parameters)
        # m.operating_conditions = results['operating_conditions']
        # _, y = m.run(results['t'])
        y = run_reactor(cls_.empirical_model, parameters, results)
        if y.ndim == 2:
            # for multivariables case take only the first solution
            y = y[:, 0]
        err += cls_.error_run(y, results["y"])
        # del m
    error._log.debug("Parameters:%s - err:%s", parameters, err)
    return (err,)


def run_reactor(model, parameters, results):
    """Run reactor."""
    m = reactor.Reactor(model, parameters)
    m.operating_conditions = results["operating_conditions"]
    _, y = m.run(results["t"])
    return y


# @binary.bin2float(0, 1, 16)


def error_binary(cls_, individual):
    """Return error for binary representation."""

    @binary.bin2float(0, 1, 16)
    def f(individual, cls_):
        return error(cls_, individual)

    return f(individual, cls_)


@logged
class Evolution(object):
    r"""
    Evolution manager.

    Evolution manager based on DEAP. The :math:`(\mu+\lambda)`
    algorithm is implemented from the **DEAP** library.

    Mate of two individuals is governed by a blend crossover function,
    while a Gaussian distribution govern the mutation.
    Finally a selection tournamenet allows to select the :math:`(\mu)`
    individuals of the next generation from the joined population
    :math:`(\mu+\lambda)`.
    """

    def __init__(
        self, npop=40, ngen=30, cxpb=0.6, mutpb=0.2, mu=None, lambda_=None, skip=1
    ):
        """
        Init the evolution manager.

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
        skip: int
            Skip rows in the results

        """
        # GA parameters
        self.__log.debug("Init Evolution")
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
        self.__log.debug("Set npop=%s", self._npop)
        self.__log.debug("Set ngen=%s", self._ngen)
        self.__log.debug("Set cxpb=%s", self._cxpb)
        self.__log.debug("Set mutpb=%s", self._mutpb)

        self._ntargets = 0
        self.ref_results = {}

        self._empirical_model = empirical_model.SFOR
        self._parameters_min = None
        self._parameters_max = None

        self._skip = skip

    def set_target(self, t, y, operating_conditions):
        """
        Set the target conditions.

        This operation has to be done as many times as necessary

        Parameters
        ----------
        t: array
            Time vector [t0, t1, ..., tN]
        y: array
            Yield vector [y0, y1, ..., yN]
        operating_conditions: array, list
            Time, temperature list [[t0, T0], ..., [tM, TM]]

        Note
        ----
        The length `N` and `M` of the arrays can be different.

        """
        if not len(t) == len(y):
            raise ValueError("Length of t and y should be the same")
        self.ref_results["run{}".format(self.n_targets)] = {
            "t": np.array(t)[:: self._skip],
            "y": np.array(y)[:: self._skip],
            "operating_conditions": operating_conditions,
        }
        self._ntargets += 1
        # self.__log.debug('Set target run(%s)', self._ntargets)

    @property
    def n_targets(self):
        """Number of target solutions."""
        return self._ntargets

    @property
    def empirical_model(self):
        """Empirical model."""
        return self._empirical_model

    @empirical_model.setter
    def empirical_model(self, model):
        # check attributes using the EmpiricalModel attributes
        if not issubclass(model, empirical_model.EmpiricalModel):
            raise PKPModelError("model has to be child of EmpiricalModel!")
        self.__log.debug("Set empirical_model to %s", model)
        self._empirical_model = model

    @staticmethod
    def error_run(y, y_t):
        """Calculate the error."""
        return np.mean((y - y_t) ** 2)

    def evolve(self, n_p=1, verbose=True):
        r"""
        Evolve the population.

        The :math:`(\mu+\lambda)` evolutionary algorithm is used.

        Parameters
        ----------
        n_p: int, default=1
            Number of precess to generate the new population.
        verbose: bool, default=True
            Print extra message

        Returns
        -------
        best: dict
            Best parameters dict

        """
        toolbox = self.toolbox

        # Process Pool of 4 workers
        if n_p > 1:
            # pool = multiprocessing.Pool(processes=n_p)
            # pool = Pool(2)
            pool = ProcessPool(nodes=n_p)
            toolbox.register("map", pool.map)
            # toolbox.register("map", futures.map)
            # self.__log.warning('n_p not supported at the moment!\n'
            #                   'Only serial run supported!')

        pop = toolbox.population(n=self._npop)
        hof = tools.HallOfFame(1)

        stats = self._set_stats()

        # pop, log = algorithms.eaSimple(pop, toolbox, cxpb=CXPB,
        #                               mutpb=MUTPB, ngen=NGEN,
        #                               stats=stats, halloffame=hof,
        #                               verbose=True)
        pop, log = algorithms.eaMuPlusLambda(
            pop,
            toolbox,
            mu=self._mu,
            lambda_=self._lambda,
            cxpb=self._cxpb,
            mutpb=self._mutpb,
            ngen=self._ngen,
            stats=stats,
            halloffame=hof,
            verbose=verbose,
        )
        # if n_p > 1:
        #    pool.close()

        self.pop = pop
        self.log = log
        self.stats = stats

        fitnesses = np.array([p.fitness.values for p in pop])
        best = pop[fitnesses.argmin()]

        self.__log.debug("best scaled %s", best)
        best = self.unscale_parameters_final(best)
        self.__log.debug("best non-scaled %s", best)
        return {p: v for p, v in zip(self.empirical_model.parameters_names(), best)}

    def _set_stats(self):
        stats = tools.Statistics(lambda ind: ind.fitness.values)
        stats.register("avg", np.mean)
        stats.register("std", np.std)
        stats.register("min", np.min)
        stats.register("max", np.max)
        return stats

    def register(self):
        """
        Register settings for the Evolution algorithm using DEAP.

        Note
        ----
        Check if this can be done inside a function
        """
        if hasattr(creator, "FitnessMin"):
            del creator.FitnessMin
        creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
        if hasattr(creator, "Individual"):
            del creator.Individual
        creator.create(
            "Individual", array.array, typecode="d", fitness=creator.FitnessMin
        )

        toolbox = base.Toolbox()
        # Attribute generator
        # toolbox.register("attr_float", random.randrange, -100, 100)

        # Structure initializers
        toolbox = self._individual(toolbox=toolbox)

        # toolbox.register('evaluate', self.error)

        self.toolbox = toolbox

    def _individual(self, toolbox):
        """
        Set individual enconding.

        This function can be used for defining settings for specific evolution
        strategies.
        """
        # Random number generation
        # random.random generates float between 0 and 1
        toolbox.register("attr_float", random.random)
        # individual uses n chromosomes (as many model parameters)
        # and attr_float
        toolbox.register(
            "individual",
            tools.initRepeat,
            creator.Individual,
            toolbox.attr_float,
            n=len(self.empirical_model.parameters_names()),
        )
        # define the fit function
        toolbox.register("evaluate", error, self)
        # define the population as list of individuals
        toolbox.register("population", tools.initRepeat, list, toolbox.individual)

        # define the mate algorithm
        # cxTwoPoint is good for intege.parr/binary chromosomes
        # toolbox.register('mate', tools.cxTwoPoint)
        # blend crossover extending of 0.1 respect to the parameters
        # range. This can produce values out of the range 0-1.
        toolbox.register("mate", tools.cxBlend, alpha=0.25)
        # define the mutate algorithm
        toolbox.register("mutate", tools.mutGaussian, mu=0, sigma=1, indpb=0.2)
        # define the select algorithm
        toolbox.register("select", tools.selTournament, tournsize=3)
        toolbox.decorate("mate", check_bounds(0, 1))
        toolbox.decorate("mutate", check_bounds(0, 1))
        return toolbox

    def parameters_range(self, parameters_min, parameters_max):
        """
        Define the range of variation of the model parameters.

        Parameters
        ----------
        parameters_min: list
            List of minimum values of the parameters
        parameters_max: list
            List of maximum values of the parameters

        """
        self.__log.debug("par min %s len %s", parameters_min, len(parameters_min))
        self.__log.debug("par min %s len %s", parameters_max, len(parameters_max))
        len_model = len(self.empirical_model.parameters_names())
        self.__log.debug("len_model %s is %s", self.empirical_model, len_model)
        if len(parameters_min) != len_model or len(parameters_max) != len_model:
            raise PKPParametersError(
                "Define parameters min and" " max with length",
                len_model,
                self.empirical_model.parameters_names(),
            )
        self._parameters_min = parameters_min
        self._parameters_max = parameters_max

    def unscale_parameters(self, norm_parameters):
        """
        Unscale parameters.

        The scaled parameters defined between 0 and 1 are converted to
        non-scaled parameters required to set the empirical model.

        Note
        ----
        first define min and max parameters using
        :method:`parameters_range`.

        Parameters
        ----------
        norm_parameters: array
            Normalized parameters array.

        """
        if self._parameters_min is None or self._parameters_max is None:
            raise AssertionError("Define at first the range of parameters")
        return self.empirical_model.unscale_parameters(
            norm_parameters, self._parameters_min, self._parameters_max
        )

    def unscale_parameters_final(self, norm_parameters):
        """
        Unscale the parameters for the final step.

        Only for final step of evolution
        """
        return self.unscale_parameters(norm_parameters)


class EvolutionBinary(Evolution):
    """Evolution class using binary representation."""

    n_decoding = 16

    def _individual(self, toolbox):
        """Set individual enconding using binary."""
        toolbox.register("attr_int", random.randint, 0, 1)
        toolbox.register(
            "individual",
            tools.initRepeat,
            creator.Individual,
            toolbox.attr_int,
            n=self.n_decoding * len(self.empirical_model.parameters_names()),
        )
        # toolbox.register('evaluate', error_binary, self)
        toolbox.register("evaluate", error_binary, self)
        toolbox.register("population", tools.initRepeat, list, toolbox.individual)

        toolbox.register("mate", tools.cxTwoPoint)
        toolbox.register("mutate", tools.mutFlipBit, indpb=0.05)
        toolbox.register("select", tools.selTournament, tournsize=3)
        return toolbox

    def unscale_parameters_final(self, individual):
        """First convert to float from binary and then unscale parameters."""

        @binary.bin2float(0, 1, 16)
        def f(individual, cls):
            return cls.unscale_parameters(individual)

        return f(individual, self)
