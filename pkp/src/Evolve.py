from pyevolve import G1DList, GSimpleGA, Selectors#, Statistics
from pyevolve import Initializators, Mutators, Consts, DBAdapters
import numpy as np
from pkp.src.Models import ModelError

class GenericOpt(object):
    """ Class which uses the pyevolve to search the global optimum.

        To initialize use the Kinetic model (e.g. Kobayashi) and the Fit one
        run object list, which supports the fitting process with the informations.
        TODO implement limits 
    """
    __NrGenerations=200
    __NrPopulation=30
    optimizer = "None"

    def __init__(self, inputs):
        print 'Genetic Optimsation initialized'
        self.inputs     = inputs
        self.weightMass = inputs['weightMass']
        self.weightRate = inputs['weightRate']

    # def updateParams(self):
    #     """ Updates the non-scaled parameter vector with the values of
    #         the scaled vector """
    #     self.unscaledParams = (
    #             self.scaledParams*(self.deltaParams) + self.minParams
    #         )
    #     # if self.kinModel._modelName == 'ArrheniusNoB':
    #     #     self.unscaledParams[0:1] = (self.scaledParams[0:1]*(
    #     #                np.log10(self.maxParams[0:1]) - np.log10(self.minParams[0:1]))
    #     #              + np.log10(self.minParams[0:1]))
    #     #     self.unscaledParams[0:1] = 10.**self.unscaledParams[0:1]
    #     # elif self.kinModel._modelName == 'Kobayashi':
    #     #     self.unscaledParams[0:4] = (self.scaledParams[0:1]
    #     #                     * (np.log10(self.maxParams[0:4])
    #     #                     - np.log10(self.minParams[0:4]))
    #     #                     + np.log10(self.minParams[0:4]))
    #     #     self.unscaledParams[0:4] = 10.**self.unscaledParams[0:4]
    #
    @property
    def deltaParams(self):
        """ computes maxParams-minParams """
        return self.maxParams - self.minParams

    def estimate(self, results, model, species):
        """ Generates the result."""
        # scales parameter using the initialParameters
        # bounds params between 0 and 1
        # self.initParams = model.parameter
        # self.minParams = model.lowerLimit
        # self.maxParams = model.upperLimit
        # initialize the two Parameter vectors 
        # (non-scaled and scaled)
        # self.unscaledParams = self.initParams
        # self.scaledParams = (
        #         (self.initParams - self.minParams)/(self.deltaParams)
        # )
        genome = G1DList.G1DList(len(model.parameter))
        genome.setParams(rangemin=0,rangemax=1)
        genome.initializator.set(Initializators.G1DListInitializatorReal)
        genome.mutator.set(Mutators.G1DListMutatorRealRange)
        # The evaluator function (objective function)
        error_func = (ModelError.ls_input_func 
            if self.optimizer == 'leastsq'
            else ModelError.min_input_func
        )
        model_error = ModelError(
            results, model, species,
            error_func, self.weightMass, self.weightRate
        )
        genome.evaluator.set(model_error.input_func)
        # Genetic Algorithm Instance
        ga = GSimpleGA.GSimpleGA(genome)
        ga.setMinimax(Consts.minimaxType["minimize"])
        # set the population size
        ga.setPopulationSize(self.__NrPopulation)
        # set the number of generation
        ga.setGenerations(self.__NrGenerations)
        # Set the Roulette Wheel selector method,
        # the number of generations and the termination criteria
        ga.selector.set(Selectors.GRouletteWheel)
        ga.terminationCriteria.set(GSimpleGA.ConvergenceCriteria)
        ga.setMutationRate(0.4)
        ga.setCrossoverRate(1.0)
        # parallel processing
        # ga.setMultiProcessing(True)
        # Sets the DB Adapter, the resetDB flag will make the Adapter recreate
        # the database and erase all data every run, you should use this flag
        # just in the first time, after the pyevolve.db was created, you can
        # omit it.
        # sqlite_adapter = DBAdapters.DBSQLite(identify="koba", resetDB=True)
        # ga.setDBAdapter(sqlite_adapter)
        # Do the evolution, with stats dump, frequency of 20 generations
        ga.evolve(freq_stats=10)
        # Gets the best individual
        best = ga.bestIndividual() # update or find best model
        print best[0:]
        best = best[0:]
        # selects the best indiviual
        # self.scaledParams = best[0:len(self.scaledParams)]
        # update the non-scaled parameterrs
        # self.__UpdateParam()
        return model.updateParameter(best).recalcMass()
