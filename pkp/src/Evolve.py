from pyevolve import G1DList, GSimpleGA, Selectors#, Statistics
from pyevolve import Initializators, Mutators, Consts, DBAdapters
import numpy as np
import GlobalOptParam           #contains the Information of the Number Of Runs for the Global Optimum search

class GenericOpt(object):
    """ Class which uses the pyevolve module to search for the global optimum.
        To initialize use the Kinetic model (e.g. Kobayashi) and the Fit one
        run object list, which supports the fitting process with the informations.
    """

    def __init__(self, KineticModel, Fit_one_runObj, Species):
        self.kinModel=KineticModel
        self.FitInfo=Fit_one_runObj
        self.Species=Species
        self.__NrGenerations=100
        self.__NrPopulation=30
        #the weights for optimization
        self.__a0=1.
        self.__a1=1.

    def __UpdateParam(self):
        """ Updates the non-scaled parameter vector with the values of
            the scaled vector"""
        self.__ParameterNonSc = (self.__ParameterSc*(self.__ParamMax - self.__ParamMin)
                                    + self.__ParamMin)
        if self.kinModel._modelName == 'ArrheniusNoB':
            self.__ParameterNonSc[0:1] = (self.__ParameterSc[0:1]*(
                       np.log10(self.__ParamMax[0:1]) - np.log10(self.__ParamMin[0:1]))
                     + np.log10(self.__ParamMin[0:1]))
            self.__ParameterNonSc[0:1] = 10.**self.__ParameterNonSc[0:1]
        elif self.kinModel._modelName == 'Kobayashi':
            self.__ParameterNonSc[0:4] = (self.__ParameterSc[0:1]
                            * (np.log10(self.__ParamMax[0:4])
                            - np.log10(self.__ParamMin[0:4]))
                            + np.log10(self.__ParamMin[0:4]))
            self.__ParameterNonSc[0:4] = 10.**self.__ParameterNonSc[0:4]

    def estimate(self):
        """Generates the result."""
        #evaluation Function:
        def kobaopt(ParameterVec):
            # rescale parameters #
            self.setScaledParameter( np.array(ParameterVec[0:len(ParameterVec)]) )
            self.__UpdateParam()
            self.kinModel.setParamVector(
                    self.__ParameterNonSc[0:len(self.__ParameterSc)]
                    )
            error = 0
            for runnedCaseNr in range(len(self.FitInfo)):
                # run empirical models using time history from phenomenlog.
                t=self.FitInfo[runnedCaseNr].Time()
                T=self.FitInfo[runnedCaseNr].Interpolate('Temp')
                yieldcalc = self.kinModel.calcMass(
                            self.FitInfo[runnedCaseNr],t,T,self.Species
                        )
                deltaYield2 = np.power(
                      max((self.FitInfo[runnedCaseNr].Yield(self.Species)))
                    - min((self.FitInfo[runnedCaseNr].Yield(self.Species))),2.0)
                nTime = len(t)
                errori = (yieldcalc-self.FitInfo[runnedCaseNr].Yield(self.Species))**2.
                #error += self.__a0 * np.sum(errori)
                error += self.__a0 * np.sum(errori) / nTime / deltaYield2
                if self.__a1 != 0:
                    ratecalc = self.kinModel.deriveC(
                        self.FitInfo[runnedCaseNr], yieldcalc)
                    deltaRate2 = np.power(
                          max(self.FitInfo[runnedCaseNr].Rate(self.Species))
                        - min((self.FitInfo[runnedCaseNr].Rate(self.Species))), 2.0)
                    errori = np.power(
                        ratecalc-self.FitInfo[runnedCaseNr].Rate(self.Species), 2.0)
                    error += self.__a1 * np.sum(errori) / nTime / deltaRate2
            error /= len(self.FitInfo)
            #error *= GlobalOptParam.ScaleFactor
            return error

        # scales parameter using the initialParameters
        self.__ParameterSc = ( (self.__ParamInit - self.__ParamMin)
                               /(self.__ParamMax - self.__ParamMin))
        #
        genome = G1DList.G1DList(len(self.__ParameterSc))
        genome.setParams(rangemin=0, rangemax=1)
        genome.initializator.set(Initializators.G1DListInitializatorReal)
        genome.mutator.set(Mutators.G1DListMutatorRealRange)
        # The evaluator function (objective function)
        genome.evaluator.set(kobaopt)
        # Genetic Algorithm Instance
        ga = GSimpleGA.GSimpleGA(genome)
        ga.setMinimax(Consts.minimaxType["minimize"])
        # set the population size
        ga.setPopulationSize(self.__NrPopulation)
        # set the number of generation
        ga.setGenerations(self.__NrGenerations)
        # Set the Roulette Wheel selector method, the number of generations and the termination criteria
        ga.selector.set(Selectors.GRouletteWheel)
        ga.terminationCriteria.set(GSimpleGA.ConvergenceCriteria)
        ga.setMutationRate(0.4)
        ga.setCrossoverRate(1.0)
        #parallel processing
        #ga.setMultiProcessing(True)
        # Sets the DB Adapter, the resetDB flag will make the Adapter recreate
        # the database and erase all data every run, you should use this flag
        # just in the first time, after the pyevolve.db was created, you can
        # omit it.
        sqlite_adapter = DBAdapters.DBSQLite(identify="koba", resetDB=True)
        ga.setDBAdapter(sqlite_adapter)
        # Do the evolution, with stats dump, frequency of 20 generations
        ga.evolve(freq_stats=1)
        # Gets the best individual
        best=ga.bestIndividual()
        #selects the bestiniviual
        self.__ParameterSc=best[0:len(self.__ParameterSc)]
        #update the non-scaled parameterrs
        self.__UpdateParam()
        return self.__ParameterNonSc
