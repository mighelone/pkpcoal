from pyevolve import G1DList, GSimpleGA, Selectors  # , Statistics
from pyevolve import Initializators, Mutators, Consts, DBAdapters
import numpy as np
# contains the Information of the Number Of Runs for the Global Optimum
# search
import GlobalOptParam


class GenericOpt(object):
    """
    Class which uses the pyevolve module to search for the global
    optimum.To initialize use the Kinetic model (e.g. Kobayashi) and
    the Fit one run object list, which supports the fitting process
    with the informations.
    """

    def __init__(self, KineticModel, Fit_one_runObj, Species):
        self.kinModel = KineticModel
        self.FitInfo = Fit_one_runObj
        self.Species = Species
        self._NrGenerations = 100
        self._NrPopulation = 30
        # the weights for optimization
        self._a0 = 1.
        self._a1 = 1.

    def setNrGenerations(self, NrOfGenerations):
        """Defines the number of generations for the generic
        optimization."""
        self._NrGenerations = NrOfGenerations

    def setNrPopulation(self, NrOfPopulation):
        """Defines the size of the population for the generic
        optimization."""
        self._NrPopulation = NrOfPopulation

    def setParamRanges(self, InitialGuess, minimum, maximum):
        """Sets the range where to evolve."""
        self.__ParamInit = np.array(
            InitialGuess)
        # [0.45,100157.0037,5849.84292004,0.55,85011813.05,12648.7543183]
        # kobap
        self._ParamMin = np.array(minimum)
        self._ParamMax = np.array(maximum)
        # initialize the two Parameter vecotrs (non-scaled and scaled)
        self._ParameterNonSc = self.__ParamInit
        self._ParameterSc = np.zeros(np.shape(self.__ParamInit))
        # scales between every parameter between 0 and 1
        self._ParameterSc = (self.__ParamInit - self._ParamMin) / \
                            (self._ParamMax - self._ParamMin)  # x

    def setWeights(self, WeightY, WeightR):
        """
        Sets the yield and the weight rate for
        the optimization equation.
        """
        self._a0 = WeightY
        self._a1 = WeightR

    def _UpdateParam(self):
        """
        Updates the non - scaled Parameter. Updates the non - scaled
        parameter vector with the values of the sclaed vector"""
        self._ParameterNonSc = self._ParameterSc * \
                               (self._ParamMax - self._ParamMin) + self._ParamMin
        if self.kinModel._modelName == 'ArrheniusNoB':
            self._ParameterNonSc[0:1] = self._ParameterSc[
                                        0:1] * (np.log10(self._ParamMax[0:1]) -
                                                np.log10(self._ParamMin[0:1])) + \
                                        np.log10(self._ParamMin[0:1])
            self._ParameterNonSc[0:1] = 10. ** self._ParameterNonSc[0:1]
        elif self.kinModel._modelName == 'Kobayashi':
            self._ParameterNonSc[0:4] = self._ParameterSc[
                                        0:1] * (np.log10(self._ParamMax[0:4]) -
                                                np.log10(self._ParamMin[0:4])) + \
                                        np.log10(self._ParamMin[0:4])
            self._ParameterNonSc[0:4] = 10. ** self._ParameterNonSc[0:4]

    def setScaledParameter(self, Parameter):
        """Sets Sclaed Parameter equal to the input parameter."""
        self._ParameterSc = Parameter

    def mkResults(self):
        """Generates the result."""
        # evaluation Function:
        def kobaopt(ParameterVec):
            # rescale parameters #
            self.setScaledParameter(
                np.array(ParameterVec[0:len(ParameterVec)]))
            self._UpdateParam()
            self.kinModel.setParamVector(
                self._ParameterNonSc[0:len(self._ParameterSc)])
            error = 0
            for runnedCaseNr in range(len(self.FitInfo)):
                # run empirical models using time history from
                # phenomenlog.
                t = self.FitInfo[runnedCaseNr].Time()
                T = self.FitInfo[runnedCaseNr].Interpolate('Temp')
                yieldcalc = self.kinModel.calcMass(
                    self.FitInfo[runnedCaseNr], t, T, self.Species)
                deltaYield2 = (
                    max((self.FitInfo[runnedCaseNr].Yield(self.Species))
                        ) -
                    min((self.FitInfo[runnedCaseNr].Yield(self.Species))
                        ))**2.
                nTime = len(t)
                errori = (
                    yieldcalc - self.FitInfo[runnedCaseNr].Yield(self.Species))**2.
                # error += self.__a0 * np.sum(errori)
                error += self._a0 * \
                         np.sum(errori) / nTime / deltaYield2
                if self._a1 != 0:
                    ratecalc = self.kinModel.deriveC(
                        self.FitInfo[runnedCaseNr], yieldcalc)
                    deltaRate2 = (
                        max(self.FitInfo[runnedCaseNr].Rate(
                            self.Species)) -
                        min((self.FitInfo[runnedCaseNr].Rate(
                            self.Species))))**2.
                    errori = (
                        ratecalc - self.FitInfo[runnedCaseNr].Rate(
                            self.Species))**2.
                    error += self._a1 * \
                             np.sum(errori) / nTime / deltaRate2

            error /= len(self.FitInfo)
            # error *= GlobalOptParam.ScaleFactor
            return error

        # scales parameter using the initialParameters
        self._ParameterSc = (self.__ParamInit - self._ParamMin) / \
                            (self._ParamMax - self._ParamMin)
        #
        genome = G1DList.G1DList(len(self._ParameterSc))
        genome.setParams(rangemin=0, rangemax=1)
        genome.initializator.set(
            Initializators.G1DListInitializatorReal)
        genome.mutator.set(Mutators.G1DListMutatorRealRange)
        # The evaluator function (objective function)
        genome.evaluator.set(kobaopt)
        # Genetic Algorithm Instance
        ga = GSimpleGA.GSimpleGA(genome)
        ga.setMinimax(Consts.minimaxType["minimize"])
        # set the population size
        ga.setPopulationSize(self._NrPopulation)
        # set the number of generation
        ga.setGenerations(self._NrGenerations)
        # Set the Roulette Wheel selector method, the number of
        # generations and the termination criteria
        ga.selector.set(Selectors.GRouletteWheel)
        ga.terminationCriteria.set(GSimpleGA.ConvergenceCriteria)
        ga.setMutationRate(0.40)
        ga.setCrossoverRate(1.0)

        # parallel processing
        # ga.setMultiProcessing(True)
        # Sets the DB Adapter, the resetDB flag will make the Adapter
        # recreate
        # the database and erase all data every run, you should use
        # this flag
        # just in the first time, after the pyevolve.db was created,
        # you can omit it.
        sqlite_adapter = DBAdapters.DBSQLite(
            identify="koba", resetDB=True)
        ga.setDBAdapter(sqlite_adapter)
        # Do the evolution, with stats dump, frequency of 20 generations
        ga.evolve(freq_stats=1)
        # Gets the best individual
        best = ga.bestIndividual()
        print best
        print 'Selected parameters'
        # selects the bestiniviual
        self._ParameterSc = best[0:len(self._ParameterSc)]
        print self._ParameterNonSc
        # update the non-scaled parameterrs
        self._UpdateParam()
        return self._ParameterNonSc
