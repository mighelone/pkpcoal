from pyevolve import G1DList, GSimpleGA, Selectors#, Statistics
from pyevolve import Initializators, Mutators, Consts, DBAdapters
import numpy as np

class GenericOpt(object):
    """Class which uses the pyevolve module to search for the global optimum.To initialize use the Kinetic model (e.g. Kobayashi) and the Fit one run object list, which supports the fitting process with the informations."""
    def __init__(self,KineticModel,Fit_one_runObj,Species):
        self.kinModel=KineticModel
        self.FitInfo=Fit_one_runObj
        self.Species=Species
        self.__NrGenerations=100
        self.__NrPopulation=30
        
    def setNrGenerations(self,NrOfGenerations):
        """Defines the number of generations for the generic optimization."""
        self.__NrGenerations=NrOfGenerations

    def setNrPopulation(self,NrOfPopulation):
        """Defines the size of the population for the generic optimization."""
        self.__NrPopulation=NrOfPopulation
        
    def setParamRanges(self,InitialGuess,minimum,maximum):
        """Sets the range where to evolve."""
        self.__ParamInit = np.array(InitialGuess) #[0.45,100157.0037,5849.84292004,0.55,85011813.05,12648.7543183] #kobap
        self.__ParamMin = np.array(minimum)
        self.__ParamMax = np.array(maximum)
        #initialize the two Parameter vecotrs (non-scaled and scaled)
        self.__ParameterNonSc = self.__ParamInit
        self.__ParameterSc = np.zeros(np.shape(self.__ParamInit))
        #scales between every parameter between 0 and 1
        self.__ParameterSc = (self.__ParamInit-self.__ParamMin)/(self.__ParamMax-self.__ParamMin) #x
        #the weights for optimization
        self.__a0=1.
        self.__a1=1.
        
    def setWeights(self,WeightY,WeightR):
        "Sets the yield and the weight rate for the optimization equation."""
        self.__a0=WeightY
        self.__a1=WeightR
        
    def __UpdateParam(self):
        """Updates the non-scaled Parameter. Updates the non-scaled parameter vector with the values of the sclaed vector"""
        self.__ParameterNonSc = self.__ParameterSc*(self.__ParamMax-self.__ParamMin) + self.__ParamMin
        
    def setScaledParameter(self,Parameter):
        """Sets Sclaed Parameter equal to the input parameter."""
        self.__ParameterSc=Parameter
        
    def mkResults(self):
        """Generates the result."""
        #evaluation Function:
        def kobaopt(ParameterVec):
            # rescale parameters #
            self.setScaledParameter( np.array(ParameterVec[0:len(ParameterVec)]) )
            self.__UpdateParam()
            self.kinModel.setParamVector(self.__ParameterNonSc[0:len(self.__ParameterSc)])
            error = 0
            for runnedCaseNr in range(len(self.FitInfo)):
                    # run models using CPD time history
                    t=self.FitInfo[runnedCaseNr].Time()
                    T=self.FitInfo[runnedCaseNr].Interpolate('Temp')
#                    w0=self.__a0/(( max((self.FitInfo[runnedCaseNr].Yield(self.Species))) -min((self.FitInfo[runnedCaseNr].Yield(self.Species))) )**2)
#                    w1=self.__a1/(max( ((self.FitInfo[runnedCaseNr].Rate(self.Species)))**2 ))
#                    yieldcalc = self.kinModel.calcMass(self.FitInfo[runnedCaseNr],t,T,self.Species)
#                    ratecalc = self.kinModel.deriveC(self.FitInfo[runnedCaseNr],yieldcalc)
#                    error += w0*np.sum((yieldcalc-self.FitInfo[runnedCaseNr].Yield(self.Species))**2) + w1*np.sum((ratecalc-self.FitInfo[runnedCaseNr].Rate(self.Species))**2)
                    yieldcalc = self.kinModel.calcMass(self.FitInfo[runnedCaseNr],t,T,self.Species)
                    error += np.sum((yieldcalc-self.FitInfo[runnedCaseNr].Yield(self.Species))**2)
            return 100.*error
        # scales parameter using the initialParameters
        self.__ParameterSc = (self.__ParamInit-self.__ParamMin)/(self.__ParamMax-self.__ParamMin)
        #
        genome = G1DList.G1DList(len(self.__ParameterSc))
        genome.setParams(rangemin=0, rangemax=1, bestRawScore=0.00000, roundDecimal=2)
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

	ga.setMutationRate(0.20)

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
        print "Best individual score: %.8f" % best.getRawScore()
        return self.__ParameterNonSc
