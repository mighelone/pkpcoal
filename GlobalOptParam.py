### PART OF MANY POINT OPTIMIZATION
#The number of runs and the range to get the global optimum can be defined in this file
#the general shape is:
# IndexToOptimize=[2,3]       to optimize Parameter number 2 and 3 in the list
# Boundaries=[[0,1],[10,11]]  for set the range for parameter 2 from 0 to 1 and for Parameter 3 from 10 to 11
# NrOfRuns=[0,0,3,4]          to vary the parameter 2 four (3+1) times and parameter 3 five times
#
#
#Arrhenius: only the final yield is optimized [A,b,E,yield]
ArrhIndexToOptimize = [3]
# the range is set from minimum yields to maximum yields in Pyrolysis.py
ArrhNrOfRuns = [0,0,0,3]
#
#
#Arrhenius with no T**b: only the final yield is optimized [A,E,yield]
ArrhNobIndexToOptimize = [2]
# the range is set from minimum yields to maximum yields in Pyrolysis.py
ArrhNoBNrOfRuns = [0,0,3]
#
#
#Kobayashi: the parameter are [A_1, E_1, A_2, E_2]
KobIndexToOptimize = [1,3,4,5]
KobBoundaries = [[6e3,12e3],[16e3,22e3],[0.3,0.6],[0.6,1.]]
KobNrOfRuns = [0,2,0,2,2,2]
#
#
#Kobayashi with optimization parameters [A_1,A_2,E_1,alpha1]
#gets PAVM_daf from Pyrolysis.py:
KobAlphaIndexToOptimize = [2,3]
def KobAlphaBoundaries(VolatileMatterinProxAnal):
    return [[6e3,12e3],[VolatileMatterinProxAnal/100.-0.15,VolatileMatterinProxAnal/100.+0.15]]
KobAlphaNrOfRuns = [0,0,3,3]
#
#
#DAEM: the parameter are [A, E_0, sigma]
DAEMIndexToOptimize = [1,2,3]
DAEMBoundaries = [[14e3,25e3],[1e3,7e3]] # one less pair; the range is set from minimum yields to maximum yields in Pyrolysis.py
DAEMNrOfRuns = [0,3,3,3]
#
#
#
### PART FOR GENERIC ALGORITH OPTIMIZATION
# number of generations
NrOfGeneration = 10
# size of population
NrOfPopulation = 100
#minimum for Arrhenius equation:[A,E]
EvAArrhMin=[1e3,-2,1e2]
#maximum for Arrhenius equation:[A,E]
EvAArrhMax=[1e15,2,1e6]
EvAArrhInit=[1e10,0,2e4]
#
#minimum For Kobayashi:
EvAKobMin=[100.,1000,1000,5000,0.,0.5]
EvAKobMax=[50000.,20000,100000000,50000,0.5,1.]
EvAKobInit=[100157.0037,5849.84292004,85011813.05,12648.7543183,0.45,0.65]
#
#For Kobayashi:
EvADAEMMin=[1e8, 12e3, 1e3] # one less pair; the range is set from minimum yields to maximum yields in Pyrolysis.py
EvADAEMMax=[1e14,25e3, 7e3] # one less pair; the range is set from minimum yields to maximum yields in Pyrolysis.py
EvADAEMInit=[2e10,20e3,5e3]
#
#Nr Of activation Enrgies to solve the dE Integral for
NrOFActivtionEnergies=50
#
#
#the weights for yields and rates for optimization
EvAWeightY=1.
EvAWeightR=1.
