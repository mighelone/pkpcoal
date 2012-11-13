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
KobIndexToOptimize = [0,1,2,3]
KobBoundaries = [[100,1e5],[100,8e3],[10000,1e8],[8e3,80e3]]
KobNrOfRuns = [3,3,3,3]
#
#
#Kobayashi with optimization parameters [A_1,A_2,E_1,alpha1]
#gets PAVM_daf from Pyrolysis.py:
KobAlphaIndexToOptimize = [0,1,2,3]
def KobAlphaBoundaries(VolatileMatterinProxAnal):
    return [[1000,1e8],[1000,1e8],[500,10000],[VolatileMatterinProxAnal/100.-0.25,VolatileMatterinProxAnal/100.+0.15]]
KobAlphaNrOfRuns = [1,1,1,1]
#
#
