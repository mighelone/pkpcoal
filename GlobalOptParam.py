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
