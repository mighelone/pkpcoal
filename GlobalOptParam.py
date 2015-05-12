### PART FOR GENERIC ALGORITH OPTIMIZATION
# number of generations
NrOfGeneration = 50 #100
# size of population
NrOfPopulation = 50
#minimum for constant rate equation:[k,t_start]
EvACRMin=[5,1e-5] #first argument is k, second is t_start, third final yield
#maximum for constant rate equation:[k,t_start]
EvACRMax=[1000,0.5]
EvACRInit=[100,5e-3]
#minimum for Arrhenius equation:[A,n,E]
EvAArrhMin=[1e3,-2,1e2]
#maximum for Arrhenius equation:[A,n,E]
EvAArrhMax=[1e15,2,1e6]
EvAArrhInit=[1e10,0,15e3]
#
#minimum For Kobayashi: A1 E1 A2 E2 y1 y2
EvAKobMin=[1e3,1000,5e3,5000,0.1,0.3]
EvAKobMax=[1e8,50000,1e8,500000,0.7,1.]
EvAKobInit=[100157.0037,5849.84292004,85011813.05,12648.7543183,0.45,0.65]
#
#For Kobayashi:
EvADAEMMin=[0.5e13,150e6,15e6] # one less pair; the range is set from minimum yields to maximum yields in Pyrolysis.py
EvADAEMMax=[5.e13,300e6,50e6] # one less pair; the range is set from minimum yields to maximum yields in Pyrolysis.py
EvADAEMInit=[1.67e13,212e6,29.4e6]
#
#Nr Of activation Enrgies to solve the dE Integral for
NrOFActivtionEnergies=50


### PART OF THE GRADIENT BASED OPTIMIZATION ####

# initial guess vector of Gradient based optimization is the ...Init vector defined for the Evolutionary algorithm

selectedGradBasedOpt = 'fmin' # see  http://docs.scipy.org/doc/scipy/reference/optimize.html#module-scipy.optimize
#selectedGradBasedOpt = 'brute'
#selectedGradBasedOpt = 'leastsq'
# possible inputs: 'fmin' or 'leastsq' (recommended) ; 'fmin_cg','fmin_bfgs','fmin_ncg','fmin_slsqp'
Tolerance = 1e0   # relative error to abort
MaxIter   = 1e3     # integer

### USED FOR BOTH ####

ScaleFactor = 1e4   # large values, gives better fittings

# set this parameter True for having a gradient-based optimization after the GA optimization
optimizGrad = True

#if optimizGrad == True:
#    NrOfGeneration = 20
#    NrOfPopulation = 100
#else:
#    NrOfGeneration = 100
#    NrOfPopulation = 100
#
