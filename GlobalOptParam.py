### PART FOR GENERIC ALGORITH OPTIMIZATION
# number of generations
NrOfGeneration = 100
# size of population
NrOfPopulation = 100
#minimum for constant rate equation:[k,t_start]
EvACRMin=[5,1e-5] #first argument is k, second is t_start, third final yield
#maximum for constant rate equation:[k,t_start]
EvACRMax=[1000,0.5]
EvACRInit=[100,5e-3] 
#minimum for Arrhenius equation:[A,n,E]
EvAArrhMin=[1e3,-2,1e2]
#maximum for Arrhenius equation:[A,n,E]
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

### PART OF THE GRADIENT BASED OPTIMIZATION ####

selectedGradBasedOpt = 'fmin' # see  http://docs.scipy.org/doc/scipy/reference/optimize.html#module-scipy.optimize
# possible inputs: 'fmin' or 'leastsq' (recommended) ; 'fmin_cg','fmin_bfgs','fmin_ncg','fmin_slsqp'
Tolerance = 1e-20 # relative error to abort
MaxIter   = 1e4 # integer
ScaleFactor = 1e8
