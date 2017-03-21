.. _input-label:

Input file
==========

All the configurations to run **PKP** are defined in a comprehensive input file,
using **YAML** syntax.

The input file consists of different sections, which reports the main information for running the detailed models for the given coal/biomass and for calibrating the empirical models based on the results of the first.

Coal settings
-------------

The first section defines the settings for the coal/biomass used::

  Coal:
     name: Pittsburgh
     proximate_analysis:
         FC: 56.46
         VM: 35.89
         Moist: 0.47
         Ash: 6.95
     ultimate_analysis:
         C: 75.23
         H: 5.16
         N: 1.43
         O: 9.83
         S: 2
     # Higher Heating Value, as received, in MJ/kg:
     HHV : 25.7 # MJ/kg
     # Tar Molecule weight, MTar:
     rho_dry : 1310 # kg/m3

In the coal section, the proximate and ultimate analyses are reported. The numbers will be normalized inside PKP, therefore if the sum is different from one, the values will be opportunely scaled. The Higher Heating Value of the coal can be specified, or if an empty field is left, it will be calculated by means of the Dulong formula.

Operating conditions settings
-----------------------------

In this section the operating conditions can be defined::

  operating_conditions:
    pressure    : 1.0 #atmosphere
    runs: 2
    run0        : [ [ 0, 500], [ 0.005, 1500], [ 0.02, 1500] ]
    run1        : [ [ 0, 500], [ 0.003, 1700], [ 0.02, 1700] ]
    run2        : [ [ 0, 500], [ 0.01, 1300], [ 0.02, 1300] ]
    run3        : [ [ 0, 500], [ 0.01, 1500], [ 0.02, 1500] ]

At first the pressure in atmosphere is reported. The fields `run0` to `run3` describes the time, temperature conditions of the different runs. Each block is composed by several points:
`[[t0, T0], [t1, T1], ..., [tn, Tn]`.
All the reported runs (in this example from 0 to 3) will be executed using the detailed models chosen. Of these runs only the first 2 (as reported by `runs`) will be used for the calibration.
This allows to check results of the calibration using different conditions.


Detailed model section
----------------------

Different sections can be defined for each detailed model implemented in `PKP`. At the moment, the following detailed models are supported:

- CPD
- CPDfortran
- Polimi
- BioPolimi

The section starts with the name of the model::

  CPD:
    active : true
    dt: 1e-4
    increment: 1
    dt_max: 1e-2
    nmr_parameters:
    fit0:
      ...
    fit1:
      ...

  Polimi:
    active: false

The most important option is `active`, which allows to switch on/off the model from be executed.
In this example, `CPD` will be executed using the operating conditions above, while `Polimi` will not be executed.
The next options are specific options to passed to the model. For example: `dt` and `dt_max` are the initial and maximum time step used by the ODE solver. `nmr_parameters` allows to specifiy the NMR parameters for CPD, instead of using the ones obtained from the internal correlation.
The `increment` parameters reduces the number of output time steps, in order to speed-up calibrations. This value should be chosen accurately to avoid to expensive computations.

Empirical model calibration section
-----------------------------------

This is a subsection of the detailed model (i.e. `CPD` or `Polimi`)::

  CPD:
    fit0:
       active: true
       model: SFOR
       species: volatiles
       parameters_min: [1e4, 50e6, 0.5]
       parameters_max: [1e9, 200e6, 0.6]
       parameters_init: [1e5, 100e6, 0.5] 
       method: evolve
       # from here parameters of evolve
       npop: 40
       ngen: 40
       mu: 40
       lambda_: 40
       cxpb: 0.7
       mutpb: 0.3
    fit1:
       active: false
       model: C2SM
       species: volatiles
       parameters_min: [1e3, 20e6, 0.3, 1e6, 100e6, 0.6]
       parameters_max: [1e6, 100e6, 0.5, 1e9, 200e6, 1]
       parameters_init: [1e5, 50e6, 0.4, 1e8, 150e6, 0.7]
       method: evolve+min
       # from here parameters of evolve
       npop: 40
       ngen: 40
       mu: 40
       lambda_: 40
       cxpb: 0.7
       mutpb: 0.3

Below the detailed model section, an arbitrary number of fits can be defined.
For each section (here `fit0` and `fit1`) the `model` name is defined:
The following models are supported: `SFOR`, `C2SM`, `DAEM`, `Biagini`.
Similarly to the detailed model section, an `active` parameter can be set to switch on/off the calibration.
The `species` define which species from the detailed model output will be fitted using the empirical model approach. `volatiles` is used to calibrate the whole volatile yields.
The `method` field defines which optimization strategy to use. The options are:

- `evolve`: use evolutionary algorithm (EA) based on mu+lambda algorithm.
- `min`: use the BFGS minimization method
- `evolve+min`: first search the optimum using the evolution algorithm, than starting from the best value performs a finer search using BFGS. This method is recommended for obtaining the best results.
  
The fields `parameters_min` and `parameters_max` define the minimum and maximum fields for each parameters of the empirical model. They are used by the evolutionary algorithm.
The field `parameters_init` defines the initial parameter set for the BFGS algorithm. Note, that it works only for the option `method: min`, otherwise in `method: evolve+min` the best results of the EA will be used.

The next fields are the parameters of the EA:

- `npop`: initial population size
- `ngen`: number of generation
- `mu`: size of the population during the generations
- `lambda_` number of new individuals generated at each generation
- `cxpb`: probability that two individuals mate generating two new individuals
- `mutpb`: probability of mutation. The sum of the two probabilities has to be: `cxpb+mutpb <=1`
    





