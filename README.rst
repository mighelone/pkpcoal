Pyrolysis Kinetic Preprocessor
==============================

**PKP** (Pyrolysis Kinetic Preprocessor) is software for preprocessing
kinetics of coal and biomass pyrolysis.

--------------

Introduction
------------

**PKP** was developed at TU Bergakademie Freiberg and Universitaet
Stuttgart. It interfaces different *detailed pyrolysis software*:

-  CPD
-  FG-DVC (not active)
-  Flashchain (not active)
-  Polimi

These softwares can be used for calibrating empirical pyrolysis models,
which can be easily included and used for CFD simulations. The following
empirical models are available:

-  SFOM (Single First Order reaction Model)
-  C2SM (Competing 2 Step Model)
-  DAEM (Distributed Activation Energy Model)

More information about PKP can be found in  [Vascellari2013]_

Installation
------------

For a local installation:

.. code:: bash

    # python setup.py install --user

or

::

    # python setup.py install

for a global installation.

Use
---

.. code:: bash

    runPKP

    usage: runPKP [-h] [-n NP] [-o RESULTS_DIR] [-d] yml_file

    PKP Runner

    positional arguments:
      yml_file        YAML input file

    optional arguments:
      -h, --help      show this help message and exit
      -n NP           Number of processor
      -o RESULTS_DIR  Results directory
      -d              Print debug messages


.. _input-file-label:

Input file
~~~~~~~~~~


Input files are written using **YAML** syntax:

.. code:: yaml

    # PKP input file
    # use empty for None
    #Proximate Analysis (in percent, as received):
    Coal:
      proximate_analysis:
        FC: 45.1
        VM: 50.6
        Moist: 19.0
        Ash: 4.3
      ultimate_analysis:
        C: 69
        H: 5
        N: 0.8
        O: 24.7
        S: 0.5
      #Higher Heating Value, as received, in MJ/kg:
      HHV :   0.0 
      #Tar Molecule weight, MTar:
      rho_dry : 1310 # kg/m3

    CPD:
      active : true
      fit: 
        fit0:
          model: SFOR
          species: volatiles
          parameters_min: [1e4, 50e6, 0.4]
          parameters_max: [1e9, 200e6, 0.8]
          parameters_init: [1e5, 100e6, 0.5] # not required by evolve
          method: evolve
          # from here parameters of evolve
          npop: 40
          ngen: 100
          mu: 40
          lambda_: 40
          cxpb: 0.5
          mutpb: 0.5
            
      dt: 1e-5
      increment: 5
      dt_max: 1e-5
      nmr_parameters: 
      solver:

    Polimi:
      active: true
      backend: dopri5
      mechanism:
      fit:
        fit0:
         model: SFOR
          species: volatiles
          parameters_min: [1e5, 50e6, 0.6]
          parameters_max: [1e8, 200e6, 0.7]
          parameters_init: [1e5, 100e6, 0.65] # not required by evolve
          method: evolve
          # from here parameters of evolve
          npop: 60
          ngen: 40
          mu: 60
          lambda_: 40
          cxpb: 0.6
          mutpb: 0.2
        
    BioPolimi:
      active: false
      fit:
      backend: dopri5
      mechanism:

    operating_conditions:
        pressure    : 1.0 #atmosphere
        runs: 3
        run0        : [ [ 0, 500], [ 0.005, 1500], [ 0.02, 1500] ]
        run1        : [ [ 0, 500], [ 0.003, 1300], [ 0.02, 1300] ]
        run2        : [ [ 0, 500], [ 0.01, 1300], [ 0.02, 1300] ]
        run3        : [ [ 0, 500], [ 0.1, 1400], [ 0.5, 2000] ]
        run4        : [ [ 0, 500], [ 0.1, 1400], [ 0.5, 2000] ]

