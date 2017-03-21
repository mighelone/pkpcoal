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
-  BioPolimi

These softwares can be used for calibrating empirical pyrolysis models,
which can be easily included and used for CFD simulations. The following
empirical models are available:

-  SFOM (Single First Order reaction Model)
-  C2SM (Competing 2 Step Model)
-  DAEM (Distributed Activation Energy Model)
- Biagini-Tognotti model

More information about PKP can be found in  [Vascellari2013]_

Installation
---------

`PKP` requires an updated python distribution. It is reccomanded to
use _Anaconda: https://www.continuum.io/downloads or to create a
virtual environment usig `venv`.

Installation using Anaconda
~~~~~~~~~~~~~~~~~~~~~~

Follow the instruction to install _Anaconda:
https://www.continuum.io/downloads.
Most of the packages should be already installed.

.. code:: bash
	  # conda install matplotlib numpy scipy pandas


Remaining packages can be installed using `pip`. For a local
installation: 

.. code:: bash
	  # pip install . --user

or:

.. code:: bash
	  # pip install . 

The missing packages will be automatically installed by `pip`.

It is possible to speed-up the calculation using just-in-time (jit)
compilation using _Numba: http://numba.pydata.org/. Numba is quite
complicate to install, and it is reccomanded to use `Anaconda`:

.. code:: bash
	  # conda install numba

Installation using virtual environmnt
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Virtual environment is used to create local python installation,
specific for a project. It allows to have an update python
distribution without touching the python system installation.

In python 2.7 `virtualenv` needs to be installed:

.. code:: bash
	  # pip install virtualenv
	  # virtualenv -p PYTHON $HOME/pkp
	  # source $HOME/pkp/bin/activate

Here a new python environment, called pkp is created in $HOME and it
is activated. The option `-p` defines the original python interpreter
used to create the new environment.

In python 3, it is a default module:

.. code:: bash
	  # python -m venv -p $HOME/pkp
	  # source $HOME/pkp/bin/activate

The required packages will be installed with:

.. code:: bash
	  # pip install . --user

or:

.. code:: bash
	  # pip install .

Installation of Cantera
~~~~~~~~~~~~~~~~~~

`Polimi` and `bioPolimi` models. require the use of_Cantera:
http://www.cantera.org/docs/sphinx/html/index.html .
Follow the instruction to build `cantera` in your system.

With anaconda, cantera can be installed:

.. code:: bash
	  # conda install -c cantera cantera=2.3.0


This command according to your system is not always working properly.


Use
---

The `runPKP` script executes `PKP`. The script can be found in the src
directory, or in the bin directory in the python installation path:
`$PYTHON_PATH/bin/runPKP`.
If a local installation was used (`--user` option with `pip`), the
script is in `$HOME/.local/bin/runPKP`.

..code:: bash
	  # runPKP input.yml -o Results

The command will execute `PKP` with `input.yml` as configuration
file. The Results will be saved in the directoryb `Results`.

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


Input files are written using **YAML** syntax. In the first part, the
information about the coal are introduced:

.. code:: yaml
    # PKP input file
    # use empty for None
    #Proximate Analysis (in percent, as received):
    # The values are normalized during the execution of PKP
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
      # with HHV=0.0 it is calculated with the Dulong formula
      HHV :   0.0 
      rho_dry : 1310 # kg/m3

In the next section the parameters for each detailed model are
defined, including the specific settings for the model

.. code:: yaml
    # name of the detailed model
    CPD:
      # define if the model is active
      active : true
      dt: 1e-5
      increment: 5
      dt_max: 1e-5
      nmr_parameters: 

The next part defines the `fit` of the CPD model (note that `fit` is a
subsession of `CPD`. Different fit can be used for the same model.

.. code:: yaml
      fit: 
        fit0:
          model: SFOR
          species: volatiles
          parameters_min: [1e4, 50e6, 0.4]
          parameters_max: [1e9, 200e6, 0.8]
          parameters_init: [1e5, 100e6, 0.5] # not required by evolve
          method: evolve+min
          # from here parameters of evolve
          npop: 40
          ngen: 100
          mu: 40
          lambda_: 40
          cxpb: 0.5
          mutpb: 0.5
       fit: 
        fit1:
          model: DAEM
          species: volatiles
          parameters_min: [1e4, 50e6, 5e6, 0.4]
          parameters_max: [1e9, 200e6, 20e6, 0.8]
          parameters_init: [1e5, 100e6, 10e6, 0.5] # not required by evolve
          method: evolve+min
          # from here parameters of evolve
          npop: 40
          ngen: 100
          mu: 40
          lambda_: 40
          cxpb: 0.5
          mutpb: 0.5


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
