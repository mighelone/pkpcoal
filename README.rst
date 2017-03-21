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
use Anaconda_, or to create a virtual environment usig venv_.

.. _Anaconda: https://www.continuum.io/downloads
.. _venv: http://docs.python-guide.org/en/latest/dev/virtualenvs/

Installation using Anaconda
~~~~~~~~~~~~~~~~~~~~~~

Follow the instruction to install Anaconda_. Most of the packages should be already installed.
If necessary install these packages using the Anaconda package manager.

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
compilation using Numba_. Since its installation is quite
complicate, it is recommended to use `Anaconda`:

.. _Numba: http://numba.pydata.org/

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

`Polimi` and `bioPolimi` models. require the use of Cantera_, while `CPD` works without
Cantera_.
Follow the instruction in the website to build Cantera_ in your system.

With anaconda, cantera can be installed:

.. code:: bash
	  # conda install -c cantera cantera=2.3.0

.. _Cantera: http://www.cantera.org/docs/sphinx/html/index.html

This command according to your system is not always working properly,
especially because of problems existing betweene difference `gcc` versions.


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
file. The Results will be saved in the directory `Results`.
If no directory is specified the results will be stored in `Results`.

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
