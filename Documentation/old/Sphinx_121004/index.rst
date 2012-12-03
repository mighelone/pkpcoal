.. Xml_Documentation documentation master file, created by
   sphinx-quickstart on Mon Jul 16 15:52:40 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

The PKP Code Documentation
==========================

This is the documentation of the PKP program. The general structur of the classes is the following:
- There are pyrolysis program (like CPD and FG-DVC) specific ones to write the configuration files, launch their program and read the specific output files and suppport the other classes with the pyrolysis program specific data like species-column. Every of the program has also it's own class calculating the species and energy balance.
-The class independent Fitting procedures.

The manual is also structured the same way. Firstly the specific classes for CPD and FG-DVC and afterwards the fitting classes.

Required packages for PKP are (actual all packages except of scipy and numpy were already automatically installed with python):

- scipy
- numpy
- os
- StringIO
- pylab
- platform


Contents:

.. toctree::

  CPDClasses.rst
  FGDVCClasses.rst
  FittingClasses.rst
  MainProgramCode.rst
  Appendix.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

