===============
The CPD Classes
===============

The CPD classes are contained in the files:
- CPD_Fit_lin_regr.py  (Writes the CPD input file and launches of the program)
- CPD_Compos_and_Energy.py  (The species and energy balance)
- Fit_one_run.py  ('CPD_Result' reads the CPD output and contains the output specific information.)

.. _SS-ReadGen:
.. _my-reference-label:

The Class writing the CPD instruction File and launching CPD
============================================================

.. py:currentmodule:: CPD_Fit_lin_regr.SetterAndLauncher


.. autoclass:: CPD_Fit_lin_regr.SetterAndLauncher
    :members:


The Class making the Species and the Energy Balance
===================================================

.. py:currentmodule:: CPD_Compos_and_Energy.SpeciesBalance

.. autoclass:: Compos_and_Energy.SpeciesBalance
    :members:

.. py:currentmodule:: CPD_Compos_and_Energy.CPD_SpeciesBalance

.. autoclass:: Compos_and_Energy.CPD_SpeciesBalance
    :members:
    :private-members:


The Reading Class containg the CPD specific output information
==============================================================

.. py:currentmodule:: Fit_one_run.CPD_Result

.. autoclass:: Fit_one_run.CPD_Result
    :members:


