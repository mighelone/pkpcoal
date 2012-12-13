========================
The CPD specific Classes
========================

The CPD classes are contained in the files:

- CPD_SetAndLaunch.py  (Writes the CPD input file and launches the program)
- Compos_and_Energy.py  (The species and energy balance)
- CPD_Result.py (Reads the specific output format of CPD) 

.. _SS-ReadGen:
.. _my-reference-label:

The Class writing the CPD instruction File and launching CPD
============================================================

.. py:currentmodule:: CPD_SetAndLaunch.SetterAndLauncher


.. autoclass:: CPD_SetAndLaunch.SetterAndLauncher
    :members:


The Class making the Species and the Energy Balance
===================================================

.. py:currentmodule:: Compos_and_Energy.SpeciesBalance

.. autoclass:: Compos_and_Energy.SpeciesBalance
    :members:

.. py:currentmodule:: Compos_and_Energy.CPD_SpeciesBalance

.. autoclass:: Compos_and_Energy.CPD_SpeciesBalance
    :members:
    :private-members:


The Reading Class containg the CPD specific output information
==============================================================

.. py:currentmodule:: CPD_Result.CPD_Result

.. autoclass:: CPD_Result.CPD_Result
    :members:


