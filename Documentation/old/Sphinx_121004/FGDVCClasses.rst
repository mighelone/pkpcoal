==================
The FG-DVC Classes
==================

The FG-FVC classes are contained in the files:
- FGDVC_Fit_lin_regr.py  (Writes the FG-DVC 'instruct.ini' and launches FG-DVC)
- FGDVC_Compos_and_Energy.py  (The species and energy balance)
- Fit_one_run.py  ('FGDVC_Result' reads the FG-DVC output and contains the output specific information.)

.. _SS-ReadGen:
.. _my-reference-label:

The Class writing the FG-DVC instruction file and launches FG-DVC
=================================================================

.. py:currentmodule:: FGDVC_Fit_lin_regr.SetterAndLauncher


.. autoclass:: FGDVC_Fit_lin_regr.SetterAndLauncher
    :members:


The Class generating the Coal File
----------------------------------

.. py:currentmodule:: ReadInputFiles.WriteFGDVCCoalFile


.. autoclass:: ReadInputFiles.WriteFGDVCCoalFile
    :members:


The Class making the Species and the Energy Balance
===================================================

.. py:currentmodule:: Compos_and_Energy.SpeciesBalance

.. autoclass:: Compos_and_Energy.SpeciesBalance
    :members:

.. py:currentmodule:: Compos_and_Energy.FGDVC_SpeciesBalance

.. autoclass:: Compos_and_Energy.FGDVC_SpeciesBalance
    :members:
    :private-members:


The Reading Class containg the CPD specific output information
==============================================================

.. py:currentmodule:: Fit_one_run.FGDVC_Result

.. autoclass:: Fit_one_run.FGDVC_Result
    :members:


