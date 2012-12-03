===========================
The FG-DVC specific Classes
===========================

The FG-FVC classes are contained in the files:

- FGDVC_SetAndLaunch.py  (Writes the FG-DVC 'instruct.ini' and launches FG-DVC)
- Compos_and_Energy.py  (The species and energy balance)
- FGDVC_Result.py  ('FGDVC_Result' reads the FG-DVC output and contains the output specific information.)

.. _SS-ReadGen:
.. _my-reference-label:

The Class writing the FG-DVC instruction file and launches FG-DVC
=================================================================

.. py:currentmodule:: FGDVC_SetAndLaunch.SetterAndLauncher


.. autoclass:: FGDVC_SetAndLaunch.SetterAndLauncher
    :members:


The Class generating the Coal File
----------------------------------

.. py:currentmodule:: InformationFiles.WriteFGDVCCoalFile


.. autoclass:: InformationFiles.WriteFGDVCCoalFile
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

.. py:currentmodule:: FGDVC_Result.FGDVC_Result

.. autoclass:: FGDVC_Result.FGDVC_Result
    :members:


