===============
The Appendix
===============

As an alternative to the Least Square optimizaion, a linear regression method (see: 'Holstein, A., Bassilakis, R., Wojtowicz, M.A., and Serio, M.A. “Kinetics of methane and tar
evolution during coal pyrolysis”. In: Proceedings of the Combustion Institute 30 (2005), pp. 2177–2185') was tested. But as this shows not that good results as the Least Square optimization, refering to computational time, the precision of the calculated parameter and the reliability (the results were quite dependent on the time step and the number and range of the runned heating rates), this idea and the realized classes are still in the Python-files but still not implemented in the main file, the 'Pyrolysis.py'. Another problem was that the oszillating yield curve of CPD cannot be used to generate reasonable results.
Here an overview of the classes and their Methods:

The linear Regression Class for FG-DVC
======================================

.. py:currentmodule:: FGDVC_Fit_lin_regr.Process


.. autoclass:: FGDVC_Fit_lin_regr.Process
    :members:


The linear Regression Class for CPD
===================================

.. py:currentmodule:: CPD_Fit_lin_regr.ProcessCPD


.. autoclass:: CPD_Fit_lin_regr.ProcessCPD
    :members:



