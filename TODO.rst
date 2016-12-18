TODO list
=========

New features
------------

This is the list of the new features to implement:

* Gradient based optimization with scipy **DONE**
* Check if the multiobjective optimization can work better than single
  OF
* Implement CPD **IN PROGRESS**
  * Check unit conversions in datafile

Small additions/modifications
-----------------------------

This is a list of small fix:

* Use cxBlend as `mate` operator in Evolution (see `kursawefct.py`
  example given in **DEAP** `examples/ga`) **DONE**
* Use decorator to limit range of parameters (see always `kursawefct`) **DONE**
* Select `Evolution` or `EvolutionBinary` class for fitting
* Move :python:`fit_results['best'] =
  dict(zip(ga.empirical_model.parameters_names, best))` from
  `__init__.py:452` to `Evolution` class.
* Add references of papers used for **PKP**
* *DONE* Use `array` instead of list for `individuals` in **DEAP** 
* For **SFOR** and **DAEM** which use constant :math:`y_0` force to use
  a narrow range of this parameter (i.e. from min to max value in the
  species)
* Set a warning when CPD is out of the correlation
* Set a warning when optimization reach the limits of the parameters
  range
* Split report files (maybe)
* Increase value in `cxBlend` to 0.25, as suggested by
  http://www.geatbx.com/docu/algindex-03.html#P550_28854
* Select limited number of points from reference solution, or perform
  a linear regression.
* Keep -x permission for CPD exe files
