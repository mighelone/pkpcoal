TODO list
=========

New features
------------

This is the list of the new features to implement:

* Gradient based optimization with scipy
* Check if the multiobjective optimization can work better than single
  OF

Small additions/modifications
-----------------------------

This is a list of small fix:

* Use cxBlend as `mate` operator in Evolution (see `kursawefct.py`
  example given in **DEAP** `examples/ga`)
* Use decorator to limit range of parameters (see always `kursawefct`)
* Select `Evolution` or `EvolutionBinary` class for fitting
* Move :python:`fit_results['best'] =
  dict(zip(ga.empirical_model.parameters_names, best))` from
  `__init__.py:452` to `Evolution` class.
* Add references of papers used for **PKP**
* Use `array` instead of list for `individuals` in **DEAP**
