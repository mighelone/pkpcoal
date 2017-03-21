.. _runPKP-label:

Run PKP
=======

This example shows how to run PKP and check results.
The input file ``input_DAEM.yml`` is used as input file.

The best practice is to first run PKP without performing the calibration. In order to avoid to change the input file, it is possible to do from the command-line::

  # runPKP input_DAEM.yml --run-only -o Results_DAEM

  -----------------------------------------

  .______    __  ___ .______   
  |   _  \  |  |/  / |   _  \  
  |  |_)  | |  '  /  |  |_)  | 
  |   ___/  |    <   |   ___/  
  |  |      |  .  \  |  |      
  | _|      |__|\__\ | _|      
  
  
  -----------------------------------------
  Pyrolysis Kinetic Preprocessor (PKP) 
  -----------------------------------------
  (c) Numerical Thermo Fluid-Dynamics      
      TU Bergakademie Freiberg             
      Michele Vascellari                   
      Michele.Vascellari@vtc.tu-freiberg.de
      -----------------------------------------
      Run PKP version 2.5+3.g93ff014.dirty
      -----------------------------------------
      pkp.runner.PKPRunner:run:Run model CPD
      pkp.runner.PKPRunner:run_model:Run 0 with CPD model
      pkp.runner.PKPRunner:run_model:Run 1 with CPD model
      pkp.runner.PKPRunner:run_model:Run 2 with CPD model
      pkp.runner.PKPRunner:run_model:Run 3 with CPD model
      pkp.runner.PKPRunner:run:Start fit of CPD model

The results of the CPD simulations can be viewed on the Results-DAEM directory. The results are saved in CSV format under the names ``Pittsburg-CPD-run0.csv`` and plotted in ``Pittsburg-CPD-run0.png``. For example the yields of the first run are showed below:

.. image:: Results-DAEM/Pittsburg-CPD-run0.png

Since the final yield of the overall volatiles is between 0.5 and 0.6, therefore the parameter ``Y`` is set in the input file to these range.
The position of the used coal is reported in the Van Kravelen diagram, showing the reference coals used for prescribing the gas composition:

.. image:: Results-DAEM/Pittsburg-CPD-run0-van_kravelen.png

In this case the gas composition is taken from coals 3, 4, 5.

Now, **PKP** runs the optimization::

  # runPKP input_DAEM.yml -o Results_DAEM
   ...
   ...
   pkp.runner.PKPRunner:fit_single:fit0 Evolution to fit CPD with DAEM
   gen	nevals	avg     	std     	min      	max     
   0  	30    	0.157593	0.136142	0.0143846	0.513655
   1  	40    	0.0569867	0.0504176	0.0066357	0.254772
   ...
   19 	40    	0.00161547	1.65685e-07	0.00161539	0.00161607
   20 	40    	0.00161541	2.0342e-08 	0.00161538	0.00161547
   pkp.runner.PKPRunner:evolve:Best population: {'sigma': (16770044.146681631, 'J/kmol'), 'y0': (0.55777098794995306, '-'), 'A0': (2379707227.2597728, '1/s'), 'E0': (127775589.0628166, 'J/kmol')}
   pkp.runner.PKPRunner:fit_single:fit0 Minimization to fit CPD with DAEM
   Optimization terminated successfully.
         Current function value: 0.000905
         Iterations: 31
         Function evaluations: 210
         Gradient evaluations: 35
	 pkp.runner.PKPRunner:minimization:Minimized value: {'sigma': (29925317.263066139, 'J/kmol'), 'y0': (0.55962327063440376, '-'), 'A0': (78768837802240.156, '1/s'), 'E0': (216567786.08743358, 'J/kmol')}


The evolution algorithm evolves for 20 generations, obtaining a best solution with a penalty function of ``0.00161538``.
The evolution history can be find in ``Pittsburg-fit0-CPD-DAEM-evolution.png``:

.. image:: Results-DAEM/Pittsburg-fit0-CPD-DAEM-evolution.png

Using the ``method: evolve_min`` a small number of generations is recommended only to define a good starting point for the following minimization step using BFGS.
The penalty function obtained from the BFGS is ``0.000905``, which shows a small increase from the evolutionary algorithm.

The results of the calibration can be observed from ``Pittsburg-fit0-CPD-DAEM-yields.png``:

.. image:: Results-DAEM/Pittsburg-fit0-CPD-DAEM-yields.png

The plot shows the comparison between the results of the CPD model (solid lines) and of the calibrated DAEM model (dashed lines) for the 4 runs. Note that only run0, 1 and 2 are used for the calibration, while run 4 is only reported for checking results of the calibration outside the range of validation.

Finally, the sum-up of the calibration precedure is reported in ``Pittsburg-fitreport.yml``.
