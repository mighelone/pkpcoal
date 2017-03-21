Run PKP
=======

This example shows how to run PKP and check results.
The input file `input_DAEM.yml` is used as input file.

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



