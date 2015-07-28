"""
A selection of launcher functions for different coal preprocessors

Launcher functions are  called with the inputs dictionary to pass
all definitions from the input/inputs file

The launcher function return a list of results objects
"""

from Models import BalancedComposition
from PreProc import MultiRun
import CPD

def Launch_CPD(inputs, output_folder):
    """ Execute CPD for each given temperature profile and
        return a list of CPD results objects
     """
    def InitAndLaunch(*pargs):
        """ initialises and execute cpd calculation """
        cpd = CPD.CPD(*pargs)
        return cpd.Run()

    operatingConditions = inputs['Coal']['OperatingConditions']
    #pressure = operatingConditions.pop('pressure')
    pressure = operatingConditions['pressure']
    from CoalThermoPhysics import Coal
    coal = Coal(inputs["Coal"])
    # Check if number of runs is demanded explicitly
    # if not iterate over all run keys
    runs_exp = operatingConditions.get('runs')
    if runs_exp:
        runs = ['run' + str(i) for i in range(runs_exp)]
    else:
        runs = [runNr for runNr in operatingConditions.keys()
                    if 'run' in runNr ]
    return MultiRun({run: InitAndLaunch(
                coal,
                pressure,
                inputs['CPD']['deltaT'],
                run,
                output_folder)
          for i, run in enumerate(runs)})
