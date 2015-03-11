"""
A selection of launcher functions for different coal preprocessors

Launcher functions are  called with the inputs dictionary to pass
all definitions from the input/inputs file

The launcher function return a list of results objects
"""

from Models import BalancedComposition
import CPD

def Launch_CPD(inputs, output_folder):
    """ Execute CPD for each given temperature profile and
        return a list of CPD results objects
     """
    def InitAndLaunch(*pargs,**kwargs):
        """ initialises and execute cpd calculation """
        print 'Running CPD: ' + kwargs.get('runNr')
        cpd = CPD.CPD(*pargs)
        return cpd.Run()

    operatingConditions = inputs['OperatingConditions']
    #pressure = operatingConditions.pop('pressure')
    pressure = operatingConditions['pressure']
    from CoalThermoPhysics import Coal
    coal = Coal(inputs["Coal"])
    return {'run'+str(run):(operatingConditions["run"+str(run)], InitAndLaunch(
                coal,
                operatingConditions["run"+str(run)],
                pressure,
                inputs['CPD']['deltaT'],
                run,
                output_folder,
                runNr=str(run))
        ) for run in range(operatingConditions['runs'])}

