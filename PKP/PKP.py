#!/usr/bin/python
"""
PKP

Usage:
      PKP.py -h | --help
      PKP.py generate (--json-input=<string> | --file-input=<loc>) --results-folder=<loc>
      PKP.py fit      (--json-input=<string> | --file-input=<loc>) --fit-target=<string> --results-folder=<loc>

Options:
    -h  --help              Shows this screen
    --json-input=<string>   Foo
    --file-input=<loc>      Bar
    --results-folder=<loc>  Baz
    --fit-target=<string>   Bla

"""

import sys
import os
import platform
from docopt import docopt

#PKP imports
import src.CPD as CPD  # writes CPD-instruct File, launches CPD

import matplotlib
import numpy as np
import pylab as plt
#
sys.path.append('src')
#
#Directories:
#gets the current directory:
workingDir=os.getcwd()+'/'


class BaseProcess(object):

    def __init__(self, inputs_dict, output="Folder"):
        self.inputs = inputs_dict
        self.output = output

    def returnResults(self):
        pass

class Generate(BaseProcess):
    """ Controls the whole process of generating input files, fitting etc. """

    SpeciesToConsider = []

    def executeSolver(self):
        """ execute all solver that are activated in the inputs file
            and return list of results objects
        """
        print self.inputs
        from src import PreprocLauncher as Launcher
        def selector():
            if self.inputs['CPD']['active']:
                return Launcher.Launch_CPD
            # TODO GO Reimplement other solvers
            # if Case.FG_select==True:
            #     Case.CheckFGdt()
            #     Case.MakeResults_FG()
            # if Case.PMSKD_select==True:
            #     Case.RunPMSKD()
            # if Case.PCCL_select==True:
            #     Case.MakeResults_PCCL()
        return selector()(self.inputs) #


class Fit(BaseProcess):

    def startFittingProcedure(self, results, selectPyrolModel = False):
        """ starts the fitting procedure to calculate modeled rates and yields
            according to preprocessor results

            Parameters:
                results: an array of preprocessor results objects
                selectPyrolModel: selects a specific pyrolysis model
                        overriding selections from inputs file

            Returns an array of fitted pyrolysis model objects
        """
        import src.PyrolModelLauncher as pml
        solver = results[0].solver
        fit = (self.inputs[solver].get('fit','NONE')
                if not selectPyrolModel else selectPyrolModel)
        if fit not in pml.__dict__:
            print "Cannot find " + fit
            return
        return getattr(pml, fit)(self.inputs, results)


    def plotResults(self, preProcResults, fittedModels):
        """ Creates plots of preProcResults against fittedModels
        """
        # NOTE: we implement a simple pyplot version for now
        #       finally a way to select different plotting backends
        #       would be desirable
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots()
        colors = ['c', 'm', 'y', 'k']
        marker = ['.','x','+']
        for color, preProc in zip(colors, preProcResults[0].iterspecies()):
            name, data =  preProc
            for i, run in enumerate(preProcResults):
                axs.scatter(
                         x=run['time(ms)']*1e-3,
                         y=run[name],
                         color = color,
                         label = name+ "_run_" + str(i),
                         marker = marker[i],
                    )

            model_data = fittedModels[name]
            axs.plot(
                 preProcResults[-1]['time(ms)']*1e-3,
                 model_data.mass,
                 color = color,
                 label = name,
                 linewidth = 2
            )

        plt.legend()
        plt.show(fig)


def ReadInputFiles(inputs_folder):
    """ Read params from input file and generate Input objects """
    import yaml
    # TODO GO print warning if it cant find file
    inputs_folder = (inputs_folder if inputs_folder else workingDir + 'inputs/')
    full_fn = inputs_folder + 'inputs.inp'
    try:
        file_stream = open(full_fn, 'r')
        inp = yaml.load(file_stream)
        return inp
    except Exception as e:
        print e

def generate(folder=False, json_string=False):
    """ a factory method for the  Generate class
        Returns result object from the solver
    """
    inputs = (json_string if json_string
                else ReadInputFiles(folder))
    gen = Generate(inputs)
    return gen.executeSolver()
    # fittedModels = Case.startFittingProcedure(results)
    # print 'calculated Species: ',Case.SpeciesToConsider
    # self.plotResults(preProcResults, fittedModels)

def fit(folder=False, results=False, selectPyrolModel=None, json_string=False):
    """ a factory method for the Fit class
        accepts a results file or a results  Object
    """
    inputs = (json_string if json_string
                else ReadInputFiles(folder))
    fit = Fit(inputs)
    fit.pyrolModel = "constantRate"
    return fit.startFittingProcedure(results)


if __name__ == "__main__":
    arguments = docopt(__doc__)
    print arguments
    if arguments['generate']:
        generate(json_string=arguments['--json-input'], folder=arguments['--file-input'])
    main()
