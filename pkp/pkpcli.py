import sys
import os
import platform

#PKP imports
import numpy as np
#
sys.path.append('src')

class BaseProcess(object):

    def __init__(self, inputs_dict, output_folder=False):
        self.inputs = inputs_dict
        self.output_folder = output_folder

    def returnResults(self):
        pass

class Generate(BaseProcess):
    """ Controls the whole process of generating input files, fitting etc. """

    SpeciesToConsider = []

    def executeSolver(self):
        """ execute all solver that are activated in the inputs file
            and return list of results objects
        """
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
        return selector()(self.inputs, self.output_folder)


class Fit(BaseProcess):

    def startFittingProcedure(self, results, selectPyrolModel=False):
        """ starts the fitting procedure to calculate modeled rates and yields
            according to preprocessor results

            Parameters:
                results: an array of preprocessor results objects
                selectPyrolModel: selects a specific pyrolysis model
                        overriding selections from inputs file

            Returns an array of fitted pyrolysis model objects
        """
        # NOTE PyrolModelLauncher objects
        # start fitting procedure immedatily
        import src.PyrolModelLauncher as pml

        if not selectPyrolModel:
            if self.inputs.get('FIT'):
                fit = self.inputs['FIT']['Model']
            else:
                fit = 'NONE'
        else:
            fit = selectPyrolModel

        if fit == 'NONE':
            print "No pyrolysis model for fitting procedure selected."
            return
        elif fit not in pml.__dict__:
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

def ReadInputFile(input_file):
    """ Read params from input file and generate Input objects """
    import yaml
    try:
        file_stream = open(input_file, 'r')
        inp = yaml.load(file_stream)
        return inp
    except Exception as e:
        print e
        sys.exit(1)

def generate(input_file=False, json_string=False, output_folder=False):
    """ a factory method for the Generate class
        Returns result object from the solver
    """
    inputs = (json_string if json_string
                else ReadInputFile(input_file))
    gen = Generate(inputs, output_folder)
    return gen.executeSolver()
    # fittedModels = Case.startFittingProcedure(results)
    # print 'calculated Species: ',Case.SpeciesToConsider
    # self.plotResults(preProcResults, fittedModels)

def fit(input_file=False, selectPyrolModel=None, json_string=False):
    """ a factory method for the Fit class
        accepts a results file or a results  Object
    """
    inputs = (json_string if json_string else ReadInputFile(input_file))
    return Fit(inputs)
