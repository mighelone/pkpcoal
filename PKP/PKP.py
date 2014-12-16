import sys
import os
import platform
import shutil
#PKP imports
import src.CPD_SetAndLaunch as CPD_SetAndLaunch   # writes CPD-instruct File, launches CPD
# import FGDVC_SetAndLaunch # writes FG-DVC-instruct File, launches FG-DVC and fittes using eq. (68 ) (BachelorThesis)
# import PCCL_SetAndLaunch  # writes PCCL instruction file and alunches the exe
# import FGDVC_Result       # contains the information of the FG-DVC output file
# import CPD_Result         # contains the information of the CPD output file
# import PCCL_Result        # contains the information of the PCCL output file
# #import src.Fitter             # The optimizer class
# import Models             # the models like Arrhenius or Kobayashi
# import FitInfo            # supports the Fitting with the yield information
# import Compos_and_Energy  # Species balance and energy balance for CPD and FG-DVC
# import InformationFiles   # reads the user input files, writes FG-DVC coalsd.exe coal generation file
#import src.GlobalOptParam     # contains the Information of the Number Of Runs for the Global Optimum search
# import Evolve             # contains the generic algortihm optimizer
#import coalPolimi

import matplotlib
import numpy as np
import pylab as plt
#
matplotlib.use('Qt4Agg')
matplotlib.rcParams['backend.qt4']='PySide'
os.environ['QT_API'] = 'pyside'
sys.path.append('src')
#
#Directories:
#gets the current directory:
workingDir=os.getcwd()+'/'
#FG-DVC Library coals folder name:
FG_LibCoalDir='input'
#FG-DVC coal generation (coalsd.exe) folder name:
FG_GenCoalDir='coals'
#FG-DVC fgdvcd.exe folder name:
FG_ExeCoalDir='FGDVC'
#Input file name for coalsd.exe:
FG_CoalGenFileName='CoalGen_FGDVC.txt'
#File name for generated coal file:
FG_CoalName='GenCoal'
#


class MainProcess(object):
    """ Controls the whole process of generating input files, fitting etc. """

    def __init__(self, read_inp_files=True, inputs_folder=False):
        self.SpeciesToConsider = [] #for GUI
        self.ProgramModelDict = {} #for GUI
        if read_inp_files:
            self.inputs = self.ReadInputFiles(inputs_folder)

 
    def ReadInputFiles(self, inputs_folder=False):
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
        return selector()(self.inputs) #


    def postProcessResults(self, results):
        import src.PyrolModelLauncher as pml
        solver = results[0].solver
        if self.inputs[solver]['fit'] not in pml.__dict__:
            print "Cannot find " + self.inputs[solver]['fit']
            return 
        return getattr(pml, self.inputs[solver]['fit'])(self.inputs, solver, results)
        


    def OptGenAlgBased(self, Fit, ParameterVecInit, ParameterVecMin, ParameterVecMax, Species):
        """ Starts a genetic algorithm and afterwards a gradient Based optimization. 
            Sets the Final Fit result as the ParamVector in the Kinetic Model. 
            Input are the Fit (Result Objects of the Detailed Models), 
            the Parameter to initialize, the two Parameter vectors defining 
            the range of the results and the Species index. 
        """
        GenAlg=Evolve.GenericOpt(self.KinModel,Fit,Species)
        GenAlg.setWeights(self.WeightY,self.WeightR)
        GenAlg.setParamRanges(ParameterVecInit,ParameterVecMin,ParameterVecMax)
        GenAlg.setNrPopulation(GlobalOptParam.NrOfPopulation)
        GenAlg.setNrGenerations(GlobalOptParam.NrOfGeneration)
        self.KinModel.setParamVector(GenAlg.mkResults())
        # afterwards grad based optimization
        if GlobalOptParam.optimizGrad == True:
            self.OptGradBased(Fit,ParameterVecInit,False,Species)


    def SpeciesEnergy(self,PyrolProgram,File,FittingModel):
        """ Carries out the species and Energy balance. """
        ##SPECIES AND ENERGY BALANCE:
        for runNr in range(self.NrOfRuns):
            if PyrolProgram=='CPD':
                print 'CPD energy and mass balance...'
                #TODO replace by dictionaries
                Compos_and_Energy.CPD_SpeciesBalance(
                    File[runNr],
                    self.UAC,
                    self.UAH,
                    self.UAN,
                    self.UAO,
                    self.UAS,
                    self.PAVM_asrec,
                    self.PAFC_asrec,
                    self.PAmoist,
                    self.PAash,
                    self.HHV,
                    self.MTar,
                    self.densityDryCoal,
                    runNr)
                
            if PyrolProgram=='FGDVC':    
                print 'FG-DVC energy and mass balance...'
                Compos_and_Energy.FGPC_SpeciesBalance(
                    File[runNr],
                    self.UAC,
                    self.UAH,
                    self.UAN,
                    self.UAO,
                    self.UAS,
                    self.PAVM_asrec,
                    self.PAFC_asrec,
                    self.PAmoist,
                    self.PAash,
                    self.HHV,
                    self.MTar,
                    self.densityDryCoal,
                    runNr,
                    'FGDVC')
#
            if PyrolProgram=='PCCL':
                print 'PCCL-Flashchain energy and mass balance...TO FINISH'
                Compos_and_Energy.FGPC_SpeciesBalance(
                    File[runNr],
                    self.UAC,
                    self.UAH,
                    self.UAN,
                    self.UAO,
                    self.UAS,
                    self.PAVM_asrec,
                    self.PAFC_asrec,
                    self.PAmoist,
                    self.PAash,
                    self.HHV,
                    self.MTar,
                    self.densityDryCoal,
                    runNr,
                    'PCCL')

        # new implementation of species energy for Kobayashi model Michele Vascellari
        if FittingModel == 'Kobayashi':
            print FittingModel
            self.extrapolateYieldKoba(PyrolProgram,File)

    def extrapolateYieldKoba(self,PyrolProgram, File):
        '''
        extrapolate the results of Detailed model to the alpha1/alpha2 parameters
        in the Kobayashi model
        It is required for species/energy calculation
        TODO only CPD is working at the moment!!!
        '''
        yields = []
        lenFile = len(File)
        Yields2Cols = File[0].DictYields2Cols()
        iSol= Yields2Cols['Solid']
        for i in range(lenFile):
            yields.append(File[i].FinalYields())
        lenSpecies = len(yields[0])
        if PyrolProgram == 'CPD':
            volatileYield = []
            for i in range(lenFile):
                volatileYield.append(1.-yields[i][iSol])
            volatileYield = np.array(volatileYield)
            maxVol = volatileYield.max()
            iMax = volatileYield.argmax()
            minVol = volatileYield.min()
            iMin = volatileYield.argmin()
            alpha1 = self.Solution[4]
            alpha2 = self.Solution[5]
            yields1 = yields[iMin] + (yields[iMax]-yields[iMin])/(maxVol-minVol)*(alpha1-minVol)
            # check results
            iTar= Yields2Cols['Tar']
            if (yields1[iTar]<0):
                ySolid =yields1[iSol]
                yTar = 0.5*(1-ySolid)
                index = [iSol,iTar,Yields2Cols['CO'],Yields2Cols['CO2'],Yields2Cols['H2O'],Yields2Cols['CH4'],Yields2Cols['Other']]
                sumy = yields1[index].sum()
                sumy = sumy - yields1[iTar]-yields1[iSol]
                factor = 1.0 + (yields1[iTar]-yTar)/sumy
                yields1[index] = yields1[index]*factor
                yields1[iTar]=yTar
                yields1[iSol]=ySolid

            yields2 = yields[iMin] + (yields[iMax]-yields[iMin])/(maxVol-minVol)*(alpha2-minVol)
            #print yields1
            #print yields2
            #print yields[iMin]
            #print yields[iMax]

            #print yields1[iTar],yields[iMin][iTar],yields[iMax][iTar],yields2[iTar]
            #print alpha1,minVol,maxVol,alpha2
            #print yields[iMin][iTar] + (yields[iMax][iTar]-yields[iMin][iTar])/(maxVol-minVol)*(alpha1-minVol)
            fit1 = CPD_Result.CPD_ResultFake(yields1)
            Compos_and_Energy.CPD_SpeciesBalance(fit1,self.UAC,self.UAH,self.UAN,self.UAO,self.UAS,self.PAVM_asrec,self.PAFC_asrec,self.PAmoist,self.PAash,
                                                 self.HHV,self.MTar,self.densityDryCoal,'Koba1')
            fit2 = CPD_Result.CPD_ResultFake(yields2)
            Compos_and_Energy.CPD_SpeciesBalance(fit2,self.UAC,self.UAH,self.UAN,self.UAO,self.UAS,self.PAVM_asrec,self.PAFC_asrec,self.PAmoist,self.PAash,
                                                 self.HHV,self.MTar,self.densityDryCoal,'Koba2')



    def fitResults(self, results):
        #TODO GO CPDFit comes from Fit_one_run.py
        if self.CPD_FittingKineticParameter_Select=='constantRate':
            self.MakeResults_CR('CPD', cpdResults)
            #currentDict={'CPD':'constantRate'}
        

        # TODO GO CPDFit comes from Fit_one_run.py
        # if self.CPD_FittingKineticParameter_Select=='constantRate':
        #     self.MakeResults_CR('CPD',CPDFile,CPDFit)
        #     currentDict={'CPD':'constantRate'}

        # elif self.CPD_FittingKineticParameter_Select=='Arrhenius':
        #     self.MakeResults_Arrh('CPD',CPDFile,CPDFit)
        #     currentDict={'CPD':'Arrhenius'}

        # elif self.CPD_FittingKineticParameter_Select=='ArrheniusNoB':
        #     self.MakeResults_ArrhNoB('CPD',CPDFile,CPDFit)
        #     currentDict={'CPD':'ArrheniusNoB'}

        # elif self.CPD_FittingKineticParameter_Select=='Kobayashi':
        #     self.MakeResults_Kob('CPD',CPDFile,CPDFit)
        #     currentDict={'CPD':'Kobayashi'}

        # elif self.CPD_FittingKineticParameter_Select=='DAEM':
        #     self.MakeResults_DEAM('CPD',CPDFile,CPDFit)
        #     currentDict={'CPD':'DAEM'}

        # elif self.CPD_FittingKineticParameter_Select==None:
        #     currentDict={'CPD':'None'}
        # else:
        #     print 'unspecified CPD_FittingKineticParameter_Select'
        #     currentDict={}
        # #
        # self.ProgramModelDict.update(currentDict)
        # #
        # self.SpeciesEnergy('CPD',CPDFile,self.CPD_FittingKineticParameter_Select)



#Main Part starting

def main():
    Case = MainProcess(inputs_folder=workingDir+"/inputs/")
    results = Case.executeSolver()
    Case.postProcessResults(results)
    print 'calculated Species: ',Case.SpeciesToConsider
    

if __name__ == "__main__":
    main()
