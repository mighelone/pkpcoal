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
#Which operating Sytem?
oSystem=platform.system()
#if oSystem == 'Darwin':
#    oSystem = 'Linux'
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

class BalancedComposition(object):
    """ Class for compostion that ensures componenents sum up to a certain value """

    def __init__(self, inp, target=100.00):
        """ From a given input dictionary and a target sum a scaled dictionary is created """
        self.target = target
        scaling_factor = target/sum(inp.values())
        self.elems = {key:value*scaling_factor for key,value in inp.iteritems()} 


    def __getitem__(self,item):
        """ """
        try:
            return self.elems[item]
        except:
            print """Warning trying to access {} which was not set in
            input file, assuming 0.0. 
            """ 
            return 0.0

    def remove_elems_rebalance(self, elems_):
        """ To compute daf composition elements can be removed 

            Usage:  .remove_elems_rebalance(['Moisture','Ash'])
        """
        return BalancedComposition({ 
            key:elem for key,elem in self.elems.iteritems() 
                if key not in elems_
        })

    def scale(self, factor):
        return BalancedComposition(self.elems, target=self.target*factor)

    def __repr__(self):
        return str(self.elems)

class MainProcess(object):
    """ Controls the whole process of generating input files, fitting etc. """

    def __init__(self, read_inp_files=True, inputs_folder=False):
        self.SpeciesToConsider = [] #for GUI
        self.ProgramModelDict = {} #for GUI
        if read_inp_files:
            self.inputs = self.ReadInputFiles(inputs_folder)
        coal = self.inputs["Coal"]
        self.proxi_ana = BalancedComposition(coal["Proximate Analysis"])
        self.ultim_ana = BalancedComposition(coal["Ultimate Analysis"])
        self.daf = self.proxi_ana.remove_elems_rebalance(['Moisture','Ash'])

 
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
        def selector(inputs):
            if inputs['CPD']['active']:
                return self.MakeResults_CPD
            # TODO GO Reimplement other solvers
            # if Case.FG_select==True:
            #     Case.CheckFGdt()
            #     Case.MakeResults_FG()
            # if Case.PMSKD_select==True:
            #     Case.RunPMSKD()
            # if Case.PCCL_select==True:
            #     Case.MakeResults_PCCL()
        return selector(self.inputs)() #


    def postProcessResults(self, results):
        for run in results:
            solver = run.solver
            getattr(self, self.inputs[solver]['fit'])(solver, run)
        

    # def CheckFGdt(self):
    #     """ Aborts, if FG-DVC is selected and the timestep is lower than 1.e-3 
    #     """
    #     # TODO GO move to FGDVC class
    #     # TODO check if this function gets callesd ? Can't we use a limiter here ?
    #     if ((self.FG_select==True) and (self.FG_dt<1e-4)):
    #         print """Please select for FG-DVC a time step greather equal 1e-4 in 'OperCond.inp'. 
    #         FG-DVC would not be able to read the time history file for a dt<1e-4. """
    #         sys.exit()

    def OptGradBased(self, Fit, ParameterVecInit, FinalYield, Species):
        """ Starts a gradient Based Optimization. 

        returns the final Fit. 
        Input are the Fit (Result Objects of the Detailed Models), the Parameter 
        to initialize and the Final Yield (all dependent on the kinetic model). 
        For Kobayashi Model set Final Yield to False (independent), for all other set a value. 
        It will be excluded from the optimization. Species is the Species Index. """
        # LS = Fitter.LeastSquarsEstimator(
        #         optimizer = self.inputs['Optimisation']['GradBasedOpt'],
        #         maxIter = self.inputs['Optimisation']['maxIter'],
        #         fitTolerance = self.inputs['Optimisation']['Tolerance'],
        #         weights = self.WeightY,self.WeightR,
        #         finalYield = FinalYield,
        # )
        # result = LS.estimateT(
        #     Fit,self.KinModel,ParameterVecInit,Species)
        # print 'Final error=', LS.Deviation()
        # return result
        pass

    def OptGenAlgBased(self,Fit,ParameterVecInit,ParameterVecMin,ParameterVecMax,Species):
        """ Starts a genetic algorithm and afterwards a gradient Based optimization. 
            Sets the Final Fit result as the ParamVector in the Kinetic Model. 
            Input are the Fit (Result Objects of the Detailed Models), the Parameter to initialize,
            the two Parameter vectors defining the range of the results and the Species index. """
        GenAlg=Evolve.GenericOpt(self.KinModel,Fit,Species)
        GenAlg.setWeights(self.WeightY,self.WeightR)
        GenAlg.setParamRanges(ParameterVecInit,ParameterVecMin,ParameterVecMax)
        GenAlg.setNrPopulation(GlobalOptParam.NrOfPopulation)
        GenAlg.setNrGenerations(GlobalOptParam.NrOfGeneration)
        self.KinModel.setParamVector(GenAlg.mkResults())
        # afterwards grad based optimization
        if GlobalOptParam.optimizGrad == True:
            self.OptGradBased(Fit,ParameterVecInit,False,Species)

    def constantRate(self, PyrolProgram, Results):
        """ Generates the results for constant Rate.

            Parameters:
                File:        results file from the prepocessor i.e CPD_Results.dat
                PyrolPrgram: name of the preprocessor
                Fit is the fit object
         """
        # TODO GO delay writing of the file until end of computation
        #         to spearate responsibilities
        # outf = PyrolProgram+'-Results_const_rate.txt'
        # outfile = open(outf, 'w')
        # outfile.write("Species     k [1/s]     t_start [s]   FinalYield\n\n")
        # init model
        KinModel = Models.ConstantRateModel(
            self.inputs['Optimisation']['ConstantRate'])
        # if PyrolProgram =='PCCL': # TODO GO reimplement
        #         self.KinModel.setDt4Intergrate(self.FG_dt)
        # uses for optimization gradient based (LS) optimizer if NrOfRuns is one ; 
        # for more runs use Evolutionary algorithm (GenAlg; global optimum)
        if self.inputs['OperatingConditions']['runs'] == 1: 

            # Iterate species in Fit Object
            # for Spec in range(2, len(Fit[0].SpeciesNames()),1):
            for Spec in Results.species():
                if Spec not in self.SpeciesToConsider:
                    self.SpeciesToConsider.append(Spec)
                #
                self.OptGradBased(Fit, KinModel.initialParameter, Fit[0].Yield(Spec)[-1], Spec)
                # TODO GO cant the plotting be delayed until computation is done
                # self.KinModel.plot(Fit, Spec)
                # TODO GO why are parameters summed ?
                if np.sum(KinModel.parameter) != np.sum(KinModel.initialParamter):
                    # if true nothing was optimized, no result to show
                    outfile.write(
                        str(Fit[0].SpeciesName(Spec)) 
                        + '\t'+'%8.4f  %8.4f  %8.4f  ' %( 
                            self.Solution[0], 
                            self.Solution[1], 
                            self.Solution[2])+'\n')

        # else:
        #     if len(ParamInit) == 2:
        #         ParamInit.append(0.0)
        #     ParamMin = GlobalOptParam.EvACRMin
        #     if len(ParamMin) == 2:
        #         ParamMin.append(0.0)
        #     ParamMax = GlobalOptParam.EvACRMax
        #     if len(ParamMax) == 2:
        #         ParamMax.append(0.0)
        #     for Spec in range(2,len(Fit[0].SpeciesNames()),1):
        #         # max Yield, min Yield
        #         m_final_predictionAll = []
        #         for i in range(len(Fit)):
        #             m_final_predictionAll.append(Fit[i].Yield(Spec)[-1])
        #         ParamInit[2]= (max(m_final_predictionAll)
        #                     + min(m_final_predictionAll))/2.
        #         ParamMin[2] = (min(m_final_predictionAll))
        #         ParamMax[2] = (max(m_final_predictionAll))
        #         #
        #         self.OptGenAlgBased(Fit,ParamInit,ParamMin,ParamMax,Spec)
        #         self.Solution = self.KinModel.ParamVector()
        #         self.KinModel.plot(Fit,Spec)
        #         #if True nothing was optimized
        #         if np.sum(self.Solution) != np.sum(ParamInit):
        #             outfile.write(str(Fit[0].SpeciesName(Spec)) 
        #                 + '\t'+'%8.4f  %8.4f  %8.4f  ' %(
        #                     self.Solution[0],
        #                     self.Solution[1],
        #                     self.Solution[2])+'\n')
        # outfile.close()
        # if oSystem=='Linux' or oSystem == 'Darwin':
        #     shutil.move(outf, 'Result/' + outf)
        # elif oSystem=='Windows':
        #     shutil.move(outf, 'Result\\' + outf)
        # else:
        #     print "The name of the operating system couldn't be found."


    def MakeResults_Arrh(self,PyrolProgram,File,Fit):
        """Generates the results for Arrhenius Rate."""
        #TODO Giant spaghetti factory
        outfile = open(PyrolProgram+'-Results_ArrheniusRate.txt', 'w')
        outfile.write("Species                         A [1/s]         b               E_a [K]    FinalYield\n\n")
        #makes Species list which contains alls species to fit:
        SpeciesList=[] # List containing the Species Index to fit
        if self.ArrhSpec=='Total':
            SpeciesList.append(Fit[0].SpeciesIndex('Total'))
            if 'Total' not in self.SpeciesToConsider:
                self.SpeciesToConsider.append('Total')
        elif self.ArrhSpec=='MainSpecies':
            SpeciesList.append(Fit[0].SpeciesIndex('Total'))
            SpeciesList.append(Fit[0].SpeciesIndex('Tar'))
            SpeciesList.append(Fit[0].SpeciesIndex('Gas'))
            if 'Total' not in self.SpeciesToConsider:
                self.SpeciesToConsider.append('Total')
            if 'Tar' not in self.SpeciesToConsider:
                self.SpeciesToConsider.append('Tar')
            if 'Gas' not in self.SpeciesToConsider:
                self.SpeciesToConsider.append('Gas')
        elif self.ArrhSpec=='allSpecies':
            for i in range(2,len(Fit[0].SpeciesNames()),1):
                if Fit[0].SpeciesName(i) not in self.SpeciesToConsider:
                    self.SpeciesToConsider.append(Fit[0].SpeciesName(i))
                SpeciesList.append(i)
        ##The single species:
        for Species in SpeciesList:
            #
            m_final_prediction=Fit[0].Yield(Species)[-1]
            PredictionV0=[0.86e15,0.01,27700,m_final_prediction]  #for Standard Arrhenius
            #
            self.KinModel=Models.ArrheniusModel(PredictionV0)
            if PyrolProgram=='PCCL':
                self.KinModel.setDt4Intergrate(self.FG_dt)
            #
            print Fit[0].SpeciesName(Species)
            ParamInit = GlobalOptParam.EvAArrhInit
            if len(ParamInit) == 4:
                ParamInit.pop(-1)
            #
            if self.NrOfRuns == 1:
                self.OptGradBased(Fit,ParamInit,Fit[0].Yield(Species)[-1],Species)
            else:
                # init the Parameter for global optimization
                m_final_predictionAll=[] # final yields
                for i in range(len(Fit)):
                    m_final_predictionAll.append(Fit[i].Yield(Species)[-1])
                ParamMin = GlobalOptParam.EvAArrhMin
                ParamMax = GlobalOptParam.EvAArrhMax
                if len(ParamMin) == 3:
                    ParamMin.append(0.0)
                if len(ParamMax) == 3:
                    ParamMax.append(0.0)
                if len(ParamInit) == 3:
                    ParamInit.append(0.0)
                ParamInit[3] = (max(m_final_predictionAll)+min(m_final_predictionAll))/2.
                ParamMin[3] = (min(m_final_predictionAll))
                ParamMax[3] = (max(m_final_predictionAll))
                #
                self.OptGenAlgBased(Fit,ParamInit,ParamMin,ParamMax,Species)
            #
            self.KinModel.plot(Fit,Species)
            self.Solution=self.KinModel.ParamVector()
            #To avoid, a species with no yield is added to the parameter file
            if np.sum(self.KinModel.ParamVector())!=np.sum(PredictionV0): 
                outfile.write(str(Fit[0].SpeciesName(Species))+'\t'+'%.6e  %6.4f  %11.4f  %7.4f  ' %(self.Solution[0],self.Solution[1],self.Solution[2],self.Solution[3])+'\n')
        outfile.close()
        if oSystem=='Linux' or oSystem == 'Darwin':
            shutil.move(PyrolProgram+'-Results_ArrheniusRate.txt','Result/'+PyrolProgram+'-Results_Arrhenius.txt')
        elif oSystem=='Windows':
            shutil.move(PyrolProgram+'-Results_ArrheniusRate.txt','Result\\'+PyrolProgram+'-Results_Arrhenius.txt')
        else:
            print "The name of the operating system couldn't be found."


    def MakeResults_ArrhNoB(self,PyrolProgram,File,Fit):
        """Generates the results for Arrhenius Rate with no correction term T**b."""
        outfile = open(PyrolProgram+'-Results_ArrheniusNoBRate.txt', 'w')
        outfile.write("Species               A [1/s]                E_a [K]      FinalYield\n\n")
        #makes Species list which contains alls species to fit:
        SpeciesList=[]
        if self.ArrhSpec=='Total':
            SpeciesList.append(Fit[0].SpeciesIndex('Total'))
            if 'Total' not in self.SpeciesToConsider:
                self.SpeciesToConsider.append('Total')
        elif self.ArrhSpec=='MainSpecies':
            SpeciesList.append(Fit[0].SpeciesIndex('Total'))
            SpeciesList.append(Fit[0].SpeciesIndex('Tar'))
            SpeciesList.append(Fit[0].SpeciesIndex('Gas'))
            if 'Total' not in self.SpeciesToConsider:
                self.SpeciesToConsider.append('Total')
            if 'Tar' not in self.SpeciesToConsider:
                self.SpeciesToConsider.append('Tar')
            if 'Gas' not in self.SpeciesToConsider:
                self.SpeciesToConsider.append('Gas')
        elif self.ArrhSpec=='allSpecies':
            for i in range(2,len(Fit[0].SpeciesNames()),1):
                if Fit[0].SpeciesName(i) not in self.SpeciesToConsider:
                    self.SpeciesToConsider.append(Fit[0].SpeciesName(i))
                SpeciesList.append(i)
        ##The single species:
        for Species in SpeciesList:
            m_final_prediction=Fit[0].Yield(Species)[-1]
            PredictionV0=[0.86e11,17700,m_final_prediction]  #for Standard Arrhenius
            self.KinModel=Models.ArrheniusModelNoB(PredictionV0)
            if PyrolProgram=='PCCL':
                self.KinModel.setDt4Intergrate(self.FG_dt)
            #
            print Fit[0].SpeciesName(Species)
            # gets final yields for all runs
            m_final_predictionAll=[]
            for i in range(len(Fit)):
                m_final_predictionAll.append(Fit[i].Yield(Species)[-1])
            ParamInit = GlobalOptParam.EvAArrhInit
            if len(ParamInit) == 4:
                ParamInit.pop(-1)
            if len(ParamInit) == 3:
                ParamInit.pop(1)
            # optimization procedure
            if self.NrOfRuns == 1:
                self.OptGradBased(Fit,ParamInit,Fit[0].Yield(Species)[-1],Species)
            else:
                ParamMin = GlobalOptParam.EvAArrhMin
                ParamMax = GlobalOptParam.EvAArrhMax
                if len(ParamMin) == 3:
                    ParamMin.pop(1)
                    ParamMin.append(0.0)
                if len(ParamMax) == 3:
                    ParamMax.pop(1)
                    ParamMax.append(0.0)
                if len(ParamInit) == 3: # if final yield was appended
                    ParamMax.pop(-1)
                if len(ParamInit) == 2:
                    ParamInit.append(0.0)
                ParamInit[2] = (max(m_final_predictionAll)+min(m_final_predictionAll))/2.
                ParamMin[2] = (min(m_final_predictionAll))
                ParamMax[2] = (max(m_final_predictionAll))
                #
                self.OptGenAlgBased(Fit,ParamInit,ParamMin,ParamMax,Species)
                #
            self.Solution=self.KinModel.ParamVector()
            self.KinModel.plot(Fit,Species)
            if np.sum(self.KinModel.ParamVector())!=np.sum(PredictionV0): #To avoid, a species with no yield is added to the parameter file
                outfile.write(str(Fit[0].SpeciesName(Species))+'\t'+'%.6e  %11.4f  %7.4f  ' %(self.Solution[0],self.Solution[1],self.Solution[2])+'\n')
        outfile.close()
        if oSystem=='Linux'  or oSystem == 'Darwin':
            shutil.move(PyrolProgram+'-Results_ArrheniusNoBRate.txt','Result/'+PyrolProgram+'-Results_ArrheniusNoB.txt')
        elif oSystem=='Windows':
            shutil.move(PyrolProgram+'-Results_ArrheniusNoBRate.txt','Result\\'+PyrolProgram+'-Results_ArrheniusNoB.txt')
        else:
            print "The name of the operating system couldn't be found."
    #
    #
    def MakeResults_Kob(self,PyrolProgram,File,Fit):
        """Generates the results for Kobayashi Rate."""
        PredictionVKob0=[7e5,8e7/8314.33,2.3e8,1.6e8/8314.33,0.4,0.9]
        outfile = open(PyrolProgram+'-Results_KobayashiRate.txt', 'w')
        outfile.write("Species              A1 [1/s]              E_a1 [kJ/mol]             A2 [1/s]             E_a2 [kJ/mol]      alpha1     alpha2\n\n")
        self.KinModel=Models.Kobayashi(GlobalOptParam.EvAKobInit) #(PredictionVKob0)
        if PyrolProgram=='PCCL':
                self.KinModel.setDt4Intergrate(self.FG_dt)
        #######
        ##The single species:
        if 'Total' not in self.SpeciesToConsider:
            self.SpeciesToConsider.append('Total')
        for Species in [Fit[0].SpeciesIndex('Total')]:
            # optimization procedure
            print Fit[0].SpeciesName(Species)
            if self.NrOfRuns == 1:
                self.OptGradBased(Fit,self.KinModel.ParamVector(),False,Species)
            else:
                self.OptGenAlgBased(Fit,self.KinModel.ParamVector(),GlobalOptParam.EvAKobMin,GlobalOptParam.EvAKobMax,Species)
            #
            self.Solution=self.KinModel.ParamVector()
            #
            self.KinModel.plot(Fit,Species)
            outfile.write(str(Fit[0].SpeciesName(Species))+'\t'+'%.6e  %11.4f  %.6e  %11.4f  %6.4f  %6.4f  ' %(self.Solution[0],self.Solution[1]*8314.33/1e6,self.Solution[2],self.Solution[3]*8314.33/1e6,self.Solution[4],self.Solution[5])+'\n')
        outfile.close()
        if oSystem=='Linux'  or oSystem == 'Darwin':
            shutil.move(PyrolProgram+'-Results_KobayashiRate.txt','Result/'+PyrolProgram+'-Results_Kobayashi.txt')
        elif oSystem=='Windows':
            shutil.move(PyrolProgram+'-Results_KobayashiRate.txt','Result\\'+PyrolProgram+'-Results_Kobayashi.txt')
        else:
            print "The name of the operating system couldn't be found."
#
#
    def MakeResults_DEAM(self,PyrolProgram,File,Fit):
        """Generates the results for DAEM model."""
#    elif (CPD_FittingKineticParameter_Select=='DAEM' and PyrolProgram=='CPD') or (FG_FittingKineticParameter_Select=='DAEM' and PyrolProgram=='FGDVC'):
        PredictionDAEM=[2e10,20e3,5e3,0.5]
        outfile = open(PyrolProgram+'-Results_DAEM.txt', 'w')
        outfile.write("Species                  A1 [1/s]               E_a1 [K]          sigma [K]   Final Yield\n\n")
        self.KinModel=Models.DAEM(PredictionDAEM)
        self.KinModel.setNrOfActivationEnergies(GlobalOptParam.NrOFActivtionEnergies)
        if PyrolProgram=='PCCL':
                self.KinModel.setDt4Intergrate(self.FG_dt)
        #######
        ParamInit = GlobalOptParam.EvADAEMInit
        if len(ParamInit) == 4:
                ParamInit.pop(-1)
        ParamMin = GlobalOptParam.EvADAEMMin
        if len(ParamMin) == 3:
            ParamMin.append(0.0)
        ParamMax = GlobalOptParam.EvADAEMMax
        if len(ParamMax) == 3:
            ParamMax.append(0.0)
        ##The single species:
        if 'Total' not in self.SpeciesToConsider:
            self.SpeciesToConsider.append('Total')
        for Species in [Fit[0].SpeciesIndex('Total')]:
            print Fit[0].SpeciesName(Species)
            m_final_predictionAll=[]
            for i in range(len(Fit)):
                m_final_predictionAll.append(Fit[i].Yield(Species)[-1])
            # optimization procedure
            if self.NrOfRuns == 1:
                self.OptGradBased(Fit,ParamInit,Fit[0].Yield(Species)[-1],Species)
            else:
                if len(ParamInit) == 3:
                    ParamInit.append(0.0)
                ParamInit[3]= ((min(m_final_predictionAll)+max(m_final_predictionAll))/2.)
                ParamMin[3]= (min(m_final_predictionAll))
                ParamMax[3]= (max(m_final_predictionAll))
                self.OptGenAlgBased(Fit,ParamInit,ParamMin,ParamMax,Species)
            #
            self.Solution=self.KinModel.ParamVector()
            #
            self.KinModel.plot(Fit,Species)
            outfile.write(str(Fit[0].SpeciesName(Species))+'\t'+'%.6e  %11.4f  %11.4f  %6.4f  ' %(self.Solution[0],self.Solution[1],self.Solution[2],self.Solution[3])+'\n')
        outfile.close()
        if oSystem=='Linux' or oSystem == 'Darwin':
            shutil.move(PyrolProgram+'-Results_DAEM.txt','Result/'+PyrolProgram+'-Results_DAEM.txt')
        elif oSystem=='Windows':
            shutil.move(PyrolProgram+'-Results_DAEM.txt','Result\\'+PyrolProgram+'-Results_DAEM.txt')
        else:
            print "The name of the operating system couldn't be found."
#
#
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


    def MakeResults_CPD(self):
        """ Execute CPD for each given temperature profile and 
            return a list of CPD results objects
         """
        def InitAndLaunch(*pargs):
            """ initialises and execute cpd calculation """
            print 'Running CPD: ' + str(run)
            cpd = CPD_SetAndLaunch.CPD(*pargs)
            return cpd.Run() 

        operatingConditions = self.inputs['OperatingConditions']
        pressure = operatingConditions.pop('pressure') 
        return [InitAndLaunch(self.ultim_ana, 
                              self.daf,
                              operatingConditions["run"+str(run+1)],
                              pressure,
                              self.inputs['CPD']['deltaT'],
                              run)
                for run in range(operatingConditions['runs'])]

    def fitCPDResults(self, cpdResults):
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

    ####FG-DVC####
    def MakeResults_FG(self):
        """generates the result for FG-DVC"""
        #writes Time-Temperature file
        FG_TimeTemp1=np.zeros(np.shape(self.CPD_TimeTemp1),order='F')
        FG_TimeTemp2=np.zeros(np.shape(self.CPD_TimeTemp2),order='F')
        FG_TimeTemp3=np.zeros(np.shape(self.CPD_TimeTemp3),order='F')
        FG_TimeTemp4=np.zeros(np.shape(self.CPD_TimeTemp4),order='F')
        FG_TimeTemp5=np.zeros(np.shape(self.CPD_TimeTemp5),order='F')
        FG_TimeTemp1[:,0]=self.CPD_TimeTemp1[:,0]*1.e-3
        FG_TimeTemp2[:,0]=self.CPD_TimeTemp2[:,0]*1.e-3
        FG_TimeTemp3[:,0]=self.CPD_TimeTemp3[:,0]*1.e-3
        FG_TimeTemp4[:,0]=self.CPD_TimeTemp4[:,0]*1.e-3
        FG_TimeTemp5[:,0]=self.CPD_TimeTemp5[:,0]*1.e-3
        FG_TimeTemp1[:,1]=self.CPD_TimeTemp1[:,1]
        FG_TimeTemp2[:,1]=self.CPD_TimeTemp2[:,1]
        FG_TimeTemp3[:,1]=self.CPD_TimeTemp3[:,1]
        FG_TimeTemp4[:,1]=self.CPD_TimeTemp4[:,1]
        FG_TimeTemp5[:,1]=self.CPD_TimeTemp5[:,1]
        #initialize the launching object
        FGDVC=FGDVC_SetAndLaunch.SetterAndLauncher()
        #set and writes Coal Files:
        if self.FG_CoalSelection==0:
            #deletes old generated file
            os.system('cd '+self.FG_MainDir+FG_GenCoalDir+' & del '+FG_CoalName+'_com.dat, '+FG_CoalName+'_kin.dat, '+FG_CoalName+'_pol.dat')
            #generates coalsd.exe input file
            MakeCoalGenFile=InformationFiles.WriteFGDVCCoalFile(FG_CoalGenFileName)
            MakeCoalGenFile.setCoalComp(self.UAC,self.UAH,self.UAO,self.UAN,self.UAS,0)
            MakeCoalGenFile.write(self.FG_MainDir+FG_GenCoalDir+'\\',FG_CoalName,option=0)
            #makes new file
            try:
                os.system('cd '+self.FG_MainDir+FG_GenCoalDir+' & '+'coalsd.exe < '+FG_CoalGenFileName+' > coalsd_pkp.log')
            except OSError:
                print 'Problems with coalsd.exe'
            os.system('copy '+self.FG_MainDir+FG_GenCoalDir+'\coalsd_pkp.log . >> log.txt')
            #tests weather the coal file was genearated:
            if os.path.exists(self.FG_MainDir+'\\'+FG_GenCoalDir+'\\'+FG_CoalName+'_com.dat')==False:
                print 30*'*','\n','The coal is may outside the libraries coals. Select manually the closest library coal.',30*'*','\n'
                MakeCoalGenFile.write(self.FG_MainDir+FG_GenCoalDir+'\\',FG_CoalName,option=10)
                os.system('cd '+self.FG_MainDir+FG_GenCoalDir+' & '+'coalsd.exe < '+FG_CoalGenFileName+' > coalsd_pkp.log')
            #sets generated file for instruct.ini
            FGDVC.set1CoalLocation(self.FG_MainDir+FG_GenCoalDir+'\\'+FG_CoalName+'_com.dat')
            FGDVC.set2KinLocation(self.FG_MainDir+FG_GenCoalDir+'\\'+FG_CoalName+'_kin.dat')
            FGDVC.set3PolyLocation(self.FG_MainDir+FG_GenCoalDir+'\\'+FG_CoalName+'_pol.dat')
        elif self.FG_CoalSelection>0 and self.FG_CoalSelection<9:
            #sets library file for instruct.ini
            FGDVC.set1CoalLocation(self.FG_MainDir+FG_LibCoalDir+'\\coal.ar'+str(self.FG_CoalSelection))
            FGDVC.set2KinLocation(self.FG_MainDir+FG_LibCoalDir+'\\kin.ar'+str(self.FG_CoalSelection))
            FGDVC.set3PolyLocation(self.FG_MainDir+FG_LibCoalDir+'\\polymr.ar'+str(self.FG_CoalSelection))
        else:
            print "select Choose Coal: 0 interpolate between library coals and generate own coal. Set 1 to 8 for a library coal.' in FGDVC.inp equal a value between 0 and 8"
        #sets FG-DVC instruct.ini parameter
        FGDVC.set5Pressure(self.FG_pressure)
        if self.FG_TarCacking==0.0:            #case: no tar cracking
            FGDVC.set6Theorie(13,0.0)
        elif self.FG_TarCacking<0.0:           #case: full tar cracking
            FGDVC.set6Theorie(15,0.0)
        else:                             #case: partial tar cracking
            FGDVC.set6Theorie(13,float(self.FG_TarCacking))
        #
        FGFile=[]
        FGFit=[]
        OpCondInp=InformationFiles.OperCondInput('OperCond.inp')
        for runNr in range(self.NrOfRuns):
            if runNr==0:
                OpCondInp.writeFGDVCtTHist(FG_TimeTemp1,self.FG_dt,self.FG_T_t_History)
            elif runNr==1:
                OpCondInp.writeFGDVCtTHist(FG_TimeTemp2,self.FG_dt,self.FG_T_t_History)
            elif runNr==2:
                OpCondInp.writeFGDVCtTHist(FG_TimeTemp3,self.FG_dt,self.FG_T_t_History)
            elif runNr==3:
                OpCondInp.writeFGDVCtTHist(FG_TimeTemp4,self.FG_dt,self.FG_T_t_History)
            elif runNr==4:
                OpCondInp.writeFGDVCtTHist(FG_TimeTemp5,self.FG_dt,self.FG_T_t_History)
            FGDVC.set7File(self.FG_T_t_History)
            FGDVC.set9AshMoisture(0.0,0.0)
            FGDVC.setTRamp_or_TFile('File') #case: models temperature history with the file
            #writes the instruct.ini and launches FG-DVC (no graphical user interface, only main file fgdvcd.exe)
            FGDVC.writeInstructFile(self.FG_MainDir+'\\'+FG_ExeCoalDir+'\\')
            FGDVC.Run('cd '+self.FG_MainDir+FG_ExeCoalDir+' & '+'fgdvcd.exe')
            #
            ###calibrate kinetic parameter:
            #read result:
            CurrentFGFile=FGDVC_Result.FGDVC_Result(self.FG_DirOut)
            # creates object, required for fitting procedures
            CurrentFGFit=FitInfo.Fit_one_run(CurrentFGFile)
            FGFile.append(CurrentFGFile)
            FGFit.append(CurrentFGFit)
            #copies file, keeping the name:
            if oSystem=='Linux' or oSystem == 'Darwin':
                shutil.copyfile(self.FG_DirOut+'gasyield.txt', 'Result/gasyield_'+str(runNr)+'.txt')
                shutil.copyfile(self.FG_DirOut+'gasrate.txt', 'Result/gasrate_'+str(runNr)+'.txt')
            elif oSystem=='Windows':
                shutil.copyfile(self.FG_DirOut+'gasyield.txt', 'Result\\gasyield_'+str(runNr)+'.txt')
                shutil.copyfile(self.FG_DirOut+'gasrate.txt', 'Result\\gasrate_'+str(runNr)+'.txt')
        #####
        M=Models.Model()
        for Species in FGFit[0].SpeciesNames():
            M.mkSimpleResultFiles(FGFit,Species)
            if (Species not in self.SpeciesToConsider) and (Species!='Temp') and (Species!='Time'):
                self.SpeciesToConsider.append(Species)
        if self.FG_FittingKineticParameter_Select=='constantRate':
            self.MakeResults_CR('FGDVC',FGFile,FGFit)
            currentDict={'FGDVC':'constantRate'}
        elif self.FG_FittingKineticParameter_Select=='Arrhenius':
            self.MakeResults_Arrh('FGDVC',FGFile,FGFit)
            currentDict={'FGDVC':'Arrhenius'}
        elif self.FG_FittingKineticParameter_Select=='ArrheniusNoB':
            self.MakeResults_ArrhNoB('FGDVC',FGFile,FGFit)
            currentDict={'FGDVC':'ArrheniusNoB'}
        elif self.FG_FittingKineticParameter_Select=='Kobayashi':
            self.MakeResults_Kob('FGDVC',FGFile,FGFit)
            currentDict={'FGDVC':'Kobayashi'}
        elif self.FG_FittingKineticParameter_Select=='DAEM':
            self.MakeResults_DEAM('FGDVC',FGFile,FGFit)
            currentDict={'FGDVC':'DAEM'}
        elif self.FG_FittingKineticParameter_Select==None:
            currentDict={'FGDVC':'None'}
        else:
            print 'uspecified FG_FittingKineticParameter_Select'
            currentDict={}
        #
        self.ProgramModelDict.update(currentDict)
        #
        self.SpeciesEnergy('FGDVC',FGFile,self.FG_FittingKineticParameter_Select)
            #
    
    ####Pc Coal Lab####
    def MakeResults_PCCL(self):
        """generates the result for PC Coal Lab"""
        #writes Time-Temperature file
        PCCL_TimeTemp1=np.zeros(np.shape(self.CPD_TimeTemp1),order='F')
        PCCL_TimeTemp2=np.zeros(np.shape(self.CPD_TimeTemp2),order='F')
        PCCL_TimeTemp3=np.zeros(np.shape(self.CPD_TimeTemp3),order='F')
        PCCL_TimeTemp4=np.zeros(np.shape(self.CPD_TimeTemp4),order='F')
        PCCL_TimeTemp5=np.zeros(np.shape(self.CPD_TimeTemp5),order='F')
        PCCL_TimeTemp1[:,0]=self.CPD_TimeTemp1[:,0]*1.e-3
        PCCL_TimeTemp2[:,0]=self.CPD_TimeTemp2[:,0]*1.e-3
        PCCL_TimeTemp3[:,0]=self.CPD_TimeTemp3[:,0]*1.e-3
        PCCL_TimeTemp4[:,0]=self.CPD_TimeTemp4[:,0]*1.e-3
        PCCL_TimeTemp5[:,0]=self.CPD_TimeTemp5[:,0]*1.e-3
        PCCL_TimeTemp1[:,1]=self.CPD_TimeTemp1[:,1]
        PCCL_TimeTemp2[:,1]=self.CPD_TimeTemp2[:,1]
        PCCL_TimeTemp3[:,1]=self.CPD_TimeTemp3[:,1]
        PCCL_TimeTemp4[:,1]=self.CPD_TimeTemp4[:,1]
        PCCL_TimeTemp5[:,1]=self.CPD_TimeTemp5[:,1]
        #initialize the launching object
        PCCL=PCCL_SetAndLaunch.SetterAndLauncher()
        #set and writes Coal Files:
        PCCL.SetUACoalParameter(self.UAC,self.UAH,self.UAN,self.UAO,self.UAS)
        PCCL.SetPACoalParameter(self.PAVM_asrec,self.PAFC_asrec,self.PAmoist,self.PAash)
        if type(self.PCCL_CoalCalFactor)==float:
            PCCL.SetCoalCalibrationFactor(self.PCCL_CoalCalFactor)
        PCCL.SetPressure(self.FG_pressure)
        PCCL.SetParticleSize(self.PCCL_ParticleSize)
        #
        PCCL.writeCoalFiles(self.PCCL_Path)
        #
        PCCL.THist(PCCL_TimeTemp1[0,1],PCCL_TimeTemp1[1,0],PCCL_TimeTemp1[-1,1],PCCL_TimeTemp1[-1,0])
        PCCL.writeInstructFiles(self.PCCL_Path,mkNewFile=True)
        PCCL.THist(PCCL_TimeTemp2[0,1],PCCL_TimeTemp2[1,0],PCCL_TimeTemp2[-1,1],PCCL_TimeTemp2[-1,0])
        PCCL.writeInstructFiles(self.PCCL_Path,mkNewFile=False)
        PCCL.THist(PCCL_TimeTemp3[0,1],PCCL_TimeTemp3[1,0],PCCL_TimeTemp3[-1,1],PCCL_TimeTemp3[-1,0])
        PCCL.writeInstructFiles(self.PCCL_Path,mkNewFile=False)
        PCCL.THist(PCCL_TimeTemp4[0,1],PCCL_TimeTemp4[1,0],PCCL_TimeTemp4[-1,1],PCCL_TimeTemp4[-1,0])
        PCCL.writeInstructFiles(self.PCCL_Path,mkNewFile=False)
        PCCL.THist(PCCL_TimeTemp5[0,1],PCCL_TimeTemp5[1,0],PCCL_TimeTemp5[-1,1],PCCL_TimeTemp5[-1,0])
        PCCL.writeInstructFiles(self.PCCL_Path,mkNewFile=False)
        PCCL.writeInstructFilesFinish()
        #
        PCCL.Run(self.PCCL_Path,self.PCCL_Exe)
        #
        PCCLFile=[]
        PCCLFit=[]
        for runNr in range(1,self.NrOfRuns+1,1):
            #read result:
            CurrentPCCLFile=PCCL_Result.PCCL_Result(self.PCCL_Path,runNr)
            # creates object, required for fitting procedures
            CurrentPCCLFit=FitInfo.Fit_one_run(CurrentPCCLFile)
            PCCLFile.append(CurrentPCCLFile)
            PCCLFit.append(CurrentPCCLFit)
            #copies file:
            if oSystem=='Windows':
                shutil.copyfile(self.PCCL_Path+'FDC1WT'+str(runNr)+'.RPT', 'Result/PCCL_gasyield_wt_'+str(runNr)+'.txt')
                shutil.copyfile(self.PCCL_Path+'FDC1NG'+str(runNr)+'.RPT', 'Result/PCCL_gasyield_ng_'+str(runNr)+'.txt')
                shutil.copyfile(self.PCCL_Path+'FDC1HC'+str(runNr)+'.RPT', 'Result/PCCL_gasyield_hc_'+str(runNr)+'.txt')
        #####
        M=Models.Model()
        for Species in PCCLFit[0].SpeciesNames():
            M.mkSimpleResultFiles(PCCLFit,Species)
            if (Species not in self.SpeciesToConsider) and (Species!='Temp') and (Species!='Time'):
                self.SpeciesToConsider.append(Species)
        if self.PCCL_FittingKineticParameter_Select=='constantRate':
            self.MakeResults_CR('PCCL',PCCLFile,PCCLFit)
            currentDict={'PCCL':'constantRate'}
        elif self.PCCL_FittingKineticParameter_Select=='Arrhenius':
            self.MakeResults_Arrh('PCCL',PCCLFile,PCCLFit)
            currentDict={'PCCL':'Arrhenius'}
        elif self.PCCL_FittingKineticParameter_Select=='ArrheniusNoB':
            self.MakeResults_ArrhNoB('PCCL',PCCLFile,PCCLFit)
            currentDict={'PCCL':'ArrheniusNoB'}
        elif self.PCCL_FittingKineticParameter_Select=='Kobayashi':
            self.MakeResults_Kob('PCCL',PCCLFile,PCCLFit)
            currentDict={'PCCL':'Kobayashi'}
        elif self.PCCL_FittingKineticParameter_Select=='DAEM':
            self.MakeResults_DEAM('PCCL',PCCLFile,PCCLFit)
            currentDict={'PCCL':'DAEM'}
        elif self.PCCL_FittingKineticParameter_Select==None:
            currentDict={'PCCL':'None'}
        else:
            print 'uspecified PCCL_FittingKineticParameter_Select'
            currentDict={}
        #
        self.ProgramModelDict.update(currentDict)
        #
        self.SpeciesEnergy('PCCL',PCCLFile)
            #


    def RunPMSKD(self):
        '''
        run PMSKD
        '''
        # create object

        try:
            coal = coalPolimi.coalPolimi(name = 'COAL', c=self.UAC,h=self.UAH,o=self.UAO,n=self.UAN,s=self.UAS,file=self.PMSKD_mechfile)
        except coalPolimi.compositionError:
            print 'Composition outside of triangle of definition'
            sys.exit()
        # organize TimeTemp
            PMSKDFile=[]
        PMSKDFit=[]
        for runNr in range(self.NrOfRuns):
            print 'Running PMSKD n. '+str(runNr)
            #print self.timeHR[runNr]
            #print self.temperatureHR[runNr]
            #set heating rate
            coal.setHeatingRate(self.timeHR[runNr],self.temperatureHR[runNr])
            #coal.setTimeStep(self.PMSKD_npoint)
            coal.solvePyrolysis()
            #plt.figure(runNr)
            #plt.plot(coal.getTemperature(),coal.getVolatile())
            #read result:
            #CurrentPMSKDFile=Coal
            # creates object, required for fitting procedures
            CurrentPMSKDFit=FitInfo.Fit_one_run(coal)
            #PMSKDFile.append(CurrentFGFile)
            PMSKDFit.append(CurrentPMSKDFit)
            #print coal.Yields_all()

            coal.reset()
            #print coal.timeHR
            #print coal.temperatureHR

        if self.PMSKD_FittingKineticParameter_Select=='constantRate':
            self.MakeResults_CR('PMSKD','',PMSKDFit)
            currentDict={'PMSKD':'constantRate'}
        elif self.PMSKD_FittingKineticParameter_Select=='Arrhenius':
            self.MakeResults_Arrh('PMSKD','',PMSKDFit)
            currentDict={'PMSKD':'Arrhenius'}
        elif self.PMSKD_FittingKineticParameter_Select=='ArrheniusNoB':
            self.MakeResults_ArrhNoB('PMSKD','',PMSKDFit)
            currentDict={'PMSKD':'ArrheniusNoB'}
        elif self.PMSKD_FittingKineticParameter_Select=='Kobayashi':
            self.MakeResults_Kob('PMSKD','',PMSKDFit)
            currentDict={'PMSKD':'Kobayashi'}
        elif self.PMSKD_FittingKineticParameter_Select=='DAEM':
            self.MakeResults_DEAM('PMSKD','',PMSKDFit)
            currentDict={'PMSKD':'DAEM'}
        elif self.PMSKD_FittingKineticParameter_Select==None:
            currentDict={'PMSKD':'None'}
            for Species in PMSKDFit[0].SpeciesNames():
                M=Models.Model()
                M.mkSimpleResultFiles(PMSKDFit,Species)
                if ((Species not in self.SpeciesToConsider) 
                    and (Species!='Temp') 
                    and (Species!='Time')):
                    self.SpeciesToConsider.append(Species)
        else:
            print 'undefined PMSKD_FittingKineticParameter_Select'
            currentDict={}
            #
        self.ProgramModelDict.update(currentDict)
        #
        #self.SpeciesEnergy('PMSKD',FGFile)

#Main Part starting

def main():
    Case = MainProcess(inputs_folder=workingDir+"/inputs/")
    Case.executeSolver()
    print 'calculated Species: ',Case.SpeciesToConsider

if __name__ == "__main__":
    main()
