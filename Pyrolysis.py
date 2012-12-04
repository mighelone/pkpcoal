import sys
sys.path.append('src')
#
import CPD_SetAndLaunch         #writes CPD-instruct File, launches CPD
import FGDVC_SetAndLaunch       #writes FG-DVC-instruct File, launches FG-DVC and fittes using eq. (68 ) (BachelorThesis)
import FGDVC_Result             #contains the information of the FG-DVC output file
import CPD_Result               #contains the information of the FG-DVC output file
import Fitter                   #The optimizer class
import Models                   #the models like Arrhenius or Kobayashi
import FitInfo                  #supports the Fitting with the yield information
import Compos_and_Energy        #Species balance and energy balance for CPD and FG-DVC
import InformationFiles         #reads the user input files, writes FG-DVC coalsd.exe coal generation file
import GlobalOptParam           #contains the Information of the Number Of Runs for the Global Optimum search
import Evolve                   #contains the generic algortihm optimizer
import os
import numpy as np
import platform
import shutil
#
#
#Use Global Optimizer? select 'Evolve' for a generic algorithm or 'ManyPoints' to use many starting points combined with a local optimization
UseGlobalOpt= 'Evolve'
#UseGlobalOpt= GlobalOptParam.GlobalOptimizeMethod
#Which operating Sytem?
oSystem=platform.system()
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
#
#
class MainProcess(object):
    """Controls the whole process of generating input files, fitting etc."""
    def __init__(self):
        self.SpeciesToConsider=[] #for GUI
        self.ProgramModelDict={} #for GUI
    #
    def ReadInputFiles(self):
        """get parameters from input files"""
        #
        #Coal File:
        #
        print 'Reading Coal.inp ...'
        CoalInput=InformationFiles.ReadFile(workingDir+'Coal.inp')
        self.PAFC_asrec=CoalInput.getValue(InformationFiles.M_PA[0])
        self.PAVM_asrec=CoalInput.getValue(InformationFiles.M_PA[1])
        self.PAmoist = CoalInput.getValue(InformationFiles.M_PA[2])
        self.PAash = CoalInput.getValue(InformationFiles.M_PA[3])
        # scale proximate analysis
        sumPA = (self.PAFC_asrec+self.PAVM_asrec + self.PAmoist + self.PAash)/100.
        self.PAFC_asrec/=sumPA
        self.PAVM_asrec/=sumPA
        self.PAmoist/=sumPA
        self.PAash/=sumPA
        #
        #gets daf values, as CPD needs daf as input:
        self.PAFC_daf, self.PAVM_daf = self.DAF(self.PAFC_asrec,self.PAVM_asrec)
        self.UAC=CoalInput.getValue(InformationFiles.M_UA[0])
        self.UAH=CoalInput.getValue(InformationFiles.M_UA[1])
        self.UAN=CoalInput.getValue(InformationFiles.M_UA[2])
        self.UAO=CoalInput.getValue(InformationFiles.M_UA[3])
        self.UAS=CoalInput.getValue(InformationFiles.M_UA[4])
        # scale ultimate analysis
        sumUA = self.UAC+self.UAH+self.UAN+self.UAO+self.UAS
        self.UAC=self.UAC/sumUA*100
        self.UAH=self.UAH/sumUA*100
        self.UAO=self.UAO/sumUA*100
        self.UAN=self.UAN/sumUA*100
        self.UAS=self.UAS/sumUA*100
        self.HHV=CoalInput.getValue(InformationFiles.M_HHV)
        self.MTar=CoalInput.getValue(InformationFiles.M_MTar)
        self.WeightY=CoalInput.getValue(InformationFiles.M_Weight[0])
        self.WeightR=CoalInput.getValue(InformationFiles.M_Weight[1])
        #
        #CPD Properties:
        #
        print 'Reading CPD.inp ...'
        CPDInput=InformationFiles.ReadFile(workingDir+'CPD.inp')
        self.CPDselect=CPDInput.UsePyrolProgr(InformationFiles.MC_sel)
        self.CPD_FittingKineticParameter_Select=CPDInput.Fitting(InformationFiles.M_selFit)
        self.CPD_ArrhSpec=CPDInput.getText(InformationFiles.M_selArrhSpec)
        self.CPDdt=[0,1,2] #0:initila dt, 1: print increment, 2: dt max
        self.CPDdt[0]=(CPDInput.getValue(InformationFiles.MC_dt[0]))
        self.CPDdt[1]=(CPDInput.getValue(InformationFiles.MC_dt[1]))
        #
        #
        #FG-DVC Properties:
        #
        print 'Reading FGDVC.inp ...'
        FGDVCInput=InformationFiles.ReadFile(workingDir+'FGDVC.inp')
        self.FG_select=FGDVCInput.UsePyrolProgr(InformationFiles.MF_sel)
        self.FG_FittingKineticParameter_Select=FGDVCInput.Fitting(InformationFiles.M_selFit)
        self.FG_ArrhSpec=FGDVCInput.getText(InformationFiles.M_selArrhSpec)
        self.FG_CoalSelection=int(FGDVCInput.getValue(InformationFiles.MF_CoalSel))
        self.FG_MainDir=FGDVCInput.getText(InformationFiles.MF_dir[0])
        self.FG_DirOut=FGDVCInput.getText(InformationFiles.MF_dir[1])
        self.FG_TarCacking=FGDVCInput.getValue(InformationFiles.MF_TarCr)
        #
        #
        #Operating Condition File:
        #
        print 'Reading OperCond.inp ...'
        OpCondInp=InformationFiles.OperCondInput('OperCond.inp')
        self.CPD_pressure=OpCondInp.getValue(InformationFiles.M_Pressure)
        self.FG_pressure=OpCondInp.getValue(InformationFiles.M_Pressure)
        #Number of FG-DVC/CPD/PCCL runs:
        self.NrOfRuns=int(OpCondInp.getValue(InformationFiles.M_NrRuns))
        self.CPD_TimeTemp1=OpCondInp.getTimePoints(InformationFiles.M_TimePoints1[0],InformationFiles.M_TimePoints1[1])
        self.CPD_TimeTemp2=OpCondInp.getTimePoints(InformationFiles.M_TimePoints2[0],InformationFiles.M_TimePoints2[1])
        self.CPD_TimeTemp3=OpCondInp.getTimePoints(InformationFiles.M_TimePoints3[0],InformationFiles.M_TimePoints3[1])
        self.CPD_TimeTemp4=OpCondInp.getTimePoints(InformationFiles.M_TimePoints4[0],InformationFiles.M_TimePoints4[1])
        self.CPD_TimeTemp5=OpCondInp.getTimePoints(InformationFiles.M_TimePoints5[0],InformationFiles.M_TimePoints5[1])
        self.CPDdt[2]=OpCondInp.getValue(InformationFiles.M_dt)
        self.FG_dt=OpCondInp.getValue(InformationFiles.M_dt)
        self.FG_T_t_History=self.FG_MainDir+'tTHistory.txt'
        #makes for CPD time in milliseconds:
        self.CPD_TimeTemp1[:,0]=self.CPD_TimeTemp1[:,0]*1.e3
        self.CPD_TimeTemp2[:,0]=self.CPD_TimeTemp2[:,0]*1.e3
        self.CPD_TimeTemp3[:,0]=self.CPD_TimeTemp3[:,0]*1.e3
        self.CPD_TimeTemp4[:,0]=self.CPD_TimeTemp4[:,0]*1.e3
        self.CPD_TimeTemp5[:,0]=self.CPD_TimeTemp5[:,0]*1.e3
        self.CPD_t_max1=self.CPD_TimeTemp1[-1,0]*1.e-3 #tmax in s, not ms
        self.CPD_t_max2=self.CPD_TimeTemp2[-1,0]*1.e-3 #tmax in s, not ms
        self.CPD_t_max3=self.CPD_TimeTemp3[-1,0]*1.e-3 #tmax in s, not ms
        self.CPD_t_max4=self.CPD_TimeTemp4[-1,0]*1.e-3 #tmax in s, not ms
        self.CPD_t_max5=self.CPD_TimeTemp5[-1,0]*1.e-3 #tmax in s, not ms
        #
        #
        #
    def DAF(self,PAFC_asRecieved,PAVM_asRecieved):
        """calculates PAFC, PAVM  from the as recieved state to the daf state of coal"""
        fractionFC=PAFC_asRecieved/(PAFC_asRecieved+PAVM_asRecieved)
        fractionVM=PAVM_asRecieved/(PAFC_asRecieved+PAVM_asRecieved)
        return 100.*fractionFC, 100.*fractionVM
        #
    def CheckFGdt(self):
        """Aborts, if FG-DVC is selected and the timestep is lower than 1.e-3 (which is FG-DVC not able to read):"""
        if ((self.FG_select==True) and (self.FG_dt<1e-4)):
            print "Please select for FG-DVC a time step greather equal 1e-4 in 'OperCond.inp'. FG-DVC would not be able to read the time history file for a dt<1e-4."
            sys.exit()
            #
            #
    def MakeResults_CR(self,PyrolProgram,File,Fit):
        """Generates the results for constant Rate."""
    #Array for the comparison of the sum of the individual species(t) with the (1-solid(t)); this part ([0]) has to be modified, when the assumption of same length of the array (e.g. when different final times are possible) is refused
#    SumSingleYieldsCalc=np.zeros([len(Fit[0].Yield('Time')),NrOfRuns]) #Array initialized for the sum of the single yields (calculated)
#    SumSingleYields=np.zeros([len(Fit[0].Yield('Time')),NrOfRuns])  #Array initialized for the sum of the single yields (CPD output)
#    SolidYieldsCalc=np.zeros([len(Fit[0].Yield('Time')),NrOfRuns])     #Array initialized for the yields of the Solids (calculated)
#    SolidYields=np.zeros([len(Fit[0].Yield('Time')),NrOfRuns])     #Array initialized for the yields of the Solids (CPD output)
    ##CONSTANT RATE
#    if (CPD_FittingKineticParameter_Select=='constantRate' and PyrolProgram=='CPD') or (FG_FittingKineticParameter_Select=='constantRate' and PyrolProgram=='FGDVC'): #CR means ConstantRate
        LS=Fitter.LeastSquarsEstimator() #LS means LeastSquares
        LS.setOptimizer('fmin')
        LS.setTolerance(1.e-18)
        LS.setWeights(self.WeightY,self.WeightR)
        outfile = open(PyrolProgram+'-Results_const_rate.txt', 'w')
        outfile.write("Species\tk [1/s]\t\tt_start [s]\t\tFinalYield\n\n")
        for Spec in range(2,len(Fit[0].SpeciesNames()),1):
            if Fit[0].SpeciesName(Spec) not in self.SpeciesToConsider:
                self.SpeciesToConsider.append(Fit[0].SpeciesName(Spec))
            #
            m_final_prediction=Fit[0].Yield(Spec)[-1]
            PredictionVector=[50,0.01,m_final_prediction] #first argument is k, second is t_start
            CR=Models.ConstantRateModel(PredictionVector)
            #
            CR.setParamVector(PredictionVector)
            CR.setParamVector(LS.estimate_T(Fit,CR,PredictionVector,Spec))
            CR.plot(Fit,Spec)
            Solution=CR.ParamVector()
            if np.sum(Solution)!=np.sum(PredictionVector):
                outfile.write(str(Fit[0].SpeciesName(Spec))+'\t'+str(Solution[0])+'\t'+str(Solution[1])+'\t'+str(Solution[2])+'\n')
#            #for the comparison of the species sum with (1-Solid)
#            for runNr in range(NrOfRuns):
#                if Fit[runNr].SpeciesName(Spec)=='Solid':
#                    SolidYieldsCalc[:,runNr]+=CRCPD.calcMass(Fit[runNr],Fit[runNr].Time(),Fit[runNr].Interpolate('Temp'),Spec)
#                    SolidYields[:,runNr]+=Fit[runNr].Yield(Spec)
#                elif Fit[runNr].SpeciesName(Spec)!='Solid' and Fit[runNr].SpeciesName(Spec)!='Temp' and Fit[runNr].SpeciesName(Spec)!='Time' and Fit[runNr].SpeciesName(Spec)!='Gas' and Fit[runNr].SpeciesName(Spec)!='Total':
#                    SumSingleYieldsCalc[:,runNr]+=CRCPD.calcMass(Fit[runNr],Fit[runNr].Time(),Fit[runNr].Interpolate('Temp'),Spec)
#                    SumSingleYields[:,runNr]+=Fit[runNr].Yield(Spec)
        outfile.close()
        if oSystem=='Linux':
            shutil.move(PyrolProgram+'-Results_const_rate.txt','Result/'+PyrolProgram+'-Results_constantRate.txt')
        elif oSystem=='Windows':
            shutil.move(PyrolProgram+'-Results_const_rate.txt','Result\\'+PyrolProgram+'-Results_constantRate.txt')
        else:
            print "The name of the operating system couldn't be found."
#        Fit[0].plt_InputVectors(Fit[0].Time(),1.-SolidYieldsCalc,1.-SolidYields,SumSingleYieldsCalc,SumSingleYields,'1-Solid; fitted','1-Solid; CPD output','Sum Yields; fitted','Sum Yields; CPD output')
#
#
    def MakeResults_Arrh(self,PyrolProgram,File,Fit):
        """Generates the results for Arrhenius Rate."""
    ##ARRHENIUS RATE
#    elif (CPD_FittingKineticParameter_Select=='Arrhenius' and PyrolProgram=='CPD') or (FG_FittingKineticParameter_Select=='Arrhenius' and PyrolProgram=='FGDVC'): #Arr means Arrhenius
        LS=Fitter.LeastSquarsEstimator()
        LS.setOptimizer('fmin')#('leastsq')   # 'leastsq' often faster, but if this does not work: 'fmin' is more reliable
        LS.setTolerance(1.e-7)
        LS.setWeights(self.WeightY,self.WeightR)
        outfile = open(PyrolProgram+'-Results_ArrheniusRate.txt', 'w')
        outfile.write("Species\tA [1/s]\t\tb\t\tE_a [K]\t\tFinalYield\n\n")
        #select one of the follwoing notations: 
        #Arr=Models.ArrheniusModel(PredictionV0)
        #Arr=Models.ArrheniusModelAlternativeNotation1(PredictionV1)
        #######
        #makes Species list which contains alls species to fit:
        SpeciesList=[]
        if self.CPD_ArrhSpec=='Total' and PyrolProgram=='CPD':
            SpeciesList.append(Fit[0].SpeciesIndex('Total'))
            if 'Total' not in self.SpeciesToConsider:
                self.SpeciesToConsider.append('Total')
        elif self.CPD_ArrhSpec=='MainSpecies' and PyrolProgram=='CPD':
            SpeciesList.append(Fit[0].SpeciesIndex('Total'))
            SpeciesList.append(Fit[0].SpeciesIndex('Tar'))
            SpeciesList.append(Fit[0].SpeciesIndex('Gas'))
            if 'Total' not in self.SpeciesToConsider:
                self.SpeciesToConsider.append('Total')
            if 'Tar' not in self.SpeciesToConsider:
                self.SpeciesToConsider.append('Tar')
            if 'Gas' not in self.SpeciesToConsider:
                self.SpeciesToConsider.append('Gas')
        elif self.CPD_ArrhSpec=='allSpecies' and PyrolProgram=='CPD':
            for i in range(2,len(Fit[0].SpeciesNames()),1):
                if Fit[0].SpeciesName(i) not in self.SpeciesToConsider:
                    self.SpeciesToConsider.append(Fit[0].SpeciesName(i))
                SpeciesList.append(i)
        elif self.FG_ArrhSpec=='Total' and PyrolProgram=='FGDVC':
            SpeciesList.append(Fit[0].SpeciesIndex('Total'))
            if 'Total' not in self.SpeciesToConsider:
                self.SpeciesToConsider.append('Total')
        elif self.FG_ArrhSpec=='MainSpecies' and PyrolProgram=='FGDVC':
            SpeciesList.append(Fit[0].SpeciesIndex('Total'))
            SpeciesList.append(Fit[0].SpeciesIndex('Tar'))
            SpeciesList.append(Fit[0].SpeciesIndex('Gas'))
            if 'Total' not in self.SpeciesToConsider:
                self.SpeciesToConsider.append('Total')
            if 'Tar' not in self.SpeciesToConsider:
                self.SpeciesToConsider.append('Tar')
            if 'Gas' not in self.SpeciesToConsider:
                self.SpeciesToConsider.append('Gas')
        elif self.FG_ArrhSpec=='allSpecies' and PyrolProgram=='FGDVC':
            for i in range(2,len(Fit[0].SpeciesNames()),1):
                if Fit[0].SpeciesName(i) not in self.SpeciesToConsider:
                    self.SpeciesToConsider.append(Fit[0].SpeciesName(i))
                SpeciesList.append(i)
        ##The single species:
        for Species in SpeciesList:
            #
            m_final_prediction=Fit[0].Yield(Species)[-1]
            PredictionV0=[0.86e15,0.01,27700,m_final_prediction]  #for Standard Arrhenius
#            PredictionV1=[10.,-20.,m_final_prediction]         #for Arrhenius notation #1
#            PredictionV2=[10.,-18.,m_final_prediction]           #for Arrhenius notation #2
            Arr=Models.ArrheniusModel(PredictionV0)
            #
            print Fit[0].SpeciesName(Species)
            if UseGlobalOpt=='ManyPoints':
                #GlobalOptimize:
                GlobalMin=Fitter.GlobalOptimizer(LS,Arr,Fit)
                m_final_predictionAll=[]
                for i in range(len(Fit)):
                    m_final_predictionAll.append(Fit[i].Yield(Species)[-1])
                LS.setTolerance(1e-2)
                ParamGlobalMin=GlobalMin.GenerateOptima(Species,GlobalOptParam.ArrhIndexToOptimize,[[min(m_final_predictionAll),max(m_final_predictionAll)]],GlobalOptParam.ArrhNrOfRuns)
                #
                LS.setTolerance(1.e-7)
                print 'Final optimization Run:'
                Arr.setParamVector(LS.estimate_T(Fit,Arr,ParamGlobalMin,Species))
            if UseGlobalOpt=='Evolve':
                m_final_predictionAll=[]
                for i in range(len(Fit)):
                    m_final_predictionAll.append(Fit[i].Yield(Species)[-1])
                GenAlg=Evolve.GenericOpt(Arr,Fit,Species)
                GenAlg.setWeights(GlobalOptParam.EvAWeightY,GlobalOptParam.EvAWeightY)
                GAArrhMinA = GlobalOptParam.EvAArrhMin[0]
                GAArrhMinB = GlobalOptParam.EvAArrhMin[1]
                GAArrhMinE = GlobalOptParam.EvAArrhMin[2]
                GAArrhMaxA = GlobalOptParam.EvAArrhMax[0]
                GAArrhMaxB = GlobalOptParam.EvAArrhMax[1]
                GAArrhMaxE = GlobalOptParam.EvAArrhMax[2]
                GAArrhInit=GlobalOptParam.EvAArrhInit
                if len(GAArrhInit)==3:
                    GAArrhInit.append((max(m_final_predictionAll)+min(m_final_predictionAll))/2.)
                else:
                    GAArrhInit[3]=(max(m_final_predictionAll)+min(m_final_predictionAll))/2.
                GenAlg.setParamRanges(GAArrhInit,[GAArrhMinA,GAArrhMinB,GAArrhMinE,min(m_final_predictionAll)],[GAArrhMaxA,GAArrhMaxB,GAArrhMaxE,max(m_final_predictionAll)])
                GenAlg.setNrPopulation(GlobalOptParam.NrOfPopulation)
                GenAlg.setNrGenerations(GlobalOptParam.NrOfGeneration)
                Arr.setParamVector(GenAlg.mkResults())
                #
                #use afterwards local optimization
                #Arr.setParamVector(LS.estimate_T(Fit,Arr,Arr.ParamVector(),Species))
            if UseGlobalOpt==False:
                Arr.setParamVector(LS.estimate_T(Fit,Arr,Arr.ParamVector(),Species))
            Arr.plot(Fit,Species)
            Solution=Arr.ParamVector()
            if np.sum(Arr.ParamVector())!=np.sum(PredictionV0): #To avoid, a species with no yield is added to the parameter file
                outfile.write(str(Fit[0].SpeciesName(Species))+'\t'+str(Solution[0])+'\t'+str(Solution[1])+'\t'+str(Solution[2])+'\t\t'+str(Solution[3])+'\n')
            #for the comparison of the species sum with (1-Solid)
#            for runNr in range(NrOfRuns):
#                if Fit[runNr].SpeciesName(Species)=='Solid':
#                    SolidYieldsCalc[:,runNr]+=ArrPCPD.calcMass(Fit[runNr],Fit[runNr].Time()[:len(Fit[0].Yield('Time'))],Fit[runNr].Interpolate('Temp'),Species)
#                    SolidYields[:,runNr]+=Fit[runNr].Yield(Species)[:len(Fit[0].Yield('Time'))],Fit[runNr]
#                elif Fit[runNr].SpeciesName(Species)!='Solid' and Fit[runNr].SpeciesName(Species)!='Temp' and Fit[runNr].SpeciesName(Species)!='Time' and Fit[runNr].SpeciesName(Species)!='Gas' and Fit[runNr].SpeciesName(Species)!='Total':
#                    SumSingleYieldsCalc[:,runNr]+=ArrPCPD.calcMass(Fit[runNr],Fit[runNr].Time()[:len(Fit[0].Yield('Time'))],Fit[runNr].Interpolate('Temp'),Species)
#                    SumSingleYields[:,runNr]+=Fit[runNr].Yield(Species)[:len(Fit[0].Yield('Time'))],Fit[runNr]
        outfile.close()
        if oSystem=='Linux':
            shutil.move(PyrolProgram+'-Results_ArrheniusRate.txt','Result/'+PyrolProgram+'-Results_Arrhenius.txt')
        elif oSystem=='Windows':
            shutil.move(PyrolProgram+'-Results_ArrheniusRate.txt','Result\\'+PyrolProgram+'-Results_Arrhenius.txt')
        else:
            print "The name of the operating system couldn't be found."
#        Fit[0].plt_InputVectors(Fit[0].Time(),1.-SolidYieldsCalc,1.-SolidYields,SumSingleYieldsCalc,SumSingleYields,'1-Solid; fitted','1-Solid; CPD output','Sum Yields; fitted','Sum Yields; CPD output')
#
#
    def MakeResults_ArrhNoB(self,PyrolProgram,File,Fit):
        """Generates the results for Arrhenius Rate with no correction term T**b."""
#    elif (CPD_FittingKineticParameter_Select=='ArrheniusNoB' and PyrolProgram=='CPD') or (FG_FittingKineticParameter_Select=='ArrheniusNoB' and PyrolProgram=='FGDVC'): #Arr means Arrhenius
        LS=Fitter.LeastSquarsEstimator()
        LS.setOptimizer('fmin')#('leastsq')   # 'leastsq' often faster, but if this does not work: 'fmin' is more reliable
        LS.setTolerance(1.e-7)
        LS.setWeights(self.WeightY,self.WeightR)
        outfile = open(PyrolProgram+'-Results_ArrheniusNoBRate.txt', 'w')
        outfile.write("Species\tA [1/s]\t\tE_a [K]\t\tFinalYield\n\n")
        #######
        ##The single species:
        #makes Species list which contains alls species to fit:
        SpeciesList=[]
        if self.CPD_ArrhSpec=='Total' and PyrolProgram=='CPD':
            SpeciesList.append(Fit[0].SpeciesIndex('Total'))
            if 'Total' not in self.SpeciesToConsider:
                self.SpeciesToConsider.append('Total')
        elif self.CPD_ArrhSpec=='MainSpecies' and PyrolProgram=='CPD':
            SpeciesList.append(Fit[0].SpeciesIndex('Total'))
            SpeciesList.append(Fit[0].SpeciesIndex('Tar'))
            SpeciesList.append(Fit[0].SpeciesIndex('Gas'))
            if 'Total' not in self.SpeciesToConsider:
                self.SpeciesToConsider.append('Total')
            if 'Tar' not in self.SpeciesToConsider:
                self.SpeciesToConsider.append('Tar')
            if 'Gas' not in self.SpeciesToConsider:
                self.SpeciesToConsider.append('Gas')
        elif self.CPD_ArrhSpec=='allSpecies' and PyrolProgram=='CPD':
            for i in range(2,len(Fit[0].SpeciesNames()),1):
                if Fit[0].SpeciesName(i) not in self.SpeciesToConsider:
                    self.SpeciesToConsider.append(Fit[0].SpeciesName(i))
                SpeciesList.append(i)
        elif self.FG_ArrhSpec=='Total' and PyrolProgram=='FGDVC':
            SpeciesList.append(Fit[0].SpeciesIndex('Total'))
            if 'Total' not in self.SpeciesToConsider:
                self.SpeciesToConsider.append('Total')
        elif self.FG_ArrhSpec=='MainSpecies' and PyrolProgram=='FGDVC':
            SpeciesList.append(Fit[0].SpeciesIndex('Total'))
            SpeciesList.append(Fit[0].SpeciesIndex('Tar'))
            SpeciesList.append(Fit[0].SpeciesIndex('Gas'))
            if 'Total' not in self.SpeciesToConsider:
                self.SpeciesToConsider.append('Total')
            if 'Tar' not in self.SpeciesToConsider:
                self.SpeciesToConsider.append('Tar')
            if 'Gas' not in self.SpeciesToConsider:
                self.SpeciesToConsider.append('Gas')
        elif self.FG_ArrhSpec=='allSpecies' and PyrolProgram=='FGDVC':
            for i in range(2,len(Fit[0].SpeciesNames()),1):
                if Fit[0].SpeciesName(i) not in self.SpeciesToConsider:
                    self.SpeciesToConsider.append(Fit[0].SpeciesName(i))
                SpeciesList.append(i)
        ##The single species:
        for Species in SpeciesList:
            m_final_prediction=Fit[0].Yield(Species)[-1]
            PredictionV0=[0.86e15,27700,m_final_prediction]  #for Standard Arrhenius
            Arr=Models.ArrheniusModelNoB(PredictionV0)
            #
            print Fit[0].SpeciesName(Species)
            if UseGlobalOpt=='ManyPoints':
                #GlobalOptimize:
                GlobalMin=Fitter.GlobalOptimizer(LS,Arr,Fit)
                m_final_predictionAll=[]
                for i in range(len(Fit)):
                    m_final_predictionAll.append(Fit[i].Yield(Species)[-1])
                LS.setTolerance(1e-2)
                ParamGlobalMin=GlobalMin.GenerateOptima(Species,GlobalOptParam.ArrhNobIndexToOptimize,[[min(m_final_predictionAll),max(m_final_predictionAll)]],GlobalOptParam.ArrhNoBNrOfRuns)
                #
                LS.setTolerance(1.e-7)
                print 'Final optimization Run:'
                Arr.setParamVector(LS.estimate_T(Fit,Arr,ParamGlobalMin,Species))
            if UseGlobalOpt=='Evolve':
                m_final_predictionAll=[]
                for i in range(len(Fit)):
                    m_final_predictionAll.append(Fit[i].Yield(Species)[-1])
                GenAlg=Evolve.GenericOpt(Arr,Fit,Species)
                GenAlg.setWeights(GlobalOptParam.EvAWeightY,GlobalOptParam.EvAWeightY)
                GAArrhMinA = GlobalOptParam.EvAArrhMin[0]
                GAArrhMinE = GlobalOptParam.EvAArrhMin[2]
                GAArrhMaxA = GlobalOptParam.EvAArrhMax[0]
                GAArrhMaxE = GlobalOptParam.EvAArrhMax[2]
                GAArrhInit=GlobalOptParam.EvAArrhInit
                if len(GAArrhInit)==3:
                    GAArrhInit.append((max(m_final_predictionAll)+min(m_final_predictionAll))/2.)
                else:
                    GAArrhInit[3]=(max(m_final_predictionAll)+min(m_final_predictionAll))/2.
                GenAlg.setParamRanges(GAArrhInit.pop(1),[GAArrhMinA,GAArrhMinE,min(m_final_predictionAll)],[GAArrhMaxA,GAArrhMaxE,max(m_final_predictionAll)])
                GenAlg.setNrPopulation(GlobalOptParam.NrOfPopulation)
                GenAlg.setNrGenerations(GlobalOptParam.NrOfGeneration)
                Arr.setParamVector(GenAlg.mkResults())
                #
                #use afterwards local optimization
                #Arr.setParamVector(LS.estimate_T(Fit,Arr,Arr.ParamVector(),Species))
            if UseGlobalOpt==False:
                Arr.setParamVector(LS.estimate_T(Fit,Arr,Arr.ParamVector(),Species))
            Solution=Arr.ParamVector()
            Arr.plot(Fit,Species)
            if np.sum(Arr.ParamVector())!=np.sum(PredictionV0): #To avoid, a species with no yield is added to the parameter file
                outfile.write(str(Fit[0].SpeciesName(Species))+'\t'+str(Solution[0])+'\t'+str(Solution[1])+'\t'+str(Solution[2])+'\n')
        outfile.close()
        if oSystem=='Linux':
            shutil.move(PyrolProgram+'-Results_ArrheniusNoBRate.txt','Result/'+PyrolProgram+'-Results_ArrheniusNoB.txt')
        elif oSystem=='Windows':
            shutil.move(PyrolProgram+'-Results_ArrheniusNoBRate.txt','Result\\'+PyrolProgram+'-Results_ArrheniusNoB.txt')
        else:
            print "The name of the operating system couldn't be found."
    #
    #
    def MakeResults_Kob(self,PyrolProgram,File,Fit):
        """Generates the results for Kobayashi Rate."""
    ##KOBAYASHI RATE
#    elif (CPD_FittingKineticParameter_Select=='Kobayashi' and PyrolProgram=='CPD') or (FG_FittingKineticParameter_Select=='Kobayashi' and PyrolProgram=='FGDVC'): #Kob means Kobayashi
#        PredictionVKob2=[10,-16,8,-20]           #for Arrhenius notation #2 [b11,b21,b12,b22] with the second indice as the reaction
#        PredictionVKob0=[2e5,1.046e8/8314.33,1.3e7,1.674e8/8314.33,PAVM_daf/100.]
# 	PredictionVKob0=[7e5,8e7/8314.33,2.3e8,1.6e8/8314.33]#,PAVM_daf/100.]
        PredictionVKob0=[7e5,8e7/8314.33,2.3e8,1.6e8/8314.33,0.4,0.9]
        LS=Fitter.LeastSquarsEstimator()
        LS.setOptimizer('fmin')#('leastsq')   # 'leastsq' often faster, but if this does not work: 'fmin' is more reliable
        LS.setTolerance(1.e-9)
        LS.setMaxIter(2000)
        LS.setWeights(1.0,1.0)
        outfile = open(PyrolProgram+'-Results_KobayashiRate.txt', 'w')
        outfile.write("Species\t\tA1 [1/s]\t\tE_a1 [K]\t\tA2 [1/s]\t\t\tE_a2 [K]\talpha1 \t\t\talpha2 \n\n")
        Kob=Models.Kobayashi(PredictionVKob0)
        #######
        ##The single species:
        if 'Total' not in self.SpeciesToConsider:
            self.SpeciesToConsider.append('Total')
        for Species in [Fit[0].SpeciesIndex('Total')]:
            print Fit[0].SpeciesName(Species)
            if UseGlobalOpt=='ManyPoints':
                #GlobalOptimize:
                GlobalMin=Fitter.GlobalOptimizer(LS,Kob,Fit)
                LS.setTolerance(1e-2)
                ParamGlobalMin=GlobalMin.GenerateOptima(Species,GlobalOptParam.KobIndexToOptimize,GlobalOptParam.KobBoundaries,GlobalOptParam.KobNrOfRuns)
                #
                LS.setTolerance(1.e-7)
                print 'Final optimization Run:'
                Kob.setParamVector(LS.estimate_T(Fit,Kob,ParamGlobalMin,Species))
            if UseGlobalOpt=='Evolve':
                GenAlg=Evolve.GenericOpt(Kob,Fit,Species)
                GenAlg.setWeights(GlobalOptParam.EvAWeightY,GlobalOptParam.EvAWeightY)
                GenAlg.setParamRanges(GlobalOptParam.EvAKobInit,GlobalOptParam.EvAKobMin,GlobalOptParam.EvAKobMax)
                GenAlg.setNrPopulation(GlobalOptParam.NrOfPopulation)
                GenAlg.setNrGenerations(GlobalOptParam.NrOfGeneration)
                Kob.setParamVector(GenAlg.mkResults())
                #
                #use afterwards local optimization
                #Kob.setParamVector(LS.estimate_T(Fit,Kob,Kob.ParamVector(),Species))
            if UseGlobalOpt==False:
                Kob.setParamVector(LS.estimate_T(Fit,Kob,Kob.ParamVector(),Species))
            Solution=Kob.ParamVector()
            #
            Kob.plot(Fit,Species)
            outfile.write(str(Fit[0].SpeciesName(Species))+'\t\t'+str(Solution[0])+'\t\t'+str(Solution[1])+'\t\t'+str(Solution[2])+'\t\t'+str(Solution[3])+'\t\t'+str(Solution[4])+'\t\t'+str(Solution[5])+'\n')
        outfile.close()
        if oSystem=='Linux':
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
        LS=Fitter.LeastSquarsEstimator()
        LS.setOptimizer('fmin')#('leastsq')   # 'leastsq' often faster, but if this does not work: 'fmin' is more reliable
        LS.setTolerance(1.e-9)
        LS.setMaxIter(2000)
        LS.setWeights(1.0,1.0)
        outfile = open(PyrolProgram+'-Results_DAEM.txt', 'w')
        outfile.write("Species\t\tA1 [1/s]\t\tE_a1 [K]\t\tsigma [K]\t\t\tFinal Yield\n\n")
        DAEM=Models.DAEM(PredictionDAEM)
        DAEM.setNrOfActivationEnergies(GlobalOptParam.NrOFActivtionEnergies)
        #######
        ##The single species:
        if 'Total' not in self.SpeciesToConsider:
            self.SpeciesToConsider.append('Total')
        for Species in [Fit[0].SpeciesIndex('Total')]:
            print Fit[0].SpeciesName(Species)
            m_final_predictionAll=[]
            for i in range(len(Fit)):
                m_final_predictionAll.append(Fit[i].Yield(Species)[-1])
            if UseGlobalOpt=='ManyPoints':
                #GlobalOptimize:
                GlobalMin=Fitter.GlobalOptimizer(LS,DAEM,Fit)
                LS.setTolerance(1e-2)
                DAEMBdr=GlobalOptParam.DAEMBoundaries
                DAEMBdr.append([min(m_final_predictionAll),max(m_final_predictionAll)])
                ParamGlobalMin=GlobalMin.GenerateOptima(Species,GlobalOptParam.DAEMIndexToOptimize,DAEMBdr,GlobalOptParam.DAEMNrOfRuns)
                #
                LS.setTolerance(1.e-7)
                print 'Final optimization Run:'
                DAEM.setParamVector(LS.estimate_T(Fit,DAEM,ParamGlobalMin,Species))
            if UseGlobalOpt=='Evolve':
                GenAlg=Evolve.GenericOpt(DAEM,Fit,Species)
                GenAlg.setWeights(GlobalOptParam.EvAWeightY,GlobalOptParam.EvAWeightY)
                EvADAEMInit=GlobalOptParam.EvADAEMInit
                EvADAEMMin=GlobalOptParam.EvADAEMMin
                EvADAEMMax=GlobalOptParam.EvADAEMMax
                EvADAEMInit.append((min(m_final_predictionAll)+max(m_final_predictionAll))/2.)
                EvADAEMMin.append(min(m_final_predictionAll))
                EvADAEMMax.append(max(m_final_predictionAll))
                GenAlg.setParamRanges(EvADAEMInit,EvADAEMMin,EvADAEMMax)
                GenAlg.setNrPopulation(GlobalOptParam.NrOfPopulation)
                GenAlg.setNrGenerations(GlobalOptParam.NrOfGeneration)
                DAEM.setParamVector(GenAlg.mkResults())
                #
                #use afterwards local optimization
                #DAEM.setParamVector(LS.estimate_T(Fit,DAEM,DAEM.ParamVector(),Species))
            if UseGlobalOpt==False:
                DAEM.setParamVector(LS.estimate_T(Fit,DAEM,DAEM.ParamVector(),Species))
            Solution=DAEM.ParamVector()
            #
            DAEM.plot(Fit,Species)
            outfile.write(str(Fit[0].SpeciesName(Species))+'\t\t'+str(Solution[0])+'\t'+str(Solution[1])+'\t\t'+str(Solution[2])+'\t\t\t'+str(Solution[3])+'\n')
        outfile.close()
        if oSystem=='Linux':
            shutil.move(PyrolProgram+'-Results_DAEM.txt','Result/'+PyrolProgram+'-Results_DAEM.txt')
        elif oSystem=='Windows':
            shutil.move(PyrolProgram+'-Results_DAEM.txt','Result\\'+PyrolProgram+'-Results_DAEM.txt')
        else:
            print "The name of the operating system couldn't be found."
#
#
    def SpeciesEnergy(self,PyrolProgram,File):
        """Carries out the species and Energy balance."""
        ##SPECIES AND ENERGY BALANCE:
        for runNr in range(self.NrOfRuns):
            if PyrolProgram=='CPD':
                print 'CPD energy and mass balance...'
                Compos_and_Energy.CPD_SpeciesBalance(File[runNr],self.UAC,self.UAH,self.UAN,self.UAO,self.UAS,self.PAVM_asrec,self.PAFC_asrec,self.PAmoist,self.PAash,self.HHV,self.MTar,runNr)
            if PyrolProgram=='FGDVC':    
                print 'FG-DVC energy and mass balance...'
                Compos_and_Energy.FGDVC_SpeciesBalance(File[runNr],self.UAC,self.UAH,self.UAN,self.UAO,self.UAS,self.PAVM_asrec,self.PAFC_asrec,self.PAmoist,self.PAash,self.HHV,self.MTar,runNr)
                #    SpecCPD=Compos_and_Energy.CPD_SpeciesBalance(File[0],UAC,UAH,UAN,UAO,PAVM_asrec,PAFC_asrec,HHV,MTar,0)
#
#
#
####CPD#####
    def MakeResults_CPD(self):
        """generates the result for CPD"""
        CPDFile=[]
        CPDFit=[]
        for runNr in range(self.NrOfRuns):
            #launches CPD
            CPD=CPD_SetAndLaunch.SetterAndLauncher()
            CPD.SetCoalParameter(self.UAC,self.UAH,self.UAN,self.UAO,self.PAVM_daf)
            CPD.CalcCoalParam()
            if runNr==0:
                CPD.SetOperateCond(self.CPD_pressure,self.CPD_TimeTemp1)
                CPD.SetNumericalParam(self.CPDdt,self.CPD_t_max1)
            elif runNr==1:
                CPD.SetOperateCond(self.CPD_pressure,self.CPD_TimeTemp2)
                CPD.SetNumericalParam(self.CPDdt,self.CPD_t_max2)
            elif runNr==2:
                CPD.SetOperateCond(self.CPD_pressure,self.CPD_TimeTemp3)
                CPD.SetNumericalParam(self.CPDdt,self.CPD_t_max3)
            elif runNr==3:
                CPD.SetOperateCond(self.CPD_pressure,self.CPD_TimeTemp4)
                CPD.SetNumericalParam(self.CPDdt,self.CPD_t_max4)
            elif runNr==4:
                CPD.SetOperateCond(self.CPD_pressure,self.CPD_TimeTemp5)
                CPD.SetNumericalParam(self.CPDdt,self.CPD_t_max5)
            CPD.writeInstructFile(workingDir)
            print 'Running CPD ...',runNr
            if oSystem=='Linux':
                CPD.Run('./'+'cpdnlg','IN.dat','CPD_'+str(runNr)+'_output.log')   #first Arg: CPD-executeable, second: Input data containing CPD input file and the output files
            elif oSystem=='Windows':
                CPD.Run('cpdnlg.exe','IN.dat','CPD_'+str(runNr)+'_output.log')   #first Arg: CPD-executeable, second: Input data containing CPD input file and the output files
            else:
                print "The name of the operating system couldn't be found."
            #
            ###calibration of the kinetic parameter:
            #read result:
            CurrentCPDFile=CPD_Result.CPD_Result(workingDir)
            # creates object, required for fitting procedures
            CurrentCPDFit=FitInfo.Fit_one_run(CurrentCPDFile)
            CPDFile.append(CurrentCPDFile)
            CPDFit.append(CurrentCPDFit)
            #
            if oSystem=='Linux':
                shutil.move('CPD_Result1.dat', 'Result/'+'CPD_'+str(runNr)+'_Result1.dat')
                shutil.move('CPD_Result2.dat', 'Result/'+'CPD_'+str(runNr)+'_Result2.dat')
                shutil.move('CPD_Result3.dat', 'Result/'+'CPD_'+str(runNr)+'_Result3.dat')
                shutil.move('CPD_Result4.dat', 'Result/'+'CPD_'+str(runNr)+'_Result4.dat')
                shutil.move('CPD_'+str(runNr)+'_output.log', 'Result/'+'CPD_'+str(runNr)+'_output.log')
            elif oSystem=='Windows':
                shutil.move('CPD_Result1.dat', 'Result\\'+'CPD_'+str(runNr)+'_Result1.dat')
                shutil.move('CPD_Result2.dat', 'Result\\'+'CPD_'+str(runNr)+'_Result2.dat')
                shutil.move('CPD_Result3.dat', 'Result\\'+'CPD_'+str(runNr)+'_Result3.dat')
                shutil.move('CPD_Result4.dat', 'Result\\'+'CPD_'+str(runNr)+'_Result4.dat')
                shutil.move('CPD_'+str(runNr)+'_output.log', 'Result\\'+'CPD_'+str(runNr)+'_output.log')
            else:
                print "The name of the operating system couldn't be found."
        #####
        if self.CPD_FittingKineticParameter_Select=='constantRate':
            self.MakeResults_CR('CPD',CPDFile,CPDFit)
            currentDict={'CPD':'constantRate'}
        elif self.CPD_FittingKineticParameter_Select=='Arrhenius':
            self.MakeResults_Arrh('CPD',CPDFile,CPDFit)
            currentDict={'CPD':'Arrhenius'}
        elif self.CPD_FittingKineticParameter_Select=='ArrheniusNoB':
            self.MakeResults_ArrhNoB('CPD',CPDFile,CPDFit)
            currentDict={'CPD':'ArrheniusNoB'}
        elif self.CPD_FittingKineticParameter_Select=='Kobayashi':
            self.MakeResults_Kob('CPD',CPDFile,CPDFit)
            currentDict={'CPD':'Kobayashi'}
        elif self.CPD_FittingKineticParameter_Select=='DAEM':
            self.MakeResults_DEAM('CPD',CPDFile,CPDFit)
            currentDict={'CPD':'DAEM'}
        else:
            print 'uspecified CPD_FittingKineticParameter_Select'
            currentDict={}
        #
        self.ProgramModelDict.update(currentDict)
        #
        self.SpeciesEnergy('CPD',CPDFile)
            #
            #
    ####FG-DVC####
    def MakeResults_FG(self):
        """generates the result for FG-DVC"""
        #writes Time-Temperature file
        FG_TimeTemp1=self.CPD_TimeTemp1
        FG_TimeTemp2=self.CPD_TimeTemp2
        FG_TimeTemp3=self.CPD_TimeTemp3
        FG_TimeTemp4=self.CPD_TimeTemp4
        FG_TimeTemp5=self.CPD_TimeTemp5
        FG_TimeTemp1[:,0]=self.CPD_TimeTemp1[:,0]*1.e-3
        FG_TimeTemp2[:,0]=self.CPD_TimeTemp2[:,0]*1.e-3
        FG_TimeTemp3[:,0]=self.CPD_TimeTemp3[:,0]*1.e-3
        FG_TimeTemp4[:,0]=self.CPD_TimeTemp4[:,0]*1.e-3
        FG_TimeTemp5[:,0]=self.CPD_TimeTemp5[:,0]*1.e-3
        #initialize the launching object
        FGDVC=FGDVC_SetAndLaunch.SetterAndLauncher()
        #set and writes Coal Files:
        if self.FG_CoalSelection==0:
            #deletes old generated file
            os.system('cd '+self.FG_MainDir+FG_GenCoalDir+' & del '+FG_CoalName+'_com.dat, '+FG_CoalName+'_kin.dat, '+FG_CoalName+'_pol.dat')
            #generates coalsd.exe input file
            MakeCoalGenFile=InformationFiles.WriteFGDVCCoalFile(FG_CoalGenFileName)
            MakeCoalGenFile.setCoalComp(self.UAC,self.UAH,self.UAO,self.UAN,self.UAS,0)
            MakeCoalGenFile.write(self.FG_MainDir+FG_GenCoalDir+'\\',FG_CoalName)
            #makes new file
            try:
                os.system('cd '+self.FG_MainDir+FG_GenCoalDir+' & '+'coalsd.exe < '+FG_CoalGenFileName+' > coalsd_pkp.log')
            except OSError:
                print 'Problems with coalsd.exe'
            os.system('copy '+self.FG_MainDir+FG_GenCoalDir+'\coalsd_pkp.log . >> log.txt')
            #tests weather the coal file was genearated:
            if os.path.exists(self.FG_MainDir+'\\'+FG_GenCoalDir+'\\'+FG_CoalName+'_com.dat')==False:
                print 30*'*','\n','The coal is may outside the libraries coals. Select manually the closest library coal.',30*'*','\n'
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
            if oSystem=='Linux':
                shutil.copyfile(self.FG_DirOut+'gasyield.txt', 'Result/gasyield_'+str(runNr)+'.txt')
                shutil.copyfile(self.FG_DirOut+'gasrate.txt', 'Result/gasrate_'+str(runNr)+'.txt')
            elif oSystem=='Windows':
                shutil.copyfile(self.FG_DirOut+'gasyield.txt', 'Result\\gasyield_'+str(runNr)+'.txt')
                shutil.copyfile(self.FG_DirOut+'gasrate.txt', 'Result\\gasrate_'+str(runNr)+'.txt')
        #####
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
        else:
            print 'uspecified FG_FittingKineticParameter_Select'
            currentDict={}
        #
        self.ProgramModelDict.update(currentDict)
        #
        self.SpeciesEnergy('FGDVC',FGFile)
            #


#Main Part starting
if __name__ == "__main__":
    Case=MainProcess()
    Case.ReadInputFiles()
    if Case.CPDselect==True:
        Case.MakeResults_CPD()
    if Case.FG_select==True:
        Case.CheckFGdt()
        Case.MakeResults_FG()
    print 'calculated Species: ',Case.SpeciesToConsider
        
