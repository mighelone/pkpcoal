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
def DAF(PAFC_asRecieved,PAVM_asRecieved):
    """calculates PAFC, PAVM  from the as recieved state to the daf state of coal"""
    fractionFC=PAFC_asRecieved/(PAFC_asRecieved+PAVM_asRecieved)
    fractionVM=PAVM_asRecieved/(PAFC_asRecieved+PAVM_asRecieved)
    return 100.*fractionFC, 100.*fractionVM
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
#get parameters from input files:
#
#Coal File:
#
print 'Reading Coal.inp ...'
CoalInput=InformationFiles.ReadFile(workingDir+'Coal.inp')
PAFC_asrec=CoalInput.getValue(InformationFiles.M_PA[0])
PAVM_asrec=CoalInput.getValue(InformationFiles.M_PA[1])
PAmoist = CoalInput.getValue(InformationFiles.M_PA[2])
PAash = CoalInput.getValue(InformationFiles.M_PA[3])
# scale proximate analysis
sumPA = (PAFC_asrec+PAVM_asrec + PAmoist + PAash)/100.
PAFC_asrec/=sumPA
PAVM_asrec/=sumPA
PAmoist/=sumPA
PAash/=sumPA



#gets daf values, as CPD needs daf as input:
PAFC_daf, PAVM_daf = DAF(PAFC_asrec,PAVM_asrec)
UAC=CoalInput.getValue(InformationFiles.M_UA[0])
UAH=CoalInput.getValue(InformationFiles.M_UA[1])
UAN=CoalInput.getValue(InformationFiles.M_UA[2])
UAO=CoalInput.getValue(InformationFiles.M_UA[3])
UAS=CoalInput.getValue(InformationFiles.M_UA[4])
# scale ultimate analysis
sumUA = UAC+UAH+UAN+UAO+UAS
UAC=UAC/sumUA*100
UAH=UAH/sumUA*100
UAO=UAO/sumUA*100
UAN=UAN/sumUA*100
UAS=UAS/sumUA*100
HHV=CoalInput.getValue(InformationFiles.M_HHV)
MTar=CoalInput.getValue(InformationFiles.M_MTar)
WeightY=CoalInput.getValue(InformationFiles.M_Weight[0])
WeightR=CoalInput.getValue(InformationFiles.M_Weight[1])
#
#
#CPD Properties:
#
print 'Reading CPD.inp ...'
CPDInput=InformationFiles.ReadFile(workingDir+'CPD.inp')
CPDselect=CPDInput.UsePyrolProgr(InformationFiles.MC_sel)
CPD_FittingKineticParameter_Select=CPDInput.Fitting(InformationFiles.M_selFit)
CPD_ArrhSpec=CPDInput.getText(InformationFiles.M_selArrhSpec)
CPDdt=[0,1,2] #0:initila dt, 1: print increment, 2: dt max
CPDdt[0]=(CPDInput.getValue(InformationFiles.MC_dt[0]))
CPDdt[1]=(CPDInput.getValue(InformationFiles.MC_dt[1]))
#
#
#FG-DVC Properties:
#
print 'Reading FGDVC.inp ...'
FGDVCInput=InformationFiles.ReadFile(workingDir+'FGDVC.inp')
FG_select=FGDVCInput.UsePyrolProgr(InformationFiles.MF_sel)
FG_FittingKineticParameter_Select=FGDVCInput.Fitting(InformationFiles.M_selFit)
FG_ArrhSpec=FGDVCInput.getText(InformationFiles.M_selArrhSpec)
FG_CoalSelection=int(FGDVCInput.getValue(InformationFiles.MF_CoalSel))
FG_MainDir=FGDVCInput.getText(InformationFiles.MF_dir[0])
FG_DirOut=FGDVCInput.getText(InformationFiles.MF_dir[1])
FG_TarCacking=FGDVCInput.getValue(InformationFiles.MF_TarCr)
#
#
#Operating Condition File:
#
print 'Reading OperCond.inp ...'
OpCondInp=InformationFiles.OperCondInput('OperCond.inp')
CPD_pressure=OpCondInp.getValue(InformationFiles.M_Pressure)
FG_pressure=OpCondInp.getValue(InformationFiles.M_Pressure)
#Number of FG-DVC/CPD/PCCL runs:
NrOfRuns=int(OpCondInp.getValue(InformationFiles.M_NrRuns))
CPD_TimeTemp1=OpCondInp.getTimePoints(InformationFiles.M_TimePoints1[0],InformationFiles.M_TimePoints1[1])
CPD_TimeTemp2=OpCondInp.getTimePoints(InformationFiles.M_TimePoints2[0],InformationFiles.M_TimePoints2[1])
CPD_TimeTemp3=OpCondInp.getTimePoints(InformationFiles.M_TimePoints3[0],InformationFiles.M_TimePoints3[1])
CPD_TimeTemp4=OpCondInp.getTimePoints(InformationFiles.M_TimePoints4[0],InformationFiles.M_TimePoints4[1])
CPD_TimeTemp5=OpCondInp.getTimePoints(InformationFiles.M_TimePoints5[0],InformationFiles.M_TimePoints5[1])
CPDdt[2]=OpCondInp.getValue(InformationFiles.M_dt)
FG_dt=OpCondInp.getValue(InformationFiles.M_dt)
FG_T_t_History=FG_MainDir+'tTHistory.txt'
#makes for CPD time in milliseconds:
CPD_TimeTemp1[:,0]=CPD_TimeTemp1[:,0]*1.e3
CPD_TimeTemp2[:,0]=CPD_TimeTemp2[:,0]*1.e3
CPD_TimeTemp3[:,0]=CPD_TimeTemp3[:,0]*1.e3
CPD_TimeTemp4[:,0]=CPD_TimeTemp4[:,0]*1.e3
CPD_TimeTemp5[:,0]=CPD_TimeTemp5[:,0]*1.e3
CPD_t_max1=CPD_TimeTemp1[-1,0]*1.e-3 #tmax in s, not ms
CPD_t_max2=CPD_TimeTemp2[-1,0]*1.e-3 #tmax in s, not ms
CPD_t_max3=CPD_TimeTemp3[-1,0]*1.e-3 #tmax in s, not ms
CPD_t_max4=CPD_TimeTemp4[-1,0]*1.e-3 #tmax in s, not ms
CPD_t_max5=CPD_TimeTemp5[-1,0]*1.e-3 #tmax in s, not ms
#
#
#
#Aborts, if FG-DVC is selected and the timestep is lower than 1.e-3 (which is FG-DVC not able to read):
if ((FG_select==True) and (FG_dt<1e-4)):
    print "Please select for FG-DVC a time step greather equal 1e-4 in 'OperCond.inp'. FG-DVC would not be able to read the time history file for a dt<1e-4."
    sys.exit()
#
#
def MakeResults(PyrolProgram,File,Fit):
    #Array for the comparison of the sum of the individual species(t) with the (1-solid(t)); this part ([0]) has to be modified, when the assumption of same length of the array (e.g. when different final times are possible) is refused
    SumSingleYieldsCalc=np.zeros([len(Fit[0].Yield('Time')),NrOfRuns]) #Array initialized for the sum of the single yields (calculated)
    SumSingleYields=np.zeros([len(Fit[0].Yield('Time')),NrOfRuns])  #Array initialized for the sum of the single yields (CPD output)
    SolidYieldsCalc=np.zeros([len(Fit[0].Yield('Time')),NrOfRuns])     #Array initialized for the yields of the Solids (calculated)
    SolidYields=np.zeros([len(Fit[0].Yield('Time')),NrOfRuns])     #Array initialized for the yields of the Solids (CPD output)
    ##CONSTANT RATE
    if (CPD_FittingKineticParameter_Select=='constantRate' and PyrolProgram=='CPD') or (FG_FittingKineticParameter_Select=='constantRate' and PyrolProgram=='FGDVC'): #CR means ConstantRate
        LS=Fitter.LeastSquarsEstimator() #LS means LeastSquares
        LS.setOptimizer('fmin')
        LS.setTolerance(1.e-18)
        LS.setWeights(WeightY,WeightR)
        outfile = open(PyrolProgram+'-Results_const_rate.txt', 'w')
        outfile.write("Species\tk [1/s]\t\tt_start [s]\t\tFinalYield\n\n")
        for Spec in range(2,len(Fit[0].SpeciesNames()),1):
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
                outfile.write(str(Fit[0].SpeciesName(Spec))+'\t'+str(Solution[0])+'\t'+str(Solution[1])+str(Solution[0])+'\n')
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
            shutil.move(PyrolProgram+'-Results_const_rate.txt','Result/'+PyrolProgram+'-Results_const_rate.txt')
        elif oSystem=='Windows':
            shutil.move(PyrolProgram+'-Results_const_rate.txt','Result\\'+PyrolProgram+'-Results_const_rate.txt')
        else:
            print "The name of the operating system couldn't be found."
#        Fit[0].plt_InputVectors(Fit[0].Time(),1.-SolidYieldsCalc,1.-SolidYields,SumSingleYieldsCalc,SumSingleYields,'1-Solid; fitted','1-Solid; CPD output','Sum Yields; fitted','Sum Yields; CPD output')
    ##ARRHENIUS RATE
    elif (CPD_FittingKineticParameter_Select=='Arrhenius' and PyrolProgram=='CPD') or (FG_FittingKineticParameter_Select=='Arrhenius' and PyrolProgram=='FGDVC'): #Arr means Arrhenius
        LS=Fitter.LeastSquarsEstimator()
        LS.setOptimizer('fmin')#('leastsq')   # 'leastsq' often faster, but if this does not work: 'fmin' is more reliable
        LS.setTolerance(1.e-7)
        LS.setWeights(WeightY,WeightR)
        outfile = open(PyrolProgram+'-Results_ArrheniusRate.txt', 'w')
        outfile.write("Species\tA [1/s]\t\tb\t\tE_a [K]\t\tFinalYield\n\n")
        #select one of the follwoing notations: 
        #Arr=Models.ArrheniusModel(PredictionV0)
        #Arr=Models.ArrheniusModelAlternativeNotation1(PredictionV1)
        #######
        #makes Species list which contains alls species to fit:
        SpeciesList=[]
        if CPD_ArrhSpec=='Total' and PyrolProgram=='CPD':
            SpeciesList.append(Fit[0].SpeciesIndex('Total'))
        elif CPD_ArrhSpec=='allSpecies' and PyrolProgram=='CPD':
            SpeciesList.append(Fit[0].SpeciesIndex('Total'))
            SpeciesList.append(Fit[0].SpeciesIndex('Tar'))
        elif CPD_ArrhSpec=='allSpecies' and PyrolProgram=='CPD':
            for i in range(2,len(Fit[0].SpeciesNames()),1):
                SpeciesList.append(i)
        elif FG_ArrhSpec=='Total' and PyrolProgram=='FGDVC':
            SpeciesList.append(Fit[0].SpeciesIndex('Total'))
        elif FG_ArrhSpec=='allSpecies' and PyrolProgram=='FGDVC':
            SpeciesList.append(Fit[0].SpeciesIndex('Total'))
            SpeciesList.append(Fit[0].SpeciesIndex('Tar'))
        elif FG_ArrhSpec=='allSpecies' and PyrolProgram=='FGDVC':
            for i in range(2,len(Fit[0].SpeciesNames()),1):
                SpeciesList.append(i)
        ##The single species:
        for Species in SpeciesList:
            #
            m_final_prediction=Fit[0].Yield(Species)[-1]
            PredictionV0=[0.86e15,0.01,27700,m_final_prediction]  #for Standard Arrhenius
            PredictionV1=[10.,-20.,m_final_prediction]         #for Arrhenius notation #1
            PredictionV2=[10.,-18.,m_final_prediction]           #for Arrhenius notation #2
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
                GAArrhInit.append((max(m_final_predictionAll)+min(m_final_predictionAll))/2.)
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
            shutil.move(PyrolProgram+'-Results_ArrheniusRate.txt','Result/'+PyrolProgram+'-Results_ArrheniusRate.txt')
        elif oSystem=='Windows':
            shutil.move(PyrolProgram+'-Results_ArrheniusRate.txt','Result\\'+PyrolProgram+'-Results_ArrheniusRate.txt')
        else:
            print "The name of the operating system couldn't be found."
#        Fit[0].plt_InputVectors(Fit[0].Time(),1.-SolidYieldsCalc,1.-SolidYields,SumSingleYieldsCalc,SumSingleYields,'1-Solid; fitted','1-Solid; CPD output','Sum Yields; fitted','Sum Yields; CPD output')
    elif (CPD_FittingKineticParameter_Select=='ArrheniusNoB' and PyrolProgram=='CPD') or (FG_FittingKineticParameter_Select=='ArrheniusNoB' and PyrolProgram=='FGDVC'): #Arr means Arrhenius
        LS=Fitter.LeastSquarsEstimator()
        LS.setOptimizer('fmin')#('leastsq')   # 'leastsq' often faster, but if this does not work: 'fmin' is more reliable
        LS.setTolerance(1.e-7)
        LS.setWeights(WeightY,WeightR)
        outfile = open(PyrolProgram+'-Results_ArrheniusNoBRate.txt', 'w')
        outfile.write("Species\tA [1/s]\t\tb\t\tE_a [K]\t\tFinalYield\n\n")
        #select one of the follwoing notations: 
        #Arr=Models.ArrheniusModel(PredictionV0)
        #Arr=Models.ArrheniusModelAlternativeNotation1(PredictionV1)
        #######
        ##The single species:
        #makes Species list which contains alls species to fit:
        SpeciesList=[]
        if CPD_ArrhSpec=='Total' and PyrolProgram=='CPD':
            SpeciesList.append(Fit[0].SpeciesIndex('Total'))
        elif CPD_ArrhSpec=='MainSpecies' and PyrolProgram=='CPD':
            SpeciesList.append(Fit[0].SpeciesIndex('Total'))
            SpeciesList.append(Fit[0].SpeciesIndex('Tar'))
            SpeciesList.append(Fit[0].SpeciesIndex('Gas'))
        elif CPD_ArrhSpec=='allSpecies' and PyrolProgram=='CPD':
            for i in range(2,len(Fit[0].SpeciesNames()),1):
                SpeciesList.append(i)
        elif FG_ArrhSpec=='Total' and PyrolProgram=='FGDVC':
            SpeciesList.append(Fit[0].SpeciesIndex('Total'))
        elif FG_ArrhSpec=='MainSpecies' and PyrolProgram=='FGDVC':
            SpeciesList.append(Fit[0].SpeciesIndex('Total'))
            SpeciesList.append(Fit[0].SpeciesIndex('Tar'))
            SpeciesList.append(Fit[0].SpeciesIndex('Gas'))
        elif FG_ArrhSpec=='allSpecies' and PyrolProgram=='FGDVC':
            for i in range(2,len(Fit[0].SpeciesNames()),1):
                SpeciesList.append(i)
        ##The single species:
        for Species in SpeciesList:
            m_final_prediction=Fit[0].Yield(Species)[-1]
            PredictionV0=[0.86e15,0,27700,m_final_prediction]  #for Standard Arrhenius
            PredictionV2=[10.,-18.,m_final_prediction]           #for Arrhenius notation #2
            Arr=Models.ArrheniusModelAlternativeNotation2(PredictionV2)
            ArrPlot=Models.ArrheniusModel([0,0,0,0]) #use Original Arrhenius Model to Plot
            #
            Arr.setMinMaxTemp(Fit[0].Yield('Temp')[0],Fit[0].Yield('Temp')[-1])
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
                GAArrhInit.append((max(m_final_predictionAll)+min(m_final_predictionAll))/2.)
                GenAlg.setParamRanges(Arr.ConvertKinFactorsToOwnNotation(GAArrhInit),Arr.ConvertKinFactorsToOwnNotation([GAArrhMinA,0,GAArrhMinE,min(m_final_predictionAll)]),Arr.ConvertKinFactorsToOwnNotation([GAArrhMaxA,0,GAArrhMaxE,max(m_final_predictionAll)]))
                Arr.setParamVector(GenAlg.mkResults())
                #
                #use afterwards local optimization
                #Arr.setParamVector(LS.estimate_T(Fit,Arr,Arr.ParamVector(),Species))
            if UseGlobalOpt==False:
                Arr.setParamVector(LS.estimate_T(Fit,Arr,Arr.ParamVector(),Species))
            Solution=Arr.ConvertKinFactors(Arr.ParamVector())
            ArrPlot.setParamVector(Solution)
            ArrPlot.plot(Fit,Species)
            if np.sum(Arr.ParamVector())!=np.sum(PredictionV2): #To avoid, a species with no yield is added to the parameter file
                outfile.write(str(Fit[0].SpeciesName(Species))+'\t'+str(Solution[0])+'\t'+str(Solution[1])+'\t'+str(Solution[2])+'\t\t'+str(Solution[3])+'\n')
        outfile.close()
        if oSystem=='Linux':
            shutil.move(PyrolProgram+'-Results_ArrheniusNoBRate.txt','Result/'+PyrolProgram+'-Results_ArrheniusNoBRate.txt');print 'mv'
        elif oSystem=='Windows':
            shutil.move(PyrolProgram+'-Results_ArrheniusNoBRate.txt','Result\\'+PyrolProgram+'-Results_ArrheniusNoBRate.txt')
        else:
            print "The name of the operating system couldn't be found."
    ##KOBAYASHI RATE
    elif (CPD_FittingKineticParameter_Select=='Kobayashi' and PyrolProgram=='CPD') or (FG_FittingKineticParameter_Select=='Kobayashi' and PyrolProgram=='FGDVC'): #Kob means Kobayashi
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
            shutil.move(PyrolProgram+'-Results_KobayashiRate.txt','Result/'+PyrolProgram+'-Results_KobayashiRate.txt')
        elif oSystem=='Windows':
            shutil.move(PyrolProgram+'-Results_KobayashiRate.txt','Result\\'+PyrolProgram+'-Results_KobayashiRate.txt')
        else:
            print "The name of the operating system couldn't be found."
    elif (CPD_FittingKineticParameter_Select=='DAEM' and PyrolProgram=='CPD') or (FG_FittingKineticParameter_Select=='DAEM' and PyrolProgram=='FGDVC'):
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
    ##SPECIES AND ENERGY BALANCE:
    for runNr in range(NrOfRuns):
        if PyrolProgram=='CPD':
	    print 'CPD energy and mass balance...'
            Compos_and_Energy.CPD_SpeciesBalance(File[runNr],UAC,UAH,UAN,UAO,UAS,PAVM_asrec,PAFC_asrec,PAmoist,PAash,HHV,MTar,runNr)
        if PyrolProgram=='FGDVC':    
            Compos_and_Energy.FGDVC_SpeciesBalance(FGFile[runNr],UAC,UAH,UAN,UAO,UAS,PAVM_asrec,PAFC_asrec,PAmoist,PAash,HHV,MTar,runNr)
#    SpecCPD=Compos_and_Energy.CPD_SpeciesBalance(File[0],UAC,UAH,UAN,UAO,PAVM_asrec,PAFC_asrec,HHV,MTar,0)
#
#
#
####CPD####
if CPDselect==True:
    CPDFile=[]
    CPDFit=[]
    for runNr in range(NrOfRuns):
        #launches CPD
        CPD=CPD_SetAndLaunch.SetterAndLauncher()
        CPD.SetCoalParameter(UAC,UAH,UAN,UAO,PAVM_daf)
        CPD.CalcCoalParam()
        if runNr==0:
            CPD.SetOperateCond(CPD_pressure,CPD_TimeTemp1)
            CPD.SetNumericalParam(CPDdt,CPD_t_max1)
        elif runNr==1:
            CPD.SetOperateCond(CPD_pressure,CPD_TimeTemp2)
            CPD.SetNumericalParam(CPDdt,CPD_t_max2)
        elif runNr==2:
            CPD.SetOperateCond(CPD_pressure,CPD_TimeTemp3)
            CPD.SetNumericalParam(CPDdt,CPD_t_max3)
        elif runNr==3:
            CPD.SetOperateCond(CPD_pressure,CPD_TimeTemp4)
            CPD.SetNumericalParam(CPDdt,CPD_t_max4)
        elif runNr==4:
            CPD.SetOperateCond(CPD_pressure,CPD_TimeTemp5)
            CPD.SetNumericalParam(CPDdt,CPD_t_max5)
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
    MakeResults('CPD',CPDFile,CPDFit)
####FG-DVC####
if FG_select==True:
    #writes Time-Temperature file
    FG_TimeTemp1=CPD_TimeTemp1
    FG_TimeTemp2=CPD_TimeTemp2
    FG_TimeTemp3=CPD_TimeTemp3
    FG_TimeTemp4=CPD_TimeTemp4
    FG_TimeTemp5=CPD_TimeTemp5
    FG_TimeTemp1[:,0]=CPD_TimeTemp1[:,0]*1.e-3
    FG_TimeTemp2[:,0]=CPD_TimeTemp2[:,0]*1.e-3
    FG_TimeTemp3[:,0]=CPD_TimeTemp3[:,0]*1.e-3
    FG_TimeTemp4[:,0]=CPD_TimeTemp4[:,0]*1.e-3
    FG_TimeTemp5[:,0]=CPD_TimeTemp5[:,0]*1.e-3
    #initialize the launching object
    FGDVC=FGDVC_SetAndLaunch.SetterAndLauncher()
    #set and writes Coal Files:
    if FG_CoalSelection==0:
        #deletes old generated file
        os.system('cd '+FG_MainDir+FG_GenCoalDir+' & del '+FG_CoalName+'_com.dat, '+FG_CoalName+'_kin.dat, '+FG_CoalName+'_pol.dat')
        #generates coalsd.exe input file
        MakeCoalGenFile=InformationFiles.WriteFGDVCCoalFile(FG_CoalGenFileName)
        MakeCoalGenFile.setCoalComp(UAC,UAH,UAO,UAN,UAS,0)
        MakeCoalGenFile.write(FG_MainDir+FG_GenCoalDir+'\\',FG_CoalName)
        #makes new file
        try:
            os.system('cd '+FG_MainDir+FG_GenCoalDir+' & '+'coalsd.exe < '+FG_CoalGenFileName+' > coalsd_pkp.log')
        except OSError:
            print 'Problems with coalsd.exe'
	os.system('copy '+FG_MainDir+FG_GenCoalDir+'\coalsd_pkp.log . >> log.txt')
        #tests weather the coal file was genearated:
        if os.path.exists(FG_MainDir+'\\'+FG_GenCoalDir+'\\'+FG_CoalName+'_com.dat')==False:
            print 30*'*','\n','The coal is may outside the libraries coals. Select manually the closest library coal.',30*'*','\n'
        #sets generated file for instruct.ini
        FGDVC.set1CoalLocation(FG_MainDir+FG_GenCoalDir+'\\'+FG_CoalName+'_com.dat')
        FGDVC.set2KinLocation(FG_MainDir+FG_GenCoalDir+'\\'+FG_CoalName+'_kin.dat')
        FGDVC.set3PolyLocation(FG_MainDir+FG_GenCoalDir+'\\'+FG_CoalName+'_pol.dat')
    elif FG_CoalSelection>0 and FG_CoalSelection<9:
        #sets library file for instruct.ini
        FGDVC.set1CoalLocation(FG_MainDir+FG_LibCoalDir+'\\coal.ar'+str(FG_CoalSelection))
        FGDVC.set2KinLocation(FG_MainDir+FG_LibCoalDir+'\\kin.ar'+str(FG_CoalSelection))
        FGDVC.set3PolyLocation(FG_MainDir+FG_LibCoalDir+'\\polymr.ar'+str(FG_CoalSelection))
    else:
        print "select Choose Coal: 0 interpolate between library coals and generate own coal. Set 1 to 8 for a library coal.' in FGDVC.inp equal a value between 0 and 8"
    #sets FG-DVC instruct.ini parameter
    FGDVC.set5Pressure(FG_pressure)
    if FG_TarCacking==0.0:            #case: no tar cracking
        FGDVC.set6Theorie(13,0.0)
    elif FG_TarCacking<0.0:           #case: full tar cracking
        FGDVC.set6Theorie(15,0.0)
    else:                             #case: partial tar cracking
        FGDVC.set6Theorie(13,float(FG_TarCacking))
    #
    FGFile=[]
    FGFit=[]
    for runNr in range(NrOfRuns):
        if runNr==0:
            OpCondInp.writeFGDVCtTHist(FG_TimeTemp1,FG_dt,FG_T_t_History)
        elif runNr==1:
            OpCondInp.writeFGDVCtTHist(FG_TimeTemp2,FG_dt,FG_T_t_History)
        elif runNr==2:
            OpCondInp.writeFGDVCtTHist(FG_TimeTemp3,FG_dt,FG_T_t_History)
        elif runNr==3:
            OpCondInp.writeFGDVCtTHist(FG_TimeTemp4,FG_dt,FG_T_t_History)
        elif runNr==4:
            OpCondInp.writeFGDVCtTHist(FG_TimeTemp5,FG_dt,FG_T_t_History)
        FGDVC.set7File(FG_T_t_History)
        FGDVC.set9AshMoisture(0.0,0.0)
        FGDVC.setTRamp_or_TFile('File') #case: models temperature history with the file
        #writes the instruct.ini and launches FG-DVC (no graphical user interface, only main file fgdvcd.exe)
        FGDVC.writeInstructFile(FG_MainDir+'\\'+FG_ExeCoalDir+'\\')
        FGDVC.Run('cd '+FG_MainDir+FG_ExeCoalDir+' & '+'fgdvcd.exe')
        #
        ###calibrate kinetic parameter:
        #read result:
        CurrentFGFile=FGDVC_Result.FGDVC_Result(FG_DirOut)
        # creates object, required for fitting procedures
        CurrentFGFit=FitInfo.Fit_one_run(CurrentFGFile)
        FGFile.append(CurrentFGFile)
        FGFit.append(CurrentFGFit)
        #copies file, keeping the name:
        shutil.copyfile(FG_DirOut+'gasyield.txt', 'gasyield_'+str(runNr)+'.txt')
        shutil.copyfile(FG_DirOut+'gasrate.txt', 'gasrate_'+str(runNr)+'.txt')        
    #####
    MakeResults('FGDVC',FGFile,FGFit)
        
