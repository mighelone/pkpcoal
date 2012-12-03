import CPD_Fit_lin_regr        #writes CPD-instruct File, launches CPD
import FGDVC_Fit_lin_regr      #writes FG-DVC-instruct File, launches FG-DVC and fittes using eq. (68 ) (BachelorThesis)
import Fit_one_run             #fittes the kinetic parameter for CPD output using eq. (46) (BachelorThesis)
import Compos_and_Energy       #Species balance and energy balance for CPD and FG-DVC
import ReadInputFiles          #reads the user input files, writes FG-DVC coalsd.exe coal generation file
import os
import numpy as np
import platform
#
def DAF(PAFC_asRecieved,PAVM_asRecieved):
    """calculates PAFC, PAVM  from the as recieved state to the daf state of coal"""
    fractionFC=PAFC_asRecieved/(PAFC_asRecieved+PAVM_asRecieved)
    fractionVM=PAVM_asRecieved/(PAFC_asRecieved+PAVM_asRecieved)
    return 100.*fractionFC, 100.*fractionVM
#
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
#get parameters from input files:
#
#Coal File:
CoalInput=ReadInputFiles.ReadFile(workingDir+'Coal.inp')
PAFC_asrec=CoalInput.getValue('Fixed Carbon:')
PAVM_asrec=CoalInput.getValue('Volatile Matter:')
#gets daf values, as CPD needs daf as input:
PAFC_daf, PAVM_daf = DAF(PAFC_asrec,PAVM_asrec)
UAC=CoalInput.getValue('UA Carbon:')
UAH=CoalInput.getValue('UA Hydrogen:')
UAN=CoalInput.getValue('UA Nitrogen:')
UAO=CoalInput.getValue('UA Oxygen:')
HHV=CoalInput.getValue('Higher Heating Value, as recieved, in J/kg:')
MTar=CoalInput.getValue('Tar Molecule weight, MTar:')
WeightY=CoalInput.getValue('Weight-Parameter yields for fitting the kinetics:')
WeightR=CoalInput.getValue('Weight-Parameter rates for fitting the kinetics:')
#
#CPD Properties:
CPDInput=ReadInputFiles.ReadFile(workingDir+'CPD.inp')
CPDselect=CPDInput.UsePyrolProgr('useCPD?:')
CPD_FittingKineticParameter_Select=CPDInput.Fitting("selected fitting Approximation: 'constantRate', 'Arrhenius', 'Kobayashi' or 'None'; selectedFit:")
CPDdt=[0,1,2] #0:initila dt, 1: print increment, 2: dt max
CPDdt[0]=(CPDInput.getValue('initial time step in s:'))
CPDdt[1]=(CPDInput.getValue('print increment, writeValue:'))
#
#FG-DVC Properties:
FGDVCInput=ReadInputFiles.ReadFile(workingDir+'FGDVC.inp')
FG_select=FGDVCInput.UsePyrolProgr('use FG-DVC?:')
FG_FittingKineticParameter_Select=FGDVCInput.Fitting("selected fitting Approximation: 'constantRate', 'Arrhenius', 'Kobayashi' or 'None'; selectedFit:")
FG_CoalSelection=int(FGDVCInput.getValue('Choose Coal: 0 interpolate between library coals and generate own coal. Set 1 to 8 for a library coal.'))
FG_MainDir=FGDVCInput.getText('main directory FG-DVC:')
FG_DirOut=FGDVCInput.getText('directory fgdvc-output:')
FG_TarCacking=FGDVCInput.getValue('Model tar cracking? If no, set tar residence time equal 0. For a partial tar cracking enter the tar residence time in s. For full tar cracking write -1.')
#
#Operating Condition File:
OpCondInp=ReadInputFiles.OperCondInput('OperCond.inp')
CPD_pressure=OpCondInp.getValue('pressure in atm:')
FG_pressure=OpCondInp.getValue('pressure in atm:')
#FG_Tstart=FGDVCInput.getValue('The starting temperature in K:')
#FG_Tend=FGDVCInput.getValue('The final pyrolysis temperature in K:')
#FG_HeatingRate=FGDVCInput.getValue('The Heating rate in K/s:')
CPD_TimeTemp=OpCondInp.getTimePoints('Time History: first column time in seconds, second column: Temperature in K. Last point must contain the final time.','End Time History')
CPDdt[2]=OpCondInp.getValue('FG-DVC: constant (numerical) time step; CPD: maximum time step')
FG_dt=OpCondInp.getValue('FG-DVC: constant (numerical) time step; CPD: maximum time step')
FG_T_t_History=FG_MainDir+'tTHistory.txt'
#makes for CPD time in milliseconds:
CPD_TimeTemp[:,0]=CPD_TimeTemp[:,0]*1.e3
CPD_t_max=CPD_TimeTemp[-1,0]*1.e-3 #tmax in s, not ms
#
#
#
#
####CPD####
if CPDselect==True:
    #launches CPD
    CPD=CPD_Fit_lin_regr.SetterAndLauncher()
    CPD.SetCoalParameter(UAC,UAH,UAN,UAO,PAVM_daf)
    CPD.CalcCoalParam()
    CPD.SetOperateCond(CPD_pressure,CPD_TimeTemp)
    CPD.SetNumericalParam(CPDdt,CPD_t_max)
    CPD.writeInstructFile(workingDir)
    if oSystem=='Linux':
        CPD.Run('./'+'CPD.out','IN.dat')   #first Arg: CPD-executeable, second: Input data containing CPD input file and the output files
    elif oSystem=='Windows':
        CPD.Run('CPD.out','IN.dat')   #first Arg: CPD-executeable, second: Input data containing CPD input file and the output files
    else:
        print "The name of the opearting system couldn't be found."
    #
    ###calibration of the kinetic parameter:
    #read result:
    CPDFile=Fit_one_run.CPD_Result(workingDir)
    # creates object, required for fitting procedures
    CPDFit=Fit_one_run.Fit_one_run(CPDFile)
    #Array for the comparison of the sum of the individual species(t) with the (1-solid(t))
    SumSingleYieldsCalc=np.zeros(len(CPDFit.Yield('Time'))) #Array initialized for the sum of the single yields (calculated)
    SumSingleYieldsCPD=np.zeros(len(CPDFit.Yield('Time')))  #Array initialized for the sum of the single yields (CPD output)
    SolidYieldsCalc=np.zeros(len(CPDFit.Yield('Time')))     #Array initialized for the yields of the Solids (calculated)
    SolidYieldsCPD=np.zeros(len(CPDFit.Yield('Time')))     #Array initialized for the yields of the Solids (CPD output)
    ##CONSTANT RATE
    if CPD_FittingKineticParameter_Select=='constantRate': #CR means ConstantRate
        PredictionVector=[50,0.01] #first argument is k, second is t_start
        LSCPD=Fit_one_run.LeastSquarsEstimator() #LS means LeastSquares
        LSCPD.setOptimizer('fmin')
        LSCPD.setTolerance(1.e-18)
        LSCPD.setWeights(WeightY,WeightR)
        CRCPD=Fit_one_run.ConstantRateModel(PredictionVector)
        outfile = open('CPD-Results_const_rate.txt', 'w')
        outfile.write("Species\t\tk [1/s]\t\tt_start [s]\n\n")
        for Spec in range(2,len(CPDFit.SpeciesNames()),1):
            CRCPD.setParamVector(PredictionVector)
            CRCPD.setParamVector(LSCPD.estimate_T(CPDFit,CRCPD,PredictionVector,Spec))
            CRCPD.plot(CPDFit,Spec)
            Solution=CRCPD.ParamVector()
            outfile.write(str(CPDFit.SpeciesName(Spec))+'\t'+str(Solution[0])+'\t'+str(Solution[1])+'\n')
            #for the comparison of the species sum with (1-Solid)
            if CPDFit.SpeciesName(Spec)=='Solid':
                SolidYieldsCalc+=CRCPD.calcMass(CPDFit,CPDFit.Time(),CPDFit.Interpolate('Temp'),Spec)
                SolidYieldsCPD+=CPDFit.Yield(Spec)
            elif CPDFit.SpeciesName(Spec)!='Solid' and CPDFit.SpeciesName(Spec)!='Temp' and CPDFit.SpeciesName(Spec)!='Time' and CPDFit.SpeciesName(Spec)!='Gas' and CPDFit.SpeciesName(Spec)!='Total':
                SumSingleYieldsCalc+=CRCPD.calcMass(CPDFit,CPDFit.Time(),CPDFit.Interpolate('Temp'),Spec)
                SumSingleYieldsCPD+=CPDFit.Yield(Spec)
        outfile.close()
        CPDFit.plt_InputVectors(CPDFit.Time(),1.-SolidYieldsCalc,1.-SolidYieldsCPD,SumSingleYieldsCalc,SumSingleYieldsCPD,'1-Solid; fitted','1-Solid; CPD output','Sum Yields; fitted','Sum Yields; CPD output')
    ##ARRHENIUS RATE
    if CPD_FittingKineticParameter_Select=='Arrhenius': #Arr means Arrhenius
        PredictionV0=[0.86e15,0,27700]  #for Standard Arrhenius
        PredictionV1=[10.,-20.]         #for Arrhenius notation #1
        PredictionV2=[10,-18]           #for Arrhenius notation #2
        LSCPD=Fit_one_run.LeastSquarsEstimator()
        LSCPD.setOptimizer('fmin')#('leastsq')   # 'leastsq' often faster, but if this does not work: 'fmin' is more reliable
        LSCPD.setTolerance(1.e-10)
        LSCPD.setWeights(WeightY,WeightR)
        outfile = open('CPD-Results_ArrheniusRate.txt', 'w')
        outfile.write("Species\t\tA [1/s]\t\tb\t\tE_a [K]\n\n")
        #select one of the follwoing notations: 
        #Arr=Fit_one_run.ArrheniusModel(PredictionV0)
        #Arr=Fit_one_run.ArrheniusModelAlternativeNotation1(PredictionV1)
        ArrCPD=Fit_one_run.ArrheniusModelAlternativeNotation2(CPDFit,PredictionV2)
        #######
        #uses a separate Arrhenius model to plot, to ensure that the result converted into standart notation (!) output vector is right
        ArrPCPD=Fit_one_run.ArrheniusModel([0,0,0])
        ##The single species:
        for Species in range(2,len(CPDFit.SpeciesNames()),1):
            print CPDFit.SpeciesName(Species)
            ArrCPD.setParamVector(LSCPD.estimate_T(CPDFit,ArrCPD,ArrCPD.ParamVector(),Species))
            Solution=ArrCPD.ConvertKinFactors(ArrCPD.ParamVector())
            #Solution=Arr.ParamVector()
            ArrPCPD.setParamVector(Solution)
            ArrPCPD.plot(CPDFit,Species)
            outfile.write(str(CPDFit.SpeciesName(Species))+'\t'+str(Solution[0])+'\t'+str(Solution[1])+'\t'+str(Solution[2])+'\n')
            #for the comparison of the species sum with (1-Solid)
            if CPDFit.SpeciesName(Species)=='Solid':
                SolidYieldsCalc+=ArrPCPD.calcMass(CPDFit,CPDFit.Time(),CPDFit.Interpolate('Temp'),Species)
                SolidYieldsCPD+=CPDFit.Yield(Species)
            elif CPDFit.SpeciesName(Species)!='Solid' and CPDFit.SpeciesName(Species)!='Temp' and CPDFit.SpeciesName(Species)!='Time' and CPDFit.SpeciesName(Species)!='Gas' and CPDFit.SpeciesName(Species)!='Total':
                SumSingleYieldsCalc+=ArrPCPD.calcMass(CPDFit,CPDFit.Time(),CPDFit.Interpolate('Temp'),Species)
                SumSingleYieldsCPD+=CPDFit.Yield(Species)
        outfile.close()
        CPDFit.plt_InputVectors(CPDFit.Time(),1.-SolidYieldsCalc,1.-SolidYieldsCPD,SumSingleYieldsCalc,SumSingleYieldsCPD,'1-Solid; fitted','1-Solid; CPD output','Sum Yields; fitted','Sum Yields; CPD output')
    ##KOBAYASHI RATE
    if CPD_FittingKineticParameter_Select=='Kobayashi': #Kob means Kobayashi
        PredictionVKob2=[10,-16,8,-20,0.5,1.0]           #for Arrhenius notation #2 [b11,b21,b12,b22] with the second indice as the reaction
        LSCPD=Fit_one_run.LeastSquarsEstimator()
        LSCPD.setOptimizer('fmin')#('leastsq')   # 'leastsq' often faster, but if this does not work: 'fmin' is more reliable
        LSCPD.setTolerance(1.e-7)
        LSCPD.setWeights(WeightY,WeightR)
        outfile = open('CPD-Results_KobayashiRate.txt', 'w')
        outfile.write("Species\t\t\tA1 [1/s]\t\t\t\tE_a1 [K]\t\tA2 [1/s]\t\t\t\tE_a2 [K]\t\t\t\talpha1 \t\t\talpha2 \n\n")
        KobCPD=Fit_one_run.KobayashiA2(CPDFit,PredictionVKob2)
        #######
        #uses a separate Arrhenius model to plot, to ensure that the result converted into standart notation (!) output vector is right
        KobPCPD=Fit_one_run.Kobayashi(CPDFit,[0,0,0,0,0,0])
        ##The single species:
        for Species in range(2,len(CPDFit.SpeciesNames()),1):
            print CPDFit.SpeciesName(Species)
            KobCPD.setParamVector(LSCPD.estimate_T(CPDFit,KobCPD,KobCPD.ParamVector(),Species))
            Solution=KobCPD.ConvertKinFactors(KobCPD.ParamVector())
            #Solution=Arr.ParamVector()
            KobPCPD.setParamVector(Solution)
            KobPCPD.plot(CPDFit,Species)
            outfile.write(str(CPDFit.SpeciesName(Species))+'\t\t'+str(Solution[0])+'\t\t'+str(Solution[1])+'\t\t'+str(Solution[2])+'\t\t'+str(Solution[3])+'\t\t'+str(Solution[4])+'\t\t'+str(Solution[5])+'\n')
            #for the comparison of the species sum with (1-Solid)
            if CPDFit.SpeciesName(Species)=='Solid':
                SolidYieldsCalc+=KobPCPD.calcMass(CPDFit,CPDFit.Time(),CPDFit.Interpolate('Temp'),Species)
                SolidYieldsCPD+=CPDFit.Yield(Species)
            elif CPDFit.SpeciesName(Species)!='Solid' and CPDFit.SpeciesName(Species)!='Temp' and CPDFit.SpeciesName(Species)!='Time' and CPDFit.SpeciesName(Species)!='Gas' and CPDFit.SpeciesName(Species)!='Total':
                SumSingleYieldsCalc+=KobPCPD.calcMass(CPDFit,CPDFit.Time(),CPDFit.Interpolate('Temp'),Species)
                SumSingleYieldsCPD+=CPDFit.Yield(Species)
        outfile.close()
        CPDFit.plt_InputVectors(CPDFit.Time(),1.-SolidYieldsCalc,1.-SolidYieldsCPD,SumSingleYieldsCalc,SumSingleYieldsCPD,'1-Solid; fitted','1-Solid; CPD output','Sum Yields; fitted','Sum Yields; CPD output')
    ##SPECIES AND ENERGY BALANCE:
    SpecCPD=Compos_and_Energy.CPD_SpeciesBalance(CPDFile,UAC,UAH,UAN,UAO,PAVM_asrec,PAFC_asrec,HHV,MTar)
#
#
#
####FG-DVC####
if FG_select==True:
    #writes Time-Temperature file
    FG_TimeTemp=CPD_TimeTemp
    FG_TimeTemp[:,0]=CPD_TimeTemp[:,0]*1.e-3
    OpCondInp.writeFGDVCtTHist(FG_TimeTemp,FG_dt,FG_T_t_History)
    #initialize the launching object
    FGDVC=FGDVC_Fit_lin_regr.SetterAndLauncher()
    #set and writes Coal Files:
    if FG_CoalSelection==0:
        #deletes old generated file
        os.system('cd '+FG_MainDir+FG_GenCoalDir+' & del '+FG_CoalName+'_com.dat, '+FG_CoalName+'_kin.dat, '+FG_CoalName+'_pol.dat')
        #generates coalsd.exe input file
        MakeCoalGenFile=ReadInputFiles.WriteFGDVCCoalFile(FG_CoalGenFileName)
        MakeCoalGenFile.setCoalComp(UAC,UAH,UAO,UAN,(100.-UAC-UAH-UAO-UAN),0)
        MakeCoalGenFile.write(FG_MainDir+FG_GenCoalDir+'\\',FG_CoalName)
        #makes new file
        os.system('cd '+FG_MainDir+FG_GenCoalDir+' & '+'coalsd.exe < '+FG_CoalGenFileName)
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
    FGDVC.set7File(FG_T_t_History)
    FGDVC.set9AshMoisture(0.0,0.0)
    FGDVC.setTRamp_or_TFile('File') #case: models temperature history with the file
    #writes the instruct.ini and launches FG-DVC (no graphical user interface, only main file fgdvcd.exe)
    FGDVC.writeInstructFile(FG_MainDir+'\\'+FG_ExeCoalDir+'\\')
    FGDVC.Run('cd '+FG_MainDir+FG_ExeCoalDir+' & '+'fgdvcd.exe')
    #
    ###calibrate kinetic parameter:
    #read result:
    FGFile=Fit_one_run.FGDVC_Result(FG_DirOut)
    # creates object, required for fitting procedures
    FGFit=Fit_one_run.Fit_one_run(FGFile)
    #Array for the comparison of the sum of the individual species(t) with the (1-solid(t))
    SumSingleYieldsCalc=np.zeros(len(FGFit.Yield('Time'))) #Array initialized for the sum of the single yields (calculated)
    SumSingleYieldsFG=np.zeros(len(FGFit.Yield('Time')))   #Array initialized for the sum of the single yields (CPD output)
    SolidYieldsCalc=np.zeros(len(FGFit.Yield('Time')))     #Array initialized for the yields of the Solids (calculated)
    SolidYieldsFG=np.zeros(len(FGFit.Yield('Time')))       #Array initialized for the yields of the Solids (CPD output)
    ##CONSTANT RATE
    if FG_FittingKineticParameter_Select=='constantRate':
        PredictionVector=[50,0.01] #first argument is k, second is t_start
        LSFG=Fit_one_run.LeastSquarsEstimator()
        LSFG.setOptimizer('fmin')
        LSFG.setTolerance(1.e-18)
        LSFG.setWeights(WeightY,WeightR)
        CRFG=Fit_one_run.ConstantRateModel(PredictionVector)
        outfile = open('FGDVC-Results_const_rate.txt', 'w')
        outfile.write("Species\t\tk [1/s]\t\tt_start [s]\n\n")
        #
        for Spec in range(2,len(FGFit.SpeciesNames()),1):
            CRFG.setParamVector(PredictionVector)
            CRFG.setParamVector(LSFG.estimate_T(FGFit,CRFG,PredictionVector,Spec))
            CRFG.plot(FGFit,Spec)
            Solution=CRFG.ParamVector()
            outfile.write(str(FGFit.SpeciesName(Spec))+'\t'+str(Solution[0])+'\t'+str(Solution[1])+'\n')
            #for the comparison of the species sum with (1-Solid)
            if FGFit.SpeciesName(Spec)=='Solid':
                SolidYieldsCalc+=CRFG.calcMass(FGFit,FGFit.Time(),FGFit.Interpolate('Temp'),Spec)
                SolidYieldsFG+=FGFit.Yield(Spec)
            elif FGFit.SpeciesName(Spec)!='Solid' and FGFit.SpeciesName(Spec)!='Temp' and FGFit.SpeciesName(Spec)!='Time':
                SumSingleYieldsCalc+=CRFG.calcMass(FGFit,FGFit.Time(),FGFit.Interpolate('Temp'),Spec)
                SumSingleYieldsFG+=FGFit.Yield(Spec)
        outfile.close()
        FGFit.plt_InputVectors(FGFit.Time(),100.-SolidYieldsCalc,100.-SolidYieldsFG,SumSingleYieldsCalc,SumSingleYieldsFG,'1-Solid; fitted','1-Solid; FG-DVC output','Sum Yields; fitted','Sum Yields; FG-DVC output')
    ##ARRHENIUS RATE
    if FG_FittingKineticParameter_Select=='Arrhenius':
        PredictionV0=[0.86e15,0,27700]  #for Standard Arrhenius
        PredictionV1=[10.,-20.]         #for Arrhenius notation #1
        PredictionV2=[10,-18]           #for Arrhenius notation #2
        LSFG=Fit_one_run.LeastSquarsEstimator()
        LSFG.setOptimizer('fmin')#('leastsq')   # 'leastsq' often faster, but if this does not work: 'fmin' is more reliable
        LSFG.setTolerance(1.e-10)
        LSFG.setWeights(WeightY,WeightR)
        outfile = open('FGDVC-Results_ArrheniusRate.txt', 'w')
        outfile.write("Species\t\tA [1/s]\t\tb\t\tE_a [K]\n\n")
        #select one of the follwoing notations: 
        #Arr=Fit_one_run.ArrheniusModel(PredictionV0)
        #Arr=Fit_one_run.ArrheniusModelAlternativeNotation1(PredictionV1)
        ArrFG=Fit_one_run.ArrheniusModelAlternativeNotation2(FGFit,PredictionV2)
        #######
        #uses a separate Arrhenius model to plot, to ensure that the converted (!) output vector is right
        ArrPFG=Fit_one_run.ArrheniusModel([0,0,0])
        ##The single species:
        for Species in range(2,len(FGFit.SpeciesNames()),1):
            print FGFit.SpeciesName(Species)
            ArrFG.setParamVector(LSFG.estimate_T(FGFit,ArrFG,ArrFG.ParamVector(),Species))
            Solution=ArrFG.ConvertKinFactors(ArrFG.ParamVector())
            #Solution=Arr.ParamVector()
            ArrPFG.setParamVector(Solution)
            ArrPFG.plot(FGFit,Species)
            outfile.write(str(FGFit.SpeciesName(Species))+'\t'+str(Solution[0])+'\t'+str(Solution[1])+'\t'+str(Solution[2])+'\n')
            #for the comparison of the species sum with (1-Solid)
            if FGFit.SpeciesName(Species)=='Solid':
                SolidYieldsCalc+=ArrPFG.calcMass(FGFit,FGFit.Time(),FGFit.Interpolate('Temp'),Species)
                SolidYieldsFG+=FGFit.Yield(Species)
            elif FGFit.SpeciesName(Species)!='Solid' and FGFit.SpeciesName(Species)!='Temp' and FGFit.SpeciesName(Species)!='Time':
                SumSingleYieldsCalc+=ArrPFG.calcMass(FGFit,FGFit.Time(),FGFit.Interpolate('Temp'),Species)
                SumSingleYieldsFG+=FGFit.Yield(Species)
        outfile.close()
        FGFit.plt_InputVectors(FGFit.Time(),100.-SolidYieldsCalc,100.-SolidYieldsFG,SumSingleYieldsCalc,SumSingleYieldsFG,'1-Solid; fitted','1-Solid; FG-DVC output','Sum Yields; fitted','Sum Yields; FG-DVC output')
    ##KOBAYASHI RATE
    if FG_FittingKineticParameter_Select=='Kobayashi': #Kob means Kobayashi
        PredictionVKob2=[2e5,1.046e8/8314.33,1.3e7,1.674e8/8314.33,0,0]           #for Arrhenius notation #2 [b11,b21,b12,b22] with the second indice as the reaction
        LSFG=Fit_one_run.LeastSquarsEstimator()
        LSFG.setOptimizer('fmin')#('leastsq')   # 'leastsq' often faster, but if this does not work: 'fmin' is more reliable
        LSFG.setTolerance(1.e-7)
        LSFG.setWeights(WeightY,WeightR)
        outfile = open('FGDVC-Results_KobayashiRate.txt', 'w')
        outfile.write("Species\t\t\tA1 [1/s]\t\t\t\tE_a1 [K]\t\tA2 [1/s]\t\t\t\tE_a2 [K]\t\t\t\talpha1 \t\t\talpha2 \n\n")
        KobFG=Fit_one_run.KobayashiA2(FGFit,PredictionVKob2)
        #######
        #uses a separate Arrhenius model to plot, to ensure that the result converted into standart notation (!) output vector is right
        KobPFG=Fit_one_run.Kobayashi(FGFit,[2e5,1.046e8/8314.33,1.3e7,1.674e8/8314.33,0,0])
        ##The single species:
        for Species in range(2,len(FGFit.SpeciesNames()),1):
            print FGFit.SpeciesName(Species)
            KobFG.setParamVector(LSFG.estimate_T(FGFit,KobFG,KobFG.ParamVector(),Species))
            Solution=KobFG.ConvertKinFactors(KobFG.ParamVector())
            #Solution=Arr.ParamVector()
            KobPFG.setParamVector(Solution)
            KobPFG.plot(FGFit,Species)
            outfile.write(str(CPDFit.SpeciesName(Species))+'\t\t'+str(Solution[0])+'\t\t'+str(Solution[1])+'\t\t'+str(Solution[2])+'\t\t'+str(Solution[3])+'\t\t'+str(Solution[4])+'\t\t'+str(Solution[5])+'\n')
            #for the comparison of the species sum with (1-Solid)
            if FGFit.SpeciesName(Species)=='Solid':
                SolidYieldsCalc+=KobPFG.calcMass(FGFit,FGFit.Time(),FGFit.Interpolate('Temp'),Species)
                SolidYieldsFG+=FGFit.Yield(Species)
            elif FGFit.SpeciesName(Species)!='Solid' and FGFit.SpeciesName(Species)!='Temp' and FGFit.SpeciesName(Species)!='Time':
                SumSingleYieldsCalc+=KobPFG.calcMass(FGFit,FGFit.Time(),FGFit.Interpolate('Temp'),Species)
                SumSingleYieldsFG+=FGFit.Yield(Species)
        outfile.close()
        FGFit.plt_InputVectors(FGFit.Time(),1.-SolidYieldsCalc,1.-SolidYieldsFG,SumSingleYieldsCalc,SumSingleYieldsFG,'1-Solid; fitted','1-Solid; FG-DVC output','Sum Yields; fitted','Sum Yields; FG-DVC output')
    #############################
    ##SPECIES AND ENERGY BALANCE:
    SpecFG=Compos_and_Energy.FGDVC_SpeciesBalance(FGFile,UAC,UAH,UAN,UAO,PAVM_asrec,PAFC_asrec,HHV,MTar)
#
