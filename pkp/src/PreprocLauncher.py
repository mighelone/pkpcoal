"""
A selection of launcher functions for different coal preprocessors

Launcher functions are  called with the inputs dictionary to pass 
all definitions from the input/inputs file

The launcher function return a list of results objects
"""

from Models import BalancedComposition
import CPD

def Launch_CPD(inputs):
    """ Execute CPD for each given temperature profile and 
        return a list of CPD results objects
     """
    def InitAndLaunch(*pargs):
        """ initialises and execute cpd calculation """
        print 'Running CPD: ' + str(run)
        cpd = CPD.CPD(*pargs)
        return cpd.Run() 

    operatingConditions = inputs['OperatingConditions']
    #pressure = operatingConditions.pop('pressure') 
    pressure = operatingConditions['pressure']
    coal = inputs["Coal"]
    proxi_ana = BalancedComposition(coal["Proximate Analysis"])
    ultim_ana = BalancedComposition(coal["Ultimate Analysis"])
    daf = proxi_ana.remove_elems_rebalance(['Moisture','Ash'])
    return [InitAndLaunch(
                ultim_ana, 
                daf,
                operatingConditions["run"+str(run+1)],
                pressure,
                inputs['CPD']['deltaT'],
                run
            ) for run in range(operatingConditions['runs'])]



# def MakeResults_FG(self):
#     """generates the result for FG-DVC"""
#
#     # def CheckFGdt(self):
#     #     """ Aborts, if FG-DVC is selected and the timestep is lower than 1.e-3 
#     #     """
#     #     # TODO GO move to FGDVC class
#     #     # TODO check if this function gets callesd ? Can't we use a limiter here ?
#     #     if ((self.FG_select==True) and (self.FG_dt<1e-4)):
#     #         print """Please select for FG-DVC a time step greather equal 1e-4 in 'OperCond.inp'. 
#     #         FG-DVC would not be able to read the time history file for a dt<1e-4. """
#     #         sys.exit()
#     #writes Time-Temperature file
#
#     FG_TimeTemp1=np.zeros(np.shape(self.CPD_TimeTemp1),order='F')
#     FG_TimeTemp2=np.zeros(np.shape(self.CPD_TimeTemp2),order='F')
#     FG_TimeTemp3=np.zeros(np.shape(self.CPD_TimeTemp3),order='F')
#     FG_TimeTemp4=np.zeros(np.shape(self.CPD_TimeTemp4),order='F')
#     FG_TimeTemp5=np.zeros(np.shape(self.CPD_TimeTemp5),order='F')
#     FG_TimeTemp1[:,0]=self.CPD_TimeTemp1[:,0]*1.e-3
#     FG_TimeTemp2[:,0]=self.CPD_TimeTemp2[:,0]*1.e-3
#     FG_TimeTemp3[:,0]=self.CPD_TimeTemp3[:,0]*1.e-3
#     FG_TimeTemp4[:,0]=self.CPD_TimeTemp4[:,0]*1.e-3
#     FG_TimeTemp5[:,0]=self.CPD_TimeTemp5[:,0]*1.e-3
#     FG_TimeTemp1[:,1]=self.CPD_TimeTemp1[:,1]
#     FG_TimeTemp2[:,1]=self.CPD_TimeTemp2[:,1]
#     FG_TimeTemp3[:,1]=self.CPD_TimeTemp3[:,1]
#     FG_TimeTemp4[:,1]=self.CPD_TimeTemp4[:,1]
#     FG_TimeTemp5[:,1]=self.CPD_TimeTemp5[:,1]
#     #initialize the launching object
#     FGDVC=FGDVC_SetAndLaunch.SetterAndLauncher()
#     #set and writes Coal Files:
#     if self.FG_CoalSelection==0:
#         #deletes old generated file
#         os.system('cd '+self.FG_MainDir+FG_GenCoalDir+' & del '+FG_CoalName+'_com.dat, '+FG_CoalName+'_kin.dat, '+FG_CoalName+'_pol.dat')
#         #generates coalsd.exe input file
#         MakeCoalGenFile=InformationFiles.WriteFGDVCCoalFile(FG_CoalGenFileName)
#         MakeCoalGenFile.setCoalComp(self.UAC,self.UAH,self.UAO,self.UAN,self.UAS,0)
#         MakeCoalGenFile.write(self.FG_MainDir+FG_GenCoalDir+'\\',FG_CoalName,option=0)
#         #makes new file
#         try:
#             os.system('cd '+self.FG_MainDir+FG_GenCoalDir+' & '+'coalsd.exe < '+FG_CoalGenFileName+' > coalsd_pkp.log')
#         except OSError:
#             print 'Problems with coalsd.exe'
#         os.system('copy '+self.FG_MainDir+FG_GenCoalDir+'\coalsd_pkp.log . >> log.txt')
#         #tests weather the coal file was genearated:
#         if os.path.exists(self.FG_MainDir+'\\'+FG_GenCoalDir+'\\'+FG_CoalName+'_com.dat')==False:
#             print 30*'*','\n','The coal is may outside the libraries coals. Select manually the closest library coal.',30*'*','\n'
#             MakeCoalGenFile.write(self.FG_MainDir+FG_GenCoalDir+'\\',FG_CoalName,option=10)
#             os.system('cd '+self.FG_MainDir+FG_GenCoalDir+' & '+'coalsd.exe < '+FG_CoalGenFileName+' > coalsd_pkp.log')
#         #sets generated file for instruct.ini
#         FGDVC.set1CoalLocation(self.FG_MainDir+FG_GenCoalDir+'\\'+FG_CoalName+'_com.dat')
#         FGDVC.set2KinLocation(self.FG_MainDir+FG_GenCoalDir+'\\'+FG_CoalName+'_kin.dat')
#         FGDVC.set3PolyLocation(self.FG_MainDir+FG_GenCoalDir+'\\'+FG_CoalName+'_pol.dat')
#     elif self.FG_CoalSelection>0 and self.FG_CoalSelection<9:
#         #sets library file for instruct.ini
#         FGDVC.set1CoalLocation(self.FG_MainDir+FG_LibCoalDir+'\\coal.ar'+str(self.FG_CoalSelection))
#         FGDVC.set2KinLocation(self.FG_MainDir+FG_LibCoalDir+'\\kin.ar'+str(self.FG_CoalSelection))
#         FGDVC.set3PolyLocation(self.FG_MainDir+FG_LibCoalDir+'\\polymr.ar'+str(self.FG_CoalSelection))
#     else:
#         print "select Choose Coal: 0 interpolate between library coals and generate own coal. Set 1 to 8 for a library coal.' in FGDVC.inp equal a value between 0 and 8"
#     #sets FG-DVC instruct.ini parameter
#     FGDVC.set5Pressure(self.FG_pressure)
#     if self.FG_TarCacking==0.0:            #case: no tar cracking
#         FGDVC.set6Theorie(13,0.0)
#     elif self.FG_TarCacking<0.0:           #case: full tar cracking
#         FGDVC.set6Theorie(15,0.0)
#     else:                             #case: partial tar cracking
#         FGDVC.set6Theorie(13,float(self.FG_TarCacking))
#     #
#     FGFile=[]
#     FGFit=[]
#     OpCondInp=InformationFiles.OperCondInput('OperCond.inp')
#     for runNr in range(self.NrOfRuns):
#         if runNr==0:
#             OpCondInp.writeFGDVCtTHist(FG_TimeTemp1,self.FG_dt,self.FG_T_t_History)
#         elif runNr==1:
#             OpCondInp.writeFGDVCtTHist(FG_TimeTemp2,self.FG_dt,self.FG_T_t_History)
#         elif runNr==2:
#             OpCondInp.writeFGDVCtTHist(FG_TimeTemp3,self.FG_dt,self.FG_T_t_History)
#         elif runNr==3:
#             OpCondInp.writeFGDVCtTHist(FG_TimeTemp4,self.FG_dt,self.FG_T_t_History)
#         elif runNr==4:
#             OpCondInp.writeFGDVCtTHist(FG_TimeTemp5,self.FG_dt,self.FG_T_t_History)
#         FGDVC.set7File(self.FG_T_t_History)
#         FGDVC.set9AshMoisture(0.0,0.0)
#         FGDVC.setTRamp_or_TFile('File') #case: models temperature history with the file
#         #writes the instruct.ini and launches FG-DVC (no graphical user interface, only main file fgdvcd.exe)
#         FGDVC.writeInstructFile(self.FG_MainDir+'\\'+FG_ExeCoalDir+'\\')
#         FGDVC.Run('cd '+self.FG_MainDir+FG_ExeCoalDir+' & '+'fgdvcd.exe')
#         #
#         ###calibrate kinetic parameter:
#         #read result:
#         CurrentFGFile=FGDVC_Result.FGDVC_Result(self.FG_DirOut)
#         # creates object, required for fitting procedures
#         CurrentFGFit=FitInfo.Fit_one_run(CurrentFGFile)
#         FGFile.append(CurrentFGFile)
#         FGFit.append(CurrentFGFit)
#         #copies file, keeping the name:
#         if oSystem=='Linux' or oSystem == 'Darwin':
#             shutil.copyfile(self.FG_DirOut+'gasyield.txt', 'Result/gasyield_'+str(runNr)+'.txt')
#             shutil.copyfile(self.FG_DirOut+'gasrate.txt', 'Result/gasrate_'+str(runNr)+'.txt')
#         elif oSystem=='Windows':
#             shutil.copyfile(self.FG_DirOut+'gasyield.txt', 'Result\\gasyield_'+str(runNr)+'.txt')
#             shutil.copyfile(self.FG_DirOut+'gasrate.txt', 'Result\\gasrate_'+str(runNr)+'.txt')
#     #####
#     M=Models.Model()
#     for Species in FGFit[0].SpeciesNames():
#         M.mkSimpleResultFiles(FGFit,Species)
#         if (Species not in self.SpeciesToConsider) and (Species!='Temp') and (Species!='Time'):
#             self.SpeciesToConsider.append(Species)
#     if self.FG_FittingKineticParameter_Select=='constantRate':
#         self.MakeResults_CR('FGDVC',FGFile,FGFit)
#         currentDict={'FGDVC':'constantRate'}
#     elif self.FG_FittingKineticParameter_Select=='Arrhenius':
#         self.MakeResults_Arrh('FGDVC',FGFile,FGFit)
#         currentDict={'FGDVC':'Arrhenius'}
#     elif self.FG_FittingKineticParameter_Select=='ArrheniusNoB':
#         self.MakeResults_ArrhNoB('FGDVC',FGFile,FGFit)
#         currentDict={'FGDVC':'ArrheniusNoB'}
#     elif self.FG_FittingKineticParameter_Select=='Kobayashi':
#         self.MakeResults_Kob('FGDVC',FGFile,FGFit)
#         currentDict={'FGDVC':'Kobayashi'}
#     elif self.FG_FittingKineticParameter_Select=='DAEM':
#         self.MakeResults_DEAM('FGDVC',FGFile,FGFit)
#         currentDict={'FGDVC':'DAEM'}
#     elif self.FG_FittingKineticParameter_Select==None:
#         currentDict={'FGDVC':'None'}
#     else:
#         print 'uspecified FG_FittingKineticParameter_Select'
#         currentDict={}
#     #
#     self.ProgramModelDict.update(currentDict)
#     #
#     self.SpeciesEnergy('FGDVC',FGFile,self.FG_FittingKineticParameter_Select)
#         #
#
# ####Pc Coal Lab####
# def MakeResults_PCCL(self):
#     """generates the result for PC Coal Lab"""
#     #writes Time-Temperature file
#     PCCL_TimeTemp1=np.zeros(np.shape(self.CPD_TimeTemp1),order='F')
#     PCCL_TimeTemp2=np.zeros(np.shape(self.CPD_TimeTemp2),order='F')
#     PCCL_TimeTemp3=np.zeros(np.shape(self.CPD_TimeTemp3),order='F')
#     PCCL_TimeTemp4=np.zeros(np.shape(self.CPD_TimeTemp4),order='F')
#     PCCL_TimeTemp5=np.zeros(np.shape(self.CPD_TimeTemp5),order='F')
#     PCCL_TimeTemp1[:,0]=self.CPD_TimeTemp1[:,0]*1.e-3
#     PCCL_TimeTemp2[:,0]=self.CPD_TimeTemp2[:,0]*1.e-3
#     PCCL_TimeTemp3[:,0]=self.CPD_TimeTemp3[:,0]*1.e-3
#     PCCL_TimeTemp4[:,0]=self.CPD_TimeTemp4[:,0]*1.e-3
#     PCCL_TimeTemp5[:,0]=self.CPD_TimeTemp5[:,0]*1.e-3
#     PCCL_TimeTemp1[:,1]=self.CPD_TimeTemp1[:,1]
#     PCCL_TimeTemp2[:,1]=self.CPD_TimeTemp2[:,1]
#     PCCL_TimeTemp3[:,1]=self.CPD_TimeTemp3[:,1]
#     PCCL_TimeTemp4[:,1]=self.CPD_TimeTemp4[:,1]
#     PCCL_TimeTemp5[:,1]=self.CPD_TimeTemp5[:,1]
#     #initialize the launching object
#     PCCL=PCCL_SetAndLaunch.SetterAndLauncher()
#     #set and writes Coal Files:
#     PCCL.SetUACoalParameter(self.UAC,self.UAH,self.UAN,self.UAO,self.UAS)
#     PCCL.SetPACoalParameter(self.PAVM_asrec,self.PAFC_asrec,self.PAmoist,self.PAash)
#     if type(self.PCCL_CoalCalFactor)==float:
#         PCCL.SetCoalCalibrationFactor(self.PCCL_CoalCalFactor)
#     PCCL.SetPressure(self.FG_pressure)
#     PCCL.SetParticleSize(self.PCCL_ParticleSize)
#     #
#     PCCL.writeCoalFiles(self.PCCL_Path)
#     #
#     PCCL.THist(PCCL_TimeTemp1[0,1],PCCL_TimeTemp1[1,0],PCCL_TimeTemp1[-1,1],PCCL_TimeTemp1[-1,0])
#     PCCL.writeInstructFiles(self.PCCL_Path,mkNewFile=True)
#     PCCL.THist(PCCL_TimeTemp2[0,1],PCCL_TimeTemp2[1,0],PCCL_TimeTemp2[-1,1],PCCL_TimeTemp2[-1,0])
#     PCCL.writeInstructFiles(self.PCCL_Path,mkNewFile=False)
#     PCCL.THist(PCCL_TimeTemp3[0,1],PCCL_TimeTemp3[1,0],PCCL_TimeTemp3[-1,1],PCCL_TimeTemp3[-1,0])
#     PCCL.writeInstructFiles(self.PCCL_Path,mkNewFile=False)
#     PCCL.THist(PCCL_TimeTemp4[0,1],PCCL_TimeTemp4[1,0],PCCL_TimeTemp4[-1,1],PCCL_TimeTemp4[-1,0])
#     PCCL.writeInstructFiles(self.PCCL_Path,mkNewFile=False)
#     PCCL.THist(PCCL_TimeTemp5[0,1],PCCL_TimeTemp5[1,0],PCCL_TimeTemp5[-1,1],PCCL_TimeTemp5[-1,0])
#     PCCL.writeInstructFiles(self.PCCL_Path,mkNewFile=False)
#     PCCL.writeInstructFilesFinish()
#     #
#     PCCL.Run(self.PCCL_Path,self.PCCL_Exe)
#     #
#     PCCLFile=[]
#     PCCLFit=[]
#     for runNr in range(1,self.NrOfRuns+1,1):
#         #read result:
#         CurrentPCCLFile=PCCL_Result.PCCL_Result(self.PCCL_Path,runNr)
#         # creates object, required for fitting procedures
#         CurrentPCCLFit=FitInfo.Fit_one_run(CurrentPCCLFile)
#         PCCLFile.append(CurrentPCCLFile)
#         PCCLFit.append(CurrentPCCLFit)
#         #copies file:
#         if oSystem=='Windows':
#             shutil.copyfile(self.PCCL_Path+'FDC1WT'+str(runNr)+'.RPT', 'Result/PCCL_gasyield_wt_'+str(runNr)+'.txt')
#             shutil.copyfile(self.PCCL_Path+'FDC1NG'+str(runNr)+'.RPT', 'Result/PCCL_gasyield_ng_'+str(runNr)+'.txt')
#             shutil.copyfile(self.PCCL_Path+'FDC1HC'+str(runNr)+'.RPT', 'Result/PCCL_gasyield_hc_'+str(runNr)+'.txt')
#     #####
#     M=Models.Model()
#     for Species in PCCLFit[0].SpeciesNames():
#         M.mkSimpleResultFiles(PCCLFit,Species)
#         if (Species not in self.SpeciesToConsider) and (Species!='Temp') and (Species!='Time'):
#             self.SpeciesToConsider.append(Species)
#     if self.PCCL_FittingKineticParameter_Select=='constantRate':
#         self.MakeResults_CR('PCCL',PCCLFile,PCCLFit)
#         currentDict={'PCCL':'constantRate'}
#     elif self.PCCL_FittingKineticParameter_Select=='Arrhenius':
#         self.MakeResults_Arrh('PCCL',PCCLFile,PCCLFit)
#         currentDict={'PCCL':'Arrhenius'}
#     elif self.PCCL_FittingKineticParameter_Select=='ArrheniusNoB':
#         self.MakeResults_ArrhNoB('PCCL',PCCLFile,PCCLFit)
#         currentDict={'PCCL':'ArrheniusNoB'}
#     elif self.PCCL_FittingKineticParameter_Select=='Kobayashi':
#         self.MakeResults_Kob('PCCL',PCCLFile,PCCLFit)
#         currentDict={'PCCL':'Kobayashi'}
#     elif self.PCCL_FittingKineticParameter_Select=='DAEM':
#         self.MakeResults_DEAM('PCCL',PCCLFile,PCCLFit)
#         currentDict={'PCCL':'DAEM'}
#     elif self.PCCL_FittingKineticParameter_Select==None:
#         currentDict={'PCCL':'None'}
#     else:
#         print 'uspecified PCCL_FittingKineticParameter_Select'
#         currentDict={}
#     #
#     self.ProgramModelDict.update(currentDict)
#     #
#     self.SpeciesEnergy('PCCL',PCCLFile)
#
#
#
#
# # def RunPMSKD(self):
# #     '''
# #     run PMSKD
# #     '''
# #     # create object
# #
# #     try:
# #         coal = coalPolimi.coalPolimi(name = 'COAL', c=self.UAC,h=self.UAH,o=self.UAO,n=self.UAN,s=self.UAS,file=self.PMSKD_mechfile)
# #     except coalPolimi.compositionError:
# #         print 'Composition outside of triangle of definition'
# #         sys.exit()
# #     # organize TimeTemp
# #         PMSKDFile=[]
# #     PMSKDFit=[]
# #     for runNr in range(self.NrOfRuns):
# #         print 'Running PMSKD n. '+str(runNr)
# #         #print self.timeHR[runNr]
# #         #print self.temperatureHR[runNr]
# #         #set heating rate
# #         coal.setHeatingRate(self.timeHR[runNr],self.temperatureHR[runNr])
# #         #coal.setTimeStep(self.PMSKD_npoint)
# #         coal.solvePyrolysis()
# #         #plt.figure(runNr)
# #         #plt.plot(coal.getTemperature(),coal.getVolatile())
# #         #read result:
# #         #CurrentPMSKDFile=Coal
# #         # creates object, required for fitting procedures
# #         CurrentPMSKDFit=FitInfo.Fit_one_run(coal)
# #         #PMSKDFile.append(CurrentFGFile)
# #         PMSKDFit.append(CurrentPMSKDFit)
# #         #print coal.Yields_all()
# #
# #         coal.reset()
# #         #print coal.timeHR
# #         #print coal.temperatureHR
# #
# #     if self.PMSKD_FittingKineticParameter_Select=='constantRate':
# #         self.MakeResults_CR('PMSKD','',PMSKDFit)
# #         currentDict={'PMSKD':'constantRate'}
# #     elif self.PMSKD_FittingKineticParameter_Select=='Arrhenius':
# #         self.MakeResults_Arrh('PMSKD','',PMSKDFit)
# #         currentDict={'PMSKD':'Arrhenius'}
# #     elif self.PMSKD_FittingKineticParameter_Select=='ArrheniusNoB':
# #         self.MakeResults_ArrhNoB('PMSKD','',PMSKDFit)
# #         currentDict={'PMSKD':'ArrheniusNoB'}
# #     elif self.PMSKD_FittingKineticParameter_Select=='Kobayashi':
# #         self.MakeResults_Kob('PMSKD','',PMSKDFit)
# #         currentDict={'PMSKD':'Kobayashi'}
# #     elif self.PMSKD_FittingKineticParameter_Select=='DAEM':
# #         self.MakeResults_DEAM('PMSKD','',PMSKDFit)
# #         currentDict={'PMSKD':'DAEM'}
# #     elif self.PMSKD_FittingKineticParameter_Select==None:
# #         currentDict={'PMSKD':'None'}
# #         for Species in PMSKDFit[0].SpeciesNames():
# #             M=Models.Model()
# #             M.mkSimpleResultFiles(PMSKDFit,Species)
# #             if ((Species not in self.SpeciesToConsider) 
# #                 and (Species!='Temp') 
# #                 and (Species!='Time')):
# #                 self.SpeciesToConsider.append(Species)
# #     else:
# #         print 'undefined PMSKD_FittingKineticParameter_Select'
# #         currentDict={}
# #         #
# #     self.ProgramModelDict.update(currentDict)
# #     #
# #     #self.SpeciesEnergy('PMSKD',FGFile)
