import numpy as np
import os
from scipy import stats
import pylab as plt
import shutil
from scipy.integrate import odeint
import scipy.interpolate
import Fit_one_run
################################
R=1.0 #8.3144621 # Gas constant only =8.3... if E should not include R
################################

class SetterAndLauncher(object):
    """This class is able to write the 'instruct.ini' and launch the 'fgdvcd.exe'. Before writing the instruct.ini (method 'writeInstructFile') set all parameter using the corresponding methods (at least necessary: set1CoalLocation, set2KinLocation, set3PolyLocation, set5Pressure, set7Ramp). After writing the instruct file, the .Run method can be used."""
    def __init__(self):
        self.Arguments1={'L1_Coalfile':'no value setted','L2_Kinfile':'no value setted','L3_Polyfile':'no value setted'}
        self.Arguments2={'L4_RunID':'pyrolysis','L5_Pressure':'no value setted','L6_Theorie':'13','L6_ResidenceTime':'0.00'}
        self.Arguments3_Ramp={'L7_TtHistory':'3','L81_ttotal':'no value setted','L82_dt':'0.00001','L83_dT':'no value setted','L84_WallTemp':'no value setted','L85_PyrolTemp':'no value setted','L86':'1.0','L87_T0':'no value setted','L88_Hr':'no value setted'}
        self.Arguments3_File={'L7_TtHistory':'5','L8_THistFile':'no value setted'}
        self.Arguments4={'L9_Ash':'0.0','L10_Moisture':'0.0'}

    def set1CoalLocation(self,PathCoalFile):
        """Sets the coal composition file directory."""
        self.Arguments1['L1_Coalfile'] = PathCoalFile

    def set2KinLocation(self,PathKinFile):
        """Sets the coal kinetic file directory."""
        self.Arguments1['L2_Kinfile'] = PathKinFile

    def set3PolyLocation(self,PathPolyFile):
        """Sets the coal polymer file directory."""
        self.Arguments1['L3_Polyfile'] = PathPolyFile

    def set4RunID(self,pyrolysisORgeology):
        """Sets weather the pyrolysis process or the geolegy process shall be modeled. For more information see FG-DVC manual."""
        self.Arguments2['L4_RunID'] = pyrolysisORgeology

    def set5Pressure(self,PressureIn_atm):
        """Sets the operating pressure (float)."""
        self.Arguments2['L5_Pressure'] = '%.2f'%PressureIn_atm

    def set6Theorie(self,Theorie,ResidenceTime):
        """Sets the theory: 13 for no or partial tar cracking, 15 for full tar cracking. The residence input should be 0.0 for no tar cracking or a time greater zero for partial tar pressure. Full tar pressure also requires 0.0 as input, as it is the characteristic input of FG-DVC (although writing here 0.0, the full tar cracking is activated)."""
        self.Arguments2['L6_Theorie'] = '%.0f'%Theorie
        self.Arguments2['L6_ResidenceTime'] = '%.2f'%ResidenceTime

    def set7Ramp(self,timeTotal,dt,dT,finalPyrolysisTemp,initialT,HeatingRate):
        """Sets the following time relevant and temperature history relevant parameter: the total simulation time 'timeTotal', the constant numerical time step 'dt', the maximum temperture step 'dT', the final pyrolysis temperature 'finalPyrolysisTemp', the temperature at t=0 'initialT', and the constant heating rate 'HeatingRate'. All these parameter are required for the case a linear heating rate should be modeled."""
        self.Arguments3_Ramp['L81_ttotal'] = '%.3f'%timeTotal
        self.Arguments3_Ramp['L82_dt'] = '%2.1e'%dt
        self.Arguments3_Ramp['L83_dT'] = '%.1f'%dT
        self.Arguments3_Ramp['L84_WallTemp'] = '%.1f'%finalPyrolysisTemp
        self.Arguments3_Ramp['L85_PyrolTemp'] = '%.1f'%finalPyrolysisTemp
        self.Arguments3_Ramp['L87_T0'] = '%.1f'%initialT
        self.Arguments3_Ramp['L88_Hr'] = '%.1f'%HeatingRate

    def set7File(self,THistoryFileLocation):
        """For the case a temperature history shall be imported, this method should be used. 'THistoryFileLocation' is the directory of the .txt file containing two columns: first column the time in seconds, the second the temperature indegree Celsius-"""
        self.Arguments3_File['L8_THistFile'] = THistoryFileLocation

    def set9AshMoisture(self,AshContent,MoistureContent):
        """Sets the amount of ash and moisture in the coal. By initializing the SetterAndLauncher object, both of these values are setted equal zero."""
        self.Arguments4['L9_Ash']='%.2f'%AshContent
        self.Arguments4['L10_Moisture']='%.2f'%MoistureContent

    def setTRamp_or_TFile(self,selectRamp_or_File):
        """Select weather the time should be defined using a linear ramp ('selectRamp_or_File'='Ramp') or with a input file ('selectRamp_or_File'='File')."""
        self.selectedTHistory=selectRamp_or_File

    def writeInstructFile(self,Filepath):
        """Writes the File 'instruct.ini' into the directory 'Filepath', which should end with FGDVC_8-2-3/FGDVC ."""
        ini=open(Filepath+'instruct.ini','w')#Keywords:1-15,args:16-70
        ini.write('COALFILE       '+self.Arguments1['L1_Coalfile']+'\n')
        ini.write('KINFILE        '+self.Arguments1['L2_Kinfile']+'\n')
        ini.write('POLYFILE       '+self.Arguments1['L3_Polyfile']+'\n')
        ini.write('RUNID          '+self.Arguments2['L4_RunID']+'\n')
        ini.write('PRESSURE       '+self.Arguments2['L5_Pressure']+'\n')
        ini.write('THEORY         '+self.Arguments2['L6_Theorie']+'   '+self.Arguments2['L6_ResidenceTime']+'\n')
        if self.selectedTHistory=='Ramp':
            ini.write('TTHISTORY      '+'3'+'\n')
            ini.write(self.Arguments3_Ramp['L81_ttotal']+' '+self.Arguments3_Ramp['L82_dt']+' '+self.Arguments3_Ramp['L83_dT']+' '+self.Arguments3_Ramp['L84_WallTemp']+' '+self.Arguments3_Ramp['L85_PyrolTemp']+' 1.0 '+self.Arguments3_Ramp['L87_T0']+' '+self.Arguments3_Ramp['L88_Hr']+'\n')
        elif self.selectedTHistory=='File':
            ini.write('TTHISTORY      '+'5'+'\n')
            ini.write(self.Arguments3_File['L8_THistFile']+'\n')
        else:
            print('no option selected. Please setTRamp_or_TFile(selectRamp_or_File)'+'\n')
        ini.write('ASH            '+self.Arguments4['L9_Ash']+'\n')
        ini.write('MOISTURE       '+self.Arguments4['L10_Moisture']+'\n')
        ini.write('START\n')
        ini.write('STOP')
        ini.close()

    def Run(self,PathFromEXE):
        """Lauchnes fgdvcd.exe. The 'PathFromEXE' should not include the .exe and end with a backslash."""
        os.system(PathFromEXE)#+'fgdvcd.exe')         #only working possibility, but python file has to be in the same dir as the fgdvcd


################################################################################
class Process(object):  #planed just for different linear heating rates,DAF Basis
    """Calculates the kinetic parameter using several FG-DVC runs with different heating rates."""
    def __init__(self,DirectoryFromFGDVCD_EXE):
        self.Yields2Cols={'Time':0,'Temp':1,'OneSpeciesVM':2,'H2O':3,'Tar':4,'CO':5,'CO2':6,'CH4':7,'C2H4':8,'HCN':9,'NH3':10,'SO2':11,'COS':12,'CS2':13,'H2S':14,'Olefin':15,'Parafin':16,'H2':17}
        self.Cols2Yields={0:'Time',1:'Temp',2:'OneSpeciesVM',3:'H2O',4:'Tar',5:'CO',6:'CO2',7:'CH4',8:'C2H4',9:'HCN',10:'NH3',11:'SO2',12:'COS',13:'CS2',14:'H2S',15:'Olefin',16:'Parafin',17:'H2'}
        self.PathFromFGDVC_EXE=DirectoryFromFGDVCD_EXE
        os.mkdir(self.PathFromFGDVC_EXE+'runs') #makes directory to save there a copy of each FG-DVC output

    def setFiles(self,CoalFile,KinFile,PolyFile):
	"""Sets the inputted coal files as class intern parameter."""
        self.CoalFile=CoalFile
        self.KinFile=KinFile
        self.PolyFile=PolyFile

    def SpeciesName(self,SpeciesIndice):
	"""Returns the Name of species with the Input column number."""
        return self.Cols2Yields[SpeciesIndice]

    def setConditions(self,pressureATM,totaltime,dt,dT,T_End,T_Begin):
	"""Defines the OPerating conditions in FG-DVC."""
        self.pressure=(pressureATM)
        self.totaltime=(totaltime)
        self.dt=(dt)
        self.dT=(dT)
        self.T_End=(T_End)
        self.T_Begin=(T_Begin)

    def CpFile(self,genSetterAndLauncher,PathFromEXE,HeatingR):
	"""Copies the current FG-DVC result file into the folder 'runs' and renames it with their heating rate."""
        shutil.copyfile(PathFromEXE+'gasyield.txt', PathFromEXE+'runs\\'+'LinRegr'+'_'+str(int(HeatingR))+'gasyield.txt')
        shutil.copyfile(PathFromEXE+'gasrate.txt', PathFromEXE+'runs\\'+'LinRegr'+'_'+str(int(HeatingR))+'gasrate.txt')

    def ReadYields(self,DirectoryWhereFGDVCoutFilesAreLocated):
	"""Reads the current FG-DVC yield result file."""
        #makes reading object
        ResultObject=Fit_one_run.FGDVC_Result(DirectoryWhereFGDVCoutFilesAreLocated)
        #gets data from this object
        self._yields=ResultObject.Yields_all()
        self._rates=ResultObject.Rates_all()
        #imports dictionaries:
        self.Yields2Cols=ResultObject.DictYields2Cols()
        self.Cols2Yields=ResultObject.DictCols2Yields()
        self.makeDt()

    def ReadRates(self,DirectoryWhereFGDVCoutFilesAreLocated):
	"""Reads the current FG-DVC rates."""
        self._rates=(np.genfromtxt(DirectoryWhereFGDVCoutFilesAreLocated+'gasrate.txt',skip_header=2))
        self._rates[:,1]=self._rates[:,1]+273.15
        self._rates[:,3:]=self._rates[:,3:]*(1./60.)
        self._rates[:,17]=self._rates[:,17]-self._rates[:,3]-self._rates[:,5]-self._rates[:,6]-self._rates[:,7]-self._rates[:,8]-self._rates[:,9]-self._rates[:,10]-self._rates[:,11]-self._rates[:,12]-self._rates[:,13]-self._rates[:,14]-self._rates[:,4]  #turns total Volatile yields,rates into H2
        
    def makeDt(self):
	"""Generates a numpy 1D Array with the time steps."""
        Dt=np.zeros(len(self._yields[:,0]))
        t=self._yields[:,0]
        Dt[0]=t[1]-t[0]
        Dt[1:-1]=(t[2:]-t[:-2])/2.
        Dt[-1]=t[-1]-t[-2]
        self.Dt=Dt
        
    def MaxRateTemp(self,SpeciesIndice):
	"""Returns the temperature where the maximum rate occurs."""
        u=self.Rate(SpeciesIndice)
        Line=np.argmax(abs(u))
        TmaxR = self._rates[Line,1]
        return TmaxR
            
    def Yield(self,speciesCol):
	"""Returns the yield list for the inputted species."""
        return self._yields[:,speciesCol]
            
        
    def Rate(self,speciesCol):
	"""Returns the yield list for the inputted species."""
        return self._rates[:,speciesCol]
            
    def Time(self):
	"""Returns the time Array."""
        return self._yields[:,0]
        
    def Derive(self,u):
        """Derive the input vector u using a CDS."""
        self.makeDt()
        yDot=np.zeros(len(u))
        yDot[0]=(u[1]-u[0])/self.Dt[0]  #approx only first order, but doesn't matter, because at the begin nd end rate=0
        yDot[1:-1]=(u[2:]-u[:-2])/(2*self.Dt[1:-1])
        yDot[-1]=(u[-1]-u[-2])/self.Dt[-1]
        return yDot


    def generateRuns(self,genSetterAndLauncher,lowestHr,highestHr,NumberOfRuns,pltHrCurves=True,PltLegend=False):
	"""Runs the FG-DVC as often and with the heating rates defined before and generates an array containing the heating rate and the temperature where the maximum rate occurs."""
        dHr=(highestHr-lowestHr)/(NumberOfRuns-1)
        #sets Conditions for generated SetterAndLauncher:
        genSetterAndLauncher.set1CoalLocation(self.CoalFile)
        genSetterAndLauncher.set2KinLocation(self.KinFile)
        genSetterAndLauncher.set3PolyLocation(self.PolyFile)
        genSetterAndLauncher.set5Pressure(self.pressure)
        genSetterAndLauncher.setTRamp_or_TFile('Ramp')  #Ramp is necessary for equation
        genSetterAndLauncher.set9AshMoisture(0.0,0.0) #on DAF basis
        ##
        ParamArray=[] #Lines:Hr, Columns: Tmax
        for RunNr in range(NumberOfRuns):
            Hr=(lowestHr+RunNr*dHr)
            #change Hr in ini-file
            genSetterAndLauncher.set7Ramp(self.totaltime,self.dt,self.dT,self.T_End,self.T_Begin,Hr)
            genSetterAndLauncher.writeInstructFile(self.PathFromFGDVC_EXE)
            #run current Hr:
            genSetterAndLauncher.Run(self.PathFromFGDVC_EXE)
            #makes Matrix: first column: hr, the Tmax for species
            CurrentHrVec=[float(Hr)]
            # read results to save them
            self.ReadYields(self.PathFromFGDVC_EXE)
            self.ReadRates(self.PathFromFGDVC_EXE)
            self.CpFile(genSetterAndLauncher,self.PathFromFGDVC_EXE,Hr)
            for spec in range(1,len(self.Cols2Yields),1):  #1 instead of 2, array size (columns) should be equal to the FG-DVC output array(columns), so array line = [Hr,T_const,T_CharAndAsh,T_H2O,...,T_H2], where 1=T will be used to ensure that T_max is below T_const
                CurrentHrVec.append(float(self.MaxRateTemp(spec)))
            ParamArray.append(CurrentHrVec)
        self.ParamArray=np.array(ParamArray)
        #plots the Heatingrate curves for the different species
        if pltHrCurves==True:
            print 'plot results'
            #firstly OneSpeciesVM, and H2 because has no rate
            #Yield OneSpeciesVM
            for i in range(NumberOfRuns):
                self.ReadYields(self.PathFromFGDVC_EXE+'runs\\'+'LinRegr'+'_'+str(int(lowestHr+i*dHr)))
                plt.plot(self.Yield(1),self.Yield(2),label=str(int(lowestHr+i*dHr)))
            plt.grid()
            if PltLegend==True:
                plt.legend()
            plt.title('One Species VM')
            plt.xlabel('T in K')
            plt.ylabel('yield in wt%')
            plt.savefig(self.PathFromFGDVC_EXE+'runs\\'+'Yields_'+'LinRegr'+'_'+'OneSpeciesVM'+'.pdf',format='pdf')
            plt.clf();plt.cla()
                #Rate:
            for i in range(NumberOfRuns):
                self.ReadYields(self.PathFromFGDVC_EXE+'runs\\'+'LinRegr'+'_'+str(int(lowestHr+i*dHr)))
                M=self.Yield(2)
                Mdot=self.Derive(M)
                plt.plot(self.Yield(1),Mdot,label=str(int(lowestHr+i*dHr)))
            plt.grid()
            if PltLegend==True:
                plt.legend(loc=4)
            plt.title('One Species VM')
            plt.xlabel('T in K')
            plt.ylabel('rate in wt%/s')
            plt.savefig(self.PathFromFGDVC_EXE+'runs\\'+'Rates_'+'LinRegr'+'_'+'OneSpeciesVM'+'.pdf',format='pdf')
            plt.clf();plt.cla()              
            # now the other species:
            #plt Yields
            for species in range(3,len(self.Cols2Yields),1):
                    for i in range(NumberOfRuns):
                        self.ReadYields(self.PathFromFGDVC_EXE+'runs\\'+'LinRegr'+'_'+str(int(lowestHr+i*dHr)))
                        plt.plot(self.Yield(1),self.Yield(species),label=str(int(lowestHr+i*dHr)))
                    plt.grid()
                    if PltLegend==True:
                        plt.legend()
                    plt.xlabel('T in K')
                    plt.ylabel('yield in wt%')
                    plt.title(self.SpeciesName(species))
                    plt.savefig(self.PathFromFGDVC_EXE+'runs\\'+'Yields'+'_'+'LinRegr'+'_'+str(self.SpeciesName(species))+'.pdf',format='pdf')
                    plt.clf();plt.cla()
            #plt Rates
            for species in range(3,len(self.Cols2Yields),1):
                    for i in range(NumberOfRuns):
                        self.ReadRates(self.PathFromFGDVC_EXE+'runs\\'+'LinRegr'+'_'+str(int(lowestHr+i*dHr)))
                        plt.plot(self.Rate(1),self.Rate(species),label=str(int(lowestHr+i*dHr)))
                    plt.grid()
                    if PltLegend==True:
                        plt.legend(loc=4)
                    plt.xlabel('T in K')
                    plt.ylabel('rate in wt%/s')
                    plt.title(self.SpeciesName(species))
                    plt.savefig(self.PathFromFGDVC_EXE+'runs\\'+'Rates'+'_'+'LinRegr'+'_'+str(self.SpeciesName(species))+'.pdf',format='pdf')
                    plt.clf();plt.cla()
                
                    

    def calcAE(self,genSetterAndLauncher):
	"""Calculates the Arrhenius rate parameter using the results generated with the method 'GenerateRuns'."""
        ParamFile=open('runs\\ParameterForSeveralRuns.txt','w')
        ParamError=open('runs\\ParameterForSeveralRuns_std_Error.txt','w')
        #for spec in range(2,len(self.Cols2Yields)-2,1):
        ParamFile.write('Species\t\tA\t\tE\n')
        ParamError.write('Species\tStandard Error of Curve\n')
        AEArray=[]
        for spec in range(2,len(self.Cols2Yields),1):
            CurrentIAE=[]
            if self.ParamArray[-1,spec]<self.T_End:     # ensure that Tmax is on ramp, not on const T; was an assumption for eq.:
                x=np.zeros(len(self.ParamArray[:,spec]))
                x=x+1
                #print x
                #print self.ParamArray
                x=x/self.ParamArray[:,spec]  #=1/Tmax
                y=np.log( self.ParamArray[:,0]/(self.ParamArray[:,spec]*self.ParamArray[:,spec]) )   #=Hr/Tmax**2
                slope, intercept, r_value, p_value, std_error = stats.linregress(x,y)
                E_calc=-R*slope
                A_calc=np.exp(intercept)*E_calc/R
                ParamFile.write(str(self.SpeciesName(spec))+'\t'+str(A_calc)+'\t'+str(E_calc)+'\n')
                ParamError.write(str(self.SpeciesName(spec))+'\t'+str(std_error)+'\n')
                CurrentIAE.append([spec,A_calc,E_calc])
                AEArray.append(CurrentIAE)
                #plot:
                #plt.plot(genSetterAndLauncher.Rate(0),genSetterAndLauncher.Rate(spec),'x')
                #plt.plot(genSetterAndLauncher.Rate(0),slope*genSetterAndLauncher.Rate(0)+intercept,'r-')
                plt.plot(x,y,'bx')
                plt.plot(x,intercept+slope*x,'r-')
                plt.xlabel('1/T_max in 1/K')
                plt.ylabel('ln[Hr/(T_max*T_max)]')#min')
                plt.grid()
                plt.title(str(self.SpeciesName(spec)))
                plt.savefig('runs\\Regression_'+str(self.SpeciesName(spec))+'.pdf',format='pdf')
                plt.clf(),plt.cla()
            else:
                print 'Species ',self.SpeciesName(spec),' not considered, because Tmax is not on the ramp where Hr>0'
        self.AEArray=np.array(AEArray)
        ParamFile.close(); ParamError.close()
        
    def getAE(self,Species):
	"""Returns the kinetic parameter for the input species."""
        B=np.where(self.AEArray[:,0]==Species)
        B=B[0];B=B[0]
        B=self.AEArray[B,:]
        return B
        
    def TtInterpol(self):
	"""Outputs a Interpolation object for T(t)."""
        OrderOfTimeInterpol=2
        t=self.Yield(0)
        T=self.Yield(1)
        T_Interpol=scipy.interpolate.interp1d( np.array(list(t)+[10*t[-1]]),np.array(list(T)+[T[-1]]),kind=OrderOfTimeInterpol,axis=-1,copy=True,bounds_error=True,fill_value=np.nan )
        return T_Interpol
          
    def CompareResults(self,genSetterAndLauncher,HrWhereToCompare):
	"""Compares the FG-DVC output plot with the calculated plot using the generated Arrhenius parameter."""
        self.generateRuns(genSetterAndLauncher,HrWhereToCompare,HrWhereToCompare,2,False,False)
        for spec in range(2,len(self.Cols2Yields),1):#(2,len(self.Cols2Yields),1):
            t=self.Yield(0)
            u=self.Yield(spec)
            uDot=self.Rate(spec)
            v=self.calcYield(genSetterAndLauncher,spec)
            vDot=self.Derive(v)
            #plot yields
            plt.plot(t,u,'b-',label='FG-DVC'+str(HrWhereToCompare))
            plt.plot(t,v,'g-',label='Hr-Approx')
            plt.xlabel('t in s')
            plt.ylabel('yield in wt%')
            plt.grid()
            plt.title(str(self.SpeciesName(spec)))
            plt.savefig('runs\\Compare_'+str(self.SpeciesName(spec))+'_Y.pdf',format='pdf')
            plt.clf(),plt.cla()
            #plot rates
            plt.plot(t,uDot,'b-',label='FG-DVC'+str(HrWhereToCompare))
            plt.plot(t,vDot,'g-',label='Hr-Approx')
            plt.xlabel('t in s')
            plt.ylabel('rate in wt%/s')
            plt.grid()
            plt.title(str(self.SpeciesName(spec)))
            plt.savefig('runs\\Compare_'+str(self.SpeciesName(spec))+'_R.pdf',format='pdf')
            plt.clf(),plt.cla()
            
    def calcYield(self,genSetterAndLauncher,Species):
	"""Calcules the yields for the input species."""
        absTol=1.e-8
        relTol=1.e-6
        ###
        u=self.Yield(Species)
        m_s0=u[-1]
        IC=u[0]
        T=self.TtInterpol()
        time=self.Yield(0)
        def dmdt(m,t):
            Param=self.getAE(Species)
            Param=Param[0]
            #Param[0]=SpeciesNr
            #Param[1]=A
            #Param[2]=E
            if Species==2:
                dmdt_out=-(Param[1]*np.exp(-Param[2]/T(t)))*(m)
                dmdt_out=np.where(abs(dmdt_out)>1.e-300,dmdt_out,0.0)
            else:
                dmdt_out=(Param[1]*np.exp(-Param[2]/T(t)))*(m_s0-m)
            return dmdt_out
        m_out=scipy.integrate.odeint(dmdt,IC,time,atol=absTol,rtol=relTol,hmax=1.e-4)
        return m_out[:,0]
            
