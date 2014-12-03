import numpy as np
import os
from scipy import stats
import pylab as plt
import shutil
from scipy.integrate import odeint
import scipy.interpolate
import Fit_one_run 
import platform

OS = platform.system()

################################
R=1.0 #8.3144621 # Gas constant only =8.3... if E should not include R
################################

class SetterAndLauncher(object):
    """This class is able to write the CPD input file and launch the CPD program. Before writing the CPD input file (method 'writeInstructFile') set all parameter using the corresponding methods (SetCoalParameter, SetOperateCond, SetNumericalParam, CalcCoalParam). After writing the instruct file, the .Run method can be used."""
    def __init__(self):
        # in the following, rarely changed parameter are defined
        self.ab=2.602e15
        self.eb=55400
        self.ebsig=1800
        self.ac=0.9
        self.ec=0
        self.ag=3.0e15
        self.eg=69000
        self.egsig=8100
        self.Acr=3.0e15
        self.Ecr=65000
        self.arad=18.4
        self.erad=6000
        self.fstable=0.03
        self.an=5.5e7
        self.en=90000
        self.ensig=0
        self.nmax=20
        
    

    def SetCoalParameter(self,fcar,fhyd,fnit,foxy,VMdaf):
        """Set the mass fraction of carbon (UA), hydrogen  (UA), nitrogen (UA), oxygen (UA) and the daf fraction of volatile matter (PA)."""
        self.fcar=fcar
        self.fhyd=fhyd
        self.fnit=fnit
        self.foxy=foxy
        self.VMdaf=VMdaf
        
    def SetOperateCond(self,pressure,TimeTemp):  #TimeTemp is an 2D-Array: [time_i,temp_i]
        """Stes the operating condtions pressure and the time-temperature array. TimeTemp is an 2D-Array. One line of this array is [time_i,temp_i]."""
        self.pressure=pressure
        self.TP=np.array(TimeTemp)
        
      
    def SetNumericalParam(self,dt,timax):       #dt: Vector [dt,increment,max dt]
        """Sets the numerical parameter and the maximum simulation time. dt is a vector with the following information: [dt_initial,print-increment,dt_max]"""
        self.dt=dt
        self.timax=timax
        
    def CalcCoalParam(self):
        """Calculates the CPD coal parameter mdel, mw, p0, sig and sets the as an attribute of the class."""
        #uses equations from: http://www.et.byu.edu/~tom/cpd/correlation.html
        #c[0,0]=c0
        #c[1,Column]=c1
        #... see table in http://www.et.byu.edu/~tom/cpd/correlation.html
        c=np.array([[0.0,0.0,0.0,0.0],
               [421.957,1301.41,0.489809,-52.1054],
               [ -8.64692,16.3879,-0.00981566,1.63872],
               [0.0463894,-0.187493,0.000133046,-0.0107548],
               [-8.47272,-454.773,0.155483,-1.23688],
               [1.18173,51.7109,-0.0243873,0.0931937],
               [1.15366,-10.0720,0.00705248,-0.165673],
               [-0.0434024,0.0760827,0.000219163,0.00409556],
               [0.556772,1.36022,-0.0110498,0.00926097],
               [-0.00654575,-0.0313561,0.000100939,-0.0000826717]])
        # calculates c0:
        if self.fcar>85.9:
            c[0,0] = 0.1183*self.fcar - 10.16
            if c[0,0]>0.36:
                c[0,0]=0.0
        elif self.foxy>12.5:
            c[0,0] = 0.014*self.foxy - 0.175
            if c[0,0]>0.15:
                c[0,0]=0.0
        else:
            c[0,0]=0.0
        #
        self.c0=c[0,0]
        #
        #
        Y=np.zeros(int(len(c[1,:])))
        for i in range (int(len(c[1,:]))):
            Y[i] = c[1,i]+c[2,i]*self.fcar+c[3,i]*self.fcar**2+c[4,i]*self.fhyd+c[5,i]*self.fhyd**2+c[6,i]*self.foxy+c[7,i]*self.foxy**2+c[8,i]*self.VMdaf+c[9,i]*self.VMdaf**2
        self.mdel=Y[0]
        self.mw=Y[1]
        self.p0=Y[2]
        self.sig=Y[3]
        
    

    def writeInstructFile(self,Dirpath):
        """Writes the File 'CPD_input.dat' into the directory Dirpath."""
        ini=open(Dirpath+'CPD_input.dat','w')#Keywords:1-15,args:16-70
        ini_str = """
{p0}                           ! p0
{c0}                           ! c0
{sig}                          ! sig+1
{mw}                           ! mw
{mdel}                         ! mdel (7 will be subtracted internally to the CPD model
{fcar/100.}                    ! fcar (daf mass fraction of carbon in unreacted coal)
{fhyd/100.}                    ! fhyd (daf mass fraction of hydrogen in unreacted coal)
{fnit/100.}                    ! fnit (daf mass fraction of nitrogen in unreacted coal)
{foxy/100.}                    ! foxy (daf mass fraction of oxygen in unreacted coal)

{ab}                           ! ab
{eb}                           ! eb
{ebsig}                        ! ebsig
{ac}                           ! ac=rho
{ec}                           ! ec
{ag}                           ! ag
{eg}                           ! eg
{egsig}                        ! egsig
{Acr}                          ! Acr (pre-exponential factor for crosslinking rate)
{Ecr}                          ! Ecr (activation energy for crosslinking rate)
#
{arad}                         ! arad (pre-exponential factor for N attack by free radical)
{erad}                         ! erad (activation energy for N attack by free radical, cal.)
{fstable}                      ! fstable (initial frac. of MW decay with no radical N attack)
{an}                           ! an (high T slow N release pre-exponential factor)
{en}                           ! en (high T slow N release activation energy, calories)
{ensig}                        ! ensig (deviation bound for distribution of en)

{pressure}                     ! pressure (atm)
{len(self.TP)}                 ! number of time points
{TP[0,0]) {TP[0,1]}            ! time(ms),temp(K)
""".format(**self.__dict__)
        time_points_str = ""
        for i in range(1,len(self.TP),1):
            time_points_str += "{}, {} \n".format(self.TP[i,0],self.TP[i,1])

        ini_str += """
{dt[0]} {dt[1]} {dt[2]}        ! dt (s),print increment,max dt (s)\n')
{self.timax}                   ! timax (maximum residence time [s] for calculations)\n')
{str(self.nmax}                ! nmax (maximum number of mers for tar molecular wt)\n')
""".format(**self.__dict__)
        print ini_str
        ini.close()

    def Run(self, run_nr, inp_file="IN.dat"):
        """Launches the CPD executable and inputs Input_File. 

           TODO: If the CPD executable is in  another directory than the 
           current working directory the script should use absolute path 
           for the CPD_exe.
        """
        if OS == 'Linux':
            exe = './cpdnlg'
        elif OS == 'Darwin':
            exe = './cpdnlg.x'
        elif OS == 'Windows':
            exe = './cpdnlg.exe'
        else:
            print "The name of the operating system couldn't be found."
            return 
        OScommand='{} < {} > CPD_{}_output.log'.format(exe, inp_file, run_nr)
        os.system(OScommand)


################################################################################

class ProcessCPD(object):  #planed just for different linear heating rates,DAF Basis
    """Calculates the kinetic parameter using several CPD runs with different heating rates."""
    def __init__(self,DirectoryFromCPD_EXE):
        os.mkdir(DirectoryFromCPD_EXE+'runs')
        self.PathFromCPD_EXE=DirectoryFromCPD_EXE

    def SpeciesName(self,SpeciesIndice):
	"""Returns the Name of species with the Input column number."""
        return self.Cols2Yields[SpeciesIndice]

    def SpeciesIndex(self,SpeciesName):
	"""Returns the column number of the input species."""
        return self.Yields2Cols[SpeciesName]

    def CpFile(self,genSetterAndLauncher,PathFromEXE,HeatingR):
	"""Copies the current FG-DVC result file into the folder 'runs' and renames it with their heating rate."""
        shutil.copyfile(PathFromEXE+'CPD_Result1.dat', PathFromEXE+'runs/'+str(HeatingR)+'_CPD_Result1.dat')
        shutil.copyfile(PathFromEXE+'CPD_Result4.dat', PathFromEXE+'runs/'+str(HeatingR)+'_CPD_Result4.dat')

    def ReadYields(self,DirectoryWhereFGDVCoutFilesAreLocated):
	"""Reads the yields and rates of the current run."""
        #makes reading object
        ResultObject=Fit_one_run.CPD_Result(DirectoryWhereFGDVCoutFilesAreLocated)
        #gets data from this object
        self._yields=ResultObject.Yields_all()
        self._rates=ResultObject.Rates_all()
        #imports dictionaries:
        self.Yields2Cols=ResultObject.DictYields2Cols()
        self.Cols2Yields=ResultObject.DictCols2Yields()
        self.makeDt()
       
    def makeDt(self):
	"""Generates a list of the time steps."""
        Dt=np.zeros(len(self._yields[:,0]))
        t=self._yields[:,0]
        Dt[0]=t[1]-t[0]
        Dt[1:-1]=(t[2:]-t[:-2])/2.
        Dt[-1]=t[-1]-t[-2]
        self.Dt=Dt
        
    def MaxRateTemp(self,SpeciesIndice):
	"""Returns the temperature where the maximum rate occurs."""
        u=self._rates[:,SpeciesIndice]
        Line=np.argmax(abs(u))
        TmaxR = self._rates[Line,1]
        return TmaxR
            
    def Time(self):
	"""Returns the time array."""
        return self._yields[:,0]
        
    def Derive(self,u):
	"""Derives the input vector using a CDS."""
        self.makeDt()
        yDot=np.zeros(len(u))
        yDot[0]=(u[1]-u[0])/self.Dt[0]  #approx only first order, but doesn't matter, because at the begin nd end rate=0
        yDot[1:-1]=(u[2:]-u[:-2])/(2*self.Dt[1:-1])
        yDot[-1]=(u[-1]-u[-2])/self.Dt[-1]
        return yDot


    def generateRuns(self,genSetterAndLauncher,tEnd,Tstart,TEnd,lowestHr,highestHr,NumberOfRuns,CPD_exe,Input_File,PltLegend=False):
	"""Runs CPD for the range of the heating rates and writes an Array containing the heating rates and the corresponding maximum rate temperatures."""
        self. tEnd=   tEnd    
        self.T_End=TEnd
        self.Tstart=Tstart
        self.CPD_exe=CPD_exe
        self.Input_File=Input_File
        dHr=(highestHr-lowestHr)/(NumberOfRuns-1)
        ParamArray=[] #Lines:Hr, Columns: Tmax for different Species
        for RunNr in range(NumberOfRuns):
            Hr=(lowestHr+RunNr*dHr)
            #change Hr in ini-file
            t1=(TEnd-Tstart)/(Hr/1000.) #Hr/1000. -> result in ms
            if t1<tEnd:
                genSetterAndLauncher.TP=np.array([[0.0,Tstart],[t1,TEnd],[tEnd,TEnd]])
            else:
                genSetterAndLauncher.TP=np.array([[0.0,Tstart],[t1,TEnd]])
            genSetterAndLauncher.writeInstructFile(self.PathFromCPD_EXE)
            #run current Hr:
            genSetterAndLauncher.Run(CPD_exe,Input_File)
            #makes Matrix: first column: hr, the Tmax for species
            CurrentHrVec=[float(Hr)]
            # read results to save them
            self.ReadYields(self.PathFromCPD_EXE)
            self.CpFile(genSetterAndLauncher,self.PathFromCPD_EXE,Hr)
            for spec in range(2,len(self.Cols2Yields),1):  #1 instead of 2, array size (columns) should be equal to the FG-DVC output array(columns), so array line = [Hr,T_const,T_CharAndAsh,T_H2O,...,T_H2], where 1=T will be used to ensure that T_max is below T_const
                CurrentHrVec.append(float(self.MaxRateTemp(spec)))
            ParamArray.append(CurrentHrVec)
        self.ParamArray=np.array(ParamArray)
                
                    

    def calcAE(self,genSetterAndLauncher):
	"""Calculates the Arrhenius parameter using the array generated in 'generateRuns'."""
        ParamFile=open('runs/ParameterForLinRegrRuns.txt','w')
        ParamError=open('runs/ParameterForLinRegrRuns_std_Error.txt','w')
        #for spec in range(2,len(self.Cols2Yields)-2,1):
        ParamFile.write('Species\t\tA\t\tE\n')
        ParamError.write('Species\tStandard Error of Curve\n')
        AEArray=[]
        for spec in range(2,len(self.Cols2Yields)-1,1):
            CurrentIAE=[]
            if self.ParamArray[-1,spec]<self.T_End:     # ensure that Tmax is on ramp, not on const T; was an assumption for eq.:
                x=np.zeros(len(self.ParamArray[:,spec]))
                x=x+1
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
                plt.savefig('runs/Regression_'+str(self.SpeciesName(spec))+'.pdf',format='pdf')
                plt.clf(),plt.cla()
            else:
                print 'Species ',self.SpeciesName(spec),' not considered, because Tmax is not on the ramp where Hr>0'
        self.AEArray=np.array(AEArray)
        ParamFile.close(); ParamError.close()
        
    def getAE(self,Species):
	"""returns the calculated kinetic parameter."""
        B=np.where(self.AEArray[:,0]==Species)
        B=B[0];B=B[0]
        B=self.AEArray[B,:]
        return B
        
    def TtInterpol(self):
	"""Returns an interpolation object of T(t)."""
        OrderOfTimeInterpol=2
        t=self._yields[:,0]
        T=self._yields[:,1]
        T_Interpol=scipy.interpolate.interp1d( np.array(list(t)+[10*t[-1]]),np.array(list(T)+[T[-1]]),kind=OrderOfTimeInterpol,axis=-1,copy=True,bounds_error=True,fill_value=np.nan )
        return T_Interpol
          
    def CompareResults(self,genSetterAndLauncher,HrWhereToCompare):
	"""Plots the genearted curve and the CPD output in a Yield vs. time diagram, each for one species."""
        self.generateRuns(genSetterAndLauncher,self.tEnd,self.Tstart,self.T_End,HrWhereToCompare,HrWhereToCompare,2,self.CPD_exe,self.Input_File)
        for spec in range(2,len(self.Cols2Yields)-1,1):#(2,len(self.Cols2Yields),1):
            t=self._yields[:,0]
            u=self._yields[:,spec]
            uDot=self._rates[:,spec]
            v=self.calcYield(genSetterAndLauncher,spec)
            vDot=self.Derive(v)
            #plot yields
            plt.plot(t,u,'b-',label='CPD'+str(HrWhereToCompare))
            plt.plot(t,v,'g-',label='Hr-Approx')
            plt.xlabel('t in s')
            plt.ylabel('yield in wt%')
            plt.grid()
            plt.title(str(self.SpeciesName(spec)))
            plt.savefig('runs/Compare_'+str(self.SpeciesName(spec))+'_Y.pdf',format='pdf')
            plt.clf(),plt.cla()
            #plot rates
            plt.plot(t,uDot,'b-',label='CPD'+str(HrWhereToCompare))
            plt.plot(t,vDot,'g-',label='Hr-Approx')
            plt.xlabel('t in s')
            plt.ylabel('rate in wt%/s')
            plt.grid()
            plt.title(str(self.SpeciesName(spec)))
            plt.savefig('runs/Compare_'+str(self.SpeciesName(spec))+'_R.pdf',format='pdf')
            plt.clf(),plt.cla()
            
    def calcYield(self,genSetterAndLauncher,Species):
	"""Calculates the yields.""" 
        absTol=1.e-8
        relTol=1.e-6
        ###
        u=self._yields[:,Species]
        m_s0=u[-1]
        IC=u[0]
        T=self.TtInterpol()
        time=self._yields[:,0]
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
      
