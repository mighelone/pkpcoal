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

class SetAndLaunchBase(object):
    runNr=0
    printIntervall=1
    writeIntervall=1

class CPD(SetAndLaunchBase):
    """ This class is able to write the CPD input file and launch the CPD program. 
        Before writing the CPD input file (method 'writeInstructFile') set all parameter 
        using the corresponding methods (SetCoalParameter, SetOperateCond, SetNumericalParam, 
        CalcCoalParam). After writing the instruct file, the .Run method can be used."""

    cpd_constants = {
           'ab'     : 2.602e15,
           'eb'     : 55400,
           'ebsig'  : 1800,
           'ac'     : 0.9,
           'ec'     : 0,
           'ag'     : 3.0e15,
           'eg'     : 69000,
           'egsig'  : 8100,
           'Acr'    : 3.0e15,
           'Ecr'    : 65000,
           'arad'   : 18.4,
           'erad'   : 6000,
           'fstable': 0.03,
           'an'     : 5.5e7,
           'en'     : 90000,
           'ensig'  : 0,
           'nmax'   : 20,
        }

    def __init__(self, 
            ultimateAnalysis,
            proximateAnalysisDaf,
            tempProfile,
            pressure,
            deltaT,
            inputs="/home/go/documents/code/pkp.git/inputs/"
        ):
        self.tempProfile = self.timeTempProfile(tempProfile) # TODO give it a better name
        self.pressure    = pressure
        # We scale ua and daf data for cpd since input is not in percents
        self.ultim_ana   = ultimateAnalysis.scale(0.01)
        self.daf         = proximateAnalysisDaf.scale(0.01)
        self.coal_param  = CPD.CalcCoalParam(self.ultim_ana, self.daf)
        self.output_dict = {
            'num_time' : len(tempProfile),
            'pressure' : self.pressure,
            'deltaT'   : deltaT,
            'strTempProfile': self.tempProfile,
            'printIntervall': self.printIntervall,
            'writeIntervall': self.writeIntervall,
        }
        # for generating the output string 
        # all the dicts are merge into one
        self.output_dict.update(self.cpd_constants)
        self.output_dict.update(self.ultim_ana.elems)
        self.output_dict.update(self.coal_param)
        self.output_dict.update(self.daf.elems)
        self.output_dict.update(self.__dict__)
        self.inputs = inputs #TODO GO fix how input folder is defined
        

    @classmethod
    def calcC0(cls, massFracCarbon, massFracOx):
        c0 = 0.0
        #TODO GO are these really mutually exclusive?
        #       what happens if c > 0.859 and ox > 0.123
        #TODO GO double check if correct
        if massFracCarbon > 0.859:
            c0 = 11.83*massFracCarbon - 10.16
            c0 = (0.0 if c0 > 0.36 else c0)
        elif massFracOx > 0.125:
            c0 = 1.4*massFracOx - 0.175
            c0 = (0.0 if c0 > 0.15 else c0)
        return c0


    @classmethod
    def CalcCoalParam(cls, ultim_ana, daf):
        """ Calculates the CPD coal parameter mdel, mw, p0, sig 
            and sets the as an attribute of the class. """
        #uses equations from: http://www.et.byu.edu/~tom/cpd/correlation.html
        #c[0,0]=c0
        #c[1,Column]=c1
        #... see table in http://www.et.byu.edu/~tom/cpd/correlation.html
        c=np.array([
               [ 0.0    ,  0.0    ,  0.0     ,  0.0],
               [ 421.957,  1301.41,  0.489809, -52.1054],
               [-864.692,  1638.79, -0.981566,  163.872],
               [ 463.894, -1874.93,  1.33046 , -107.548],
               [-847.272, -45477.3,  15.5483 , -123.688],
               [ 11817.3,  517109., -243.873 ,  931.937],
               [ 115.366, -1007.20,  0.705248, -16.5673],
               [-434.024,  760.827,  2.19163 ,  40.9556],
               [ 55.6772,  136.022, -1.10498 ,  0.926097],
               [-65.4575, -313.561,  1.00939 , -0.826717],
        ])


        Y = np.zeros(len(c[1,:]))
        for i, yi in enumerate(Y): 
            Y[i] = (c[1,i] 
                     + c[2,i]*ultim_ana['Carbon'] 
                     + c[3,i]*ultim_ana['Carbon']**2 
                     + c[4,i]*ultim_ana['Hydrogen']
                     + c[5,i]*ultim_ana['Hydrogen']**2 
                     + c[6,i]*ultim_ana['Oxygen'] 
                     + c[7,i]*ultim_ana['Oxygen']**2
                     + c[8,i]*daf['Volatile Matter']
                     + c[9,i]*daf['Volatile Matter']**2)

        return {'c0'   : CPD.calcC0(ultim_ana['Carbon'], ultim_ana['Oxygen']),
                'mdel' : Y[0],
                'mw'   : Y[1],
                'p0'   : Y[2],
                'sig'  : Y[3]}
        
    

    def timeTempProfile(self,tempProfile):
        """ Returns a string from yaml read temp profile, basically reversing the yaml read function
            Probably there is a more direct way
        """
        #return '\n'.join(["{} {}".format(time, temp) for time, temp in tempProfile.iteritems()])
        return '\n'.join([' '.join(map(str,_)) for _ in tempProfile])

    def writeInstructFile(self, Dirpath):
        """Writes the File 'CPD_input.dat' into the directory Dirpath."""
        ini=open(self.inputs+'CPD_input.dat','w')#Keywords:1-15,args:16-70
        #TODO GO is fcar,fhyd ... from daf or ua ?
        #TODO GO where does timax and nmax come frome
        ini_str = """{p0:<20} ! p0
{c0:<20} ! c0
{sig:<20} ! sig+1
{mw:<20} ! mw
{mdel:<20} ! mdel (7 will be subtracted internally to the CPD model
{Carbon:<20} ! fcar (daf mass fraction of carbon in unreacted coal)
{Hydrogen:<20} ! fhyd (daf mass fraction of hydrogen in unreacted coal)
{Nitrogen:<20} ! fnit (daf mass fraction of nitrogen in unreacted coal)
{Oxygen:<20} ! foxy (daf mass fraction of oxygen in unreacted coal)
{ab:<20} ! ab
{eb:<20} ! eb
{ebsig:<20} ! ebsig
{ac:<20} ! ac=rho
{ec:<20} ! ec
{ag:<20} ! ag
{eg:<20} ! eg
{egsig:<20} ! egsig
{Acr:<20} ! Acr (pre-exponential factor for crosslinking rate)
{Ecr:<20} ! Ecr (activation energy for crosslinking rate)
{arad:<20} ! arad (pre-exponential factor for N attack by free radical)
{erad:<20} ! erad (activation energy for N attack by free radical, cal.)
{fstable:<20} ! fstable (initial frac. of MW decay with no radical N attack)
{an:<20} ! an (high T slow N release pre-exponential factor)
{en:<20} ! en (high T slow N release activation energy, calories)
{ensig:<20} ! ensig (deviation bound for distribution of en)
{pressure:<20} ! pressure (atm)
{num_time:<20} ! number of time points
{strTempProfile}
{deltaT} {printIntervall} {deltaT} ! dt (s), print increment, max dt (s))
0.03                    ! timax (maximum residence time [s] for calculations))
{nmax:<20} ! nmax (maximum number of mers for tar molecular wt))
""".format(**self.output_dict)
        ini.write(ini_str)
        ini.close()

    def Run(self, inp_file="IN.dat"):
        """Launches the CPD executable and inputs Input_File. 

        """
        self.runNr += 1 
        import PKP.bins
        cpdExec =  os.path.dirname(PKP.bins.__file__)
        if OS == 'Linux':
            exe = '{}/cpdnlg'.format(cpdExec)
        elif OS == 'Darwin':
            exe = '{}/cpdnlg.x'.format(cpdExec)
        elif OS == 'Windows':
            exe = '{}cpdnlg.exe'.format(cpdExec)
        else:
            print "The name of the operating system couldn't be found."
            return 
        OScommand='{} < {}{} > {}CPD_{}_output.log'.format(
                        exe, self.inputs, inp_file, self.inputs, self.runNr)
        print OScommand
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
      