#MODIFIED FOR CPD VERSION: CPD_NLG (1999)
##############Global####################################################################################################
##############
import numpy as np
from StringIO import StringIO
import pylab as plt
import scipy as sp
from scipy.optimize import fmin
from scipy.optimize import fmin_cg
from scipy.optimize import fmin_bfgs
from scipy.optimize import fmin_ncg
from scipy.optimize import leastsq
from scipy.optimize import fmin_slsqp
from scipy.integrate import odeint
import scipy.integrate
import scipy.interpolate
#plt.ion()

####general classes for reading the results and pass the results to the general fitting class
#
#for CPD:    
class CPD_Result(object):
    """Reads the CPD input and saves the values in one array. The results include the yields and the rates. The rates were calculated using a CDS. This class also contains the dictionaries for the columns in the array - the name of the species. These dictionaries are CPD-Version dependent and the only thing which has to be changed for the case of a new release of CPD with a new order of species in the result files."""
    def __init__(self,FilePath):
        self.__path=FilePath
        path1=FilePath+'CPD_Result1.dat'
        path4=FilePath+'CPD_Result4.dat'
        #assembles relevant columns from the two files 'CPD_Result1.dat' and 'CPD_Result4.dat' to one Matrix:
        yields1=(np.genfromtxt(path1,skip_header=1,skip_footer=1))  #last line is twice in this file
        yields4=(np.genfromtxt(path4,skip_header=1))
        self.__yields=np.zeros(( int(len(yields1[:,0])),11) )  #shapes new Matrix containing all necessary information; Files have the same number of lines, 12 because 12 relevant species
        self.__yields[:,0:2]=yields1[:,0:2] #  0=Time   1=Temperature
        self.__yields[:,2:6]=yields1[:,4:8] #  4=ftar    5=fgas   6=fsolid   7=ftot
        self.__yields[:,6:11]=yields4[:,2:7]#   2=fh2O     3=fco2     4=fch4     5=fco     6=fother
        self.__yields[:,0]=self.__yields[:,0]*1.E-3            #t in s instead of ms
        #calculate rates:
        self.__rates=np.zeros(np.shape(self.__yields))      #same dimension
        #cp time and temperature into rates matrix:
        self.__rates[:,0]=self.__yields[:,0];self.__rates[:,1]=self.__yields[:,1]
        for col in range(2,11,1):
            self.__rates[0,col]=(self.__yields[1,col]-self.__yields[0,col])/(self.__yields[1,0]-self.__yields[0,0])
            self.__rates[1:-1,col]=(self.__yields[2:,col]-self.__yields[:-2,col])/(self.__yields[2:,0]-self.__yields[:-2,0])
            self.__rates[-1,col]=(self.__yields[-1,col]-self.__yields[-2,col])/(self.__yields[-1,0]-self.__yields[-2,0])
        print '\nimported data-fields, size(rows,columns): ',self.__yields.shape,'\n'
        #Yields2Cols: updated for FG-DVC Version 8.2.2, one comment: 'Char_and_Ash' is actual the whole coal mass, but this part is named in the FG-DVC output file that way
        self.Yields2Cols={'Time':0,'Temp':1,'Tar':2,'Gas':3,'Solid':4,'Total':5,'H2O':6,'CO2':7,'CH4':8,'CO':9,'Other':10}
        self.Cols2Yields={0:'Time',1:'Temp',2:'Tar',3:'Gas',4:'Solid',5:'Total',6:'H2O',7:'CO2',8:'CH4',9:'CO',10:'Other'}    

    def Yields_all(self):
        """Returns the whole result matrix of the yields."""
        return self.__yields
        
    def Rates_all(self):
        """Returns the whole result matrix of the Rates."""
        return self.__rates
        
    def DictYields2Cols(self):
        """Returns the whole Dictionary Yield names to Columns of the matrix"""
        return self.Yields2Cols
        
    def DictCols2Yields(self):
        """Returns the whole Dictionary Columns of the matrix to Yield names"""
        return self.Cols2Yields

    def FinalYields(self):
        """Returns the last line of the Array, containing the yields at the time=time_End"""
        return self.__yields[-1,:]

    def Name(self):
        """returns 'CPD' as the name of the Program"""
        return 'CPD'

#for FG-DVC
class FGDVC_Result(object):
    """Reads the FG-DVC input and saves the values in one array. The results include the yields (from 'gasyields.txt') and the rates. The rates for all species except the solids (here a CDS is used) are imported from 'gasrates.txt'. The H_2 yields were calculated by subtract all other species except parafins and olefins from the total yields (see FG-DVC manual). This H_2-yieldsd curve was smoothed and derived using a CDS to generate the H_2 rates. The parafins and olefins are added into the tar. This class also contains the dictionaries for the columns in the array - the name of the species. These dictionaries are FG-DVC-Version dependent and the only thing which has to be changed for the case of a new release of FG-DVC with a new order of species in the result files (this was made for Versions 8.2.2. and 8.2.3.)."""
    def __init__(self,FilePath): #for FG-DVC Version 8.2.2 and 8.2.3
        self.__path=FilePath
        self.__yields=(np.genfromtxt(self.__path+'gasyield.txt',skip_header=2))   #,names=True??
        self.__yields[:,1]=self.__yields[:,1]+273.15            #T in K instead degree C
        self.__rates=(np.genfromtxt(self.__path+'gasrate.txt',skip_header=2))     #,names=True??
        self.__rates[:,3:]=self.__rates[:,3:]*(1./60.)          #rate in 1/s instead 1/min
        self.__rates[:,1]=self.__rates[:,1]+273.15              #T in K instead degree C
        #use absolute instead of percentage values:
        self.__yields[:,2:]=self.__yields[:,2:]/100.
        self.__rates[:,2:]=self.__rates[:,2:]/100.
        #calculates Yields and rates of H2 from 'Total'
        self.__yields[:,17]=self.__yields[:,17]-self.__yields[:,3]-self.__yields[:,5]-self.__yields[:,6]-self.__yields[:,7]-self.__yields[:,8]-self.__yields[:,9]-self.__yields[:,10]-self.__yields[:,11]-self.__yields[:,12]-self.__yields[:,13]-self.__yields[:,14]  #converts Total yields into H2 yields, H2 is not reported anywhere else
        self.__rates[:,17]=self.__rates[:,17]-self.__rates[:,3]-self.__rates[:,5]-self.__rates[:,6]-self.__rates[:,7]-self.__rates[:,8]-self.__rates[:,9]-self.__rates[:,10]-self.__rates[:,11]-self.__rates[:,12]-self.__rates[:,13]-self.__rates[:,14]-self.__rates[:,4]
        #calculates rate of Solids as in the gasrate.txt only the solid yields are reported:
        self.__rates[0,2]=(self.__yields[1,2]-self.__yields[0,2])/(self.__yields[1,0]-self.__yields[0,0])
        self.__rates[1:-1,2]=(self.__yields[2:,2]-self.__yields[:-2,2])/(self.__yields[2:,0]-self.__yields[:-2,0])
        self.__rates[-1,2]=(self.__yields[-1,2]-self.__yields[-2,2])/(self.__yields[-1,0]-self.__yields[-2,0])
        #Yields2Cols: updated for FG-DVC Version 8.2.2, one comment: 'Char_and_Ash' is actual the whole coal mass, but this part is named in the FG-DVC output file that way
        self.Yields2Cols={'Time':0,'Temp':1,'Solid':2,'H2O':3,'Tar':4,'CO':5,'CO2':6,'CH4':7,'C2H4':8,'HCN':9,'NH3':10,'SO2':11,'COS':12,'CS2':13,'H2S':14,'Olefin':15,'Parafin':16,'H2':17}
        self.Cols2Yields={0:'Time',1:'Temp',2:'Solid',3:'H2O',4:'Tar',5:'CO',6:'CO2',7:'CH4',8:'C2H4',9:'HCN',10:'NH3',11:'SO2',12:'COS',13:'CS2',14:'H2S',15:'Olefin',16:'Parafin',17:'H2'}
        #Filter for H2
        alpha=0.3  #weight of the neighbor values
        NumberOfFilterRuns=50
        for n in range(NumberOfFilterRuns): #smooths the yield curve of H2
            self.__yields[1:-1,17]=(alpha/2)*(self.__yields[2:,17]+self.__yields[:-2,17])+(1-alpha)*self.__yields[1:-1,17]
        #gets Rate H2 by derive H2 yields:
        self.__rates[0,17]=(self.__yields[1,17]-self.__yields[0,17])/(self.__yields[1,0]-self.__yields[0,0])
        self.__rates[1:-1,17]=(self.__yields[2:,17]-self.__yields[:-2,17])/(self.__yields[2:,0]-self.__yields[:-2,0])
        self.__rates[-1,17]=(self.__yields[-1,17]-self.__yields[-2,17])/(self.__yields[-1,0]-self.__yields[-2,0])
        print '\nimported data-fields, size(rows,columns): ',self.__yields.shape,'\n'
        #merge the parafins and olefins into the tar. In the rates, the parafins and olefins are already included in the FG-DVC output.
        self.__yields[:,4]=self.__yields[:,4]+self.__yields[:,15]+self.__yields[:,16]
        
    def Yields_all(self):
        """Returns the whole result matrix of the yields."""
        return self.__yields
        
    def Rates_all(self):
        """Returns the whole result matrix of the Rates."""
        return self.__rates
        
    def DictYields2Cols(self):
        """Returns the whole Dictionary Yield names to Columns of the matrix"""
        return self.Yields2Cols
        
    def DictCols2Yields(self):
        """Returns the whole Dictionary Columns of the matrix to Yield names"""
        return self.Cols2Yields

    def FinalYields(self):
        """Returns the last line of the Array, containing the yields at the time=time_End"""
        return self.__yields[-1,:]
        
    def FilePath(self):
        """Returns the FG-DVC File path"""
        return self.__path

    def Name(self):
        """returns 'FG-DVC' as the name of the Program"""
        return 'FG-DVC'

    
###################################################################################################################################
#general fitting class for the use for CPD and FG-DVC:
class Fit_one_run(object):
    """Imports from the Result objects the arrays. It provides the fitting objects with the yields and rates over time for the specific species. This class futher offers the option to plot the generated fitting results."""
    def __init__(self,ResultObject):
        #importes yields and rates arrays
        self.__yields=ResultObject.Yields_all()
        self.__rates=ResultObject.Rates_all()
        #imports dictionaries:
        self.Yields2Cols=ResultObject.DictYields2Cols()
        self.Cols2Yields=ResultObject.DictCols2Yields()
        self.__ImportedResultObject=ResultObject
        
    def plt_InputVectors(self,xVector,y1Vector,y2Vector,y3Vector,y4Vector,y1Name,y2Name,y3Name,y4Name):
        """Plots the y input Vectors vs. the x input vector."""
        plt.plot(xVector,y1Vector,label=y1Name)
        plt.plot(xVector,y2Vector,label=y2Name)
        plt.plot(xVector,y3Vector,label=y3Name)
        plt.plot(xVector,y4Vector,label=y4Name)
        plt.xlabel('t in s')
        plt.ylabel('yield in wt%')
        plt.legend()
        plt.grid()
        plt.savefig('Compare_'+self.__ImportedResultObject.Name()+'.pdf',format='pdf')
        plt.clf(),plt.cla()

    def plt_YieldVsTime(self,ColumnNumber):
        """Plots the original yield output of the pyrolysis program (as e.g. CPD) of the species marke with the columns number"""
        plt.plot(self.__yields[:,0],self.__yields[:,ColumnNumber],label=self.SpeciesName(ColumnNumber))
 #       plt.axis([0,max(self.yields[:,Yields2Cols['Time']]),0,max(self.yields[:,Yields2Cols['Total_Yields']])])
        plt.xlabel('t in s')
        plt.ylabel('yield in wt%')
        plt.legend()
        plt.grid()
        plt.savefig('Yields_'+self.SpeciesName(ColumnNumber)+'.pdf',format='pdf')
        plt.clf(),plt.cla()

    def plt_RateVsTime(self,ColumnNumber):
        """Plots the original rates output of the pyrolysis program (as e.g. CPD) of the species marke with the columns number"""
        plt.plot(self.__rates[:,0],self.__rates[:,ColumnNumber],label=self.SpeciesName(ColumnNumber))
 #       plt.axis([0,max(self.yields[:,Yields2Cols['Time']]),0,max(self.yields[:,Yields2Cols['Total_Yields']])])
        plt.xlabel('t in s')
        plt.ylabel('rate in wt%/s')
        plt.legend()
        plt.grid()
        plt.savefig('Rates_'+self.SpeciesName(ColumnNumber)+'.pdf',format='pdf')
        plt.clf(),plt.cla()

    def Time(self):
        """Returns the time vector"""
        return self.__yields[:,self.SpeciesIndex('Time')]

    def Yield(self, species):
        """Returns the Vector of the species yield(t). The species can be inputted with the Column number (integer) or the name corresponding to the dictionary saved in the result class (string)."""
        if type(species)==int:
            return self.__yields[:,(species)]
        if type(species)==str:
            return self.__yields[:,self.SpeciesIndex(species)]
        #elif species != 'Solid' and species != self.Yields2Cols['Solid']:
        #    return self.__yields[:,self.SpeciesIndex(species)]
        #elif species == 'Solid' or species == self.Yields2Cols['Solid']:
        #    return self.MassVM_s()

    def Rate(self, species):
        """Returns the Vector of the species rate(t). The species can be inputted with the Column number (integer) or the name corresponding to the dictionary saved in the result class (string)."""
        if type(species)==int:
            return self.__rates[:,(species)]
        if type(species)==str:
            return self.__rates[:,self.SpeciesIndex(species)]
        #elif species != 'Solid' and species != self.Yields2Cols('Solid'):
        #    return self.__rates[:,self.SpeciesIndex(species)]
        #elif species == 'Solid' or species == self.Yields2Cols('Solid'):
        #    return self.RateSingleSpec('Solid')

    def MassCoal(self):
        """returns the Vector with the solid mass(t)"""
        return self.Yield('Solid')

    def MassVM_s(self):
        """Returns the Vector with the mass of the volatile Matter over time"""
        wholeCoalMass= self.Yield('Solid')
        VM_s=wholeCoalMass[:]-wholeCoalMass[-1]             #makes VM_{solid} without char and ash at the time steps
        print VM_s
        return VM_s

    def NPoints(self):
        """returns number of Points for each species over time. Is equal the number of time points."""
        number=len(self.__yields[:,self.SpeciesIndex('Time')])
        return int(number)
        
    def RateSingleSpec(self,NameSpecies):
        """Returns the Rate of the species (inputted as string) by calculate it from the yields by using a CDS"""
        u=self.__yields[:,self.SpeciesIndex(NameSpecies)]
        dt=self.Dt()
        uDot=np.zeros(self.NPoints())
        uDot[0]=(u[1]-u[0])/dt[0]
        uDot[1:-1]=(u[2:]-u[:-2])/(2*dt[1:-1])
        uDot[-1]=(u[-1]-u[-2])/dt[-1]
        return uDot

    
    def SpeciesName(self,ColumnNumber):
        """Returns the species name (string) of the recieved column number (integer)"""
        return self.Cols2Yields[ColumnNumber]
    
    def SpeciesNames(self):
        """Returns a list with all species names (including time and temperature)."""
        ListOfSpeciesNames=[]
        for i in range(len(self.Yields2Cols)):
            ListOfSpeciesNames.append(self.SpeciesName(int(i)))
        return ListOfSpeciesNames

    def SpeciesIndex(self,species):
        """Returns the species column number (integer) of the recieved species name (string)"""
        return self.Yields2Cols[species]

    def Dt(self):                                #between the original points
        """Returns the vector with the time steps dt."""
        Dt=np.zeros(self.NPoints())
        t=self.Time()
        Dt[0]=t[1]-t[0]
        Dt[1:-1]=(t[2:]-t[:-2])/2.
        Dt[-1]=t[-1]-t[-2]
        return Dt

    def DtC(self):                            #for central points between the original points
        """Returns the vector with the time steps dt_C. This time steps are for points between the original points, so the lenght of this vector is the lenght of the time vector minus one."""
        Dtc=np.zeros(self.NPoints()-1)
        t=self.Time()
        Dtc=t[1:]-t[0:-1]
        return Dtc

    def Interpolate(self,Species):
        """Outputs the interpolation object Species(time)."""
        #ODE-Solver needs a continous T(t), managed with an interpolation, the T(t) interpolation is extended to a value 10*t with const. T for ODE (requires more time t the end):
        OrderOfTimeInterpolation=1
        t=self.Time()
        T= self.Yield(Species)
        self.T_interpol=sp.interpolate.interp1d(np.array(list(t)+[10*t[-1]]),np.array(list(T)+[T[-1]]), kind=OrderOfTimeInterpolation, axis=-1, copy=True, bounds_error=True,fill_value=np.nan)
        return self.T_interpol
        

    def LineNumberMaxRate(self,Species):
        """Returns the line with the maximum Rate of the inputted species."""
        u=self.Rate(Species)
        Line=np.argmax(u)
        return Line
    
    def Name(self):
        """Returns the Name of the imported Result object (e.g. 'CPD')"""
        return self.__ImportedResultObject.Name()
        
###################################################################################################################################


class TwoPointEstimator(object):
    """Solves the devolatilization reaction analytically using two arbitrary selected points and the constant rate model. Unprecise. Should only be used for tests."""
    def __init__(self):
        print '2 Point initialized'

    def estimate(self,fgdvc,Species,PointLocation):
    #first part: calculates a list of VM_s and the precursor values k and VM(0) for VM as one species##########################
        t=fgdvc.Time()
        TimePoint=int(PointLocation*(fgdvc.NPoints()))
        if species == 'Solid' or species == self.Yields2Cols('Solid'):
            VM_s=fgdvc.MassVM_s()
            VM_0=VM_s[0]
            k_VM=(np.log(VM_s[TimePoint]/VM_0))/(-t[TimePoint])
        else:
            u=fgdvc.Yield(Species)
            u_0=u[-1]
            k_VM=(np.log(1- u[TimePoint]/u_0))/(-t[TimePoint])
        return k_VM
    

class LeastSquarsEstimator(object):
    """Optimizes the Fitting curve using the Least Squares for the Yields and the Rates."""
    def __init__(self):
        print 'Least Square initialized'
        self.Pre_Tolerance=1.e-5   #Tolerance used for the self.improve_? functions
        self.Fit_Tolerance=1.e-10  #Tolerance for main run
        self.PreMaxIter=None       #Number of maximum iterations used for the self.improve_? functions
        self.MaxIter=None          #Number of maximum iterations for main run

    def improve_E(self,fgdvc,model,t,T,Parameter_Vector,Name,):
        """Additional option: Only the Activation Energy in the Arrhenius Equation is optimized. Actual not necessary."""
        def E_adjust(ActivationEnergie):
            model.setParamVector([Parameter_Vector[0],Parameter_Vector[1],ActivationEnergie])
            print model.ParamVector()
            uDot=fgdvc.Rate(Name)
            vDot=model.calcRate(fgdvc,t,T,Name)
            return sum((uDot-vDot)**2)
        IC_ActEner=Parameter_Vector[2] #InitialCondition for Optimzation
        #Optimizing
        Optimized_E1=fmin(E_adjust,IC_ActEner,ftol=self.Pre_Tolerance,maxiter=self.PreMaxIter)  #shifts the rates curves by the Activation Energy together
        #reformate, because Vec has shape: [0.1, 5.0, array([ 26468.75])]
        Optimized_AE=Optimized_E1[0]
        return Optimized_AE
       

    def improve_a(self,fgdvc,model,t,T,Parameter_Vector,Name):
        """Additional option: Only the preexponential factor in the Arrhenius Equation is optimized. Actual not necessary."""
        IC_a=Parameter_Vector
        def a_adjust(aPre):
            model.setParamVector([aPre,Parameter_Vector[1],Parameter_Vector[2]])
            print model.ParamVector()
            uDot=fgdvc.Rate(Name)
            vDot=model.calcRate(fgdvc,t,T,Name)
            return sum((uDot-vDot)**2)
        Optimized_a1=fmin(a_adjust,IC_a[0],ftol=self.Pre_Tolerance,maxiter=self.PreMaxIter)
        Optimized_a=Optimized_a1[0]
        return Optimized_a

    def estimate_T(self,fgdvc,model,Parameter_Vector,Name,preLoopNumber=0):
        """The main optimization method. Optimizes the Fitting curve using the Least Squares for the weighted Yields and the weighted Rates considering the temperatur history. Requires at input: The corresponding Fit_one_run object, the Model object, the kinetic parameter list, a name (e.g. the species). preLoopNumber is the number of running the  improve_E and improve_a routines. So the standard setting of preLoopNumber is equal zero. It may be used if there is only a very bad convergence when optimize all three parameter."""
        t=fgdvc.Time()
        Dt=fgdvc.Dt()
        DtC=fgdvc.DtC()
        T=fgdvc.Interpolate('Temp')
        u=fgdvc.Yield(Name)
        #updated Vector: CurrentVector
        CurrentVector=Parameter_Vector
        #!!preloopnumber should be zero if optimizer is 'leastsq' or the model is a constant rate model, only needed for the other optimizer
        for i in range(preLoopNumber):
            #get better value for activation Energy:
            Optimized_E=self.improve_E(fgdvc,model,t,T,CurrentVector,Name)
            model.setParamVector([CurrentVector[0],CurrentVector[1],Optimized_E])
            CurrentVector=model.ParamVector()
            #get better value for a:
            Optimized_a=self.improve_a(fgdvc,model,t,T,CurrentVector,Name)
            CurrentVector=[Optimized_a,CurrentVector[1],Optimized_E]
        #####improve all parameters
        w0=self.a0/(( max((fgdvc.Yield(Name))) -min((fgdvc.Yield(Name))) )**2)
        print 'a0', self.a0
        print 'w0 ', w0
        w1=self.a1/(max( ((fgdvc.Rate(Name)))**2 ))
        print 'a1', self.a1
        print 'w1 ',w1
        
        def LeastSquareFunction(Parameter):
            model.setParamVector(Parameter)
            print model.ParamVector()
            v=model.calcMass(fgdvc,t,T,Name)
            #uDot=(u[1:]-u[0:-1])/DtC #rausnehmen, ausserhalb der fkt?
            uDot=fgdvc.Rate(Name)#model.deriveC(fgdvc,u)
            vDot=model.deriveC(fgdvc,v)
            #vDot=(v[1:]-v[0:-1])/DtC
            if self.selectedOptimizer=='leastsq':
                Dot1=np.array(len(u),dtype='d'); Dot2=np.array(len(uDot),dtype='d')
                Dot2=w1*(((uDot-vDot)**2)*Dt) #the rate term
                Dot1=w0*(((u-v)**2)*Dt)        #the yield term
                Error=np.zeros(len(Dot1)+len(Dot2),dtype='d')
                #makes a long array, containing both, the rates and yields
                Error[0:len(Dot1)]=Dot1
                Error[len(Dot1):]=Dot2
                print np.sum(Error)
            else:
                SumYields=np.sum( ((u-v)**2)*Dt )
                #SumRates=np.sum( ((uDot-vDot)**2)*DtC )
                SumRates=np.sum( ((uDot-vDot)**2)*Dt )
                Error= w0*SumYields+w1*SumRates
                print Error
            return Error
        model.setParamVector(CurrentVector)
        if self.selectedOptimizer=='fmin':
            OptimizedVector=fmin(LeastSquareFunction,model.ParamVector(),ftol=self.Fit_Tolerance,maxiter=self.MaxIter)
            return OptimizedVector
        elif self.selectedOptimizer=='fmin_cg':
            OptimizedVector=fmin_cg(LeastSquareFunction,model.ParamVector(),gtol=self.Fit_Tolerance,maxiter=self.MaxIter)
            return OptimizedVector
        elif self.selectedOptimizer=='fmin_bfgs':
            OptimizedVector=fmin_bfgs(LeastSquareFunction,model.ParamVector(),gtol=self.Fit_Tolerance,maxiter=self.MaxIter)
            return OptimizedVector
        elif self.selectedOptimizer=='fmin_ncg':
            OptimizedVector=fmin_ncg(LeastSquareFunction,model.ParamVector(),avextol=self.Fit_Tolerance)#,maxiter=self.MaxIter)
            return OptimizedVector
        elif self.selectedOptimizer=='fmin_slsqp':
            OptimizedVector=fmin_slsqp(LeastSquareFunction,model.ParamVector(),acc=self.Fit_Tolerance)#,maxiter=self.MaxIter)
            return OptimizedVector
        elif self.selectedOptimizer=='leastsq':
            OptimizedVector=leastsq(LeastSquareFunction,model.ParamVector(),ftol=self.Fit_Tolerance)
            return OptimizedVector[0]
        else:
            print "No Optimizer was selected. Please choose: 'fmin' or 'fmin_cg' or 'fmin_bfgs' or 'leastsq' or 'fmin_slsqp'\n"
        print  'Optm Vec:   ',OptimizedVector
        


    def setWeights(self,WeightMass,WeightRates):
        """Stes the weights for the yields and the rates for the fitting procedure. See manual for equation."""
        self.a0=WeightMass
        self.a1=WeightRates
        
    def setTolerance(self, ToleranceForFminFunction):
        """Sets the tolerance as a abort criterion for the fitting procedure."""
        self.Fit_Tolerance=ToleranceForFminFunction

    def setPreTolerance(self, ToleranceForFminFunction):
        """Sets the tolerance as a abort criterion for the prefitting procedure (if preLoopNumber in estimate_T is not equal zero)."""
        self.Pre_Tolerance=ToleranceForFminFunction

    def setPreMaxIter(self,MaxiumNumberOfIterationInPreProcedure):
        """Sets the maximum number of iteration oin the optimizer as a abort criterion for the prefitting procedure (if preLoopNumber in estimate_T is not equal zero)."""
        self.PreMaxIter=int(MaxiumNumberOfIterationInPreProcedure)

    def setMaxIter(self,MaxiumNumberOfIterationInMainProcedure):
        """Sets the maximum number of iteration oin the optimizer as a abort criterion for the fitting procedure."""
        self.MaxIter=int(MaxiumNumberOfIterationInMainProcedure)

    def setOptimizer(self,ChosenOptimizer):
        """Select one optimizer of the scipy.optimizer library: 'fmin','fmin_cg','fmin_bfgs','fmin_ncg','fmin_slsqp' or 'leastsq'. According to experience 'fmin' (or also 'leastsq') generates at best the results."""
        self.selectedOptimizer=ChosenOptimizer


###################################################################################################################################    
    
class Model(object):
    """Parent class of the children ConstantRateModel, the three Arrhenius Models (notations) and the Kobayashi models."""
    def __init__(self,Name):
        print 'Model initialized'
        
    def pltYield(self,fgdvc,xValueToPlot,yValueToPlot):
        """Plots the yields (to select with yValueToPlot) over Time or Temperature (to slect with xValueToPlot)."""
        plt.plot(fgdvc.Yield(xValueToPlot),fgdvc.Yield(yValueToPlot))
        if xValueToPlot=='Time':
            plt.xlabel('t in s')
        if xValueToPlot=='Temp':
            plt.xlabel('T in K')
        if type(yValueToPlot)==int:
            SpeciesForTitle=fgdvc.SpeciesName(yValueToPlot)
        if type(yValueToPlot)==str:
            SpeciesForTitle=yValueToPlot
        plt.title(SpeciesForTitle)
        plt.ylabel('yield in wt%')
        plt.legend()
        plt.grid()
        plt.savefig('Yields_'+yValueToPlot+'VS'+xValueToPlot+'.pdf',format='pdf')
        plt.clf(),plt.cla()

    def pltRate(self,fgdvc,xValueToPlot,yValueToPlot):
        """Plots the rates (to select with yValueToPlot) over Time or Temperature (to slect with xValueToPlot)."""
        plt.plot(fgdvc.Rate(fgdvc.SpeciesIndex(xValueToPlot)),fgdvc.Rate(fgdvc.SpeciesIndex(yValueToPlot)),label=yValueToPlot)
        if xValueToPlot=='Time':
            plt.xlabel('t in s')
        if xValueToPlot=='Temp':
            plt.xlabel('T in K')
        plt.ylabel('rate in wt%/s')
        if type(yValueToPlot)==int:
            SpeciesForTitle=fgdvc.SpeciesName(yValueToPlot)
        if type(yValueToPlot)==str:
            SpeciesForTitle=yValueToPlot
        plt.title(SpeciesForTitle)
        plt.legend()
        plt.grid()
        plt.savefig('Rates_'+yValueToPlot+'VS'+xValueToPlot+'.pdf',format='pdf')
        plt.clf(),plt.cla()

    def plot(self,fgdvc,Species):
        """Plot the yield and the rates over time with two curves: one is the original data, the other the fitting curve. Also file 'PyrolysisProgramName-Species.out' (e.g. 'CPD-CO2.out') containing the time (s), yields (kg/kg), rates (kg/(kg s))."""
        #plots:
        #Yields to compare
        u=fgdvc.Yield(Species)
        v=self.calcMass(fgdvc,fgdvc.Time(),fgdvc.Interpolate('Temp'),Species)
        if type(Species)==int:
            SpeciesForTitle=fgdvc.SpeciesName(Species)
        if type(Species)==str:
            SpeciesForTitle=Species
        if SpeciesForTitle=='Solid':
            plt.plot(fgdvc.Time(),u,label=fgdvc.Name())
            plt.plot(fgdvc.Time(),v,label='fit')
            plt.plot(fgdvc.Time(),(1.-u[:]),label='Sum yields')
            plt.plot(fgdvc.Time(),(1.-v[:]),label='fit')
        else:
            plt.plot(fgdvc.Time(),u,label=fgdvc.Name())
            plt.plot(fgdvc.Time(),v,label='fit')            
        plt.title(SpeciesForTitle)
        plt.xlabel('t in s')
        plt.ylabel('yield fraction in kg/kg_coal')
        plt.legend()
        plt.grid()
        plt.savefig(fgdvc.Name()+'-Fit_result_'+SpeciesForTitle+'_Y.pdf',format='pdf')
        plt.clf(),plt.cla()
        #Rates to compare
        u=fgdvc.Rate(Species)
        plt.plot(fgdvc.Time(),u,label=fgdvc.Name())
        w=self.deriveC(fgdvc,v)
        plt.plot(fgdvc.Time(),w,label='fit')
        if type(Species)==int:
            SpeciesForTitle=fgdvc.SpeciesName(Species)
        if type(Species)==str:
            SpeciesForTitle=Species
        plt.title(SpeciesForTitle)
        plt.xlabel('t in s')
        plt.ylabel('rate in 1/s')#min')
        plt.legend()
        plt.grid()
        plt.savefig(fgdvc.Name()+'-Fit_result_'+SpeciesForTitle+'_R.pdf',format='pdf')
        plt.clf(),plt.cla()
        #writes result file
        t=fgdvc.Yield('Time')
        T=fgdvc.Yield('Temp')
        resultFile=open(fgdvc.Name()+'-Fit_result_'+SpeciesForTitle+'.out','w')
        resultFile.write('    Time       Temperature    Yields       Rates \n')
        for i in range(len(t)):
            resultFile.write('%7e  %11e %7e %8e \n' % (t[i], T[i], v[i], w[i]))
        resultFile.close()

    def deriveC(self,fgdvc,yVector):
        """Returns a CDS of the inputted yVector."""
        dt=fgdvc.Dt()
        yDot=np.zeros(fgdvc.NPoints())
        yDot[0]=(yVector[1]-yVector[0])/dt[0]
        yDot[1:-1]=(yVector[2:]-yVector[:-2])/(2*dt[1:-1])
        yDot[-1]=(yVector[-1]-yVector[-2])/dt[-1]
        return yDot
        
    def calcRate(self,fgdvc,time,T,Name):
        """Generates the Rates using the yields vector and a CDS."""
        m=self.calcMass(fgdvc,time,T,Name)
        mDot=self.deriveC(fgdvc,m)
        return mDot

    def setParamVector(self,ParameterList):
        """Sets the Vector containing the kinetic parameter of the Model (refering to the child model)."""
        self._ParamVector=ParameterList

    def ParamVector(self):
        """Returns the Vector containing the kinetic parameter of the Model (refering to the child model)."""
        return self._ParamVector

    def ErrorYield(self,fgdvc,Species):                                  #calculates avrg Error between curves
        """Returns the absolute deviation per point between the fitted and the original yield curve."""
        v=self.calcMass(fgdvc,fgdvc.Time(),fgdvc.Interpolate('Temp'),Species)
        u=fgdvc.Yield(Species)
        absE=(np.sum(u-v))/fgdvc.NPoints()
        return absE

    def ErrorRate(self,fgdvc,Species):                                  #calculates avrg Error between curves
        """Returns the absolute deviation per point between the fitted and the original rate curve."""
        v=self.calcMass(fgdvc,fgdvc.Time(),fgdvc.Interpolate('Temp'),Species)
        uDot=fgdvc.Rate(Species)
        vDot=self.deriveC(fgdvc,v)
        absE=(np.sum(uDot-vDot))/fgdvc.NPoints()
        return absE
    
    def alpha_d(self,fgdvc):
        """returns alpha_d which is the final yield of tar devided by amount of volatile Matter."""
        uTar=fgdvc.Yield('Tar')
        uVM=1.-fgdvc.Yield('Solid')
        return uTar[-1]/uVM
    
    def beta_d(self,fgdvc,species):
        """returns beta_d which is the final yield of a species devided by amount of volatile Matter."""
        u=fgdvc.Yield(species)
        uVM=1.-fgdvc.Yield('Solid')
        return u[-1]/uVM


################subclasses###########################

class ConstantRateModel(Model):
    """The model calculating the mass with m(t)=m_s0+(m_s0-m_s,e)*e**(-k*(t-t_start)) from the ODE dm/dt = -k*(m-m_s,e). The Parameter to optimize are k and t_start."""
    def __init__(self,InitialParameterVector):
        print 'Constant rate initialized'
        self._ParamVector=InitialParameterVector

    def calcMass(self,fgdvc,t,T,SpeciesToCalc):
        ParamVector=self.ParamVector()
        u=fgdvc.Yield(SpeciesToCalc)
        u_0=u[0]; u_s0=u[-1]
        if SpeciesToCalc=='Solid' or SpeciesToCalc==int(fgdvc.Yields2Cols['Solid']):
            v=u_s0+(u_0-u_s0)*np.exp(-ParamVector[0]*(t-ParamVector[1]))
        else:
            v=u_s0*(1.-np.exp(-ParamVector[0]*(t-ParamVector[1])))
        v=np.where(t>ParamVector[1],v,u_0)
        return v

                
class ArrheniusModel(Model):
    """The Arrhenius model in the standart notation: dm/dt=A*(T**b)*exp(-E/T)*(m_s-m) with the parameter a,b,E to optimize."""
    def __init__(self,InitialParameterVector):
        print 'Arrhenuis Model initialized'
        self._ParamVector=InitialParameterVector
        self.ODE_hmax=1.e-2
 
    def calcMass(self,fgdvc,time,T,Name):
        """Outputs the mass(t) using the model specific equation."""
        #numercal values:
        absoluteTolerance = 1.0e-8
        relativeTolerance = 1.0e-6
        ##################
        u=fgdvc.Yield(Name)
        m_s0=u[-1]
        def dmdt(m,t):
            #A=parameter[0]
            #b=parameter[1]
            #E=parameter[2]
            ParamVec=self.ParamVector()
            if Name == 'Solid' or Name == fgdvc.Yields2Cols['Solid']:
                dmdt_out=-ParamVec[0]*(T(t)**ParamVec[1])*np.exp(-ParamVec[2]/(T(t)))*(m-m_s0)
                dmdt_out=np.where(abs(dmdt_out)>1.e-300,dmdt_out,0.0) #sets values<0 =0.0, otherwise it will further cause problems (nan)
            else:
                dmdt_out= ParamVec[0]*(T(t)**ParamVec[1])*np.exp(-ParamVec[2]/(T(t)))*(m_s0-m)
                #print 'ParamVec[0]',ParamVec[0]
                #print '(T(t)**ParamVec[1])',(T(t)**ParamVec[1])
                #print 'np.exp(-ParamVec[2]/(T(t)))', np.exp(-ParamVec[2]/(T(t)))
                #print '(m_s0-m)',(m_s0-m)
            return dmdt_out
        InitialCondition=[u[0]]
        m_out=sp.integrate.odeint(dmdt,InitialCondition,time,atol=absoluteTolerance,rtol=relativeTolerance,hmax=self.ODE_hmax)
        return m_out[:,0]


    def ConvertKinFactors(self,ParameterVector):
        """Dummy. Function actual has to convert the parameter into the standart Arrhenius notation."""
        #does nothing, just to have the same way of use for all notations
        return [ParameterVector[0],ParameterVector[1],ParameterVector[2]]

        
        
class ArrheniusModelAlternativeNotation1(ArrheniusModel):
    """Arrhenius model with a notation having a better optimization behaviour: dm/dt=exp[k0-a*(T0/T(t)-1)]*(ms-m). See the documentation for the reference. The parameters to optimize are k0 and a."""
    def __init__(self,InitialParameterVector):
        print 'Arrhenuis Model class initialized'
        self._ParamVector=InitialParameterVector
        self.ODE_hmax=1.e-2
        
    def calcMass(self,fgdvc,time,T,Name):
        """Outputs the mass(t) using the model specific equation."""
        #numercal values:
        absoluteTolerance = 1.0e-8
        relativeTolerance = 1.0e-6
        ##################
        ParamVec=self.ParamVector()
        u=fgdvc.Yield(Name)
        m_s0=u[-1]
        twhereTmax=time[fgdvc.LineNumberMaxRate(Name)]
        self.T0=T(twhereTmax)
        def dmdt(m,t):
            #k0=parameter[0]
            #n=0
            #a=parameter[1]
            if Name == 'Solid' or Name == fgdvc.Yields2Cols['Solid']:
                dmdt_out=-np.exp( ParamVec[0]  - ParamVec[1]*( self.T0/T(t) - 1 ))*(m-m_s0) #-ParamVec[0]*( (T(t)/self.T0)**ParamVec[1] )*np.exp( -ParamVec[2]*(self.T0/(T(t))-1) )*(m)    #IC
                dmdt_out=np.where(abs(dmdt_out)>1.e-300,dmdt_out,0.0) #sets values<0 =0.0, otherwise it will further cuase problem(nan)
            else:
                dmdt_out= np.exp( ParamVec[0] - ParamVec[1]*( self.T0/T(t) - 1 ))*(m_s0-m)
            return dmdt_out
        InitialCondition=[u[0]]
        m_out=sp.integrate.odeint(dmdt,InitialCondition,time,atol=absoluteTolerance,rtol=relativeTolerance)
        return m_out[:,0]

    def ConvertKinFactors(self,ParameterVector):
        """Converts the own kinetic factors back to the standard Arrhenius kinetic factors."""
        #k=ParameterVector[0]
        #a=ParameterVector[1]
        A=np.exp( ParameterVector[0] + ParameterVector[1] )
        E=ParameterVector[1]*self.T0
        return [A,0.0,E]

    def ConvertKinFactorsToOwnNotation(self,fgdvc,ParameterVector,Species):
        """Converts the standard Arrhenius kinetic factors backk to the factors of the own notation."""
        time=fgdvc.Time()
        twhereTRmax=time[fgdvc.LineNumberMaxRate(Species)]
        T=fgdvc.Interpolate('Temp')
        self.T0=T(twhereTRmax)
        A=ParameterVector[0]
        if ParameterVector[1]!=0:
            print 'b is not considered in this Notation and set to 0'
        E=ParameterVector[2]
        a=E/self.T0
        k0=np.log(A)-a
        return [k0,a]

        
class ArrheniusModelAlternativeNotation2(ArrheniusModel):
    """Arrhenius model with a notation having a better optimization behaviour: dm/dt=exp[c*(b1*(1/T(t)-1/T_min)-b2*(1/T(t)-1/T_max))]*(ms-m); with c=(1/T_max-1/Tmin)**(-1). See the documentation for the reference. The parameters to optimize are b1 and b2."""
    def __init__(self,fgdvc,InitialParameterVector):
        print 'Arrhenius Model class initialized'
        self._ParamVector=InitialParameterVector
        T_general=fgdvc.Rate('Temp')
        self.T_min=T_general[0]
        self.T_max=T_general[-1]
        self.c=1./(1./self.T_max-1./self.T_min)
        self.ODE_hmax=1.e-2
        self.fgdvc=fgdvc
        
    def calcMass(self,fgdvc,time,T,Name):
        """Outputs the mass(t) using the model specific equation."""
        #numercal values:
        absoluteTolerance = 1.0e-8
        relativeTolerance = 1.0e-6
        ##################
        ParamVec=self.ParamVector()
        u=fgdvc.Yield(Name)
        #uDot=fgdvc.Rate(Name)
        m_s0=u[-1]
        def dmdt(m,t):
            #b1=parameter[0]
            #b2=parameter[1]
            if Name == 'Solid' or Name == fgdvc.Yields2Cols['Solid']:
                dmdt_out= -np.exp( self.c*( ParamVec[0]*(1./T(t)-1./self.T_min) - ParamVec[1]*(1./T(t)-1./self.T_max) ) )*(m-m_s0)
                dmdt_out=np.where(abs(dmdt_out)>1.e-300,dmdt_out,0.0) #sets values<0 =0.0, otherwise it will further cause problem(nan)
            else:
                dmdt_out= np.exp( self.c*( ParamVec[0]*(1./T(t)-1./self.T_min) - ParamVec[1]*(1./T(t)-1./self.T_max) ) )*(m_s0-m)
            return dmdt_out
        InitialCondition=[u[0]]
        m_out=sp.integrate.odeint(dmdt,InitialCondition,time,atol=absoluteTolerance,rtol=relativeTolerance,hmax=self.ODE_hmax)
        return m_out[:,0]

    def ConvertKinFactors(self,ParameterVector):
        """Converts the own kinetic factors back to the standard Arrhenius kinetic factors."""
        #b1=ParameterVector[0]
        #b2=ParameterVector[1]
        A=np.exp( -self.c*ParameterVector[0]/self.T_min + self.c*ParameterVector[1]/self.T_max )
        E=self.c*ParameterVector[1]-self.c*ParameterVector[0]
        return [A,0,E]

    def ConvertKinFactorsToOwnNotation(self,fgdvc,ParameterVector):
        """Converts the standard Arrhenius kinetic factors backk to the factors of the own notation."""
        T_general=fgdvc.Rate('Temp')
        self.T_min=T_general[0]    #maybe change later
        self.T_max=T_general[-1]
        self.c=1/(1/self.T_max-1/self.T_min)
        #A=ParameterVector[0]
        #b not used
        #E=ParameterVector[2]
        if ParameterVector[1]!=0:
            print 'b is not considered in this Notation and set to 0'
        b2=( (self.T_max/self.c)*np.log(ParameterVector[0])-(ParameterVector[2]*self.T_max)/(self.c*self.T_min) )*(1-self.T_max/self.T_min)**(-1)
        b1=b2-ParameterVector[2]/self.c
        return [b1,b2]

class Kobayashi(Model):
    """Calculates the devolatilization reaction using the Kobayashi model. The Arrhenius equation inside are in the standard notation."""
    def __init__(self,fgdvc,InitialParameterVector):
        print 'Kobayashi Model initialized'
        self._ParamVector=InitialParameterVector
        self.ODE_hmax=1.e-2
        self.fgdvc=fgdvc
        
    def calcMass(self,fgdvc,time,T,Name):
        """Outputs the mass(t) using the model specific equation."""
        #alpha has to be greater equal zero:
        if self._ParamVector[4]<0:
            self._ParamVector[4]=0
        if self._ParamVector[5]<0:
            self._ParamVector[5]=0
        absoluteTolerance = 1.0e-8
        relativeTolerance = 1.0e-6
        self.tList=[0.0]
        #deletes to solve trapezian rule:
        time=np.delete(time,0)
        self.Integral=0.0
        self.k1k2=[0.0]
        ParamVec=self.ParamVector()
        u=fgdvc.Yield(Name)
        #
        def dmdt(m,t):
            k1=ParamVec[0]*np.exp(-ParamVec[1]/(T(t)))
            k2=ParamVec[2]*np.exp(-ParamVec[3]/(T(t)))
            self.tList.append(t)
            self.k1k2.append(k1+k2)
            self.Integral+=0.5*(self.tList[-1]-self.tList[-2])*(self.k1k2[-1]+self.k1k2[-2])
            dmdt_out = ( (ParamVec[4]*k1+ParamVec[5]*k2)*np.exp(-self.Integral) )
            dmdt_out=np.where(abs(dmdt_out)>1.e-300,dmdt_out,0.0) #sets values<0 =0.0, otherwise it will further cause problem(nan)
            return dmdt_out
        InitialCondition=[u[0]]
        m_out=sp.integrate.odeint(dmdt,InitialCondition,time,atol=absoluteTolerance,rtol=relativeTolerance,hmax=self.ODE_hmax)
        m_out=m_out[:,0]
        m_out=np.insert(m_out,0,0.0)
        return m_out

    def ConvertKinFactors(self,ParameterVector):
        """Outputs the Arrhenius equation factors in the shape [A1,E1,A2,E2,alpha1,alpha2]. Here where the real Arrhenius model is in use only a dummy function."""
        return ParameterVector

class KobayashiA2(Model):
    """Calculates the devolatilization reaction using the Kobayashi model. The Arrhenius equation inside are in the secend alternative notation (see class ArrheniusModelAlternativeNotation2)."""
    def __init__(self,fgdvc,InitialParameterVector):
        print 'Kobayashi Model initialized'
        self._ParamVector=InitialParameterVector
        T_general=fgdvc.Rate('Temp')
        self.T_min=min(T_general)
        self.T_max=max(T_general)
        self.c=1./(1./self.T_max-1./self.T_min)
        self.ODE_hmax=1.e-2
        self.fgdvc=fgdvc

    def calcMass(self,fgdvc,time,T,Name):
        """Outputs the mass(t) using the model specific equation."""
        #alpha has to be greater equal zero:
        if self._ParamVector[4]<0:
            self._ParamVector[4]=0
        if self._ParamVector[5]<0:
            self._ParamVector[5]=0
        absoluteTolerance = 1.0e-8
        relativeTolerance = 1.0e-6
        self.tList=[0.0]
        #deletes to solve trapezian rule:
        time=np.delete(time,0)
        self.Integral=0.0
        self.k1k2=[0.0]
        ParamVec=self.ParamVector()
        u=fgdvc.Yield(Name)
        #
        def dmdt(m,t):
            k1=np.exp( self.c*( ParamVec[0]*(1./T(t)-1./self.T_min) - ParamVec[1]*(1./T(t)-1./self.T_max) ) )
            k2=np.exp( self.c*( ParamVec[2]*(1./T(t)-1./self.T_min) - ParamVec[3]*(1./T(t)-1./self.T_max) ) )
            self.tList.append(t)
            self.k1k2.append(k1+k2)
            self.Integral+=0.5*(self.tList[-1]-self.tList[-2])*(self.k1k2[-1]+self.k1k2[-2])
            dmdt_out = ( (ParamVec[4]*k1+ParamVec[5]*k2)*np.exp(-self.Integral) )
            dmdt_out=np.where(abs(dmdt_out)>1.e-300,dmdt_out,0.0) #sets values<0 =0.0, otherwise it will further cause problem(nan)
            return dmdt_out
        InitialCondition=[u[0]]
        m_out=sp.integrate.odeint(dmdt,InitialCondition,time,atol=absoluteTolerance,rtol=relativeTolerance,hmax=self.ODE_hmax)
        m_out=m_out[:,0]
        m_out=np.insert(m_out,0,0.0)
        return m_out

    def ConvertKinFactors(self,ParameterVector):
        """Converts the alternative notaion Arrhenius factors into the satndard Arrhenius factors and return them in the shape  [A1,E1], [A2,E2]"""
        A1=np.exp( -self.c*ParameterVector[0]/self.T_min + self.c*ParameterVector[1]/self.T_max )
        E1=self.c*ParameterVector[1]-self.c*ParameterVector[0]
        A2=np.exp( -self.c*ParameterVector[2]/self.T_min + self.c*ParameterVector[3]/self.T_max )
        E2=self.c*ParameterVector[3]-self.c*ParameterVector[2]
        return [A1,E1,A2,E2,ParameterVector[4],ParameterVector[5]]
