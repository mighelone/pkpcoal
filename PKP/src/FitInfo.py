import numpy as np
import pylab as plt
import scipy as sp
import scipy.interpolate
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
        if self.Name()=='CPD':
            return self.Yield('Solid')
        elif self.Name()=='FGDVC':
            return 1.-self.Yield('Total')

    def MassVM_s(self):
        """Returns the Vector with the mass of the volatile Matter over time"""
        wholeCoalMass= self.MassCoal()
        VM_s=wholeCoalMass[:]-wholeCoalMass[-1]             #makes VM_{solid} without char and ash at the time steps
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
        try:
            return self.Cols2Yields[ColumnNumber]
        except KeyError:
            return None
    
    def SpeciesNames(self):
        """Returns a list with all species names (including time and temperature)."""
        ListOfSpeciesNames=[]
        for i in range(len(self.Yields2Cols)):
            ListOfSpeciesNames.append(self.SpeciesName(int(i)))
        return ListOfSpeciesNames

    def SpeciesIndex(self,species):
        """Returns the species column number (integer) of the recieved species name (string)"""
        try:
            return self.Yields2Cols[species]
        except KeyError:
            return None

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
        """Outputs the interpolation object (e.g.Species(time))."""
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
        
