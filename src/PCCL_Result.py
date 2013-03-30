#Modified for FG-DVC 8.2.2 and 8.2.3
#for further Versions see compatibility of the output files
import numpy as np
import sys
import scipy.interpolate
#
#
class PCCL_Result(object):
    """Reads the PC Coal Lab input and saves the values in one array. The results include the yields of the species. The rates for all species are calculated using a CDS. This class also contains the dictionaries for the columns in the array - the name of the species. These dictionaries might be PC Coal Lab - Version dependent and the only thing which has to be changed for the case of a new release of FG-DVC with a new order of species in the result files (this was made for Version 4.1)."""
    def __init__(self,FilePath,NrOfRun,dt,OrderOfInterpolation=1): 
        if float(NrOfRun) > 5:
            print 'Only a maximum nr of Runs of five is possible for PC Coal Lab.'
            sys.exit()
        self.__path = FilePath
        self.__nrRun = str(NrOfRun)
        WTyields = np.genfromtxt(self.__path+'FDC1WT'+self.__nrRun+'.RPT',skip_header=28)
        NGyields = np.genfromtxt(self.__path+'FDC1NG'+self.__nrRun+'.RPT',skip_header=28)
        HCyields = np.genfromtxt(self.__path+'FDC1HC'+self.__nrRun+'.RPT',skip_header=28)
        NSpec = 18 #number of reported species
        # transform this into an own array
        self.__yields = np.zeros([len(WTyields[:,0]),NSpec])
        self.__yields[:,0]  = WTyields[:,0] # time
        self.__yields[:,1]  = WTyields[:,1]+273. # temperature, T in K instead degree C
        self.__yields[:,2]  = WTyields[:,2]/100. # total yields
        self.__yields[:,3]  = WTyields[:,3] # gas
        self.__yields[:,4]  = WTyields[:,4]/100. # tar
        self.__yields[:,5]  = WTyields[:,5]/100. # char
        self.__yields[:,6]  = NGyields[:,4]/100. # CO2
        self.__yields[:,7]  = NGyields[:,5]/100. # H2O
        self.__yields[:,8]  = NGyields[:,6]/100. # CO
        self.__yields[:,9]  = NGyields[:,7]/100. # Hydrocarbons
        self.__yields[:,10] = NGyields[:,8]/100. # HCN
        self.__yields[:,11] = HCyields[:,4]/100. # CH4
        self.__yields[:,12] = HCyields[:,5]/100. # C2H4
        self.__yields[:,13] = HCyields[:,6]/100. # C2H6
        self.__yields[:,14] = HCyields[:,7]/100. # C3H6
        self.__yields[:,15] = HCyields[:,8]/100. # C3H8
        self.__yields[:,16] = HCyields[:,9]/100. # H2
        self.__yields[:,17] = HCyields[:,10]/100. # H2S
        # Yields2Cols
        self.Yields2Cols={'Time':0,'Temp':1,'Total':2,'Gas':3,'Tar':4,'Char':5,'CO2':6,'H2O':7,'CO':8,'Hydrocarbons':9,'HCN':10,'CH4':11,'C2H4':12,'C2H6':13,'C3H6':14,'C3H8':15,'H2':16,'H2S':17}
        self.Cols2Yields={0:'Time',1:'Temp',2:'Total',3:'Gas',4:'Tar',5:'Char',6:'CO2',7:'H2O',8:'CO',9:'Hydrocarbons',10:'HCN',11:'CH4',12:'C2H4',13:'C2H6',14:'C3H6',15:'C3H8',16:'H2',17:'H2S'}
        #
        # Interpolation
        # when appending or deleting species: change the array size of self.__yields and self.__rates
#        NPoints = int(self.__yields[-1,0]/dt + 1) # =FinalTime/dt
#        self.__yields = np.zeros([NPoints,NSpec]) # the final array with interpolation points
#        #over the temperature, which is linear ramp or constant value -> first order fir is exact solution
#        Interplt=scipy.interpolate.interp1d(yieldsSmall[:,0],yieldsSmall[:,1], kind=1, axis=-1, copy=True, bounds_error=True,fill_value=np.nan) #interpol Obj
#        for j in range(NPoints): #over all Points
#            self.__yields[j,1] = Interplt(j*dt) # defines Temp
#            self.__yields[j,0] = j*dt           # defines Time
#        for i in range(2,NSpec,1): #over all Species
#            Interplt=scipy.interpolate.interp1d(yieldsSmall[:,0],yieldsSmall[:,i], kind=OrderOfInterpolation, axis=-1, copy=True, bounds_error=True,fill_value=np.nan) #interpol Obj
#            for j in range(NPoints): #over all Points
#                self.__yields[j,i] = Interplt(j*dt)
#        #
        # get the rates:
        self.__rates = np.zeros(np.shape(self.__yields))
        self.__rates[:,0] = self.__yields[:,0]  #Time
        self.__rates[:,1] = self.__yields[:,1]  #Temperature
        #calculates rate of Totals as in the gasrate.txt only the solid yields are reported:
        for i in range(2,NSpec,1):
            self.__rates[0,i]=(self.__yields[1,i]-self.__yields[0,i])/(self.__yields[1,0]-self.__yields[0,0])
            self.__rates[1:-1,i]=(self.__yields[2:,i]-self.__yields[:-2,i])/(self.__yields[2:,0]-self.__yields[:-2,0])
            self.__rates[-1,i]=(self.__yields[-1,i]-self.__yields[-2,i])/(self.__yields[-1,0]-self.__yields[-2,0])
#        #
#        #
#        print '\nimported PC Coal Lab data-fields, size(rows,columns): ',yieldsSmall.shape,'  and interpolated to: ',self.__yields.shape,'\n'
    
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
        return 'FGDVC'

if __name__ == "__main__":
    PR = PCCL_Result('C:\\Users\\MaP\\PCCL\\',1,1e-3,2)
    a = PR.Yields_all()
    print '\n\nTotal',a[:,2]
    print '\n\nTar',a[:,4]
    print '\n\nCO',a[:,8]
    print '\n\nCH4',a[:,11]
    print '\n\nH2',a[:,16]
    PR.Rates_all()
    
    