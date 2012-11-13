#MODIFIED FOR CPD VERSION: CPD_NLG (1999)
import numpy as np
#
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
