#Modified for FG-DVC 8.2.2 and 8.2.3
#for further Versions see compatibility of the output files
import numpy as np

#
class FGDVC_Result(object):
    """Reads the FG-DVC input and saves the values in one array. The results include the yields (from 'gasyields.txt') and the rates. The rates for all species except the solids (here a CDS is used) are imported from 'gasrates.txt'. The H_2 yields were calculated by subtract all other species except parafins and olefins from the total yields (see FG-DVC manual). This H_2-yield curve was smoothed and derived using a CDS to generate the H_2 rates. The parafins and olefins are added into the tar. This class also contains the dictionaries for the columns in the array - the name of the species. These dictionaries are FG-DVC-Version dependent and the only thing which has to be changed for the case of a new release of FG-DVC with a new order of species in the result files (this was made for Versions 8.2.2. and 8.2.3.)."""
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
        #calculates Yields and rates of H2 from 'Total'(original 'Total')
        self.__yields[:,17]=self.__yields[:,17]-self.__yields[:,3]-self.__yields[:,5]-self.__yields[:,6]-self.__yields[:,7]-self.__yields[:,8]-self.__yields[:,9]-self.__yields[:,10]-self.__yields[:,11]-self.__yields[:,12]-self.__yields[:,13]-self.__yields[:,14]  #converts Total yields into H2 yields, H2 is not reported anywhere else
        self.__rates[:,17]=self.__rates[:,17]-self.__rates[:,3]-self.__rates[:,5]-self.__rates[:,6]-self.__rates[:,7]-self.__rates[:,8]-self.__rates[:,9]-self.__rates[:,10]-self.__rates[:,11]-self.__rates[:,12]-self.__rates[:,13]-self.__rates[:,14]-self.__rates[:,4]
        #make 'Total' yields instead of 'Char+Ash':
        self.__yields[:,2]=1.0-self.__yields[:,2]
        #calculates rate of Totals as in the gasrate.txt only the solid yields are reported:
        self.__rates[0,2]=(self.__yields[1,2]-self.__yields[0,2])/(self.__yields[1,0]-self.__yields[0,0])
        self.__rates[1:-1,2]=(self.__yields[2:,2]-self.__yields[:-2,2])/(self.__yields[2:,0]-self.__yields[:-2,0])
        self.__rates[-1,2]=(self.__yields[-1,2]-self.__yields[-2,2])/(self.__yields[-1,0]-self.__yields[-2,0])
        #calculates the 'gas' yields (are all light gases, no tar,Parafins,Olefins):
        LightGases=np.zeros([len(self.__yields[:,0]),1],'d')
        LightGases[:,0]=self.__yields[:,2]-self.__yields[:,4]-self.__yields[:,15]-self.__yields[:,16]
         #merge both in the Array
        self.__yields=np.concatenate((self.__yields,LightGases),axis=1)
        #calculates the 'gas' rates (are all light gases, no tar):
        LightGases[:,0]=( self.__yields[1,18]-self.__yields[0,18] )/(self.__yields[1,0]-self.__yields[0,0])
        LightGases[1:-1,0]=(self.__yields[2:,18]-self.__yields[:-2,18])/(self.__yields[2:,0]-self.__yields[:-2,0])
        LightGases[-1,0]=(self.__yields[-1,18]-self.__yields[-2,18])/(self.__yields[-1,0]-self.__yields[-2,0])
         #merge both in the Array
        self.__rates=np.concatenate((self.__rates,LightGases),axis=1)
        #Yields2Cols: updated for FG-DVC Version 8.2.2, one comment: 'Char_and_Ash' is actual the whole coal mass, but this part is named in the FG-DVC output file that way
        self.Yields2Cols={'Time':0,'Temp':1,'Total':2,'H2O':3,'Tar':4,'CO':5,'CO2':6,'CH4':7,'C2H4':8,'HCN':9,'NH3':10,'SO2':11,'COS':12,'CS2':13,'H2S':14,'Olefin':15,'Parafin':16,'H2':17,'Gas':18}
        self.Cols2Yields={0:'Time',1:'Temp',2:'Total',3:'H2O',4:'Tar',5:'CO',6:'CO2',7:'CH4',8:'C2H4',9:'HCN',10:'NH3',11:'SO2',12:'COS',13:'CS2',14:'H2S',15:'Olefin',16:'Parafin',17:'H2',18:'Gas'}
        #Filter for H2
        alpha=0.3  #weight of the neighbor values
        NumberOfFilterRuns=50
        for n in range(NumberOfFilterRuns): #smooths the yield curve of H2
            self.__yields[1:-1,17]=(alpha/2.)*(self.__yields[2:,17]+self.__yields[:-2,17])+(1-alpha)*self.__yields[1:-1,17]
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
        return 'FGDVC'
