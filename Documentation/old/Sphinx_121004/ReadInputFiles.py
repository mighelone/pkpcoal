import numpy as np
import scipy
from scipy import interpolate
######
class ReadFile(object):
    """general parent class for the reading objects CPDFile and FGDVCFile"""
    def __init__(self,InputFile):
        self.InputFile=InputFile
        self.Input=np.array(self.readLines())
        
    def readLines(self): #generates list from the input file
        """reads the input File line by line"""
        lines = [line.strip() for line in open(self.InputFile)]
        return lines

    def UsePyrolProgr(self,FileNote):
        """gets the information, whether Pyrolsis Program will be in use. Enter 'Yes' or 'True' for the case it should be used"""
        NameLine=np.where(self.Input[:]==FileNote)
        NameLine=NameLine[0]; NameLine=NameLine[0]
        ValueLine=NameLine+1       #the following line contains the value
        Value=self.Input[ValueLine]
        if Value=="'Yes'" or Value=="'yes'" or Value=="'true'" or Value=="'TRUE'" or Value=="'True'" or Value=='Yes' or Value=='yes' or Value=='true' or Value=='TRUE' or Value=='True':
            return True
        else:
            return False
        
    def Fitting(self,FileNote):
        """outputs the Fitting mode for Pyrolysis Program output (string: 'constantRate','Arrhenius','Kobayashi'). Possible input: 'constantRate', 'Arrhenius' or 'Kobayashi'"""
        NameLine=np.where(self.Input[:]==FileNote)
        NameLine=NameLine[0]; NameLine=NameLine[0]
        ValueLine=NameLine+1       #the following line contains the value
        Value=self.Input[ValueLine]
        if Value=='constantRate' or Value=="'constantRate'" or Value=='CONSTANTRATE' or Value=="'CONSTANTRATE'":
            return 'constantRate'
        elif Value=='Arrhenius'  or Value=="'Arrhenius'" or Value=='ARRHENIUS'  or Value=="'ARRHENIUS'":
            return 'Arrhenius'
        elif Value=='Kobayashi' or Value=="'Kobayashi'" or Value=='KOBAYASHI' or Value=="'KOBAYASHI'":
            return 'Kobayashi'
        else:
            return None
    
    def getValue(self,FileNote):
        """output the data of the line below the FileNote as a float"""
        NameLine=np.where(self.Input[:]==FileNote)
        NameLine=NameLine[0]; NameLine=NameLine[0]
        ValueLine=NameLine+1       #the following line contains the value
        Value=float(self.Input[ValueLine])
        return Value
        
    def getText(self,FileNote):
        """output the data of the line below the FileNote as a string"""
        NameLine=np.where(self.Input[:]==FileNote)
        NameLine=NameLine[0]; NameLine=NameLine[0]
        ValueLine=NameLine+1       #the following line contains the value
        Value=self.Input[ValueLine]
        return Value
                   
class WriteFGDVCCoalFile(object):
    """writes the file, which will be inputted into the FG-DVC coal generator"""
    def __init__(self,CoalGenFile):
        self.__FileName=CoalGenFile

    def setCoalComp(self,Carbon,Hydrogen,Oxygen,Nitrogen,Sulfur,SulfurPyrite):
        """Enter the coal composition with values in percent which have to sum up to 100""" 
        self.__C=Carbon
        self.__H=Hydrogen
        self.__O=Oxygen
        self.__N=Nitrogen
        self.__S=Sulfur
        self.__Sp=SulfurPyrite
        
    def write(self,CoalsDirectory,CoalResultFileName):
        """writes the FG-DVC coal generator input file"""
        self.__File = open(CoalsDirectory+self.__FileName, 'w')
        self.__File.write(str(self.__C)+' '+str(self.__H)+' '+str(self.__O)+' '+str(self.__N)+' '+str(self.__S)+' '+str(self.__Sp)+'\n')
        self.__File.write('0\n4\n'+CoalResultFileName+'_com.dat\n'+CoalResultFileName+'_kin.dat\n'+CoalResultFileName+'_pol.dat\n'+'5')
        self.__File.close()
        
class OperCondInput(ReadFile):
    """Reads the input file for the operating conditions and also writes the temperature-history file, required by FG-DVC."""
    def __init__(self,InputFile):
        self.InputFile=InputFile
        self.Input=np.array(self.readLines())
            
    def getTimePoints(self,FileNoteBegin,FileNoteEnd):
        """reads the time points in the shape 'time, temperature' for the lines between the line with the FileNoteBegin and the line with the FileNoteEnd"""
        BeginLine=np.where(self.Input[:]==FileNoteBegin)
        BeginLine=BeginLine[0]; BeginLine=BeginLine[0]
        ###
        EndLine=np.where(self.Input[:]==FileNoteEnd)
        EndLine=EndLine[0]; EndLine=EndLine[0]
        TimeTemp=np.genfromtxt(self.InputFile,skip_header=BeginLine+1,skip_footer=self.NumberOfLines()-EndLine, delimiter=',')
        return TimeTemp
        
    def writeFGDVCtTHist(self,tTPoints,dt,OutputFilePath):
        """Writes output file for FG-DVC containing in first column time in s, in the second tempearure in degree Celsius. FG-DVC will import this file. The time-temperature array has to be a numpy.array, dt a float, OutputFilePath a string."""
        tTFile=open(OutputFilePath,'w')
        #Transforms the T in degree Celsius, as required as input in FG-DVC t-T File (see manual V8.2.3 page 12)
        tTPoints[:,1]=tTPoints[:,1]-273.15
        THist=scipy.interpolate.interp1d(tTPoints[:,0],tTPoints[:,1],kind='linear')
        t=0.0
        while t<tTPoints[-1,0]:
            tTFile.write('%.5f  %.2f \n' % (t,THist(t)))
            t+=dt
        tTFile.close()

    def NumberOfLines(self):
        return int(len(self.Input[:]) )
