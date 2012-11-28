import InformationFiles

class WriteCoalFile(object):
    """Writes the Coal.inp file using the output of the GUI."""
    def __init__(self,InfosFromGUIObject):
        self.Info=InfosFromGUIObject
        CoalFile=open('Coal.inp','w')
        CoalFile.write(self.__mkstr())
        CoalFile.close()
    
    def __mkstr(self):
        """Writes the content of the coal file"""
        FileStr='#Proximate Analysis (in percent, as recieved):\n'
        PA = InformationFiles.M_PA
        UA = InformationFiles.M_UA
        PAFC,PAVM,PAMoi,PAAsh=self.Info.PA()
        UAC,UAH,UAN,UAO,UAS=self.Info.UA()
        FileStr+=PA[0]+'\n'
        FileStr+=PAFC+'\n'
        FileStr+=PA[1]+'\n'
        FileStr+=PAVM+'\n'
        FileStr+=PA[2]+'\n'
        FileStr+=PAMoi+'\n'
        FileStr+=PA[3]+'\n'
        FileStr+=PAAsh+'\n'
        FileStr+='\n'
        FileStr+= '#Ultimate Analysis (in percent):\n'
        FileStr+=UA[0]+'\n'
        FileStr+=UAC+'\n'
        FileStr+=UA[1]+'\n'
        FileStr+=UAH+'\n'
        FileStr+=UA[2]+'\n'
        FileStr+=UAN+'\n'
        FileStr+=UA[3]+'\n'
        FileStr+=UAO+'\n'
        FileStr+=UA[4]+'\n'
        FileStr+=UAS+'\n'
        FileStr+='\n'
        #
        MRaw,MTar,HHV=self.Info.MwsHHV()
        FileStr+=InformationFiles.M_HHV+'\n'
        FileStr+=HHV+'\n'
        FileStr+=InformationFiles.M_MRaw+'\n'
        FileStr+=MRaw+'\n'
        FileStr+=InformationFiles.M_MTar+'\n'
        FileStr+=MTar+'\n'
        FileStr+='\n'
        #
        WY,WR = self.Info.WeightYR()
        FileStr+=InformationFiles.M_Weight[0]+'\n'
        FileStr+=WY+'\n'
        FileStr+=InformationFiles.M_Weight[1]+'\n'
        FileStr+=WR+'\n'
        return FileStr

        