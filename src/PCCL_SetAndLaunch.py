import sys
import numpy as np
import os
import shutil

class SetterAndLauncher(object):
    """This class is able to write the PC Coal Lab input file and launch CPD program. Before writing the CPD input file (method 'writeInstructFile') set all parameter using the corresponding methods. After writing the instruct file, the .Run method can be used."""
    def __init__(self):
        # in the following, rarely changed parameter are defined
        ### THE OPERATING CONDITIONS ###
        self.TestReactor  = 'WG ' # wireGrid
        self.TestSeq = 'FFTG'#'FFQQ' # freeform, mainly heating rate is varied (but FF tells that not all other are kept constant. see man 4.7, 4.13)
#        self.OutputRes = 'FD' # Fully dynamics
#        self.OutputRes = 'UY' # ultimate yields only
        self.OutputRes = 'UD' # full dynmaics AND ultimate yields
        self.Gas = 'N2' # must be nitrogen
        self.Tinit = 25. #in degree Clesius, F6.0 char
        self.Tultim = 1500. #in degree Clesius, F6.0 char
        self.Hr = 10000 # heating rate in K/s , F7.0 char
        self.O2 = 0.  # O2 conc must be zero for wire grid, F6.1
        self.pressure = 0.1 # pressure in MPa
        self.tHold = 0.011 # hold time, t in s, F8.3
        self.particleDiam = 100.0 # particle diameter in micrometer, F6.1
        # see page 4.8
        self.RS = ['Y','Y','Y','Y','Y','Y','Y','Y','Y','Y','N']  # all species reported, Burnout report not possible for wire grid (last entry)
        self.RM = ['N','N','N']                                  # all three [SFOR, DAEM, C2SM] are fitted
        self.SA = ['N','N','N','N','N','N','N','N','N','N']              # only the overall yields are fitted
        ### THE COAL PROPERTIES ###
        self.Fuel = 'CO '   #CO for Coal, PC for petroleeum coke and BM for BioMass
        self.CoalLable = 'TESTCOAL                                 '  #40 chars
        # P A
        self.VM = 45.0  # must be printed in format F5.1 e.g.' 41.4'
        self.FC = 45.0
        self.Moist = 5.0
        self.Ash = 5.
        # U A
        self.UAC = 75.0
        self.UAH = 5.0
        self.UAO = 19.0
        self.UAN = 1.0
        self.UAS = 0.0
        self.CalibrFactor = False  #set False or equal a value
        #
        #Counter for the first column in the testplan file
        self.CoalCounter = 1 #for Coal File
        self.OCCounter = 1 #for operating Condition / Testplan


    def SetUACoalParameter(self,fcar,fhyd,fnit,foxy,fsul=0.0):
        """Set the mass fraction of carbon, hydrogen, nitrogen, oxygen and sulfur (there is standard is zero)."""
        sumUA = fcar+fhyd+fnit+foxy+fsul
        if sumUA != 100.:
            sumUA = sumUA/100.
            self.UAC = fcar/sumUA
            self.UAH = fhyd/sumUA
            self.UAO = foxy/sumUA
            self.UAN = fnit/sumUA
            self.UAS = fsul/sumUA
        else:
            self.UAC = fcar
            self.UAH = fhyd
            self.UAO = foxy
            self.UAN = fnit
            self.UAS = fsul
    
    def SetPACoalParameter(self,volatileMatter,fixedCarbon,Moisture,Ash):
        """Sets the Proximate Analysis of the Coal using the as recieved state."""
        sumVM = volatileMatter+fixedCarbon+Moisture+Ash
        if sumVM != 100.:
            sumVM = sumVM/100.
            self.VM = volatileMatter/sumVM
            self.FC = fixedCarbon/sumVM
            self.Moist = Moisture/sumVM
            self.Ash = Ash/sumVM
        else:
            self.VM = volatileMatter
            self.FC = fixedCarbon
            self.Moist = Moisture
            self.Ash = Ash
        
    def SetCoalCalibrationFactor(self,CalibrationFactor=False):
        """Here the PC Coal Lab Coal Calibration factor can be defined (float) if available."""
        self.CalibrFactor = CalibrationFactor

    def THist(self,Tstart,t1,Tfinal,tfinal):
        """Define the heating rate ramp. Therefore enter the Start Temperature, the Temperature and time when the linear heating ends and the final Temperature and time. All are in Kelvin and are converted into Degree Celsius."""
        self.Tinit = Tstart-273.
        self.Hr = (Tfinal-Tstart)/t1
        self.Tultim = Tfinal-273.
        self.tHold = tfinal-t1
    
    def SetPressure(self,pressure):
        """Defines the pressure in atm"""
        atm2MPa = 1.01325e-1
        self.pressure = pressure*atm2MPa
        if self.pressure < 0.01:
            self.pressure = 0.01 #because otherwise with respect to the input format a pressure of 0.0 would be printed
            
    def SetParticleSize(self,DiameterinMicrons):
        """Defines the particle diameter. Input it in micrometer."""
        self.particleDiam = DiameterinMicrons
    
    def writeCoalFiles(self,Dirpath,mkNewFile=False):
        """Writes the Coal File 'Coalpc.dat' into the given path."""
        self.DirPath = Dirpath
        if mkNewFile == True:
            self.CoalCounter = 1
        if self.CoalCounter == 1:
            ini=open(Dirpath+'Coalpc.dat','w') #mk new file
        elif self.CoalCounter < 6:
            ini=open(Dirpath+'Coalpc.dat','a') #append to old file
        else:
            print 'Do not make more than five PC Coal Lab runs in one Coal File.'
        ini.write( self.Fuel )
        ini.write( self.CoalLable )
        ini.write( '%4.1f ' % self.VM )
        ini.write( '%4.1f ' % self.FC )
        ini.write( '%4.1f ' % self.Moist )
        ini.write( '%4.1f ' % self.Ash )
        ini.write( '%4.1f ' % self.UAC )
        ini.write( '%4.1f ' % self.UAH )
        ini.write( '%4.1f ' % self.UAO )
        ini.write( '%4.1f ' % self.UAN )
        ini.write( '%4.1f' % self.UAS )
        if self.CalibrFactor == False:
            ini.write( 'N\n' )
        elif type(self.CalibrFactor) == float:
            ini.write( 'Y\n'+('%6.4f' % self.CalibrFactor) + '\n' )
        ini.close()   
        self.CoalCounter += 1

    def writeInstructFiles(self,Dirpath,mkNewFile=False):
        """Writes the Testplan File 'Testplan.dat' into the given path. If New File is set to True, a new one is generated. Otherwise the next lines are added to the existing Testplan.dat"""
        self.DirPath = Dirpath
        if mkNewFile == True:
            self.OCCounter = 1
        if self.OCCounter == 1:
            ini=open(Dirpath+'Testplan.dat','w') #mk new file
        elif self.OCCounter < 6:
            ini=open(Dirpath+'Testplan.dat','a') #append to old file
        else:
            print 'Do not make more than five PC Coal Lab runs in one Testplan.'
            sys.exit()
        ini.write( '0'+('%1i' % self.OCCounter)+' '  )
        ini.write( self.TestReactor )
        ini.write( self.TestSeq + ' ' )
        ini.write( self.OutputRes + ' ' )
        ini.write( self.Gas + ' ' )
        ini.write( '%4.0f.' % self.Tinit )
        ini.write( '%5.0f.' % self.Tultim )
        ini.write( '%6.0f.' % self.Hr )
        ini.write( '%6.1f' % self.O2 )
        ini.write( '%6.2f' % self.pressure )
        ini.write( '%8.3f' % self.tHold )
        ini.write( ('%6.1f' % self.particleDiam )+'\n' )
        ini.close()
        self.OCCounter += 1
        
    def writeInstructFilesFinish(self):
        """Writes the last statements into the Instruct File to finish it."""
        ini=open(self.DirPath+'Testplan.dat','a') #append to old file
        for i in self.RS:
            ini.write( i )
        ini.write( '\n' )
        for i in self.RM:
            ini.write( i )
        ini.write( '\n' )
        for i in self.SA:
            ini.write( i )
        ini.write( '\n' )

    def Run(self,PathToExe,NameOfExe):
        """Launches PCCoV41M1to7Par.exe."""
        OScd=PathToExe
        OSExe=NameOfExe
        os.system('cd '+OScd+' & '+OSExe)

if __name__ == "__main__":
    sl = SetterAndLauncher()
    sl.SetPACoalParameter(45.,42,10.,3.)
    sl.SetUACoalParameter(70,7,5,12)
#    sl.SetCoalCalibrationFactor(0.5227)
    sl.THist(300.,0.03,1500.,0.1)
    sl.SetPressure(1.3)
    #
    sl.writeCoalFiles('C:\\Users\\MaP\\PCCL\\')
    sl.THist(400.,0.1,1200.,0.11)
    sl.writeInstructFiles('C:\\Users\\MaP\PCCL\\')
    sl.THist(400.,0.09,1200.,0.1)
    sl.writeInstructFiles('C:\\Users\\MaP\PCCL\\')
    sl.THist(400.,0.08,1200.,0.1)
    sl.writeInstructFiles('C:\\Users\\MaP\PCCL\\')
    sl.THist(400.,0.07,1200.,0.1)
    sl.writeInstructFiles('C:\\Users\\MaP\PCCL\\')
    sl.THist(400.,0.06,1200.,0.1)
    sl.writeInstructFiles('C:\\Users\\MaP\PCCL\\')
    sl.writeInstructFilesFinish()
    sl.Run('C:\\Users\\MaP\\PCCL','PCCoV41M1to7Par.exe')
