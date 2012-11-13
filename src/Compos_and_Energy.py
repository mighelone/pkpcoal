#PARAMETER:
#Molecular weights [Formeln und Tabellen, Peatec, 8.Auflage] (g/mol):
MH=2.01588*0.5
MC=12.011
MN=28.0134*0.5
MO=31.9988*0.5
MS=32.06
MRaw=100.
#Latent Heat of Water in J/kg :
rH2O=2263073
#Enthalpie of Formation in J/mol:
hfCO2=-3.935428e8
hfSO2=-2.968612e+08
hfH2O=-2.418428e8
hfCH4=-7.4829298e7
hfCO =-1.105397e8
hfO2 = -847.6404
hfN2 = 1429.881
hfH2 = 2448.595
hfChar = -101.268
######################################
import sys

class SpeciesBalance(object):
    """This is the parent Class for the CPD and FG-DVC specific Species Balances, containing general methods like the Dulong formular."""
                
    def SpeciesIndex(self,species):
        """Returns the column number of the input species."""
        return self.Yields2Cols[species]
        
    def Dulong(self):
        """Uses the Dulong formular to calculate the Higher heating value. The output is in J/kg."""
        #HHV = 32.79 MJ/kg fC + 150.4 (fH - fO/8) + 9.26 fS + 4.97 fO + 2.42 fN
        HHV = 32.79*self.UAC + 150.4*(self.UAH - self.UAO/8) + 9.26*0.0 + 4.97*self.UAO + 2.42*self.UAN
        return HHV*1.0e6      
    
    def correctUA(self):
        """Scale Ultimate Analysis to have sum=1"""
        #UAC(new)=UAC+(1-UAH-UAN-UAO-UAC)
        #self.UAC=1.-self.UAH-self.UAN-self.UAO
	sumUA = self.UAC+self.UAH+self.UAO+self.UAN+self.UAS
	self.UAC = self.UAC/sumUA
	self.UAH = self.UAH/sumUA
	self.UAO = self.UAO/sumUA
	self.UAN = self.UAN/sumUA
	self.UAS = self.UAS/sumUA
	sumUA = self.UAC+self.UAH+self.UAO+self.UAN+self.UAS

###############CPD specific class#############################

class CPD_SpeciesBalance(SpeciesBalance):
    """This class calculates the Species and the Energy balance for CPD. See the manual for the formulas and more details."""
    def __init__(self,CPD_ResultObject,UAC,UAH,UAN,UAO,UAS,PAVM,PAFC,PAmoist,PAash,HHV,MTar,RunNr):
        #imports the dictionaries for the species defined in 'CPD_Fit_one_run.py' -> CPD_Result()
        self.Yields2Cols=CPD_ResultObject.DictYields2Cols()
        self.Cols2Yields=CPD_ResultObject.DictCols2Yields()
        #imports final yield composition defined in CPD_Fit_one_run.py' -> CPD_Result()
        self.Yields=CPD_ResultObject.FinalYields()
        #0:'Time', 1:'Temp', 2:'Tar', 3:'Gas', 4:'Solid', 5:'Total', 6:'H2O', 7:'CO2', 8:'CH4', 9:'CO', 10:'Other'
        self.UAC=UAC/100.
        self.UAH=UAH/100.
        self.UAN=UAN/100.
        self.UAO=UAO/100.
        self.UAS=UAS/100.
        self.PAVM=PAVM/100.
        self.PAFC=PAFC/100.
	self.PAash=PAash/100.
	self.PAmoist=PAmoist/100.
        if HHV==0:
            HHV=self.Dulong()
        self.HHV=HHV
        self.MTar=MTar
        #corrects UA, if UA<1: Sulfur -> Carbon
        self.correctUA()
        self.CPDBalanceFile=open('CPD-BalanceResults'+str(RunNr)+'.txt','w')
        self.CPDBalanceFile.write('= Coal Properties =\n==Ultimate Analysis==\n')
        self.CPDBalanceFile.write('|C |'+str('%6.3f' %(self.UAC*100.))+'%|\n')
        self.CPDBalanceFile.write('|H |'+str('%6.3f' %(self.UAH*100.))+'%|\n')
        self.CPDBalanceFile.write('|N |'+str('%6.3f' %(self.UAN*100.))+'%|\n')
        self.CPDBalanceFile.write('|O |'+str('%6.3f' %(self.UAO*100.))+'%|\n')
        self.CPDBalanceFile.write('|S |'+str('%6.3f' %(self.UAS*100.))+'%|\n\n')
        self.CPDBalanceFile.write('== Proximate analysis ==\n')
        self.CPDBalanceFile.write('|   |  AR% | dry% | daf% |\n')
        self.CPDBalanceFile.write('|VM |'+str('%6.3f|' %PAVM)+str('%6.3f|' %(100*PAVM/(100.-PAmoist)))+str('%6.3f|' %(100*PAVM/(100.-PAmoist-PAash)))+'\n')
        self.CPDBalanceFile.write('|FC |'+str('%6.3f|' %PAFC)+str('%6.3f|' %(100*PAFC/(100.-PAmoist)))+str('%6.3f|' %(100*PAFC/(100.-PAmoist-PAash)))+'\n')
        self.CPDBalanceFile.write('|ash|'+str('%6.3f|' %PAash)+str('%6.3f|' %(100*PAash/(100.-PAmoist)))+str('%6.3f|' %0.0)+'\n')
        self.CPDBalanceFile.write('|H2O|'+str('%6.3f|' %PAmoist)+str('%6.3f|' %0.0)+str('%6.3f|' %0.0)+'\n\n')
        #if sum Yields != 1.0: scales every Yield up
        self.__correctYields()
        #print "Sum of Yields, input" ,( self.Yields[self.SpeciesIndex('Solid')]+self.Yields[self.SpeciesIndex('Tar')]+self.Yields[self.SpeciesIndex('CO')]+self.Yields[self.SpeciesIndex('CO2')]+self.Yields[self.SpeciesIndex('H2O')]+self.Yields[self.SpeciesIndex('CH4')] + self.Yields[self.SpeciesIndex('Other')])
        #print "Sum of Yields, input + Nitrogen" ,( self.Yields[self.SpeciesIndex('Solid')]+self.Yields[self.SpeciesIndex('Tar')]+self.Yields[self.SpeciesIndex('CO')]+self.Yields[self.SpeciesIndex('CO2')]+self.Yields[self.SpeciesIndex('H2O')]+self.Yields[self.SpeciesIndex('CH4')] + self.Yields[self.SpeciesIndex('Other')] +self.UAN)
        #print "UA ", (self.UAC+self.UAH+self.UAN+self.UAO)
        self.__Q_React() #calculate heat of reaction for raw coal
        self.__CheckOxygen()
        self.__CheckOthers()
        self.__TarComp()
        self.__hfRaw()
        self.__hfTar()
        self.__QPyro()
        self.__closeResultFile()


        
    def __correctYields(self):
        """Correct the amount of the yields 'Other'."""
        SumYieldsWitoutN=( self.Yields[self.SpeciesIndex('Solid')]+
            self.Yields[self.SpeciesIndex('Tar')]+
            self.Yields[self.SpeciesIndex('CO')]+
            self.Yields[self.SpeciesIndex('CO2')]+
            self.Yields[self.SpeciesIndex('H2O')]+
            self.Yields[self.SpeciesIndex('CH4')]+
            self.Yields[self.SpeciesIndex('Other')])
        self.Yields[self.SpeciesIndex('Solid')]=self.Yields[self.SpeciesIndex('Solid')]/SumYieldsWitoutN
        self.Yields[self.SpeciesIndex('Tar')]=self.Yields[self.SpeciesIndex('Tar')]/SumYieldsWitoutN
        self.Yields[self.SpeciesIndex('CO')]=self.Yields[self.SpeciesIndex('CO')]/SumYieldsWitoutN
        self.Yields[self.SpeciesIndex('CO2')]=self.Yields[self.SpeciesIndex('CO2')]/SumYieldsWitoutN
        self.Yields[self.SpeciesIndex('H2O')]=self.Yields[self.SpeciesIndex('H2O')]/SumYieldsWitoutN
        self.Yields[self.SpeciesIndex('CH4')]=self.Yields[self.SpeciesIndex('CH4')]/SumYieldsWitoutN
        self.Yields[self.SpeciesIndex('Other')]=self.Yields[self.SpeciesIndex('Other')]/SumYieldsWitoutN
        

    def __CheckOxygen(self):
        """Checks weather the amount of Oxygen in the light gases is lower than in the Ultimate Analysis. If not, the amount of these species decreases. If yes, the tar contains oxygen."""
        # eq. (1) in Michele's report
        Sum1=(1.*self.Yields[self.SpeciesIndex('H2O')])/(MO+2*MH) + (2.*self.Yields[self.SpeciesIndex('CO2')])/(MC+2*MO) + (1.*self.Yields[self.SpeciesIndex('CO')])/(MC+MO)
        Y0_cpd=MO*Sum1
        self.gamma=(self.UAO/Y0_cpd) #gamma in eq. (2) in Michele's report
        #
        if self.gamma<1.: #if gamma < 1 then is Y0_cpd > UAO
            # eq. (4) in Michele's report:
            Sum4=self.Yields[self.SpeciesIndex('H2O')]+self.Yields[self.SpeciesIndex('CO2')]+self.Yields[self.SpeciesIndex('CO')]
            self.Yields[self.SpeciesIndex('Other')]=self.Yields[self.SpeciesIndex('Other')]+(1-self.gamma)*Sum4
            # eq. (3) in Michele's report:
            self.Yields[self.SpeciesIndex('H2O')]=self.gamma*self.Yields[self.SpeciesIndex('H2O')]
            self.Yields[self.SpeciesIndex('CO2')]=self.gamma*self.Yields[self.SpeciesIndex('CO2')]
            self.Yields[self.SpeciesIndex('CO')] =self.gamma*self.Yields[self.SpeciesIndex('CO')]


    
    def __CheckOthers(self):
        """If the yield of nitrogen (is equal the UA of nitrogen) is lower than the species 'Other' the difference is set equal Methane."""
        if self.UAN < self.Yields[self.SpeciesIndex('Other')]:
            # eq. (6) in Michele's report:
            self.Yields[self.SpeciesIndex('CH4')]=self.Yields[self.SpeciesIndex('CH4')]+(self.Yields[self.SpeciesIndex('Other')]-self.UAN)
        #Result File:
        self.CPDBalanceFile.write('== Final yields ==\n')
        self.CPDBalanceFile.write('Char '+ str('%7.5f\n' %self.Yields[self.SpeciesIndex('Solid')]))
        self.CPDBalanceFile.write('Tar  '+ str('%7.5f\n' %self.Yields[self.SpeciesIndex('Tar')]))
        self.CPDBalanceFile.write('H2O  '+ str('%7.5f\n' %self.Yields[self.SpeciesIndex('H2O')]))
        self.CPDBalanceFile.write('CO2  '+ str('%7.5f\n' %self.Yields[self.SpeciesIndex('CO2')]))
        self.CPDBalanceFile.write('CH4  '+ str('%7.5f\n' %self.Yields[self.SpeciesIndex('CH4')]))
        self.CPDBalanceFile.write('CO   '+ str('%7.5f\n' %self.Yields[self.SpeciesIndex('CO')]))
        self.CPDBalanceFile.write('N2   '+ str('%7.5f\n' %self.UAN)+'\n')
        #self.CPDBalanceFile.write('Gas:  '+ str(self.Yields[self.SpeciesIndex('Gas')])+'\n\n')
        self.Yields[self.SpeciesIndex('Other')] = str(self.UAN)
        self.CPDBalanceFile.write('=== Fluent ===\n')
	char = self.Yields[self.SpeciesIndex('Solid')]*(1.-self.PAash)*100.
	vol = (1-self.Yields[self.SpeciesIndex('Solid')])*(1.-self.PAash)*100.
	ash = 100. - char - vol

	self.CPDBalanceFile.write('Volatile component fraction: '+ str('%7.3f\n' %char))
	self.CPDBalanceFile.write('Combustibile fraction:       '+ str('%7.3f\n' %vol))
	self.CPDBalanceFile.write('Ash fraction:                '+ str('%7.3f\n\n' %ash))

           
    def __TarComp(self): #saves [m,n,o] from C_m H_n O_o
        """Calculates and returns the tar composition."""
        # eq. (7) in Michele's report:
#        print "Sum of Yields" ,( self.Yields[self.SpeciesIndex('Solid')]+self.Yields[self.SpeciesIndex('Tar')]+self.Yields[self.SpeciesIndex('CO')]+self.Yields[self.SpeciesIndex('CO2')]+self.Yields[self.SpeciesIndex('H2O')]+self.Yields[self.SpeciesIndex('CH4')] +self.UAN )

#check Sum of UA
#        sumy = self.Yields[self.SpeciesIndex('Solid')]+self.Yields[self.SpeciesIndex('Tar')]+self.Yields[self.SpeciesIndex('CO')]+self.Yields[self.SpeciesIndex('CO2')]+self.Yields[self.SpeciesIndex('H2O')]+self.Yields[self.SpeciesIndex('CH4')] 
#        self.Yields[self.SpeciesIndex('Other')] = 1.0 - sumy

        #check Sum UA, sums up carbon part: Sulfur->Carbon, newer, higher amount of char (incresed by sulfur)
        #Carbon:
	sum = 0
        self.muCTar = ( self.UAC/MC - self.Yields[self.SpeciesIndex('Solid')]/MC - self.Yields[self.SpeciesIndex('CO2')]/(MC+2.*MO) - self.Yields[self.SpeciesIndex('CO')]/(MC+MO) - self.Yields[self.SpeciesIndex('CH4')]/(MC+4.*MH) )*(self.MTar/self.Yields[self.SpeciesIndex('Tar')])
        #Hydrogen:        
        self.muHTar = ( self.UAH/MH - 2.*self.Yields[self.SpeciesIndex('H2O')]/(2.*MH+MO) - 4.*self.Yields[self.SpeciesIndex('CH4')]/(MC+4.*MH) )*(self.MTar/self.Yields[self.SpeciesIndex('Tar')])
        #Oxygen:
        self.muOTar = ( self.UAO/MO - self.Yields[self.SpeciesIndex('H2O')]/(2.*MH+MO) - 2.*self.Yields[self.SpeciesIndex('CO2')]/(MC+2.*MO) - self.Yields[self.SpeciesIndex('CO')]/(MC+MO) )*(self.MTar/self.Yields[self.SpeciesIndex('Tar')])
        if abs(self.muOTar)<1.e-10:
            self.muOTar=0.0
	self.muSTar =  self.UAS/MS  *(self.MTar/self.Yields[self.SpeciesIndex('Tar')])
	self.muNTar =  (self.UAN/MN - 2.*self.Yields[self.SpeciesIndex('Other')]/(2.*MN)) *(self.MTar/self.Yields[self.SpeciesIndex('Tar')])
        if abs(self.muNTar)<1.e-10:
            self.muNTar=0.0
        #print 'TarComposition: ', 'C=',self.muCTar,',H=', self.muHTar,',O=',self.muOTar,',N=',self.muNTar,',S=',self.muSTar
        self.CPDBalanceFile.write('== Tar ==\n=== Atomic composition ===\n')
        self.CPDBalanceFile.write(str('C_%5.3f ' %self.muCTar)+str('H_%5.3f ' %self.muHTar)+str('O_%5.3f ' %self.muOTar)+str('N_%5.3f ' %self.muNTar)+str('S_%5.3f\n' %self.muSTar))
        self.CPDBalanceFile.write(str('Molecular weigth = %5.3f kg/kmol\n\n' %self.MTar))
        #self.CPDBalanceFile.write('Carbon:   '+ str(self.muCTar)+'\n')
        #self.CPDBalanceFile.write('Hydrogen: '+ str(self.muHTar)+'\n')
        #self.CPDBalanceFile.write('Oxygen:   '+ str(self.muOTar)+'\n')
        #self.CPDBalanceFile.write('Nytrogen:   '+ str(self.muNTar)+'\n')
        #self.CPDBalanceFile.write('Sulphur:   '+ str(self.muSTar)+'\n\n')
	self.CPDBalanceFile.write('=== TAR oxidation ===\n')
	nuO2 =0.5*(2.*self.muCTar+self.muHTar/2.+2.*self.muSTar-self.muOTar)
	self.CPDBalanceFile.write(str('TAR + %5.3f O2 --> ' %nuO2)+str('%5.3f CO2 + ' %self.muCTar)+str('%5.3f H2O + ' %(self.muHTar/2))+str('%5.3f N2 + ' %(self.muNTar/2))+str('%5.3f SO2\n\n' %self.muSTar))
	self.CPDBalanceFile.write('=== TAR partial oxidation ===\n')
	nuO2 =0.5*(self.muCTar+2.*self.muSTar-self.muOTar)
	self.CPDBalanceFile.write(str('TAR + %5.3f O2 -->' %nuO2)+str(' %5.3f CO + ' %self.muCTar)+str('%5.3f H2 + ' %(self.muHTar/2))+str('%5.3f N2 + ' %(self.muNTar/2))+str('%5.3f SO2\n\n' %self.muSTar))
        # check mass weight
	MTarCheck=self.muCTar*MC+self.muHTar*MH +self.muOTar*MO+self.muNTar*MN+self.muSTar*MS 
        return [self.muCTar,self.muHTar,self.muOTar,self.muSTar,self.muNTar]
        

    def __Q_React(self):
        """Calculates the Heat of Reaction of the coal cumbustion."""
        # eq. (8) in Michele's report
        HHVdaf=self.HHV/(self.PAFC+self.PAVM)
        # eq. (9) in Michele's report:
        LHVdaf = HHVdaf - ((MO+2.*MH)/(2.*MH))*self.UAH*rH2O
	LHVar = LHVdaf * (self.PAFC+self.PAVM) -self.PAmoist * rH2O
        # Q_react in J/mol; eq. (10) in Michele's report:
        self.Q_react=LHVdaf*MRaw
        self.CPDBalanceFile.write('== Heating value ==\n')
        self.CPDBalanceFile.write('|   | ar MJ/kg|daf MJ/kg|\n') 
        self.CPDBalanceFile.write('|HHV|'+str('%9.3f|' %(self.HHV/1e6))+str('%9.3f|\n' %(HHVdaf/1e6)))
        self.CPDBalanceFile.write('|LHV|'+str('%9.3f|' %(LHVar/1e6))+str('%9.3f|\n\n' %(HHVdaf/1e6)))
        
    def __hfRaw(self):
        """Calculates the heat of formation of the coal molecule and writes it into the output file."""
        # eq. (12) in Michele's report:
        muCRaw=self.UAC*(MRaw/MC)
        muHRaw=self.UAH*(MRaw/MH)
        muORaw=self.UAO*(MRaw/MO)
        muNRaw=self.UAN*(MRaw/MN)
	muSRaw=self.UAS*(MRaw/MS)
        # eq. (13) in Michele's report:
#        print 'Q_React =',self.Q_react/1e6,'MJ/kg'
#        print 'muCRaw =',muCRaw
#        print 'muHRaw =',muHRaw
#        print 'muORaw =',muORaw
#        print 'H2O =',muHRaw
#        print 'muORaw =',muORaw        
        muCO2 = muCRaw
        muH2O = 0.5*muHRaw
        muN2 = 0.5*muNRaw
	muSO2 = muSRaw
        muO2 = 0.5*(2*muCO2 + muH2O + 2*muSO2 - muORaw)
#        print 'muCO2=',muCO2
#        print 'muH2O=',muH2O
#        print 'muO2=',muO2 
#        print 'muN2=',muN2        
#        print 'mu*hfCO2=',muCO2*hfCO2
#        print 'muhfH2O=',muH2O*hfH2O
#        print 'muhfO2=',muO2*hfO2 
#        print 'muhfN2=',muN2*hfN2               
        self.hfraw = self.Q_react + muCO2*hfCO2 + muH2O*hfH2O + muN2*hfN2 +muSO2*hfSO2 - muO2*hfO2
        #self.CPDBalanceFile.write('Heat of Formation in J/mol: \n')
        #self.CPDBalanceFile.write('Raw Coal Molecule: '+ str(self.hfraw)+'\n')
        
    def __nysEq15(self):
        """Calculates the stoichiometric coefficients of the devolatilization reaction."""
        #self.hfRaw()
        # eq. (15) in Michele's report:
        self.nyTar=self.Yields[self.SpeciesIndex('Tar')]*MRaw/self.MTar
        self.nyChar=self.Yields[self.SpeciesIndex('Solid')]*MRaw/MC
        self.nyH2O=self.Yields[self.SpeciesIndex('H2O')]*MRaw/(2.*MH+MO)
        self.nyCO2=self.Yields[self.SpeciesIndex('CO2')]*MRaw/(MC+2.*MO)
        self.nyCH4=self.Yields[self.SpeciesIndex('CH4')]*MRaw/(MC+4.*MH)
        self.nyCO=self.Yields[self.SpeciesIndex('CO')]*MRaw/(MC+MO)   
        self.nyN2=self.Yields[self.SpeciesIndex('Other')]*MRaw/(2.*MN)           
#        print 'nyTar=',self.nyTar    
#        print 'nyChar=',self.nyChar    
#        print 'nyH2O=',self.nyH2O    
#        print 'nyCO2=',self.nyCO2    
#        print 'nyCH4=',self.nyCH4    
#        print 'nyCO=',self.nyCO    
#        print 'nyN2=',self.nyN2  
      #  
#        print 'Tar=',self.Yields[self.SpeciesIndex('Tar')]    
#        print 'Char=',self.Yields[self.SpeciesIndex('Solid')] 
#        print 'H2O=',self.Yields[self.SpeciesIndex('H2O') ]
#        print 'CO2=',self.Yields[self.SpeciesIndex('CO2')] 
#        print 'CO=',self.Yields[self.SpeciesIndex('CO')] 
#        print 'N2=',self.Yields[self.SpeciesIndex('Other')]         
        
    def __hfTar(self):
        """Calculates the heat of formation of tar and writes it into the CPD Composition file."""
        self.__nysEq15()
        # eq. (16) in Michele's report:
        Sum16 = self.nyH2O*hfH2O + self.nyCO2*hfCO2+self.nyCH4*hfCH4 + self.nyCO*hfCO + self.nyChar*hfChar+self.nyN2*hfN2
        hfTar = ( self.hfraw - Sum16)/self.nyTar
        self.CPDBalanceFile.write(str('=== Enthalpy of formation ===\nhf = %10.3f J/kmol\n' %hfTar))
    
    def __QPyro(self):
        """Calculates the heat of the pyrolysis process, assuming heat of tar formation is equal zero."""
        self.__nysEq15()
        # eq. (17) in Michele's report:
        Sum17 = self.nyH2O*hfH2O + self.nyCO2*hfCO2 + self.nyCH4*hfCH4 + self.nyCO*hfCO
        QPyro = -(self.hfraw - Sum17)/(MRaw)
        # eq. (18) in Michele's report:
        QPyroVM=QPyro/(1-self.Yields[self.SpeciesIndex('Solid')])
        #self.CPDBalanceFile.write('Pyrolysis Heat in J/g: \n')
        #self.CPDBalanceFile.write('Q= '+ str(QPyro)+'\n')
        #self.CPDBalanceFile.write('refered to VM= '+ str(QPyroVM)+'\n')
        
    def __closeResultFile(self):
        """Closes the CPD Composition file."""
        self.CPDBalanceFile.close()


###############FG-DVC specific class###########################

class FGDVC_SpeciesBalance(SpeciesBalance):
    """This class calculates the Species and the Energy balance for FG-DVC. See the manual for the formulars and more details."""
    def __init__(self,FGDVC_ResultObject,UAC,UAH,UAN,UAO,PAVM,PAFC,HHV,MTar,RunNr):
        #imports the dictionaries for the species defined in 'CPD_Fit_one_run.py' -> CPD_Result()
        self.Yields2Cols=FGDVC_ResultObject.DictYields2Cols()
        self.Cols2Yields=FGDVC_ResultObject.DictCols2Yields()
        #imports final yield composition defined in CPD_Fit_one_run.py' -> CPD_Result()
        self.Yields=FGDVC_ResultObject.FinalYields()
        #0:'Time', 1:'Temp', 2:'Tar', 3:'Gas', 4:'Solid', 5:'Total', 6:'H2O', 7:'CO2', 8:'CH4', 9:'CO', 10:'Other'
        self.UAC=UAC/100.
        self.UAH=UAH/100.
        self.UAN=UAN/100.
        self.UAO=UAO/100.
        self.PAVM=PAVM/100.
        self.PAFC=PAFC/100.
        if HHV==0:
            HHV=self.Dulong()
        self.HHV=HHV
        self.MTar=MTar
        #corrects UA, if UA<1: Sulfur -> Carbon
        self.correctUA()
        self.FGBalanceFile=open('FGDVC_'+str(RunNr)+'_BalanceResults.txt','w')
        self.FGBalanceFile.write('Coal Properties input:\nUltimate Analysis:\n')
        self.FGBalanceFile.write('UAC= '+str(self.UAC*100.)+'\n')
        self.FGBalanceFile.write('UAH= '+str(UAH)+'\n')
        self.FGBalanceFile.write('UAN= '+str(UAN)+'\n')
        self.FGBalanceFile.write('UAO= '+str(UAO)+'\n')
        self.FGBalanceFile.write('Proximate analysis, as recieved:\n')
        self.FGBalanceFile.write('PAVM='+str(PAVM)+'\n')
        self.FGBalanceFile.write('PAFC='+str(PAFC)+'\n')
        self.FGBalanceFile.write('Higher Heating Value, as recieved:\n')
        self.FGBalanceFile.write('HHV= '+str(HHV)+' MJ/kg\n\n\n')
        #considered yields: char, tar, CO, CO2, H2O, CH4, N2, H2, O2
        self.__correctYields()
        #the missing of the UA is included into carbon
        self.correctUA()
        #calculates the Tar composition:
        self.__TarComp()
        #writes Composition results into file:
        self.__writeSpeciesResults()
        #Energy Balance:
        self.__EnergyBalance()
        self.__writeEnergyResults()
        self.__closeFile()
        
    def __correctYields(self):
        """Further, only the yields of char, tar, CO, CO2, H2O, CH4, N2, H2, O2 are considered. Modifies the yields, merge species like Olefins parafins, HCN to tar."""
        self.Char=1.0-self.Yields[self.SpeciesIndex('Total')]
        self.CO=self.Yields[self.SpeciesIndex('CO')]
        self.CO2=self.Yields[self.SpeciesIndex('CO2')]
        self.H2O=self.Yields[self.SpeciesIndex('H2O')]
        self.CH4=self.Yields[self.SpeciesIndex('CH4')]
        #self.N2=self.UAN  ?? elemental Nitrogen or all in Tar
        self.H2=self.Yields[self.SpeciesIndex('H2')]
        self.Tar=self.Yields[self.SpeciesIndex('Tar')]
        Sum=self.Char+self.CO+self.CO2+self.H2O+self.CH4+self.H2+self.Tar#+self.N2 # ??
        self.Tar += 1.-Sum
        print 'Sum of modified Species: ', self.Char+self.CO+self.CO2+self.H2O+self.CH4+self.H2+self.Tar#+self.N2 #??
        
    def __TarComp(self):
        """Calculates the Tar composition using analyisis of the Ultimate Analysis."""
        #The species in tar (kg) per kg of coal
        #Carbon:
        SumC=self.Char+self.CO*MC/(MC+MO)+self.CO2*MC/(MC+2.*MO)+self.CH4*MC/(MC+4.*MH)
        TarC=self.UAC-SumC
        #Hydrogen:
        SumH=2.*self.H2O*MH/(2*MH+MO)+4.*self.CH4*MH/(4*MH+MC)+self.H2
        TarH=self.UAH-SumH
        #Oxygen:
        SumO=self.CO*MO/(MC+MO)+2.*self.CO2*MO/(MC+2.*MO)+self.H2O*MO/(2.*MH+MO)
        TarO=self.UAO-SumO
        #Nitrogen:
        TarN=self.UAN
        #
        #now, tar coefficients (molar, not mass)
        self.TarC=(TarC*self.MTar)/(self.Tar*MC)
        self.TarH=(TarH*self.MTar)/(self.Tar*MH)
        self.TarO=(TarO*self.MTar)/(self.Tar*MO)
        self.TarN=(TarN*self.MTar)/(self.Tar*MN)
        print "Mass Weight of Tar " ,  self.TarC*MC+self.TarH*MH +self.TarO*MO+self.TarN
    
    def __writeSpeciesResults(self):
        """Writes the Species Balance results into the result file."""
        self.FGBalanceFile.write('modified yields for further species and energy calculations:\n')
        self.FGBalanceFile.write('Char:  '+str(self.Char)+'\n')
        self.FGBalanceFile.write('CO:    '+str(self.CO)+'\n')
        self.FGBalanceFile.write('CO2:   '+str(self.CO2)+'\n')
        self.FGBalanceFile.write('H2O:   '+str(self.H2O)+'\n')
        self.FGBalanceFile.write('CH4:   '+str(self.CH4)+'\n')
        self.FGBalanceFile.write('H2:    '+str(self.H2)+'\n')
        self.FGBalanceFile.write('Tar:   '+str(self.Tar)+'\n\n')
        self.FGBalanceFile.write('The tar molecule composition:\n')
        self.FGBalanceFile.write('Carbon:   '+str(self.TarC)+'\n')
        self.FGBalanceFile.write('Hydrogen: '+str(self.TarH)+'\n')
        self.FGBalanceFile.write('Oxygen:   '+str(self.TarO)+'\n')
        self.FGBalanceFile.write('Nitrogen: '+str(self.TarN)+'\n\n\n')
        
    def __EnergyBalance(self):
        """Calculates the heat of formation of tar."""
        #Reaction enthalpie (J/mol) for: (E..Energy, M...molar)
        CharEM = hfChar+hfO2-hfCO2
        H2EM   = hfH2+0.5*hfO2-hfH2O
        CH4EM  = hfCH4+2.*hfO2-hfCO2-2.*hfH2O
        COEM   = hfCO+0.5*hfO2-hfCO2
#        print 'CharEM ', CharEM
#        print 'H2EM', H2EM
#        print 'CH4EM', CH4EM
#        print 'COEM', COEM, '\n'
        #Reaction enthalpie (J/kg), abs. including yields:
#        print 'CharEM/MC',CharEM/MC
#        print 'H2EM/(2.*MH)',H2EM/(2.*MH)
#        print '(CH4EM)/(4.*MH+MC)',(CH4EM)/(4.*MH+MC)
#        print '(COEM)/(MC+MO)',(COEM)/(MC+MO), '\n'
        CharE=(CharEM*self.Char)/MC
        H2E=(H2EM*self.H2)/(2.*MH)
        CH4E=(CH4EM*self.CH4)/(4.*MH+MC)
        COE=(COEM*self.CO)/(MC+MO)
#        print 'CharE', CharE
#        print 'H2E', H2E
#        print 'CH4E', CH4E
#        print 'COE', COE, '\n'
        #Reaction enthalpie Tar (J/kg):
        print 'HHV as rec =',self.HHV/1e6,'MJ/kg'
        HHVdaf=self.HHV/(self.PAFC+self.PAVM)
        print 'HHV daf =',HHVdaf/1e6,'MJ/kg'
        LHVdaf = HHVdaf - ((MO+2.*MH)/(2.*MH))*self.UAH*rH2O
        self.LHVdaf=LHVdaf
        print 'LHV-daf =',LHVdaf/1e6,'MJ/kmol'
        TarE=LHVdaf-CharE-H2E-CH4E-COE
#        print 'TarE', TarE
        #Reaction enthalpie Tar in J/mol
        TarEM=(TarE*self.MTar)/self.Tar
#        print 'TarEM',TarEM
        #coeffs
        nuN2=self.TarN*0.5
        nuH2O=self.TarH*0.5
        nuCO2=self.TarC
        nuO2=0.5*(2*self.TarC+0.5*self.TarH-self.TarO)
        #heat of formation for tar:
#        print 'O2', nuO2*hfO2
#        print 'CO2', nuCO2*hfCO2
#        print 'H2O', nuH2O*hfH2O
#        print 'N2', nuN2*hfN2
        self.hfTar=TarEM-nuO2*hfO2+nuCO2*hfCO2+nuH2O*hfH2O+nuN2*hfN2
    
    def __writeEnergyResults(self):
        """Writes the Energy results into the result file."""
        self.FGBalanceFile.write('The lower heating value of daf coal:\n')
        self.FGBalanceFile.write('LHV:   %.4f MJ/mol\n' % (self.LHVdaf*1.e-6))
        self.FGBalanceFile.write('The enthalpie of foramtion for the tar:\n')
        self.FGBalanceFile.write('hfTar: %.4f MJ/mol\n' % (self.hfTar*1.e-6))
        
    def __closeFile(self):
        """Closes the FG-DVC Composition file."""
        self.FGBalanceFile.close()
        
