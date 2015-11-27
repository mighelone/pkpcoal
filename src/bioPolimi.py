import numpy as np
from coal import coal
from coalPolimi import triangle, compositionError, coalPolimi
import sys
#sys.path.append('/usr/local/lib/python2.7/site-packages/')
import cantera
from scipy.integrate import odeint,ode

import warnings


class bioPolimi(coalPolimi):
    '''
    coal class customized for polimi model
    inherit methods from general coal class

    No N and S are assigned for this coal
    '''
    #def __init__(self,name='',c=0.8,h=0.05,o=0.15,n=0,s=0,file='Biomass.xml'):
    #    '''
    #    init class
    #    '''
    #    #super(bioPolimi, self).__init__(name=name, c=c, h=h, o=o, n=0, s=0, file=file)

    def _referenceBiomass(self):
        '''
        define the composition of the reference components and biomasses
        '''

        self._bio1 = coal(name='CELL',c=6*self._Mc,h=10*self._Mh,o=5*self._Mo,n=0,s=0)
        self._bio2 = coal(name='HCE',c=5*self._Mc,h=8*self._Mh,o=4*self._Mo,n=0,s=0)
        self._bio3 = coal(name='LIGC',c=15*self._Mc,h=14*self._Mh,o=4*self._Mo,n=0,s=0)
        self._bio4 = coal(name='LIGO',c=20*self._Mc,h=22*self._Mh,o=10*self._Mo,n=0,s=0)
        self._bio5 = coal(name='LIGH',c=22*self._Mc,h=28*self._Mh,o=9*self._Mo,n=0,s=0)
        self._bioS1 = coal(name='S1',c=0.448,h=0.0616,o=0.4906,n=0,s=0)
        self._bioS2 = coal(name='S2',c=0.59,h=0.0544,o=0.3556,n=0,s=0)
        self._bioS3 = coal(name='S3',c=0.619,h=0.064,o=0.3172,n=0,s=0)

    def _calculateComposition(self):
        '''
        calculate the composition of the actual biomass according to the reference biomasses
        '''
        self._referenceBiomass()
        # determine in which triangle the coal lies
        # out means outside
        self._inside = 'out'
        # triangle 012
        t123=triangle(x0=self._bioS1.getVanKravelen(),x1=self._bioS2.getVanKravelen(),x2=self._bioS3.getVanKravelen())
        if t123.isInside(self.getVanKravelen()):
            self._inside = '123'
            self.triangle = t123
            self._interpolateBiomass()
            return
        raise compositionError('Composition outside of triangle!')

    def _interpolateBiomass(self):
        '''
        interpolate coal using reference
        '''
        c0 = self._bioS1 # Biomass 1
        c1 = self._bioS2 # Biomass 2
        c2 = self._bioS3 # Biomass 2

        matrix = np.array(
            [[c0._c,c1._c,c2._c],
             [c0._h,c1._h,c2._h],
             [c0._o,c1._o,c2._o]])
        b = np.array([self._c,self._h,self._o])
        composition = np.linalg.solve(matrix,b)
        components = np.zeros(6)
        # Calculate components based on interpolated biomass S1 = 0.6 CELL + 0.4 HCE, S2 = 0.8 LIGO + 0.2 LIGC, S3 = 0.8 LIGH + 0.2 LIGC
        components[0] = composition[0]*0.6                      # CELL
        components[1] = composition[0]*0.4                      # HCE
        components[2] = (composition[1]*0.2 + composition[2]*0.2) # LIGC
        components[3] = composition[1]*0.8                      # LIGO
        components[4] = composition[2]*0.8                      # LIGH
        components[5] = 0.0 #Moisture
        # define string
        self._bioComposition = np.append(0.0,components)
        self._compositionString = 'CELL:'+str(components[0])+',HCE:'+str(components[1])+',LIGC:'+str(components[2])+',LIGO:'+str(components[3])+',LIGH:'+str(components[4])+',ACQUA:'+str(components[5])
        print self._compositionString

    def __repr__(self):
        out = coal.__repr__(self)
        out += '\nCELL:'+str(self._bioComposition[0])
        out += '\nHCE:'+ str(self._bioComposition[1])
        out += '\nLIGC:'+str(self._bioComposition[2])+'\n'
        out += '\nLIGO:'+str(self._bioComposition[3])+'\n'
        out += '\nLIGH:'+str(self._bioComposition[4])+'\n'
        return out

    def getBioComposition(self):
        return self._bioComposition

    def getRawBiomass(self):
        return self._y[:, self._solid.speciesIndex('CELL')] + \
               self._y[:, self._solid.speciesIndex('HCE')] + \
               self._y[:, self._solid.speciesIndex('LIGC')] + \
               self._y[:, self._solid.speciesIndex('LIGO')] + \
               self._y[:, self._solid.speciesIndex('LIGH')]

    def getCharCoal(self):
        char_species = ['Char',
                        'CELL',
                        'CELLA',
                        'HCE',
                        'HCE1',
                        'HCE2',
                        'ACQUA',
                        'LIGC',
                        'LIGO',
                        'LIGH',
                        'LIGOH',
                        'LIGCC',
                        'LIG',
                        'GCO2',
                        'GCO',
                        'GCOH2',
                        'GH2',
                        'GCH4',
                        'GCH3OH',
                        'GC2H4']
        return self.get_sumspecies(char_species)

    #def getCH4(self):
    #    return self._y[:,self._solid.speciesIndex('CH4')]

    #def getH2(self):
    #    return self._y[:,self._solid.speciesIndex('H2')]

    #def getH2O(self):
    #    return self._y[:,self._solid.speciesIndex('H2O')]


    #def getCO(self):
    #    return self._y[:,self._solid.speciesIndex('CO')]
    #def getCO2(self):
    #    return self._y[:,self._solid.speciesIndex('CO2')]
    #def getCH2(self):
    #    return self._y[:,self._solid.speciesIndex('CH2O')]

    #def getTAR(self):
    #    return self._y[:,self._solid.speciesIndex('XYLOSE')]

    #def getLightGases(self):
    #    return self.getCO() + \
    #           self.getCO2() + \
    #           self.getH2O() + \
    #           self.getH2() + \
    #           self.getCH4() + \
    #           self.getCH2() + \
    #           self._y[:,self._solid.speciesIndex('CH3OH')]

    #def getMetaplast(self):
    #    metaplast = self._y[:,self._solid.speciesIndex('GCO2')] +\
    #           self._y[:,self._solid.speciesIndex('GCO')] +\
    #           self._y[:,self._solid.speciesIndex('GCOH2')]
    #    return metaplast


    def getVolatile(self):
        """
        :rtype : object
        """
        volatile_species = ['HAA',
                            'HMFU',
                            'LVG',
                            'XYLOSE',
                            'Glyoxal',
                            'Phenol',
                            'pCoumaryl',
                            'C11H12O4',
                            'C3H6O2',
                            'C3H4O2',
                            'ALD3',
                            'MECHO',
                            'C2H6OH',
                            'C2H4',
                            'CH3OH',
                            'CH4',
                            'CO2',
                            'CO',
                            'H2O',
                            'H2',
                            'GCO2',
                            'GCO',
                            'GCOH2',
                            'GH2',
                            'EtOH',
                            'HCOOH']
        return self.get_sumspecies(volatile_species)


    def getTime(self):
        return self.time

    def getTemperature(self):
        return self._getInterpTemperature(self.time)


    def _calculateYields(self):
        '''
        return yields as daf fraction
        using format required by PKP
        '''
        self._yields=np.zeros((int(len(self.time)), 4))  #shapes new Matrix containing all necessary information;
        self._yields[:, 0]=self.time
        self._yields[:, 1]=self.getTemperature() # temp
        self._yields[:, 2]=self.getVolatile()#self.getTAR() # total volatile yield
        self._yields[:, 3]=self.getCharCoal()#1.-self._yields[:,2]
        #self._yields[:,3]=self.getLightGases()
        #self._yields[:,5]=self.getVolatile()
        #self._yields[:,4]=1.-self._yields[:,5]
        #self._yields[:,6]=self.getH2O()
        #self._yields[:,7]=self.getCO2()
        #self._yields[:,8]=self.getCH4()
        #self._yields[:,9]=self.getCO()
        #self._yields[:,10]=self.getH2()
        #self._yields[:,11]=self._y[:,self._solid.speciesIndex('CH2')]
        #self._yields[:,12]=self._y[:,self._solid.speciesIndex('CH3O')]
        #self._yields[:,13]=self._y[:,self._solid.speciesIndex('BTX2')]
        self.Yields2Cols={'Time':0,'Temp':1,'Total':2,'Solid':3}#,'HCE':4,'HCE1':5,'HCE2':6,'LIGC':7,'LIGH':8,'LIGO':9,'LIG':10,'LIGCC':11,'LIGOH':12,'Char':13,'HAA':14,'HMFU':15,'LVG':16,'XYLOSE':17,'Glyoxal':18,'Phenol':19,'pCoumaryl':20,'C11H12O4':21,'C3H6O2':22,'C3H4O2':23,'C3H6O':24,'CH3CHO':25,'C2H6OH':26,'C2H4':27,'CH3OH':28,'CH2O':29,'CH4':30,'CO2':31,'CO':32,'H2O':33,'H2':34,'GCO2':35,'GCO':36,'GCOH2':36,'GH2':37,'EtOH':38}
        self.Cols2Yields={0:'Time',1:'Temp',2:'Total',3:'Solid'}#2:'Tar',3:'Gas',4:'Solid',5:'Total',
                          #6:'H2O',7:'CO2',8:'CH4',9:'CO',10:'H2',
                          #11:'CH2',12:'CH3O',13:'BTX2'}

    def Yields_all(self):
        """Returns the whole result matrix of the yields."""
        self._calculateYields()
        return self._yields

    def Rates_all(self):
        """Returns the whole result matrix of the rates"""
        self.__rates=np.zeros(( int(len(self.time)),4) )  #shapes new Matrix containing all necessary information; 14!!!
        return self.__rates

    def DictYields2Cols(self):
        """Returns the whole Dictionary Yield names to Columns of the matrix"""
        return self.Yields2Cols

    def DictCols2Yields(self):
        """Returns the whole Dictionary Columns of the matrix to Yield names"""
        return self.Cols2Yields

    def FinalYields(self):
        """Returns the last line of the Array, containing the yields at the time=time_End"""
        return self._yields[-1, :]

    def Name(self):
        """returns 'CPD' as the name of the Program"""
        return 'bioPolimi'




