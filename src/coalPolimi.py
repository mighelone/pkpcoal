__author__ = 'vascella'


import numpy as np
from coal import coal
import cantera
from scipy.integrate import odeint,ode

import warnings




class triangle(object):
    '''
    calculate properties of coal
    using triangolation
    '''
    def __init__(self,x0=np.array([0,0]),x1=np.array([1,0]),x2=np.array([0,1])):
        self.x0 = x0
        self.x1 = x1
        self.x2 = x2

    def _coeff(self,x=np.array([0,0])):
        '''
        calculate coefficient of linear combination of x
        x-x0 = a*(x1-x0)+b*(x2-x0)
        '''
        v1 = self.x1-self.x0
        v2 = self.x2-self.x0
        v = x-self.x0
        matr=np.transpose(np.array([v1,v2]))
        return np.linalg.solve(matr,v)

    def isInside(self,x=np.array([0,0])):
        '''
        verify is point x is inside the triangle
        '''
        coeff = self._coeff(x)
        if (coeff[0]>=0 and coeff[1] >=0 and sum(coeff)<=1):
            return True
        else:
            return False

class compositionError(Exception):
    pass

class coalPolimi(coal):
    '''
    coal class customized for polimi model
    inherit methods from general coal class

    No N and S are assigned for this coal
    '''
    def __init__(self,name='',c=0.8,h=0.05,o=0.15,n=0,s=0,file='COAL.xml'):
        '''
        init class
        '''
        coal.__init__(self,name=name,c=c,h=h,o=o,n=0.,s=0.)
        self._calculateCoalComposition()
        self._setCanteraObject(file=file)

        # define default parameter for pyrolysis
        self.setHeatingRate()
        self.setTimeStep()

    def reset(self):
        '''
        reset object to the initial condition
        '''
        self._coalCantera.TPY = self._coalCantera.T, self._coalCantera.P, self._compositionString


    def _referenceCoals(self):
        '''
        define the composition of the reference coals
        '''

        self._coal1 = coal(name='COAL1',c=12*self._Mc,h=11*self._Mh,o=0,n=0,s=0)
        self._coal2 = coal(name='COAL2',c=14*self._Mc,h=10*self._Mh,o=1*self._Mo,n=0,s=0)
        self._coal3 = coal(name='COAL3',c=12*self._Mc,h=12*self._Mh,o=5*self._Mo,n=0,s=0)
        self._char  = coal(name='CHAR',c=1*self._Mc,h=0,o=0,n=0,s=0)

    def _calculateCoalComposition(self):
        '''
        calculate the composition of the actual coal according to the reference coals
        '''
        self._referenceCoals()
        # determine in which triangle the coal lies
        # out means outside
        self._inside = 'out'
        # triangle 012
        t012=triangle(x0=self._char.getVanKravelen(),x1=self._coal1.getVanKravelen(),x2=self._coal2.getVanKravelen())
        if t012.isInside(self.getVanKravelen()):
            self._inside = '012'
            self.triangle = t012
            self._interpolateCoal()
            return
        t023=triangle(x0=self._char.getVanKravelen(),x1=self._coal2.getVanKravelen(),x2=self._coal3.getVanKravelen())
        if t023.isInside(self.getVanKravelen()):
            self._inside = '023'
            self.triangle = t023
            self._interpolateCoal()
            return
        t123=triangle(x0=self._coal1.getVanKravelen(),x1=self._coal2.getVanKravelen(),x2=self._coal3.getVanKravelen())
        if t123.isInside(self.getVanKravelen()):
            self._inside = '123'
            self.triangle = t123
            self._interpolateCoal()
            return
        raise compositionError('Composition outside of triangle!')

    def _interpolateCoal(self):
        '''
        interpolate coal using reference
        '''
        if self._inside == '012':
            c0 = self._char
            c1 = self._coal1
            c2 = self._coal2
        elif self._inside == '023':
            c0 = self._char
            c1 = self._coal2
            c2 = self._coal3
        elif self._inside == '123':
            c0 = self._coal1
            c1 = self._coal2
            c2 = self._coal3
        else:
            print 'composition outside of triangle'
            compositionError('Composition outside of triangle!')

        matrix = np.array(
            [[c0._c,c1._c,c2._c],
             [c0._h,c1._h,c2._h],
             [c0._o,c1._o,c2._o]])

        b = np.array([self._c,self._h,self._o])
        composition = np.linalg.solve(matrix,b)
        # define string
        if self._inside == '012':
            self._coalComposition = np.append(composition,0.)
            self._compositionString = 'CHAR:'+str(composition[0])+',COAL1:'+str(composition[1])+',COAL2:'+str(composition[2])
        elif self._inside == '023':
            self._coalComposition = np.array([composition[0],0.0,composition[1],composition[2]])
            self._compositionString = 'CHAR:'+str(composition[0])+',COAL2:'+str(composition[1])+',COAL3:'+str(composition[2])
        elif self._inside == '123':
            self._coalComposition = np.append(0.0,composition)
            self._compositionString = 'COAL1:'+str(composition[0])+',COAL2:'+str(composition[1])+',COAL3:'+str(composition[2])

    def _setCanteraObject(self,file='COAL.xml'):
        '''
        set cantera solution object
        '''
        #self._coalCantera = IdealGasMix('COAL.xml')
        self._coalCantera = cantera.Solution(file)
        #self._coalCantera.set(T=300,P=OneAtm,Y=self._compositionString)
        self._coalCantera.TPY =  300, cantera.one_atm ,self._compositionString

    def setHeatingRate(self,time=np.array([0,0.1]),temperature=np.array([400,1000])):
        '''
        define heating rate during pyrolysis
        using two NP array for time and temperature
        '''
        self.timeHR = time
        self.temperatureHR = temperature

    #def setTimeStep(self,dt=1e-3):
    #    '''
    #    define time step for printing results of pyrolysis
    #    '''
    #    self.dt = dt
    #    self.time = np.linspace(min(self.timeHR),max(self.timeHR))
    def setTimeStep(self,npoint=100):
        '''
        define time step for printing results of pyrolysis
        '''
        self.time = np.linspace(min(self.timeHR),max(self.timeHR),npoint)
        self.dt = self.time[1]-self.time[0]

    def _getInterpTemperature(self,t):
        '''
        get the interpolated temperature for a given time
        using the heating rate
        '''
        # check if the time is inside the range
        temp = np.interp(t,self.timeHR,self.temperatureHR)
        return temp

    def solvePyrolysis(self):
        '''
        solve pyrolysis
        '''
        def dmidt(t,m):
            '''
            calculate the derivative of the mass of each species i
            dm_i/dt = omega_i * Mw_i
            '''
            self._updateReactor(t,m)
            #omegai = self._coalCantera.netProductionRates()
            omegai = self._coalCantera.net_production_rates
            dmdt = omegai * self._Mw / self._rhoDry
            return dmdt

        #self._Mw = self._coalCantera.molarMasses()
        #m0=self._coalCantera.massFractions()
        ##sol = ode(dydt).set_integrator('dopri5',rtol=1e-9,atol=1e-6) #, method='bdf')
        ##sol = ode(dydt).set_integrator('vode',method='bdf',rtol=1e-9,atol=1e-6)
        #sol = ode(dmidt).set_integrator('vode',method='bdf',rtol=1e-9,atol=1e-5)
        ##sol = ode(dmidt).set_integrator('vode',method='bdf',rtol=1e-4,atol=1e-2)
        #sol.set_initial_value(m0,0)
        #self._y = [m0]
        #self._r = [dmidt(0,m0)]
        #for t in self.time[1:]:
        #    sol.integrate(t)
        #    self._y=np.concatenate((self._y, [sol.y]))
        #    self._r=np.concatenate((self._y, [dmidt(sol.t,sol.y)]))
        #    #print 'coal0='+str(sol.y[iCoal3])
        #
        #del(sol)
        backend = 'dopri5'
        self._Mw = self._coalCantera.molecular_weights
        t0 = self.timeHR[0]
        #m0=self._coalCantera.massFractions()
        m0=self._coalCantera.Y
        solver = ode(dmidt).set_integrator(backend, nsteps=1)
        solver.set_initial_value(m0, t0)
        solver._integrator.iwork[2] = -1
        self.time = [t0]
        self._y = [m0]
        self._r = [dmidt(t0,m0)]
        warnings.filterwarnings("ignore", category=UserWarning)
        timeEnd = np.max(self.timeHR)
        while solver.t < timeEnd:
            solver.integrate(timeEnd, step=True)
            self.time = np.concatenate((self.time, [solver.t]))
            self._y=np.concatenate((self._y, [solver.y]))
            self._r=np.concatenate((self._r, [dmidt(solver.t,solver.y)]))

        warnings.resetwarnings()







    def _updateReactor(self,t,m):
        ''' update reactor '''
        temp = self._getInterpTemperature(t)
        #print 'temp='+str(temp)
        pressure=self._coalCantera.P
        self._coalCantera.TPY = temp, pressure, m

    def __repr__(self):
        out = coal.__repr__(self)
        out += '\nCHAR:'+str(self._coalComposition[0])
        out += '\nCOAL1:'+str(self._coalComposition[1])
        out += '\nCOAL2:'+str(self._coalComposition[2])
        out += '\nCOAL3:'+str(self._coalComposition[3])+'\n'
        return out


    def getCoalComposition(self):
        return self._coalComposition


    def getRawCoal(self):
        return self._y[:,self._coalCantera.species_index('COAL1')] + \
               self._y[:,self._coalCantera.species_index('COAL2')] + \
               self._y[:,self._coalCantera.species_index('COAL3')]

    def getCharCoal(self):
        return self._y[:,self._coalCantera.species_index('CHAR')] +\
               self._y[:,self._coalCantera.species_index('CHARH')] +\
               self._y[:,self._coalCantera.species_index('CHARG')]

    def getCH4(self):
        return self._y[:,self._coalCantera.species_index('CH4')]

    def getH2(self):
        return self._y[:,self._coalCantera.species_index('H2')]

    def getH2O(self):
        return self._y[:,self._coalCantera.species_index('H2O')]


    def getCO(self):
        return self._y[:,self._coalCantera.species_index('CO')]
    def getCO2(self):
        return self._y[:,self._coalCantera.species_index('CO2')]
    def getCH2(self):
        return self._y[:,self._coalCantera.species_index('CH2')]

    def getTAR(self):
        return self._y[:,self._coalCantera.species_index('VTAR1')]+ \
               self._y[:,self._coalCantera.species_index('VTAR2')]+ \
               self._y[:,self._coalCantera.species_index('VTAR3')]

    def getLightGases(self):
        return self.getCO() + \
               self.getCO2() + \
               self.getH2O() + \
               self.getH2() + \
               self.getCH4() + \
               self.getCH2() + \
               self._y[:,self._coalCantera.species_index('CH3O')] + \
               self._y[:,self._coalCantera.species_index('BTX2')]

    def getMetaplast(self):
        metaplast = self._y[:,self._coalCantera.species_index('GCH2')] +\
               self._y[:,self._coalCantera.species_index('TAR1')] +\
               self._y[:,self._coalCantera.species_index('GBTX2')] +\
               self._y[:,self._coalCantera.species_index('GCH4')] +\
               self._y[:,self._coalCantera.species_index('GCOH2')] +\
               self._y[:,self._coalCantera.species_index('GCO2S')] +\
               self._y[:,self._coalCantera.species_index('GH2O')] +\
               self._y[:,self._coalCantera.species_index('GCOL')] +\
               self._y[:,self._coalCantera.species_index('TAR2')] +\
               self._y[:,self._coalCantera.species_index('GCO2TS')] +\
               self._y[:,self._coalCantera.species_index('GCOAL3')] +\
               self._y[:,self._coalCantera.species_index('GCO2')] +\
               self._y[:,self._coalCantera.species_index('TAR3')] +\
               self._y[:,self._coalCantera.species_index('GCOLS')]
        return metaplast


    def getVolatile(self):
        """

        :rtype : object
        """
        return self.getTAR() + self.getLightGases()


    def getTime(self):
        return self.time

    def getTemperature(self):
        return self._getInterpTemperature(self.time)


    def _calculateYields(self):
        '''
        return yields as daf fraction
        using format required by PKP
        '''
        self.__yields=np.zeros(( int(len(self.time)),14) )  #shapes new Matrix containing all necessary information;
        self.__yields[:,0]=self.time
        self.__yields[:,1]=self.getTemperature() # temp
        self.__yields[:,2]=self.getTAR() # total volatile yield
        self.__yields[:,3]=self.getLightGases()
        self.__yields[:,5]=self.getVolatile()
        self.__yields[:,4]=1.-self.__yields[:,5]
        self.__yields[:,6]=self.getH2O()
        self.__yields[:,7]=self.getCO2()
        self.__yields[:,8]=self.getCH4()
        self.__yields[:,9]=self.getCO()
        self.__yields[:,10]=self.getH2()
        self.__yields[:,11]=self._y[:,self._coalCantera.species_index('CH2')]
        self.__yields[:,12]=self._y[:,self._coalCantera.species_index('CH3O')]
        self.__yields[:,13]=self._y[:,self._coalCantera.species_index('BTX2')]
        self.Yields2Cols={'Time':0,'Temp':1,'Tar':2,'Gas':3,'Solid':4,'Total':5,
                          'H2O':6,'CO2':7,'CH4':8,'CO':9,'H2':10,
                          'CH2':11,'CH3O':12,'BTX2':13}
        self.Cols2Yields={0:'Time',1:'Temp',2:'Tar',3:'Gas',4:'Solid',5:'Total',
                          6:'H2O',7:'CO2',8:'CH4',9:'CO',10:'H2',
                          11:'CH2',12:'CH3O',13:'BTX2'}

    def Yields_all(self):
        """Returns the whole result matrix of the yields."""
        self._calculateYields()
        return self.__yields

    def Rates_all(self):
        """Returns the whole result matrix of the rates"""
        self.__rates=np.zeros(( int(len(self.time)),14) )  #shapes new Matrix containing all necessary information;
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
        return 'PMSKD'




