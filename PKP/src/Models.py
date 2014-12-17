"""
A selection of different Pyrolysis models

A model pyrolysis model usually provides method to calculate the yield of 
individual species the rates and energy balancing methods

"""
import numpy as np
import pylab as plt
import scipy as sp
import scipy.integrate
import scipy.interpolate
import platform

class BalancedComposition(object):
    """ Class for compostion that ensures componenents sum up to a 
        certain value 
    """

    def __init__(self, inp, target=100.00):
        """ From a given input dictionary and a target sum a 
            dictionary with scaled composition is created 
        """
        self.target = target
        scaling_factor = target/sum(inp.values())
        self.elems = {key:value*scaling_factor for key,value in inp.iteritems()} 


    def __getitem__(self,item):
        """ """
        try:
            return self.elems[item]
        except:
            print """Warning trying to access {} which was not set in
            input file, assuming 0.0. 
            """ 
            return 0.0

    def remove_elems_rebalance(self, elems_):
        """ To compute daf composition elements can be removed 

            Usage:  .remove_elems_rebalance(['Moisture','Ash'])
        """
        return BalancedComposition({ 
            key:elem for key,elem in self.elems.iteritems() 
                if key not in elems_
        })

    def scale(self, factor):
        return BalancedComposition(self.elems, target=self.target*factor)

    def __repr__(self):
        return str(self.elems)


class Model(object):
    """ Parent class of the children ConstantRateModel, 
        the three Arrhenius Models (notations) and the Kobayashi models. 

        TimeVectorToInterplt allows the option to define the discrete time points, 
        where to interpolate the results. If set to False (standard), then is are 
        the outputed results equal the dt to solve the ODE. If set 
        TimeVectorToInterplt=[t0,t1,t2,t3,t4] (t: floats) then is the 
        yields result returned at method calcMass the yields at [t0,t1,t2,t3,t4], 
        linear interploated."""
    
    def __init__(self, name):
        print "Initialised {} Model".format(name)
        self.name = name
    
    def computeTimeDerivative(self, mass, deltaT=False, times=False):
        """ Return time derivatives for a given deltat array
            and a mass array
        """
        from numpy import gradient
        # rate[0]    = (mass[1]  - mass[0])/deltat[0]
        # rate[1:-1] = (mass[2:] - mass[:-2])/(2*deltat[1:-1])
        # rate[-1]   = (mass[-1] - mass[-2])/deltat[-1]
        if deltaT:
            return gradient(mass,deltat)
        else:
            return gradient(np.array([mass,times]))
        
    def calcRate(self, preProcResult, time, temp, species):
        """ computes actual release reates for a given species
            by a time array and temperature array

            the preProcResults are used for initial values
            #TODO GO can time be taken from preProcResult?
        """
        mass = self.calcMass(preProcResult, time, temp, species)
        return self.computeTimeDerivative(mass, times = time)


    @classmethod
    def modelErrori(cls, target, model):
        """ compute the deviation between modeled values and the target values
            from the pre processor
        """
        return target - model

    @classmethod
    def modelErrorSquared(cls, target, model):
        """ compute the deviation between modeled values and the target values
            from the pre processor
        """
        return np.power(target - model, 2.0)

    @classmethod
    def summedModelError(cls, target, model):
        return np.sum(Model.modelError(target, model))

    @classmethod
    def totModelErrorSquaredPerc(cls, target, model):
        return np.sum(Model.modelErrorSquared(target, model))/len(target)


    @classmethod
    def yieldDelta(cls, mass):
        return max(mass) - min(mass)

    def ErrorYield(self, preProcResult, Species):
        """ Returns the absolute deviation per point between the fitted
            pyrolysis model and the original preprocessor yield 
        TODO GO: why only take some arrays from preProcResults
        """
        modeled_mass = self.calcMass(
                            preProcResult,
                            preProcResult["Time"],
                            preProcResult["Temp"],
                            Species
                        )
        target_mass = preProcResult[Species]
        return (np.sum(target_mass - modeled_mass))/(len(target_mass))

    def ErrorRate(self, preProcResult, species):
        """ Returns the absolute deviation per point between 
            the fitted and the original rate curve.
        """
        mass = self.calcMass(
                preProcResult, 
                preProcResult.Time(),
                preProcResult.Interpolate('Temp'),
                species)

        uDot = preProcResult.Rate(Species)
        mDot = self.deriveC(preProcResult, mass)
        absE = (np.sum(uDot-vDot))/preProcResult.NPoints()
        return absE
    
    # def mkSimpleResultFiles(self,preProcResult_list,Species):
    #     """Simple result file if no fitting is carried out. Writes only the transformed results into a file."""
    #     if type(Species)==int:
    #         SpeciesForTitle=preProcResult_list[0].SpeciesName(Species)
    #     if type(Species)==str:
    #         SpeciesForTitle=Species
    #     #
    #     for runnedCaseNr in range(len(preProcResult_list)):
    #         t=preProcResult_list[runnedCaseNr].Yield('Time')
    #         T=preProcResult_list[runnedCaseNr].Yield('Temp')
    #         u=[]
    #         for Nr in range(len(preProcResult_list)):
    #             u_=preProcResult_list[Nr].Yield(Species)
    #             u.append(u_)
    #         w=self.deriveC(preProcResult_list[runnedCaseNr],u[runnedCaseNr])
    #         if oSystem=='Linux':
    #             resultFile=open('Result/'+preProcResult_list[runnedCaseNr].Name()+'-Fit_result_'+SpeciesForTitle+'_'+str(runnedCaseNr)+'.out','w')
    #         elif oSystem=='Windows':
    #             resultFile=open('Result\\'+preProcResult_list[runnedCaseNr].Name()+'-Fit_result_'+SpeciesForTitle+'_'+str(runnedCaseNr)+'.out','w')
    #         else:
    #             print 'Models: Operating Platform cannot be specified.'
    #         resultFile.write('    Time       Temperature Yields(original) Rates(original) \n')
    #         for i in range(len(u[runnedCaseNr])):
    #             resultFile.write('%7e  %11e %7e %8e \n' % (t[i], T[i], u[runnedCaseNr][i], w[i]))
    #         resultFile.close()
            
    def _mkInterpolatedRes(self,InputVecYields,Time):
        """ Generates the result vector. Outputs the result vector at the 
            corresponding time steps corresponding to the imported time at 
            method calcMass. Requiered for Pc Coal Lab (only few reported points).
        """
        # t for for interplt, t_points, y_points
        return np.interp(Time,self.constDtVec,InputVecYields)
    
    def setDt4Intergrate(self,constantDt):
     """ constantDt allows the option to define numerical time step to solve the ODE. 
        
        The outputted results ever equal the imported time list 
        (when applying method calcMass Time = [t0,t1,t2,t3,t4]. If these time steps 
        are too large, then is this defined dt used to solve the ODE and the results 
        are linear interploated that way that they correspond to the imported time 
        vector. To reset it, just set constantDt to False.
        """
     if constantDt != False:
         self.constDt = float(constantDt)
    
    def _mkDt4Integrate(self,Time):
        """ Time is the original time vector calculated by exact model, e.g. CPD. 
            This class generates the internal dt vector if the dt defined by the 
            user file is too large. A time step must be defined in Method 
            setDt4Intergrate before.
        """
        if self.constDt != False:
         self.constDtVec = np.arange(Time[0],Time[-1],self.constDt)

################childrenclasses####################

class ConstantRateModel(Model):
    """ The model calculating the mass 
        with m(t)  =  m_s0+(m_s0-m_s,e)*exp(-k*(t-t_start)) 
             dm/dt = -k*(m-m_s,e). 
        The Parameter to optimize are k and t_start. 
    """
    #TODO store initial parameter to see test if parameters have been changed
    #TODO GO parameter[2] is not specified in inputs example

    def __init__(self, parameter):
        Model.__init__(self,"ConstantRate")
        self.initialParameter = parameter
        self.k                = parameter["k"] 
        self.start_time       = parameter["tstart"]
        self.final_yield      = parameter.get('finalYield',False)
        # if set to false, the numerical time step corresponding to the outputed 
        # by the detailled model (e.g CPD) is used; define a value to use instead this 
        self.constDt = False 

    def updateParameter(self, parameter):
        self.k          = parameter[0] 
        self.start_time = parameter[1]
        ##print "Parameter update " + str(parameter)

    def calcMass(self, init_mass, time, temp=False):
        """
            Inputs: 
                time an array of time values

            the function might get an additional temp argument 
            which captured in *args or **kwargs
        """
        # we care only about time value
        # starting at start time
        time = time - self.start_time 
        # the yield still retained in the coal
        # this should converge to zero at large 
        # time
        print self.k
        retained_mass = self.final_yield * np.exp(-self.k*time)
        released_mass = self.final_yield - retained_mass 

        # if we are interested in the solid mass
        # TODO GO what is this solid thing going on here 
        if False: #species == 'Solid':
            released_mass += solid_mass*np.exp(-self.k*time)

        # why choosing between released or solid mass 
        # start_time is small then time
        # released_mass = np.where(time > self.start_time, released_mass, solid_mass) 
        if self.constDt == False: # TODO GO shouldnt interpolation be used for var dt?
            #print "modeled_mass " + str(released_mass)
            return released_mass
        else: #returns the short, interpolated list (e.g. for PCCL)
            return self._mkInterpolatedRes(released_mass, time)
    
    @property
    def parameter(self):
        return np.array([self.k,self.start_time])
                

class ArrheniusModel(Model):
    """The Arrhenius model in the standart notation: dm/dt=A*(T**b)*exp(-E/T)*(m_s-m) with the parameter a,b,E to optimize."""
    def __init__(self,InitialParameterVector):
        print 'Arrhenuis Model initialized'
        self._modelName = 'Arrhenius'
        self._ParamVector=InitialParameterVector
        self.ODE_hmax=1.e-2
        self.constDt = False # if set to false, the numerical time step corresponding to the outputted by the dtailled model (e.g CPD) is used; define a value to use instead this 
 
    def calcMass(self,preProcResult,time,T,Name):
        """Outputs the mass(t) using the model specific equation."""
        #numercal values:
        absoluteTolerance = 1.0e-8
        relativeTolerance = 1.0e-6
        ##################
        u=preProcResult.Yield(Name)
        ParamVec=self.ParamVector()
        m_s0=ParamVec[3]
        # question whether the dt from DetailledModel result file or from a constant dt should be used
        if self.constDt == False: # dt for integrate = dt from DM result file
            timeInt = time
        else: #if dt in DM results file has too large dt
            self._mkDt4Integrate(time)
            timeInt = self.constDtVec
        def dmdt(m,t):
            if Name == 'Solid':# or Name == preProcResult.Yields2Cols['Solid']:
                dmdt_out=-ParamVec[0]*(T(t)**ParamVec[1])*np.exp(-ParamVec[2]/(T(t)))*(m-m_s0)
                dmdt_out=np.where(abs(dmdt_out)>1.e-300,dmdt_out,0.0) #sets values<0 =0.0, otherwise it will further cause problems (nan)
            else:
                dmdt_out= ParamVec[0]*(T(t)**ParamVec[1])*np.exp(-ParamVec[2]/(T(t)))*(m_s0-m)
                dmdt_out=np.where(abs(dmdt_out)>1.e-300,dmdt_out,0.0) #sets values<0 =0.0, otherwise it will further cause problems (nan)
            return dmdt_out
        InitialCondition=[u[0]]
        m_out=sp.integrate.odeint(dmdt,InitialCondition,timeInt,atol=absoluteTolerance,rtol=relativeTolerance,hmax=self.ODE_hmax) 
        if self.constDt == False:
            if (ParamVec[0]<0 or ParamVec[2]<0):
                m_out[:,0]=float('inf')
                return m_out[:,0]
            else:
                return m_out[:,0]
        else: #returns the short, interpolated list (e.g. for PCCL)
            return self._mkInterpolatedRes(m_out[:,0],time)

    def ConvertKinFactors(self,ParameterVector):
        """Dummy. Function actual has to convert the parameter into the standart Arrhenius notation."""
        #does nothing, just to have the same way of use for all notations
        return ParameterVector

                
class ArrheniusModelNoB(Model):
    """The Arrhenius model in the standart notation: dm/dt=A*exp(-E/T)*(m_s-m) with the parameter a,b,E to optimize."""
    def __init__(self,InitialParameterVector):
        print 'Arrhenuis Model initialized'
        self._modelName = 'ArrheniusNoB'
        self._ParamVector=InitialParameterVector
        self.ODE_hmax=1.e-2
        self.constDt = False # if set to false, the numerical time step corresponding to the outputted by the dtailled model (e.g CPD) is used; define a value to use instead this 
 
    def calcMass(self,preProcResult,time,T,Name):
        """Outputs the mass(t) using the model specific equation."""
        #numercal values:
        absoluteTolerance = 1.0e-8
        relativeTolerance = 1.0e-6
        ##################
        u=preProcResult.Yield(Name)
        ParamVec=self.ParamVector()
        m_s0=ParamVec[2]
        # question whether the dt from DetailledModel result file or from a constant dt should be used
        if self.constDt == False: # dt for integrate = dt from DM result file
            timeInt = time
        else: #if dt in DM results file has too large dt
            self._mkDt4Integrate(time)
            timeInt = self.constDtVec
        def dmdt(m,t):
            if Name == 'Solid':# or Name == preProcResult.Yields2Cols['Solid']:
                dmdt_out=-ParamVec[0]*np.exp(-ParamVec[1]/(T(t)))*(m-m_s0)
                dmdt_out=np.where(abs(dmdt_out)>1.e-300,dmdt_out,0.0) #sets values<0 =0.0, otherwise it will further cause problems (nan)
            else:
                dmdt_out= ParamVec[0]*np.exp(-ParamVec[1]/(T(t)))*(m_s0-m)
                dmdt_out=np.where(abs(dmdt_out)>1.e-300,dmdt_out,0.0) #sets values<0 =0.0, otherwise it will further cause problems (nan)
            return dmdt_out
        InitialCondition=[u[0]]
        m_out=sp.integrate.odeint(dmdt,InitialCondition,timeInt,atol=absoluteTolerance,rtol=relativeTolerance,hmax=self.ODE_hmax) 
        if self.constDt == False:
            if (ParamVec[0]<0 or ParamVec[1]<0):
                m_out[:,0]=float('inf')
                return m_out[:,0]
            else:
                return m_out[:,0]
        else: #returns the short, interpolated list (e.g. for PCCL)
            return self._mkInterpolatedRes(m_out[:,0],time)

    def ConvertKinFactors(self,ParameterVector):
        """Dummy. Function actual has to convert the parameter into the standart Arrhenius notation."""
        #does nothing, just to have the same way of use for all notations
        return ParameterVector
        
        
class ArrheniusModelAlternativeNotation1(ArrheniusModel):
    """Arrhenius model with a notation having a better optimization behaviour: dm/dt=exp[k0-a*(T0/T(t)-1)]*(ms-m). See the documentation for the reference. The parameters to optimize are k0 and a."""
    def __init__(self,InitialParameterVector):
        print 'Arrhenuis Model class initialized'
        self._ParamVector=InitialParameterVector
        self.ODE_hmax=1.e-2
        self.constDt = False # if set to false, the numerical time step corresponding to the outputted by the dtailled model (e.g CPD) is used; define a value to use instead this 
        
    def calcMass(self,preProcResult,time,T,Name):
        """Outputs the mass(t) using the model specific equation."""
        #numercal values:
        absoluteTolerance = 1.0e-8
        relativeTolerance = 1.0e-6
        ##################
        ParamVec=self.ParamVector()
        u=preProcResult.Yield(Name)
        m_s0=ParamVec[2]
        twhereTmax=time[preProcResult.LineNumberMaxRate(Name)]
        self.T0=T(twhereTmax)
        # question whether the dt from DetailledModel result file or from a constant dt should be used
        if self.constDt == False: # dt for integrate = dt from DM result file
            timeInt = time
        else: #if dt in DM results file has too large dt
            self._mkDt4Integrate(time)
            timeInt = self.constDtVec
        def dmdt(m,t):
            if Name == 'Solid':# or Name == preProcResult.Yields2Cols['Solid']:
                dmdt_out=-np.exp( ParamVec[0]  - ParamVec[1]*( self.T0/T(t) - 1 ))*(m-m_s0) #-ParamVec[0]*( (T(t)/self.T0)**ParamVec[1] )*np.exp( -ParamVec[2]*(self.T0/(T(t))-1) )*(m)    #IC
                dmdt_out=np.where(abs(dmdt_out)>1.e-300,dmdt_out,0.0) #sets values<0 =0.0, otherwise it will further cuase problem(nan)
            else:
                dmdt_out= np.exp( ParamVec[0] - ParamVec[1]*( self.T0/T(t) - 1 ))*(m_s0-m)
            return dmdt_out
        InitialCondition=[u[0]]
        m_out=sp.integrate.odeint(dmdt,InitialCondition,timeInt,atol=absoluteTolerance,rtol=relativeTolerance)
        if self.constDt == False:
            return m_out[:,0]
        else: #returns the short, interpolated list (e.g. for PCCL)
            return self._mkInterpolatedRes(m_out[:,0],time)

    def ConvertKinFactors(self,ParameterVector):
        """Converts the own kinetic factors back to the standard Arrhenius kinetic factors."""
        #k=ParameterVector[0]
        #a=ParameterVector[1]
        A=np.exp( ParameterVector[0] + ParameterVector[1] )
        E=ParameterVector[1]*self.T0
        return [A,0.0,E,ParameterVector[2]]

    def ConvertKinFactorsToOwnNotation(self,preProcResult,ParameterVector,Species):
        """Converts the standard Arrhenius kinetic factors backk to the factors of the own notation."""
        time=preProcResult.Time()
        twhereTRmax=time[preProcResult.LineNumberMaxRate(Species)]
        T=preProcResult.Interpolate('Temp')
        self.T0=T(twhereTRmax)
        A=ParameterVector[0]
        if ParameterVector[1]!=0:
            print 'b is not considered in this Notation and set to 0'
        E=ParameterVector[2]
        a=E/self.T0
        k0=np.log(A)-a
        return [k0,a,ParameterVector[3]]

        
class ArrheniusModelAlternativeNotation2(ArrheniusModel):
    """Arrhenius model with a notation having a better optimization behaviour: dm/dt=exp[c*(b1*(1/T(t)-1/T_min)-b2*(1/T(t)-1/T_max))]*(ms-m); with c=(1/T_max-1/Tmin)**(-1). See the documentation for the reference. The parameters to optimize are b1 and b2."""
    def __init__(self,InitialParameterVector):
        print 'Arrhenius Model class initialized'
        self._ParamVector=InitialParameterVector
        self.T_min=300
        self.T_max=1500
        self.constDt = False # if set to false, the numerical time step corresponding to the outputted by the dtailled model (e.g CPD) is used; define a value to use instead this 
        
    def setMinMaxTemp(self,Tmin,Tmax):
        """Sets the temperature constants, see the equation."""
        self.T_min=Tmin
        self.T_max=Tmax
        
    def calcMass(self,preProcResult,time,T,Name):
        """Outputs the mass(t) using the model specific equation."""
        self.c=1./(1./self.T_max-1./self.T_min)
        self.ODE_hmax=1.e-2
        self.preProcResult=preProcResult
        #numercal values:
        absoluteTolerance = 1.0e-8
        relativeTolerance = 1.0e-6
        ##################
        ParamVec=self.ParamVector()
        u=preProcResult.Yield(Name)
        #uDot=preProcResult.Rate(Name)
        m_s0=ParamVec[2]
        # question whether the dt from DetailledModel result file or from a constant dt should be used
        if self.constDt == False: # dt for integrate = dt from DM result file
            timeInt = time
        else: #if dt in DM results file has too large dt
            self._mkDt4Integrate(time)
            timeInt = self.constDtVec
        def dmdt(m,t):
            if Name == 'Solid':# or Name == preProcResult.Yields2Cols['Solid']:
                dmdt_out= -np.exp( self.c*( ParamVec[0]*(1./T(t)-1./self.T_min) - ParamVec[1]*(1./T(t)-1./self.T_max) ) )*(m-m_s0)
                dmdt_out=np.where(abs(dmdt_out)>1.e-300,dmdt_out,0.0) #sets values<0 =0.0, otherwise it will further cause problem(nan)
            else:
                dmdt_out= np.exp( self.c*( ParamVec[0]*(1./T(t)-1./self.T_min) - ParamVec[1]*(1./T(t)-1./self.T_max) ) )*(m_s0-m)
                dmdt_out=np.where(abs(dmdt_out)>1.e-300,dmdt_out,0.0) #sets values<0 =0.0, otherwise it will further cause problem(nan)
            return dmdt_out
        InitialCondition=[u[0]]
        m_out=sp.integrate.odeint(dmdt,InitialCondition,timeInt,atol=absoluteTolerance,rtol=relativeTolerance,hmax=self.ODE_hmax)
        ArrStandartParam=self.ConvertKinFactors(ParamVec)
        if self.constDt == False:
            if (ArrStandartParam[0]<0 or ArrStandartParam[2]<0):
                m_out[:,0]=float('inf')
                return m_out[:,0]
            else:
                return m_out[:,0]
        else: #returns the short, interpolated list (e.g. for PCCL)
            return self._mkInterpolatedRes(m_out[:,0],time)

    def ConvertKinFactors(self,ParameterVector):
        """Converts the own kinetic factors back to the standard Arrhenius kinetic factors."""
        #b1=ParameterVector[0]
        #b2=ParameterVector[1]
        self.c=1./(1./self.T_max-1./self.T_min)
        A=np.exp( -self.c*ParameterVector[0]/self.T_min + self.c*ParameterVector[1]/self.T_max )
        E=self.c*ParameterVector[1]-self.c*ParameterVector[0]
        m_s0=ParameterVector[2]
        return [A,0,E,m_s0]

    def ConvertKinFactorsToOwnNotation(self,ParameterVector):
        """Converts the standard Arrhenius kinetic factors back to the factors of the own notation."""
        self.c=1/(1/self.T_max-1/self.T_min)
        #A=ParameterVector[0]
        #b not used
        #E=ParameterVector[2]
        if ParameterVector[1]!=0:
            print 'b is not considered in this Notation and set to 0'
        b2=( (self.T_max/self.c)*np.log(ParameterVector[0])-(ParameterVector[2]*self.T_max)/(self.c*self.T_min) )*(1-self.T_max/self.T_min)**(-1)
        b1=b2-ParameterVector[2]/self.c
        return [b1,b2,ParameterVector[3]]

class Kobayashi(Model):
    """Calculates the devolatilization reaction using the Kobayashi model. The Arrhenius equation inside are in the standard notation."""
    def __init__(self,InitialParameterVector):
        print 'Kobayashi Model initialized'
        self._modelName = 'Kobayashi'
        self._ParamVector=InitialParameterVector
        self.ODE_hmax=1.e-2
        self.constDt = False # if set to false, the numerical time step corresponding to the outputted by the dtailled model (e.g CPD) is used; define a value to use instead this 
        
    def calcMass(self,preProcResult,time,T,Name):
        """Outputs the mass(t) using the model specific equation. The input Vector is [A1,E1,A2,E2,alpha1,alpha2]"""
        # question whether the dt from DetailledModel result file or from a constant dt should be used
        if self.constDt == False: # dt for integrate = dt from DM result file
            timeInt = time
        else: #if dt in DM results file has too large dt
            self._mkDt4Integrate(time)
            timeInt = self.constDtVec
        self.preProcResult=preProcResult
        timeInt=np.delete(timeInt,0)
        self.__Integral=0.0
        tList=[0.0]
        k1k2=[0.0]
        ParamVec=self.ParamVector()
        #
        def dmdt(m,t):
            k1=ParamVec[0]*np.exp(-ParamVec[1]/(T(t)))
            k2=ParamVec[2]*np.exp(-ParamVec[3]/(T(t)))
            tList.append(t)
            k1k2.append(k1+k2)
            self.__Integral+=0.5*(tList[-1]-tList[-2])*(k1k2[-1]+k1k2[-2])
            dmdt_out = ( (ParamVec[4]*k1+ParamVec[5]*k2)*np.exp(-self.__Integral) )
            dmdt_out=np.where(abs(dmdt_out)>1.e-300,dmdt_out,0.0) #sets values<0 =0.0, otherwise it will further cause problem(nan)
            return dmdt_out
        InitialCondition=[0]
        m_out=sp.integrate.odeint(dmdt,InitialCondition,timeInt,atol=1.e-5,rtol=1.e-4,hmax=1.e-2)
        m_out=m_out[:,0]
        m_out=np.insert(m_out,0,0.0)
        if self.constDt == False:
            if (ParamVec[0]<0 or ParamVec[1]<0 or ParamVec[2]<0 or ParamVec[3]<0 or ParamVec[4]<0  or ParamVec[5]>1  ):
                m_out[:]=float('inf')
                return m_out
            else:
                return m_out
        else: #returns the short, interpolated list (e.g. for PCCL)
            return self._mkInterpolatedRes(m_out,time)
        

class KobayashiPCCL(Model):
    """Calculates the devolatilization reaction using the Kobayashi model. The Arrhenius equation inside are in the standard notation. The fitting parameter are as in PCCL A1,A2,E1,alpha1. TimeVectorToInterplt allows the option to define the discrete time points, where to interpolate the results. If set to False (standard), then is are the outputted results equal the dt to solve the ODE."""
    def __init__(self,InitialParameterVector):
        print 'Kobayashi Model initialized'
        self._modelName = 'KobayashiPCCL'
        self._ParamVector=InitialParameterVector
        self.ODE_hmax=1.e-2
        self.constDt = False # if set to false, the numerical time step corresponding to the outputted by the dtailled model (e.g CPD) is used; define a value to use instead this 
        
    def calcMass(self,preProcResult,time,T,Name):
        """Outputs the mass(t) using the model specific equation."""
        # question whether the dt from DetailledModel result file or from a constant dt should be used
        if self.constDt == False: # dt for integrate = dt from DM result file
            timeInt = time
        else: #if dt in DM results file has too large dt
            self._mkDt4Integrate(time)
            timeInt = self.constDtVec
        self.preProcResult=preProcResult
        timeInt=np.delete(timeInt,0)
        self.__Integral=0.0
        tList=[0.0]
        k1k2=[0.0]
        ParamVec=self.ParamVector()
        #
        def dmdt(m,t):
            k1=ParamVec[0]*np.exp(-ParamVec[2]/(T(t)))
            k2=ParamVec[1]*np.exp(-(ParamVec[2]+self.__E2diff)/(T(t)))
            tList.append(t)
            k1k2.append(k1+k2)
            self.__Integral+=0.5*(tList[-1]-tList[-2])*(k1k2[-1]+k1k2[-2])
            dmdt_out = ( (ParamVec[3]*k1+self.__alpha2*k2)*np.exp(-self.__Integral) )
            dmdt_out=np.where(abs(dmdt_out)>1.e-300,dmdt_out,0.0) #sets values<0 =0.0, otherwise it will further cause problem(nan)
            return dmdt_out
        InitialCondition=[0]
        m_out=sp.integrate.odeint(dmdt,InitialCondition,timeInt,atol=1.e-5,rtol=1.e-4,hmax=1.e-2)
        m_out=m_out[:,0]
        m_out=np.insert(m_out,0,0.0)
        if self.constDt == False:
            if (ParamVec[0]<0 or ParamVec[1]<0 or ParamVec[2]<0 or ParamVec[3]<0):
                m_out[:]=float('inf')
                return m_out
            else:
                return m_out
        else: #returns the short, interpolated list (e.g. for PCCL)
            return self._mkInterpolatedRes(m_out,time)
        

    def ConvertKinFactors(self,ParameterVector):
        """Outputs the Arrhenius equation factors in the shape [A1,E1,A2,E2]. Here where the real Arrhenius model is in use only a dummy function."""
        P=self.ParamVector()
        return [P[0],P[1],P[2],P[3]]
    
    def setKobWeights(self,alpha2):
        """Sets the two Kobayashi weights alpha2."""
        self.__alpha2=alpha2
        
    def KobWeights(self):
        """Returns the two Kobayashi weights alpha2."""
        return self.__alpha2
    
    def setE2Diff(self,DifferenceE1E2):
        """Sets the dE in E2=E1+dE."""
        self.__E2diff=DifferenceE1E2
        
    def E2Diff(self):
        """Returns the dE in E2=E1+dE."""
        return self.__E2diff

class KobayashiA2(Model):
    """Calculates the devolatilization reaction using the Kobayashi model. The Arrhenius equation inside are in the secend alternative notation (see class ArrheniusModelAlternativeNotation2)."""
    def __init__(self,InitialParameterVector):
        print 'Kobayashi Model initialized'
        self._ParamVector=InitialParameterVector
        self.ODE_hmax=1.e-2
        self.constDt = False # if set to false, the numerical time step corresponding to the outputted by the dtailled model (e.g CPD) is used; define a value to use instead this 

    def calcMass(self,preProcResult,time,T,Name):
        """Outputs the mass(t) using the model specific equation."""
        # question whether the dt from DetailledModel result file or from a constant dt should be used
        if self.constDt == False: # dt for integrate = dt from DM result file
            timeInt = time
        else: #if dt in DM results file has too large dt
            self._mkDt4Integrate(time)
            timeInt = self.constDtVec
        self.preProcResult=preProcResult
        T_general=preProcResult.Rate('Temp')
        self.T_min=min(T_general)
        self.T_max=max(T_general)
        self.c=1./(1./self.T_max-1./self.T_min)
        #alpha has to be greater equal zero:
        absoluteTolerance = 1.0e-8
        relativeTolerance = 1.0e-6
        self.tList=[0.0]
        #deletes to solve trapezian rule:
        timeInt=np.delete(timeInt,0)
        self.Integral=0.0
        self.k1k2=[0.0]
        ParamVec=self.ParamVector()
        u=preProcResult.Yield(Name)
        #
        def dmdt(m,t):
            k1=np.exp( self.c*( ParamVec[0]*(1./T(t)-1./self.T_min) - ParamVec[1]*(1./T(t)-1./self.T_max) ) )
            k2=np.exp( self.c*( ParamVec[2]*(1./T(t)-1./self.T_min) - ParamVec[3]*(1./T(t)-1./self.T_max) ) )
            self.tList.append(t)
            self.k1k2.append(k1+k2)
            self.Integral+=0.5*(self.tList[-1]-self.tList[-2])*(self.k1k2[-1]+self.k1k2[-2])
            dmdt_out = ( (self.__alpha1*k1+self.__alpha2*k2)*np.exp(-self.Integral) )
            dmdt_out=np.where(abs(dmdt_out)>1.e-300,dmdt_out,0.0) #sets values<0 =0.0, otherwise it will further cause problem(nan)
            return dmdt_out
        InitialCondition=[u[0]]
        m_out=sp.integrate.odeint(dmdt,InitialCondition,timeInt,atol=absoluteTolerance,rtol=relativeTolerance,hmax=self.ODE_hmax)
        m_out=m_out[:,0]
        m_out=np.insert(m_out,0,0.0)
        if self.constDt == False:
            return m_out
        else: #returns the short, interpolated list (e.g. for PCCL)
            return self._mkInterpolatedRes(m_out,time)

    def ConvertKinFactors(self,ParameterVector):
        """Converts the alternative notaion Arrhenius factors into the satndard Arrhenius factors and return them in the shape  [A1,E1], [A2,E2]"""
        A1=np.exp( -self.c*ParameterVector[0]/self.T_min + self.c*ParameterVector[1]/self.T_max )
        E1=self.c*ParameterVector[1]-self.c*ParameterVector[0]
        A2=np.exp( -self.c*ParameterVector[2]/self.T_min + self.c*ParameterVector[3]/self.T_max )
        E2=self.c*ParameterVector[3]-self.c*ParameterVector[2]
        return [A1,E1,A2,E2]

    def setKobWeights(self,alpha1,alpha2):
        """Sets the two Kobayashi weights alpha1 and alpha2."""
        self.__alpha1=alpha1
        self.__alpha2=alpha2
        
    def KobWeights(self):
        """Returns the two Kobayashi weights alpha1 and alpha2."""
        return self.__alpha1, self.__alpha2
        

class DAEM(Model):
    """Calculates the devolatilization reaction using the Distributed Activation Energy Model."""
    def __init__(self,InitialParameterVector):
        print 'DAEM initialized'
        self._modelName = 'DAEM'
        self._ParamVector=InitialParameterVector
        self.ODE_hmax=1.e-2
        self.NrOfActivationEnergies=50
        self.constDt = False # if set to false, the numerical time step corresponding to the outputted by the dtailled model (e.g CPD) is used; define a value to use instead this 
    
    def setNrOfActivationEnergies(self,NrOfE):
        """Define for how many activation energies of the range of the whole distribution the integral shall be solved (using Simpson Rule)."""
        self.NrOfActivationEnergies=NrOfE
        
    def NrOfActivationEnergies(self):
        """Returns the number of activation enrgies the integral shall be solved for (using Simpson Rule)."""
        return self.NrOfActivationEnergies
        
    def calcMass(self,preProcResult,time,T,Name):
        """Outputs the mass(t) using the model specific equation."""
        self.E_List=np.arange(int(self._ParamVector[1]-3.*self._ParamVector[2]),int(self._ParamVector[1]+3.*self._ParamVector[2]),int((6.*self._ParamVector[2])/self.NrOfActivationEnergies)) #integration range E0 +- 3sigma, see [Cai 2008]
        # question whether the dt from DetailledModel result file or from a constant dt should be used
        if self.constDt == False: # dt for integrate = dt from DM result file
            timeInt = time
        else: #if dt in DM results file has too large dt
            self._mkDt4Integrate(time)
            timeInt = self.constDtVec
        #Inner Integral Funktion
        def II_dt(t,E_i):
            return np.exp( -E_i/T(t) )
        #outer Integral for one activation energy from t0 to tfinal
        #stores all values of the inner Integrals (time,ActivationEnergy) in a 2D-Array
        InnerInts=np.zeros([len(timeInt),len(self.E_List)])
        CurrentInnerInt=np.zeros(len(timeInt))
        for Ei in range(len(self.E_List)):
            CurrentInnerInt[:]=II_dt(timeInt[:],self.E_List[Ei])
            InnerInts[1:,Ei] = sp.integrate.cumtrapz(CurrentInnerInt,timeInt[:])
        #
        def OI_dE(EIndex,tIndex):
            m = np.exp(-self._ParamVector[0]*InnerInts[tIndex,EIndex])*(1./(self._ParamVector[2]*(2.*np.pi)**0.5))*np.exp(-(self.E_List[EIndex]-self._ParamVector[1])**2/(2.*self._ParamVector[2]**2)) 
#            print 'InnerInt',InnerInt,'mass',dm_dt
            return m
        m_out=np.zeros(np.shape(timeInt))
        mE=np.zeros(np.shape(self.E_List))
        for ti in range(len(timeInt)):
            for Ei in range(len(self.E_List)):
                mE[Ei]=OI_dE(Ei,ti)
            m_out[ti]=sp.integrate.simps(mE,self.E_List)
        #descaling
        m_out = self._ParamVector[3]*(1.-m_out)
        if self.constDt == False:
            return m_out
        else: #returns the short, interpolated list (e.g. for PCCL)
            return self._mkInterpolatedRes(m_out,time)
