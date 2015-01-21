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

    def __init__(self, name, parameter, constDt=False):
        print "Initialised {} Model".format(name)
        self.name = name
        self.initialParameter = parameter
        self.constDt = constDt

    def computeTimeDerivative(self, mass, deltaT=False, times=False):
        """ Return time derivatives for a given deltat array
            and a mass array
        """
        from numpy import gradient
        if deltaT:
            return gradient(mass, deltat)
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
        pass

    def ErrorRate(self, preProcResult, species):
        """ Returns the absolute deviation per point between
            the fitted and the original rate curve.
        """
        pass

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

    def _mkDt4Integrate(self, Time):
        """ Time is the original time vector calculated by exact model, e.g. CPD.
            This class generates the internal dt vector if the dt defined by the
            user file is too large. A time step must be defined in Method
            setDt4Intergrate before.
        """
        if self.constDt != False:
         self.constDtVec = np.arange(Time[0], Time[-1], self.constDt)

class ModelError(object): 
    # TODO merge this with the pyrolysis model

    def __init__(self, runs, model, species, func, weightMass, weightRate):
        self.runs = runs 
        self.model = model
        self.func = func
        self.species = species 
        self.weightMass = weightMass 
        self.weightRate = weightRate 

    def input_func(self, parameter):
        """ The main function which returns the error, which serves as
            input for the optimiser and computes the errors per run for
            a given model, input parameter and precompt results

            Arguments:
            ---------
                    parameter: input parameter for the model e.g.:
                               pre-exp factor and tinit for const rate
                    func:
                    model:
                    runs:
                    species: name of the species to be fitted, needs
                             to be stored in runs 
        """
        # rename it and make a class function
        self.model.updateParameter(parameter)
        # collect errors of individual runs
        ret = [self.errorPerRun(run) for run in self.runs]

        # If we have a simple scalar list just sum the errors
        # else we component wise sum the error and return a vector
        # of errors per point
        self.error = (sum(ret) if type(ret[0]) != list else map(np.add, ret))
        return self.error

    def errorPerRun(self,run):
        """ Evaluate the the error per run compared to pre comp
            
            Computation of the error is based on given function func,
            since we either want a the global error or the error per point
            for least squares
         """
        times      = run['time(ms)']*1e-3
        targetMass = run[self.species]
        targetRate = run[self.species] #FIXME
        self.model.final_yield = targetMass[-1] #FIXME does this make sense?
        modeledMass = self.model.calcMass(
                init_mass = targetMass[0],
                time = times,
                temp = run.interpolate('temp'),
            )
        dt = False # FIXME
        modeledRate = self.model.computeTimeDerivative(modeledMass, times=times)
        # normalisation factor
        def norm(weight, target):
            return weight/np.power(Model.yieldDelta(target), 2.0)
        normMass = norm(self.weightMass, targetMass)
        normRate = norm(self.weightRate, targetRate)
        return self.func(targetRate, modeledRate,
                        targetMass, modeledMass,
                        normRate, normMass, dt)

    @classmethod
    def ls_input_func(cls, tr, mr, tm, mm, nr, nm, dt):
        ErrorRate = Model.modelErrorSquared(tr, mr)*dt
        ErrorMass = Model.modelErrorSquared(tm, mm)*dt
        return (ErrorMass * nm + ErrorRate * nr) * self.scaleFactor * dt

    @classmethod
    def min_input_func(cls, tr, mr, tm, mm, nr, nm, dt):
        ErrorMass = Model.totModelErrorSquaredPerc(tm, mm)
        ErrorRate = Model.totModelErrorSquaredPerc(tr, mr)
        return (ErrorMass*nm + ErrorRate*nr)/len(tm)

################childrenclasses####################

class constantRate(Model):
    """ The model calculating the mass
        with m(t)  =  m_s0+(m_s0-m_s,e)*exp(-k*(t-t_start))
             dm/dt = -k*(m-m_s,e).
        The Parameter to optimize are k and t_start.
    """
    #TODO store initial parameter to see test if parameters have been changed
    #TODO GO parameter[2] is not specified in inputs example

    def __init__(self, parameter):
        Model.__init__(self, "ConstantRate", parameter)
        self.k           = parameter["k"]
        self.start_time  = parameter["tstart"]
        self.final_yield = parameter.get('finalYield',False)
        # if set to false, the numerical time step corresponding to the outputed
        # by the detailled model (e.g CPD) is used; define a value to use instead this

    def __repr__(self):
        return  "Const Rate k {} tstart {}".format(self.k, self.start_time)


    def updateParameter(self, parameter):
        self.k          = parameter[0]
        self.start_time = parameter[1]
        return self

    def recalcMass(self):
        """ recalculate mass releas after updateParameter

            reuses init_mass, time and temp from previous
            computation, needed for genetic algorhythm
            since we only get the best parameters back
            and need to adjust the model
        """
        self.calcMass(self.init_mass, self.time)
        return self

    def calcMass(self, init_mass, time, temp=False):
        """ Computes the released mass over time

            Inputs:
                time: array of time values

        """
        # TODO GO Could it be beneficial to take k and
        #         and t_start as arguments?
        # we care only about time value
        # starting at start time
        self.init_mass = init_mass # store for recalc
        time = time - self.start_time
        # the yield still retained in the coal
        # this should converge to zero at large
        # time
        retained_mass = self.final_yield * np.exp(-self.k*time)
        released_mass = self.final_yield - retained_mass

        # if we are interested in the solid mass
        # TODO GO what is this solid thing going on here
        if False: #species == 'Solid':
            released_mass += solid_mass*np.exp(-self.k*time)

        self.mass = released_mass
        self.time = time
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
        return np.array([self.k, self.start_time])


class arrheniusRate(Model):
    """ The Arrhenius model in the standart notation:
            dm/dt=A*(T**b)*exp(-E/T)*(m_s-m)
        with the parameters a,b,E to optimize.
    """
    ODE_hmax = 1.e-2 #TODO what is this??
    absoluteTolerance = 1.0e-8
    relativeTolerance = 1.0e-6

    def __init__(self, parameter):
        Model.__init__(self, "ArrhenuisRate", parameter)
        self.A = parameter["preExp"]
        self.beta = parameter["beta"]
        self.E  = parameter["activationEnergy"]

    @property
    def parameter(self):
        return np.array([self.A, self.beta, self.E])

    def updateParameter(self, parameter):
        self.A    = parameter[0]
        self.beta = parameter[1]
        self.E    = parameter[2]
        ##print "Parameter update " + str(parameter)

    def calcMass(self, init_mass, time, temp=False):
        """Outputs the mass(t) using the model specific equation."""
        """ dm/dt=A*(T**b)*exp(-E/T)*(m_s-m)  """
        def dmdt(m, t):
            T  = temp(t) # so temp is a function that takes t and returns T
            dm = self.final_yield - m # finalYield
            if False:
                dmdt_ = (-self.A * dm  #FIXME this doesnt make sense!
                          * np.power(T, beta)
                          * np.exp(-self.E/T)
                            )
            else:
                dmdt_ =  (self.A * dm
                          * np.power(T, self.beta)
                          * np.exp(-self.E/T)) # TODO this should be Ta instead of E

            # sets values < 0 to 0.0, to avoid further problems
            return np.where(dmdt_ > 1e-64, dmdt_, 0.0)

        m_out = sp.integrate.odeint(
                func=dmdt,
                y0=[init_mass],
                t=time,
                atol=self.absoluteTolerance,
                rtol=self.relativeTolerance,
                hmax=self.ODE_hmax,
            )

        self.mass = m_out
        return m_out

    def ConvertKinFactors(self,ParameterVector):
        """ Dummy function actual has to convert the
            parameter to standard Arrhenius notation.
        """
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
