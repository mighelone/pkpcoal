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
        certain value  """

    def __init__(self, inp, target=100.00):
        """ From a given input dictionary and a target sum a
            dictionary with scaled composition is created """

        self.target = target
        self.basis = sum(inp.values())
        scaling_factor = target/self.basis
        self.elems = {key:value*scaling_factor for key,value in inp.iteritems()}

    def __getitem__(self,item):
        """ """
        try:
            return self.elems[item]
        except:
            print """Warning trying to access {} which was not set in
            input file, assuming 0.0.
            """.format(item)
            return 0.0

    def __contains__(self, key):
        return key in self.elems

    def __iter__(self):
        for name, mass in self.elems.iteritems():
            yield name, mass

    def remove_elem_mass_rebalance(self, elem, amount):
        """ select an element an subtract a percentage
            mass and rebalance the result """

        from copy import deepcopy
        input_dict = deepcopy(self.elems)
        input_dict[elem] = input_dict[elem] - amount
        return BalancedComposition(input_dict, self.target)

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

    def values(self):
        return self.elems.values()

    def iteritems(self):
        return self.elems.iteritems()

class Model(object):
    """ Parent class of the children ConstantRateModel,
        the three Arrhenius Models (notations) and the Kobayashi models.

        Parameter:
                name      = model name e.g constantRate
                parameter = initial parameter which are adapted during
                            optimisation
                species   = Name of the modeled species
                runs      = preProc results as list

        TimeVectorToInterplt allows the option to define the discrete time points,
        where to interpolate the results. If set to False (standard), then is are
        the outputed results equal the dt to solve the ODE. If set
        TimeVectorToInterplt=[t0,t1,t2,t3,t4] (t: floats) then is the
        yields result returned at method calcMass the yields at [t0,t1,t2,t3,t4],
        linear interploated. """

    def __init__(self, name, parameter, parameterBounds, inputs,
            species, calcMass, recalcMass, runs=False, constDt=False):
        print "Initialised {} Model".format(name)
        self.name = name
        self.initialParameter = parameter
        self.parameter = parameter
        self.parameterBounds = parameterBounds
        self.species = species
        self.constDt = constDt
        self.calcMass = calcMass
        self.finalYield = 1.0
        self.runs = (runs if runs else [])
        self.postGeneticOpt = True
        self.recalcMass = recalcMass
        self.inputs = inputs

    def printOptimalResults(self,opt_string,opt_value,opt_parameters):
        print('{} \nF(x) = {:e}'.format(opt_string,opt_value))
        print(''.join('{}\t{}\t{:f}\n'.format(i, self.paramNames[i], xi) for i,xi in enumerate(opt_parameters)))

    def fit(self, delta=0.05, finalOpt=True, **kwargs):
        # print 'initial parameter: ' + str(self.initialParameter)
        # do a rough estimation step first

        def check_limits(x,bounds):
            '''
            check if the optimized solution lies on the boundary limits
            :param x:
            :param bounds:
            :return:
            '''
            check = []
            for i,xi in enumerate(x):
                check.extend([xi <= bounds[i][0],xi >= bounds[i][1]])
            if any(check):
                print('WARNING: optimal solution lies in the bounds\nTry to increase the boundling limits!')
                print('var\tmin\topt\tmax')
                print(''.join('{:d}\t{:f}\t{:f}\t{:f}\n'.format(i,bounds[i][0],xi,bounds[i][1]) for i,xi in enumerate(x)))
            return None

        print 'preliminary optimisation: ' + self.species
        delta = delta
        #delta = 0.25
        optimization_mode = self.inputs['FitMethod'].split('+')
        if len(optimization_mode) == 2:
            if optimization_mode[0] == 'brute':
                print "Start brute force optimization..."
                from scipy.optimize import brute
                delta = self.inputs['bruteSettings']['delta']
                preOptimizedParameter = brute(
                    func=self.error_func,
                    ranges=self.parameterBounds,
                    finish=None,
                    Ns=int(1/delta))
                opt_value = self.error_func(preOptimizedParameter)
                self.printOptimalResults('Brute force optimization',opt_value,preOptimizedParameter)
                check_limits(preOptimizedParameter,self.parameterBounds)
                # TODO MV: to limit the minimization with the constrains from the brute search step it might be too limiting
                postOptBounds = [(optParam*(1.0-delta), optParam*(1.0+delta))
                            for optParam in preOptimizedParameter]
                self.parameter = preOptimizedParameter
            elif optimization_mode[0] == 'evolve':
                print('Start Evolve GA optimization...')
                preOptimizedParameter = self.geneticOpt()
                opt_value = self.error_func(preOptimizedParameter)
                self.printOptimalResults('GA evolve optimization',opt_value,preOptimizedParameter)
                #TODO add name of variable in the print
                check_limits(preOptimizedParameter,self.parameterBounds)
                postOptBounds = self.parameterBounds
        else:
            postOptBounds = self.parameterBounds
            preOptimizedParameter = self.parameter

        # minimize with fmin
        from scipy.optimize import minimize
        if not finalOpt:
            return self
        print 'Start fmin optimization...'
        optimizedParameter = minimize(
                fun  = self.error_func,
                x0   = preOptimizedParameter,
                bounds = postOptBounds,
                **kwargs
        )
        print('Minimize Optimization\nF(x) = {:e}'.format(optimizedParameter.fun))
        print(''.join('{:d}\t{:f}\n'.format(i,xi) for i,xi in enumerate(optimizedParameter.x)))
        self.printOptimalResults('Fmin optimization', optimizedParameter.fun, optimizedParameter.x)
        check_limits(optimizedParameter.x,postOptBounds)
        if not optimizedParameter.success:
            print ("WARNING final optimisation failed\nStatus: ",
                   optimizedParameter.status,
                   "using preliminary optimisation results")
            if optimizedParameter.fun < opt_value:
                self.parameter = optimizedParameter.x
            else:
                self.parameter = preOptimizedParameter
        else:
            self.parameter = optimizedParameter.x
        return self

    def fittedYield(self):
        # NOTE needs fit to be run before, probably
        # some checking
        optParams = self.parameter
        return self.recalcMass(optParams, time=self.runs[self.runs.keys()[0]]['time'])

    def updateParameter(self, parameter):
        self.parameter  = parameter
        return self

    def computeTimeDerivative(self, mass, deltaT=False, times=False):
        """ Return time derivatives for a given deltat array
            and a mass array """

        from numpy import gradient
        if deltaT:
            return gradient(mass, deltat)
        else:
            return False #FIXME gradient(np.array([mass,times]))

    def calcRate(self, preProcResult, time, temp, species):
        """ computes actual release reates for a given species
            by a time array and temperature array

            the preProcResults are used for initial values """
            #TODO GO can time be taken from preProcResult?
        
        mass = self.calcMass(preProcResult, time, temp, species)
        return self.computeTimeDerivative(mass, times = time)

    def geneticOpt(self):
        from pyevolve import G1DList, GSimpleGA, Selectors
        from pyevolve import Initializators, Mutators, Consts, DBAdapters
        # define genome
        genome = G1DList.G1DList(len(self.parameter))
        genome.setParams(rangemin=0.0, rangemax=1.0)
        genome.initializator.set(Initializators.G1DListInitializatorReal)
        genome.mutator.set(Mutators.G1DListMutatorRealRange)
        # The evaluator function (objective function)
        genome.evaluator.set(self.error_funcGA)
        # Genetic Algorithm Instance
        ga = GSimpleGA.GSimpleGA(genome)
        ga.setMinimax(Consts.minimaxType["minimize"])
        # set the population size
        ga.setPopulationSize(self.inputs['evolveSettings']['nPopulation']) #FIXME
        # set the number of generation
        ga.setGenerations(self.inputs['evolveSettings']['nGenerations']) #FIXME
        # Set the Roulette Wheel selector method,
        # the number of generations and the termination criteria
        ga.selector.set(Selectors.GRouletteWheel)
        ga.terminationCriteria.set(GSimpleGA.ConvergenceCriteria)
        sqlite_adapter = DBAdapters.DBSQLite(identify="GA", resetDB=True)
        ga.setDBAdapter(sqlite_adapter)
        # Do the evolution, with stats dump, frequency of 20 generations
        ga.evolve(freq_stats=1)
        # Gets the best individual
        best=ga.bestIndividual()
        #selects the bestiniviual
        return self.dimensionalParameters(best)

    def nonDimensionalParameters(self,parameters):
        '''
        convert dimensional parameters to non dimensional form
        fof GA.
        xnd = (x-min(x))/(min(x)-max(x))
        :param parameters: dimensional parameter
        :return: list non dimensional parameter
        '''
        return [(parameters[i]-self.parameterBounds[i][0])/
                (self.parameterBounds[i][1]-self.parameterBounds[i][1])
                for i,_ in enumerate(parameters)]

    def dimensionalParameters(self,nonDimensionalParameters):
        '''
        convert non-dimensional parameters to dimensional form
        after GA
        xnd = (x-min(x))/(min(x)-max(x))
        x = min(x) + (max(x) - min(x))*xnd
        :param parameters: dimensional parameter
        :return: list dimensional parameter
        '''
        return [self.parameterBounds[i][0] +
                (self.parameterBounds[i][1]-self.parameterBounds[i][0])*nonDimensionalParameters[i]
                for i,_ in enumerate(nonDimensionalParameters)]

    def error_funcGA(self, nonDimensionalParameters, func="cummulative", weightMass=1.0, weightRate=0.0):
        # converto to dimensional parameters
        parameters = self.dimensionalParameters(nonDimensionalParameters)
        return self.error_func(parameters)


    @classmethod
    def yieldDelta(cls, mass):
        return max(mass) - min(mass)

    def ErrorYield(self, preProcResult, Species):
        """ Returns the absolute deviation per point between the fitted
            pyrolysis model and the original preprocessor yield """
        # TODO GO: why only take some arrays from preProcResults
        pass

    def ErrorRate(self, preProcResult, species):
        """ Returns the absolute deviation per point between
            the fitted and the original rate curve. """
        pass

    def _mkInterpolatedRes(self,InputVecYields,Time):
        """ Generates the result vector. Outputs the result vector at the
            corresponding time steps corresponding to the imported time at
            method calcMass. Requiered for Pc Coal Lab (only few reported 
            points). """
        # t for for interplt, t_points, y_points
        return np.interp(Time,self.constDtVec,InputVecYields)

    def setDt4Intergrate(self,constantDt):
     """ constantDt allows the option to define numerical time step to solve the ODE.

        The outputted results ever equal the imported time list
        (when applying method calcMass Time = [t0,t1,t2,t3,t4]. If these time steps
        are too large, then is this defined dt used to solve the ODE and the results
        are linear interploated that way that they correspond to the imported time
        vector. To reset it, just set constantDt to False. """

     if constantDt != False:
         self.constDt = float(constantDt)

    def _mkDt4Integrate(self, Time):
        """ Time is the original time vector calculated by exact model, e.g. CPD.
            This class generates the internal dt vector if the dt defined by the
            user file is too large. A time step must be defined in Method
            setDt4Intergrate before. """

        if self.constDt != False:
         self.constDtVec = np.arange(Time[0], Time[-1], self.constDt)

        self.func = func

    def error_func(self, parameter, func="cummulative", weightMass=1.0, weightRate=0.0):
        """ Function for the optimizer, computes model error as a function
            of the input parameter

            Arguments:
            ---------
                    parameter: input parameter for the model e.g.:
                               pre-exp factor and tinit for const rate
                    func:      name of the function for evaluating the
                               error e.g. cummulative, or perpoint
        """
        # collect errors of individual runs
        func = Model.cumulative_error
        ret = [self.errorPerRun(parameter, run, func, weightMass, weightRate)
                    for run in self.runs.values()]

        # If we have a simple scalar list just sum the errors
        # else we component wise sum the error and return a vector
        # of errors per point
        self.error = (sum(ret) if type(ret[0]) != list else map(np.add, ret))
        #print parameter, self.error
        return self.error

    def errorPerRun(self, parameter, run, func, weightMass, weightRate):
        """ Evaluate the error per run compared to pre comp

            Computation of the error is based on given function func,
            since we either want a the global error or the error per point
            for least squares """

        times      = run['time']
        targetMass = run[self.species]
        targetRate = run[self.species] #FIXME
        modeledMass = self.calcMass(parameter, 0.0, times, run.interpolate('temp'))
        modeledRate = self.computeTimeDerivative(modeledMass, times=times)
        dt = False
        # normalisation factor
        def norm(weight, target):
            return weight/np.power(Model.yieldDelta(target), 2.0)
        normMass = norm(weightMass, targetMass)
        normRate = norm(weightRate, targetRate)
        return func(targetRate, modeledRate,
                    targetMass, modeledMass,
                    normRate, normMass, dt)

    @classmethod
    def perPointError(cls, tr, mr, tm, mm, nr, nm, dt):
        ErrorRate = Model.modelErrorSquared(tr, mr)*dt
        ErrorMass = Model.modelErrorSquared(tm, mm)*dt
        return (ErrorMass * nm + ErrorRate * nr) * self.scaleFactor * dt

    @classmethod
    def cumulative_error(cls, tr, mr, tm, mm, nr, nm, dt):
        """ Parameter tr: target rate
                      mr: model rate
                      tm: target mass
                      mm: model mass
                      nr: weight paramater rate
                      nm: weight parameter mass """

        ErrorMass = Model.totModelErrorAbsPerc(tm, mm)
        ErrorRate = Model.totModelErrorAbsPerc(tr, mr)
        return (ErrorMass*nm + ErrorRate*nr)

    @classmethod
    def modelErrori(cls, target, model):
        """ compute the deviation between modeled values and the target values
            from the pre processor """
        return abs(target - model)

    @classmethod
    def modelErrorSquared(cls, target, model):
        """ compute the deviation between modeled values and the target values
            from the pre processor """
        return np.power(target - model, 2.0)

    @classmethod
    def modelErrorAbs(cls, target, model):
        """ compute the deviation between modeled values and the target values
            from the pre processor """
        return np.abs(target - model)

    @classmethod
    def totModelErrorAbsPerc(cls, target, model):
        return np.sum(Model.modelErrorAbs(target, model))/len(target)

    @classmethod
    def summedModelError(cls, target, model):
        return np.sum(Model.modelError(target, model))

    @classmethod
    def totModelErrorSquaredPerc(cls, target, model):
        return np.sum(Model.modelErrorSquared(target, model))/len(target)

################ children classes ####################

class constantRate(Model):
    """ The model calculating the mass
        with m(t)  =  m_s0+(m_s0-m_s,e)*exp(-k*(t-t_start))
             dm/dt = -k*(m-m_s,e).
        The Parameter to optimize are k and t_start. """

    def __init__(self, inputs, runs, species):
        self.paramNames  = ['k', 'tstart', 'finalYield']
        parameter   = [inputs['constantRate'][paramName] for paramName in self.paramNames]
        paramBounds = [inputs['constantRate'].get(paramName+"Bounds",(None,None))
                         for paramName in self.paramNames]
        Model.__init__(self, "ConstantRate", parameter, paramBounds, inputs,
            species, self.calcMassConstRate, self.recalcMass, runs)

    # def __repr__(self):
    #     return  "Const Rate k {} tstart {}".format(self.k, self.start_time)

    def recalcMass(self, parameter, time):
        """ recalculate mass release after updateParameter

            reuses init_mass, time and temp from previous
            computation, needed for genetic algorithm
            since we only get the best parameters back
            and need to adjust the model """
        return self.calcMassConstRate(parameter, init_mass=None, time=time)

    def calcMassConstRate(self, parameter, init_mass, time, temp=False):
        """ Computes the released mass over time

            Inputs:
                parameter: array of parameter e.g [k, tstart, finalYield]
                init_mass: not implemented
                time: array of time values from the preproc

            NOTE: the function doesnt implicitly modifiy state of the
                  constantRate Model, to modify the released mass
                  use instance.mass = calcMassConstRate(...) """
        # TODO better use a classmethod?
        k = parameter[0]
        start_time = parameter[1]
        final_yield = parameter[2]

        # we care only about time value
        # starting at start time
        time_ = time - start_time

        # NOTE dont compute values for negative times
        # the yield should be zero anyways
        time_ = time_.clip(0)

        # the yield still retained in the coal
        # this should converge to zero at large
        # time
        retained_mass = final_yield * np.exp(-k*time_)
        released_mass = final_yield - retained_mass

        # if we are interested in the solid mass
        # TODO GO what is this solid thing going on here
        if False: #species == 'Solid':
            released_mass += solid_mass*np.exp(-k*time)

        # why choosing between released or solid mass
        # start_time is small then time
        # released_mass = np.where(time > self.start_time, released_mass, solid_mass)
        if self.constDt == False: # TODO GO shouldnt interpolation be used for var dt?
            #print "modeled_mass " + str(released_mass)
            return released_mass
        else: #returns the short, interpolated list (e.g. for PCCL)
            return self._mkInterpolatedRes(released_mass, time)

    @property
    def k(self):
        return self.parameter[0]

    @property
    def start_time(self):
        return self.parameter[1]

class arrheniusRate(Model):
    """ The Arrhenius model in the standart notation:
            dm/dt=A*(T**b)*exp(-E/T)*(m_s-m)
        with the parameters a,b,E to optimize. """

    def __init__(self, inputs, runs, species):
        self.paramNames  = ['preExp', 'activationEnergy']
        parameter   = [inputs['arrheniusRate'][paramName] for paramName in self.paramNames]
        paramBounds = [inputs['arrheniusRate'].get(paramName+"Bounds",(None,None))
                         for paramName in self.paramNames]
        Model.__init__(self, "ArrheniusRate", parameter, paramBounds, inputs,
            species, self.calcMassArrhenius, self.recalcMassArrhenius, runs=runs)
        self.updateParameter(self.parameter)
        # FIXME this assumes that the final yield is run independent
        sel_run = runs.keys()[0] # FIXME
        self.final_yield = runs[sel_run][species][-1] # FIXME
        self.lowerT = inputs['arrheniusRate'].get('lowerDevolTemp', False)
        self.temp = runs[sel_run].interpolate('temp')

    def updateParameter(self, parameter):
        self.A = parameter[0]
        self.E = parameter[1]
        return self
        ##print "Parameter update " + str(parameter)

    def recalcMassArrhenius(self, parameter, time):
        return self.calcMassArrhenius(parameter, init_mass=0.0, time=time, temp=self.temp)

    def calcMassArrhenius(self, parameter, init_mass, time, temp):
        """ Outputs the mass(t) using the model specific equation.
            dm/dt=A*(T**0)*exp(-E/T)*(m_s-m) """
        # TODO replace Ta by E

        A, E = parameter[0], parameter[1]
        def dmdt(m, t):
            T = temp(t) # so temp is a function that takes t and returns T

            if self.lowerT and T < self.lowerT:
               return 0.0
            dm = self.final_yield - m # finalYield
            exp = np.exp(-E/T)
            if exp == np.inf:
                print """Warning overflow in the exponential
                      term detected, change parameter bounds
                      of the activation Temperature """
                return 0.0
            if False:
                dmdt_ = (-A * dm  #FIXME this doesnt make sense!
                          * np.power(T, beta)
                          * np.exp(-E/T)
                            )
            else:
                dmdt_ = (init_mass + A*dm*exp)
            # print "dmdt, t, m, dm, T ",  dmdt_, t, m, dm, T
            # sets values < 0 to 0.0, to avoid further problems
            return float(dmdt_) #np.where(dmdt_ > 1e-64, dmdt_, 0.0)
        m_out = sp.integrate.odeint(func=dmdt, y0=0.0, t=time)
        m_out = m_out[:, 0]
        if self.constDt == False: # TODO GO shouldnt interpolation be used for var dt?
            #print "modeled_mass " + str(released_mass)
            return m_out
        else: #returns the short, interpolated list (e.g. for PCCL)
            return self._mkInterpolatedRes(m_out, time)


class C2SM(Model):
    """
    The C2SM Kobayashi model
    """
    def __init__(self,inputs, runs, species):
        self.paramNames = ['alpha1','A1','E1','alpha2','A2','E2']
        parameter = [inputs['C2SM'][paramName] for paramName in self.paramNames]
        paramBounds = [inputs['C2SM'].get(paramName+"Bounds",(None,None))
                         for paramName in self.paramNames]
        Model.__init__(self, "C2SM", parameter, paramBounds, inputs,
            species, self.calcMassC2SM, self.recalcMassC2SM, runs=runs)
        self.updateParameter(parameter)
        sel_run = runs.keys()[0] # FIXME
        #self.final_yield = runs[sel_run][species][-1] # FIXME
        #self.temp = runs[sel_run].interpolate('temp')


    def updateParameter(self, parameter):
        self.alpha1 = parameter[0]
        self.A1 = parameter[1]
        self.E1 = parameter[2]
        self.alpha2 = parameter[3]
        self.A2 = parameter[4]
        self.E2 = parameter[5]
        return self

    def recalcMassC2SM(self, parameter, time):
        return self.calcMassC2SM(parameter, init_mass=0.0, time=time, temp=self.temp)

    #def calcMassC2SM(self,parameter, init_mass, time, temp):
    #    return 0.0

    def calcMassC2SM(self,parameter, init_mass, time, temp):
         self.updateParameter(parameter=parameter)
         def dmdt(m, t):
             '''
             the function calculate the following yields
             ds/dt = - (k1 + k2) s ! raw coal fraction consumption
             dy/dt = (alpha1 * k1 + alpha2 * k2) * s ! volatile products

             y + s + c = 1
             :param m: = [s, y]
             :param t: time
             :return:
             '''
             s = m[0]
             y = m[1]
             T = temp(t) # so temp is a function that takes t and returns T
             exp1 = np.exp(-self.E1/T)
             exp2 = np.exp(-self.E2/T)
             if exp1 == np.inf or exp2 == np.inf:
                 print """Warning overflow in the exponential
                       term detected, change parameter bounds
                       of the activation Temperature """
                 return 0.0
             k1 = self.A1*exp1
             k2 = self.A2*exp2
             dsdt = -(k1+k2)*s
             dydt = (self.alpha1*k1 + self.alpha2*k2)*s
             return [dsdt, dydt]
    #
         m_out = sp.integrate.odeint(func=dmdt, y0=[1.0,0.0], t=time)
         m_out = m_out[:, 1] # return only the volatile yield y
         if self.constDt == False: # TODO GO shouldnt interpolation be used for var dt?
             #print "modeled_mass " + str(released_mass)
             return m_out
         else: #returns the short, interpolated list (e.g. for PCCL)
             return self._mkInterpolatedRes(m_out, time)

# class Kobayashi(Model):
#     """ Calculates the devolatilization reaction using the Kobayashi model.
#         The Arrhenius equation inside are in the standard notation. """
#
#     def __init__(self,InitialParameterVector):
#         print 'Kobayashi Model initialized'
#         self._modelName = 'Kobayashi'
#         self._ParamVector=InitialParameterVector
#         self.ODE_hmax=1.e-2
#         self.constDt = False # if set to false, the numerical time step corresponding to the outputted by the dtailled model (e.g CPD) is used; define a value to use instead this
#
#     def calcMass(self,preProcResult,time,T,Name):
#         """ Outputs the mass(t) using the model specific equation.
#             The input Vector is [A1, E1, A2, E2, alpha1, alpha2] """
#         # question whether the dt from DetailledModel result file or
#         # from a constant dt should be used
#         if self.constDt == False: # dt for integrate = dt from DM result file
#             timeInt = time
#         else: #if dt in DM results file has too large dt
#             self._mkDt4Integrate(time)
#             timeInt = self.constDtVec
#         self.preProcResult=preProcResult
#         timeInt=np.delete(timeInt,0)
#         self.__Integral=0.0
#         tList=[0.0]
#         k1k2=[0.0]
#         ParamVec=self.ParamVector()
#         #
#         def dmdt(m,t):
#             k1=ParamVec[0]*np.exp(-ParamVec[1]/(T(t)))
#             k2=ParamVec[2]*np.exp(-ParamVec[3]/(T(t)))
#             tList.append(t)
#             k1k2.append(k1+k2)
#             self.__Integral+=0.5*(tList[-1]-tList[-2])*(k1k2[-1]+k1k2[-2])
#             dmdt_out = ( (ParamVec[4]*k1+ParamVec[5]*k2)*np.exp(-self.__Integral) )
#             dmdt_out=np.where(abs(dmdt_out)>1.e-300,dmdt_out,0.0) #sets values<0 =0.0, otherwise it will further cause problem(nan)
#             return dmdt_out
#         InitialCondition=[0]
#         m_out=sp.integrate.odeint(dmdt,InitialCondition,timeInt,atol=1.e-5,rtol=1.e-4,hmax=1.e-2)
#         m_out=m_out[:,0]
#         m_out=np.insert(m_out,0,0.0)
#         if self.constDt == False:
#             if (ParamVec[0]<0 or ParamVec[1]<0 or ParamVec[2]<0 or
#                 ParamVec[3]<0 or ParamVec[4]<0 or ParamVec[5]>1):
#                 m_out[:]=float('inf')
#                 return m_out
#             else:
#                 return m_out
#         else: #returns the short, interpolated list (e.g. for PCCL)
#             return self._mkInterpolatedRes(m_out,time)
#
#
# class KobayashiPCCL(Model):
#     """ Calculates the devolatilization reaction using the Kobayashi model.
#         The Arrhenius equation inside are in the standard notation. The fitting
#         parameter are as in PCCL A1, A2, E1, alpha1. TimeVectorToInterplt
#         allows the option to define the discrete time points, where to
#         interpolate the results. If set to False (standard), then is are the
#         outputted results equal the dt to solve the ODE. """
#
#     def __init__(self,InitialParameterVector):
#         print 'Kobayashi Model initialized'
#         self._modelName = 'KobayashiPCCL'
#         self._ParamVector=InitialParameterVector
#         self.ODE_hmax=1.e-2
#         self.constDt = False # if set to false, the numerical time step corresponding to the outputted by the dtailled model (e.g CPD) is used; define a value to use instead this
#
#     def calcMass(self,preProcResult,time,T,Name):
#         """ Outputs the mass(t) using the model specific equation. """
#         # question whether the dt from DetailledModel result file or from a constant dt should be used
#         if self.constDt == False: # dt for integrate = dt from DM result file
#             timeInt = time
#         else: #if dt in DM results file has too large dt
#             self._mkDt4Integrate(time)
#             timeInt = self.constDtVec
#         self.preProcResult=preProcResult
#         timeInt=np.delete(timeInt,0)
#         self.__Integral=0.0
#         tList=[0.0]
#         k1k2=[0.0]
#         ParamVec=self.ParamVector()
#         #
#         def dmdt(m,t):
#             k1=ParamVec[0]*np.exp(-ParamVec[2]/(T(t)))
#             k2=ParamVec[1]*np.exp(-(ParamVec[2]+self.__E2diff)/(T(t)))
#             tList.append(t)
#             k1k2.append(k1+k2)
#             self.__Integral+=0.5*(tList[-1]-tList[-2])*(k1k2[-1]+k1k2[-2])
#             dmdt_out = ( (ParamVec[3]*k1+self.__alpha2*k2)*np.exp(-self.__Integral) )
#             dmdt_out=np.where(abs(dmdt_out)>1.e-300,dmdt_out,0.0) #sets values<0 =0.0, otherwise it will further cause problem(nan)
#             return dmdt_out
#         InitialCondition=[0]
#         m_out=sp.integrate.odeint(dmdt,InitialCondition,timeInt,atol=1.e-5,rtol=1.e-4,hmax=1.e-2)
#         m_out=m_out[:,0]
#         m_out=np.insert(m_out,0,0.0)
#         if self.constDt == False:
#             if (ParamVec[0]<0 or ParamVec[1]<0 or ParamVec[2]<0 or ParamVec[3]<0):
#                 m_out[:]=float('inf')
#                 return m_out
#             else:
#                 return m_out
#         else: #returns the short, interpolated list (e.g. for PCCL)
#             return self._mkInterpolatedRes(m_out,time)
#
#
#     def ConvertKinFactors(self,ParameterVector):
#         """ Outputs the Arrhenius equation factors in the shape
#             [A1, E1, A2, E2]. Here where the real Arrhenius model
#             is in use only a dummy function. """
#
#         P=self.ParamVector()
#         return [P[0],P[1],P[2],P[3]]
#
#     def setKobWeights(self,alpha2):
#         """ Sets the two Kobayashi weights alpha2. """
#         self.__alpha2=alpha2
#
#     def KobWeights(self):
#         """ Returns the two Kobayashi weights alpha2. """
#         return self.__alpha2
#
#     def setE2Diff(self,DifferenceE1E2):
#         """ Sets the dE in E2=E1+dE. """
#         self.__E2diff=DifferenceE1E2
#
#     def E2Diff(self):
#         """Returns the dE in E2=E1+dE."""
#         return self.__E2diff
#
# class KobayashiA2(Model):
#     """ Calculates the devolatilization reaction using the Kobayashi model.
#         The Arrhenius equation inside are in the secend alternative notation
#         (see class ArrheniusModelAlternativeNotation2). """
#
#     def __init__(self,InitialParameterVector):
#         print 'Kobayashi Model initialized'
#         self._ParamVector=InitialParameterVector
#         self.ODE_hmax=1.e-2
#         self.constDt = False # if set to false, the numerical time step corresponding to the outputted by the dtailled model (e.g CPD) is used; define a value to use instead this
#
#     def calcMass(self,preProcResult,time,T,Name):
#         """ Outputs the mass(t) using the model specific equation. """
#         # question whether the dt from DetailledModel result file or from a constant dt should be used
#         if self.constDt == False: # dt for integrate = dt from DM result file
#             timeInt = time
#         else: #if dt in DM results file has too large dt
#             self._mkDt4Integrate(time)
#             timeInt = self.constDtVec
#         self.preProcResult=preProcResult
#         T_general=preProcResult.Rate('Temp')
#         self.T_min=min(T_general)
#         self.T_max=max(T_general)
#         self.c=1./(1./self.T_max-1./self.T_min)
#         #alpha has to be greater equal zero:
#         absoluteTolerance = 1.0e-8
#         relativeTolerance = 1.0e-6
#         self.tList=[0.0]
#         #deletes to solve trapezian rule:
#         timeInt=np.delete(timeInt,0)
#         self.Integral=0.0
#         self.k1k2=[0.0]
#         ParamVec=self.ParamVector()
#         u=preProcResult.Yield(Name)
#         #
#         def dmdt(m,t):
#             k1=np.exp( self.c*( ParamVec[0]*(1./T(t)-1./self.T_min) - ParamVec[1]*(1./T(t)-1./self.T_max) ) )
#             k2=np.exp( self.c*( ParamVec[2]*(1./T(t)-1./self.T_min) - ParamVec[3]*(1./T(t)-1./self.T_max) ) )
#             self.tList.append(t)
#             self.k1k2.append(k1+k2)
#             self.Integral+=0.5*(self.tList[-1]-self.tList[-2])*(self.k1k2[-1]+self.k1k2[-2])
#             dmdt_out = ( (self.__alpha1*k1+self.__alpha2*k2)*np.exp(-self.Integral) )
#             dmdt_out=np.where(abs(dmdt_out)>1.e-300,dmdt_out,0.0) #sets values<0 =0.0, otherwise it will further cause problem(nan)
#             return dmdt_out
#         InitialCondition=[u[0]]
#         m_out=sp.integrate.odeint(dmdt,InitialCondition,timeInt,atol=absoluteTolerance,rtol=relativeTolerance,hmax=self.ODE_hmax)
#         m_out=m_out[:,0]
#         m_out=np.insert(m_out,0,0.0)
#         if self.constDt == False:
#             return m_out
#         else: #returns the short, interpolated list (e.g. for PCCL)
#             return self._mkInterpolatedRes(m_out,time)
#
#     def ConvertKinFactors(self,ParameterVector):
#         """ Converts the alternative notaion Arrhenius factors into the
#             standard Arrhenius factors and return them in the
#             shape [A1,E1], [A2,E2] """
#
#         A1=np.exp( -self.c*ParameterVector[0]/self.T_min + self.c*ParameterVector[1]/self.T_max )
#         E1=self.c*ParameterVector[1]-self.c*ParameterVector[0]
#         A2=np.exp( -self.c*ParameterVector[2]/self.T_min + self.c*ParameterVector[3]/self.T_max )
#         E2=self.c*ParameterVector[3]-self.c*ParameterVector[2]
#         return [A1,E1,A2,E2]
#
#     def setKobWeights(self,alpha1,alpha2):
#         """ Sets the two Kobayashi weights alpha1 and alpha2. """
#         self.__alpha1=alpha1
#         self.__alpha2=alpha2
#
#     def KobWeights(self):
#         """ Returns the two Kobayashi weights alpha1 and alpha2. """
#         return self.__alpha1, self.__alpha2
#
#
# class DAEM(Model):
#     """ Calculates the devolatilization reaction using the
#         Distributed Activation Energy Model. """
#     def __init__(self,InitialParameterVector):
#         print 'DAEM initialized'
#         self._modelName = 'DAEM'
#         self._ParamVector=InitialParameterVector
#         self.ODE_hmax=1.e-2
#         self.NrOfActivationEnergies=50
#         self.constDt = False # if set to false, the numerical time step corresponding to the outputted by the dtailled model (e.g CPD) is used; define a value to use instead this
#
#     def setNrOfActivationEnergies(self,NrOfE):
#         """ Define for how many activation energies of the range
#             of the whole distribution the integral shall be solved
#             (using Simpson Rule)."""
#         self.NrOfActivationEnergies=NrOfE
#
#     def NrOfActivationEnergies(self):
#         """ Returns the number of activation enrgies the integral shall
#             be solved for (using Simpson Rule). """
#         return self.NrOfActivationEnergies
#
#     def calcMass(self,preProcResult,time,T,Name):
#         """ Outputs the mass(t) using the model specific equation. """
#         self.E_List=np.arange(int(self._ParamVector[1]-3.*self._ParamVector[2]),int(self._ParamVector[1]+3.*self._ParamVector[2]),int((6.*self._ParamVector[2])/self.NrOfActivationEnergies)) #integration range E0 +- 3sigma, see [Cai 2008]
#         # question whether the dt from DetailledModel result file or from a constant dt should be used
#         if self.constDt == False: # dt for integrate = dt from DM result file
#             timeInt = time
#         else: #if dt in DM results file has too large dt
#             self._mkDt4Integrate(time)
#             timeInt = self.constDtVec
#         #Inner Integral Funktion
#         def II_dt(t,E_i):
#             return np.exp( -E_i/T(t) )
#         #outer Integral for one activation energy from t0 to tfinal
#         #stores all values of the inner Integrals (time,ActivationEnergy) in a 2D-Array
#         InnerInts=np.zeros([len(timeInt),len(self.E_List)])
#         CurrentInnerInt=np.zeros(len(timeInt))
#         for Ei in range(len(self.E_List)):
#             CurrentInnerInt[:]=II_dt(timeInt[:],self.E_List[Ei])
#             InnerInts[1:,Ei] = sp.integrate.cumtrapz(CurrentInnerInt,timeInt[:])
#         #
#         def OI_dE(EIndex,tIndex):
#             m = np.exp(-self._ParamVector[0]*InnerInts[tIndex,EIndex])*(1./(self._ParamVector[2]*(2.*np.pi)**0.5))*np.exp(-(self.E_List[EIndex]-self._ParamVector[1])**2/(2.*self._ParamVector[2]**2))
# #            print 'InnerInt',InnerInt,'mass',dm_dt
#             return m
#         m_out=np.zeros(np.shape(timeInt))
#         mE=np.zeros(np.shape(self.E_List))
#         for ti in range(len(timeInt)):
#             for Ei in range(len(self.E_List)):
#                 mE[Ei]=OI_dE(Ei,ti)
#             m_out[ti]=sp.integrate.simps(mE,self.E_List)
#         #descaling
#         m_out = self._ParamVector[3]*(1.-m_out)
#         if self.constDt == False:
#             return m_out
#         else: #returns the short, interpolated list (e.g. for PCCL)
#             return self._mkInterpolatedRes(m_out,time)
