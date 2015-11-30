"""
A selection of different Pyrolysis models

A model pyrolysis model usually provides method to calculate the yield of
individual species the rates and energy balancing methods

"""
import numpy as np
import scipy as sp
import scipy.integrate
import matplotlib.pyplot as plt


class BalancedComposition(object):
    """ Class for compostion that ensures componenents sum up to a
        certain value  """

    def __init__(self, inp, target=100.00):
        """ From a given input dictionary and a target sum a
            dictionary with scaled composition is created """

        self.target = target
        self.basis = sum(inp.values())
        scaling_factor = target / self.basis
        self.elems = {key: value * scaling_factor for key, value in inp.iteritems()}


    def __getitem__(self, item):
        """ """
        try:
            return self.elems[item]
        except:
            raise TypeError('%r item is not defined' % item)


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
        return BalancedComposition({key: elem for key, elem in self.elems.iteritems()
                                    if key not in elems_})

    def scale(self, factor):
        return BalancedComposition(self.elems, target=self.target * factor)

    def __repr__(self):
        return str(self.elems)

    def values(self):
        return self.elems.values()

    def iteritems(self):
        return self.elems.iteritems()

class ParticleHeatUp(object):

    def __init__(self, cp, rho, d, sigma, epsilon, T_wall, Nu=2.0):
        self.cp = cp
        self.rho = rho
        self.d = d
        self.sigma = sigma
        self.epsilon = epsilon
        self.Nu = 2.0
        self.T_wall = T_wall

    def computeParticleTemperature(self, time, t_gas,
            start_temp=300, dt=1e-5, t_end=0.1):
        def updateTg(T_old, T_gas, dt):
            return T_old + dTdt(T_gas, T_old)*dt

        def dTdt(gas_temp, particle_temp):
            T4diff = (particle_temp**4.0 - self.T_wall**4.0)*self.sigma*self.epsilon
            h = self.Nu*self.lambda_g(gas_temp)/self.d
            Tdiff = (particle_temp - gas_temp)
            a = -3.0/(self.cp*self.rho*self.d/2.0)
            return  a*(T4diff + h*Tdiff)

        Temps, Times = [start_temp], [0]
        for i in range(int(t_end/dt)):
            time_ = dt*(i+1) # current time
            T = np.interp(time_, time, t_gas) # interpolate temperature
            start_temp = updateTg(start_temp, T, dt) # update temperature
            Temps.append(start_temp)
            Times.append(time_)
        return Times, Temps

    def particle_response_time(self, cp, rho, d, T_gas, lambda_g):
        """ compute particle response time based on
            Formula 1.93b Page 44 from PARTICLES, DROPS AND BUBBLES:
            FLUID DYNAMICS AND NUMERICAL METHODS, E. Loth """

        return (cp*rho*d**2)/(12.0*lambda_g(T_gas))

    def rho_g(self, temp):
        """ p/rho  = R_s*T"""
        return 1.0 /(8.314*temp)

    def lambda_g(self, temp):
        """ compute the conductivity of dry air"""
        # http://www.nist.gov/data/PDFfiles/jpcrd283.pdf
        gamma = 25.98e-3 # [W/mK]
        c1,c2,c3,c4,c5,c6,c7 = [0.2395, 0.0064, 1.0, -1.9261,
                                2.0038, -1.0755, 0.2294]
        d1,d2,d3,d4,d5 = [0.4022, 0.3566, -0.1631, 0.1380,
                          -0.0201]
        T = temp/132.5
        r = self.rho_g(temp)/314.3
        return (gamma * (c1*T + c2*T**0.5 + c3
                        + c4*T**(-1.0) + c5*T**(-2.0)
                        + c6*T**(-3.0) + c7*T**(-4.0)
                        + d1*r + d2*r**2.0 + d3*r**3.0
                        + d4*r**4.0 + d5*r**4.0))

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
        self._error = None

    def printOptimalResults(self, opt_string, opt_value, opt_parameters):
        print('{} \nF(x) = {:e}'.format(opt_string, opt_value))
        print(''.join('{}\t{}\t{:f}\n'.format(i, self.paramNames[i], xi) for i, xi in enumerate(opt_parameters)))

    def fit(self, delta=0.05, finalOpt=True, **kwargs):
        # print 'initial parameter: ' + str(self.initialParameter)
        # do a rough estimation step first

        def check_limits(x, bounds):
            """ check if the optimized solution lies on the boundary limits
            :param x:
            :param bounds:
            :return: """
            check = []
            for i, xi in enumerate(x):
                check.extend([xi <= bounds[i][0], xi >= bounds[i][1]])
            if any(check):
                print('WARNING: optimal solution lies in the bounds\nTry to increase the boundling limits!')
                print('var\tmin\topt\tmax')
                print(
                ''.join('{:d}\t{:f}\t{:f}\t{:f}\n'.format(i, bounds[i][0], xi, bounds[i][1]) for i, xi in enumerate(x)))
            return None
        print 'preliminary optimisation: ' + self.species
        delta = delta
        # delta = 0.25
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
                    Ns=int(1 / delta))
                opt_value = self.error_func(preOptimizedParameter)
                self.printOptimalResults('Brute force optimization', opt_value, preOptimizedParameter)
                check_limits(preOptimizedParameter, self.parameterBounds)
                # TODO MV: to limit the minimization with the constrains from
                # the brute search step it might be too limiting
                postOptBounds = [(optParam * (1.0 - delta), optParam * (1.0 + delta))
                                 for optParam in preOptimizedParameter]
                self.parameter = preOptimizedParameter
            elif optimization_mode[0] == 'evolve':
                print('Start Evolve GA optimization...')
                preOptimizedParameter = self.geneticOpt()
                opt_value = self.error_func(preOptimizedParameter)
                self.printOptimalResults('GA evolve optimization', opt_value, preOptimizedParameter)
                # TODO add name of variable in the print
                check_limits(preOptimizedParameter, self.parameterBounds)
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
            fun=self.error_func,
            x0=preOptimizedParameter,
            bounds=postOptBounds,
            **kwargs
        )
        self.printOptimalResults('Fmin optimization', optimizedParameter.fun, optimizedParameter.x)
        check_limits(optimizedParameter.x, postOptBounds)
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

    def plot_yield(self, axis="time"):
        f, ax = plt.subplots(1, 2)
        f.set_figwidth(15)
        cols = plt.rcParams['axes.color_cycle']
        for col, (runNr, run) in zip(cols, self.runs.iteritems()):
            modelYield = self.fittedYield(run)
            modelRate = self.fittedRate(run)
            targetYield = run[self.species]
            ax[0].plot(run[axis], targetYield, ls='--', label=run.solver, color=col)
            ax[0].plot(run[axis], modelYield, label=self.name, color=col )
            ax[1].plot(run[axis], run.rate(self.species), ls='--', label=run.solver, color=col)
            ax[1].plot(run[axis], modelRate, label=self.name, color=col)
        xlabel = 'Time [s]' if axis=="time" else 'Temp [K]'
        ylabels = ['Yield [kg/kg_daf]', 'Yield Rate [kg/kg_daf/s]']
        for a, ylabel in zip(ax, ylabels):
            a.get_yaxis().get_label().set_text(ylabel)
            a.get_xaxis().get_label().set_text(xlabel)
        plt.legend()
        return f, ax

    # def fit(self, **kwargs):
    #     # print 'initial parameter: ' + str(self.initialParameter)
    #     # do a rough estimation step first
    #     from scipy.optimize import brute
    #     import time
    #     print 'preliminary optimisation: ' + self.species
    #     nThreads = kwargs.get('splits', 1)
    #     delta = kwargs.get('delta', 0.05)*nThreads # delta per Thread
    #     bounds = self.parameterBounds
    #     def preOpt(bounds):
    #         print "bounds ", bounds , " delta ", delta
    #         return brute(
    #             func=self.error_func,
    #             ranges=bounds,
    #             finish=None,
    #             full_output=True,
    #             Ns=int(1/delta))
    #     lower = lambda l, u: np.linspace(l, u-((u-l)/float(nThreads)), nThreads)
    #     upper = lambda l, u: np.linspace(l+((u-l)/float(nThreads)), u, nThreads)
    #     lowerupper = lambda l1, u1: (lower(l1, u1).tolist(), upper(l1, u1).tolist())
    #     for _ in range(kwargs.get('preOptRefinements', 1)):
    #         boundsP = [lowerupper(l, u) for l, u in bounds]
    #         boundsP = list(product(*boundsP)) if nThreads > 1 else bounds
    #         start_time = time.time()
    #         if nThreads > 1:
    #             optParams = parmap(preOpt, boundsP)
    #             smallest = np.inf
    #             choosenParams = None
    #             for i, opt in enumerate(optParams):
    #                 if opt[1] < smallest:
    #                     smallest =  opt[1]
    #                     choosenParams = opt[0]
    #         else:
    #             optParams = preOpt(boundsP)
    #             choosenParams, smallest  = optParams[0:2]
    #         print "--- run {} params {} error {} {}s seconds ---".format(
    #                 _, str(choosenParams), smallest, (time.time() - start_time))
    #         n = kwargs.get('narrowing', 0.75)
    #         bounds = [(param*(1.0-delta*n**_), param*(1.0+delta*n**_)) for param in choosenParams]
    #         self.parameter = choosenParams
    #         if kwargs.get('plot'):
    #             self.plot_yield()
    #     return self

    def fittedYield(self, run):
        # NOTE needs fit to be run before, probably
        # some checking
        optParams = self.parameter
        return self.recalcMass(optParams, run=run)

    def fittedRate(self, run):
        """ returns the release rate of the species """
        import numpy as np
        species = self.fittedYield(run)
        time = run["time"]
        dt = np.diff(time)
        dt = np.append(time[0], dt)
        return np.gradient(species, dt)

    def updateParameter(self, parameter):
        self.parameter = parameter
        return self

    def computeTimeDerivative(self, mass, deltaT=False, times=False):
        """ Return time derivatives for a given deltat array
            and a mass array """
        from numpy import gradient
        if deltaT:
            return gradient(mass, deltat)
        else:
            return False  # FIXME gradient(np.array([mass,times]))

    def calcRate(self, preProcResult, time, temp, species):
        """ computes actual release reates for a given species
            by a time array and temperature array

            the preProcResults are used for initial values """
        mass = self.calcMass(preProcResult, time, temp, species)
        return self.computeTimeDerivative(mass, times=time)

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
        ga.setPopulationSize(self.inputs['evolveSettings']['nPopulation'])  # FIXME
        # set the number of generation
        ga.setGenerations(self.inputs['evolveSettings']['nGenerations'])  # FIXME
        # Set the Roulette Wheel selector method,
        # the number of generations and the termination criteria
        ga.selector.set(Selectors.GRouletteWheel)
        ga.terminationCriteria.set(GSimpleGA.ConvergenceCriteria)
        sqlite_adapter = DBAdapters.DBSQLite(identify="GA", resetDB=True)
        ga.setDBAdapter(sqlite_adapter)
        # Do the evolution, with stats dump, frequency of 20 generations
        ga.evolve(freq_stats=1)
        # Gets the best individual
        best = ga.bestIndividual()
        # selects the bestiniviual
        return self.dimensionalParameters(best)

    def nonDimensionalParameters(self, parameters):
        """
        convert dimensional parameters to non dimensional form
        fof GA.
        xnd = (x-min(x))/(min(x)-max(x))
        :param parameters: dimensional parameter
        :return: list non dimensional parameter
        """
        return [(parameters[i] - self.parameterBounds[i][0]) /
                (self.parameterBounds[i][1] - self.parameterBounds[i][1])
                for i, _ in enumerate(parameters)]

    def dimensionalParameters(self, nonDimensionalParameters):
        """
        convert non-dimensional parameters to dimensional form
        after GA
        xnd = (x-min(x))/(min(x)-max(x))
        x = min(x) + (max(x) - min(x))*xnd
        :param nonDimensionalParameters: non dimensional parameter used in GA
        :return: list dimensional parameter
        """
        return [self.parameterBounds[i][0] +
                (self.parameterBounds[i][1] - self.parameterBounds[i][0]) * nonDimensionalParameters[i]
                for i, _ in enumerate(nonDimensionalParameters)]

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

    def _mkInterpolatedRes(self, InputVecYields, Time):
        """ Generates the result vector. Outputs the result vector at the
            corresponding time steps corresponding to the imported time at
            method calcMass. Requiered for Pc Coal Lab (only few reported
            points). """
        # t for for interplt, t_points, y_points
        return np.interp(Time, self.constDtVec, InputVecYields)

    def setDt4Intergrate(self, constantDt):
        """ constantDt allows the option to define numerical time step to solve the ODE.

        The outputted results ever equal the imported time list
        (when applying method calcMass Time = [t0,t1,t2,t3,t4]. If these time steps
        are too large, then is this defined dt used to solve the ODE and the results
        are linear interploated that way that they correspond to the imported time
        vector. To reset it, just set constantDt to False. """

        if constantDt is not False:
            self.constDt = float(constantDt)

    def _mkDt4Integrate(self, Time):
        """ Time is the original time vector calculated by exact model, e.g. CPD.
            This class generates the internal dt vector if the dt defined by the
            user file is too large. A time step must be defined in Method
            setDt4Intergrate before. """

        if self.constDt is not False:
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
        # We simply sum up the individual error per run
        self.error = sum(ret)
        return self.error

    def errorPerRun(self, parameter, run, func, weightMass, weightRate):
        """ Evaluate the error per run compared to pre comp

            Computation of the error is based on given function func,
            since we either want a the global error or the error per point
            for least squares """

        times = run['time']
        targetMass = run[self.species]
        targetRate = run[self.species]  # FIXME
        modeledMass = self.calcMass(parameter, 0.0, times, run.interpolate('temp'))
        modeledRate = self.computeTimeDerivative(modeledMass, times=times)
        dt = False
        # normalisation factor
        def norm(weight, target):
            return weight / np.power(Model.yieldDelta(target), 2.0)
        normMass = norm(weightMass, targetMass)
        normRate = norm(weightRate, targetRate)
        return func(targetRate, modeledRate,
                    targetMass, modeledMass,
                    normRate, normMass, dt)

    @classmethod
    def perPointError(cls, tr, mr, tm, mm, nr, nm, dt):
        ErrorRate = Model.modelErrorSquared(tr, mr) * dt
        ErrorMass = Model.modelErrorSquared(tm, mm) * dt
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
        return ErrorMass * nm + ErrorRate * nr

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
        return np.sum(Model.modelErrorAbs(target, model)) / len(target)

    @classmethod
    def summedModelError(cls, target, model):
        return np.sum(Model.modelError(target, model))

    @classmethod
    def totModelErrorSquaredPerc(cls, target, model):
        return np.sum(Model.modelErrorSquared(target, model)) / len(target)

################ children classes ####################

class constantRate(Model):
    """ The model calculating the mass
        with m(t)  =  m_s0+(m_s0-m_s,e)*exp(-k*(t-t_start))
             dm/dt = -k*(m-m_s,e).
        The Parameter to optimize are k and t_start. """

    def __init__(self, inputs, runs, species):
        self.paramNames = ['k', 'tstart', 'finalYield']
        parameter = [inputs['constantRate'][paramName] for paramName in self.paramNames]
        paramBounds = [inputs['constantRate'].get(paramName + "Bounds", (None, None))
                       for paramName in self.paramNames]
        Model.__init__(self, "ConstantRate", parameter, paramBounds, inputs,
                       species, self.calcMassConstRate, self.recalcMass, runs)

    # def __repr__(self):
    #     return  "Const Rate k {} tstart {}".format(self.k, self.start_time)

    def recalcMass(self, parameter, time, temp):
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
        retained_mass = final_yield * np.exp(-k * time_)
        released_mass = final_yield - retained_mass

        # if we are interested in the solid mass
        # TODO GO what is this solid thing going on here
        if False:  # species == 'Solid':
            released_mass += solid_mass * np.exp(-k * time)

        # why choosing between released or solid mass
        # start_time is small then time
        # released_mass = np.where(time > self.start_time, released_mass, solid_mass)
        if self.constDt is False:  # TODO GO shouldnt interpolation be used for var dt?
            # print "modeled_mass " + str(released_mass)
            return released_mass
        else:  # returns the short, interpolated list (e.g. for PCCL)
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
        self.paramNames = ['preExp', 'activationEnergy']
        parameter = [inputs['arrheniusRate'][paramName] for paramName in self.paramNames]
        paramBounds = [inputs['arrheniusRate'].get(paramName + "Bounds", (None, None))
                       for paramName in self.paramNames]
        Model.__init__(self, "ArrheniusRate", parameter, paramBounds, inputs,
                       species, self.calcMassArrhenius, self.recalcMassArrhenius, runs=runs)

        self.updateParameter(self.parameter)
        # FIXME this assumes that the final yield is run independent
        sel_run = runs.keys()[0]  # FIXME
        self.final_yield = runs[sel_run][species][-1]  # FIXME
        self.lowerT = inputs['arrheniusRate'].get('lowerDevolTemp', False)

    def updateParameter(self, parameter):
        self.A = parameter[0]
        self.E = parameter[1]
        return self


    def recalcMassArrhenius(self, parameter, time):
        print "!! Warning dont use recalcMass, use recalcMassPerRun"
        return self.calcMassArrhenius(parameter, init_mass=0.0, time=time, temp=self.temp)

    def recalcMassPerRunArrhenius(self, parameter, run):
        return self.calcMassArrhenius(parameter, init_mass=0.0, time=run['time'], temp=run.interpolate('temp'))

    def calcMassArrhenius(self, parameter, init_mass, time, temp):
        """ Outputs the mass(t) using the model specific equation.
            dm/dt=A*(T**0)*exp(-E/T)*(m_s-m) """
        # TODO replace Ta by E

        A, E = parameter[0], parameter[1]

        def dmdt(m, t):
            T = temp(t)  # so temp is a function that takes t and returns T

            if self.lowerT and T < self.lowerT:
                return 0.0
            dm = self.final_yield - m  # finalYield
            exp = np.exp(-E / T)
            if exp == np.inf:
                print """Warning overflow in the exponential
                      term detected, change parameter bounds
                      of the activation Temperature """
                return 0.0
            if False:
                dmdt_ = (-A * dm  # FIXME this doesnt make sense!
                         * np.power(T, beta)
                         * np.exp(-E / T)
                         )
            else:
                dmdt_ = (init_mass + A * dm * exp)
            # print "dmdt, t, m, dm, T ",  dmdt_, t, m, dm, T
            # sets values < 0 to 0.0, to avoid further problems
            return float(dmdt_)  # np.where(dmdt_ > 1e-64, dmdt_, 0.0)

        m_out = sp.integrate.odeint(func=dmdt, y0=0.0, t=time)
        m_out = m_out[:, 0]
        if self.constDt is False:  # TODO GO shouldnt interpolation be used for var dt?
            # print "modeled_mass " + str(released_mass)
            return m_out
        else:  # returns the short, interpolated list (e.g. for PCCL)
            return self._mkInterpolatedRes(m_out, time)


class C2SM(Model):
    """
    The C2SM Kobayashi model
    note that the optimization uses the log values of A1 and A2
    """

    def __init__(self, inputs, runs, species):
        self.paramNames = ['alpha1', 'A1', 'E1', 'alpha2', 'A2', 'E2']
        parameter = [inputs['C2SM'][paramName] for paramName in self.paramNames]
        paramBounds = [inputs['C2SM'].get(paramName + "Bounds", (None, None))
                       for paramName in self.paramNames]
        # modify parameters 1 and 4 to log
        for i in [1,4]:
            parameter[i] = np.log(parameter[i])
            paramBounds[i] = [np.log(par) for par in paramBounds[i]]
            self.paramNames[i] = 'log'+self.paramNames[i]
        Model.__init__(self, "C2SM", parameter, paramBounds, inputs,
                       species, self.calcMassC2SM, self.recalcMassC2SM, runs=runs)
        sel_run = runs.keys()[0]  # FIXME
        # self.final_yield = runs[sel_run][species][-1] # FIXME
        # self.temp = runs[sel_run].interpolate('temp')

    def recalcMassC2SM(self, parameter, time):
        return self.calcMassC2SM(parameter, init_mass=0.0, time=time, temp=self.temp)

    # def calcMassC2SM(self,parameter, init_mass, time, temp):
    #    return 0.0

    def calcMassC2SM(self, parameter, init_mass, time, temp):

        def dmdt(m, t):
            """
             the function calculate the following yields
             ds/dt = - (k1 + k2) s ! raw coal fraction consumption
             dy/dt = (alpha1 * k1 + alpha2 * k2) * s ! volatile products

             y + s + c = 1
             :param m: = [s, y]
             :param t: time
             :return:
             """
            s = m[0]
            y = m[1]
            alpha1 = parameter[0]
            A1 = np.exp(parameter[1])
            E1 = parameter[2]
            alpha2 = parameter[3]
            A2 = np.exp(parameter[4])
            E2 = parameter[5]
            T = temp(t)  # so temp is a function that takes t and returns T
            exp1 = np.exp(-E1 / T)
            exp2 = np.exp(-E2 / T)
            if exp1 == np.inf or exp2 == np.inf:
                print """Warning overflow in the exponential
                       term detected, change parameter bounds
                       of the activation Temperature """
                return 0.0
            k1 = A1 * exp1
            k2 = A2 * exp2
            dsdt = -(k1 + k2) * s
            dydt = (alpha1 * k1 + alpha2 * k2) * s
            return [dsdt, dydt]
            #

        m_out = sp.integrate.odeint(func=dmdt, y0=[1.0, 0.0], t=time)
        m_out = m_out[:, 1]  # return only the volatile yield y
        if self.constDt is False:  # TODO GO shouldnt interpolation be used for var dt?
            # print "modeled_mass " + str(released_mass)
            return m_out
        else:  # returns the short, interpolated list (e.g. for PCCL)
            return self._mkInterpolatedRes(m_out, time)