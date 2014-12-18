import numpy as np
from scipy.optimize import fmin
from scipy.optimize import fmin_cg
from scipy.optimize import fmin_bfgs
from scipy.optimize import fmin_ncg
from scipy.optimize import leastsq
from scipy.optimize import fmin_slsqp

#not really precise, just for tests (Analytical solution)
def OptGradBased(inputs, model, results, species):
    """ Starts a gradient Based Optimization and returns the final Fit.

    Parameter:
        inputs: dictionary containing the optimisation parameter
        model:  the pyrolsis model for which rate and yield defining
                parameters are optimised
        finalYield: the final yield of the considered species
        species: the name of the species to optimise

    Note:
        For Kobayashi Model set Final Yield to False (independent), for all
        other set a value. It will be excluded from the optimization.
    """
    # Called by PyrolModelLauncher for every species
    ls = LeastSquaresEstimator(inputs['Optimisation'])
    result = ls.estimate(results, model, species)
    print 'Final error= ' +  str(ls.deviation)
    return result

def OptGenAlgBased(inputs, model, results, species):
    """ Starts a genetic algorithm and afterwards a gradient Based optimization.
        Sets the Final Fit result as the ParamVector in the Kinetic Model.
        Input are the Fit (Result Objects of the Detailed Models),
        the Parameter to initialize, the two Parameter vectors defining
        the range of the results and the Species index.
    """
    genAlg = Evolve.GenericOpt(inputs)
    model.update(genAlg.estimate())
    # afterwards grad based optimization
    if GlobalOptParam.optimizGrad == True:
        self.OptGradBased(Fit, ParameterVecInit, False, Species)

class TwoPointEstimator(object):
    """Solves the devolatilization reaction analytically using two arbitrary selected points and the constant rate model. Unprecise. Should only be used for tests."""
    def __init__(self):
        print '2 Point initialized'

    def estimate(self,fgdvc,Species,PointLocation):
    #first part: calculates a list of VM_s and the precursor values k and VM(0) for VM as one species##########################
        t=fgdvc.Time()
        TimePoint=int(PointLocation*(fgdvc.NPoints()))
        if Species == 'Solid' or Species == self.Yields2Cols('Solid'):
            VM_s=fgdvc.MassVM_s()
            VM_0=VM_s[0]
            k_VM=(np.log(VM_s[TimePoint]/VM_0))/(-t[TimePoint])
        else:
            u=fgdvc.Yield(Species)
            u_0=u[-1]
            k_VM=(np.log(1- u[TimePoint]/u_0))/(-t[TimePoint])
        return k_VM
#
#the optimizer to use:
class LeastSquaresEstimator(object):
    """ Optimizes the Fitting curve using the Least Squares
        for Yields and the Rates.
    """

    Pre_Tolerance = 1.e-5 # Tolerance used for the self.improve_? functions
    PreMaxIter = 50       # Number of maximum iterations used for
                          # the self.improve_? functions

    def __init__(self, optimisation_parameter, finalYield=False):
        """
        Parameters:
            optimizer: one of 'fmin', 'fmin_cg', 'fmin_bfgs','fmin_ncg',
                              'fmin_slsqp', 'leastsq'
                       Experience is that 'fmin' and 'leastsq' gives best results

            fitTolerance:
            weights: weights for yields and rates for the fitting procedure
            finalYield: final yield for the ODE to optimize. False if model
                        is Kobayashi equation
                        (independend of final yield,standard Setting).
                        Must be applied for all models except the Kobayashi model"
        """
        print 'Least Square initialized'
        # Select one optimizer of the scipy.optimizer library:
        #        'fmin','fmin_cg','fmin_bfgs','fmin_ncg','fmin_slsqp' or 'leastsq'.
        # According to experience 'fmin' (or also 'leastsq') generates at best the results."""
        opt = optimisation_parameter
        self.optimizer    = opt['GradBasedOpt'] #FIXME:GO this will change for runs > 1
        self.maxIter      = opt['maxIter']
        self.scaleFactor  = opt['scaleFactor']
        self.fitTolerance = opt['Tolerance'],
        self.weightMass   = opt['weightMass']
        self.weightRate   = opt['weightRate']
        self.FinalY       = finalYield


    # def improve_E(self,fgdvc,model,t,T,Parameter_Vector,Name,):
    #     """Additional option: Only the Activation Energy in the Arrhenius Equation is optimized. Actual not necessary."""
    #     def E_adjust(ActivationEnergie):
    #         model.setParamVector([Parameter_Vector[0],Parameter_Vector[1],ActivationEnergie])
    #         print model.ParamVector()
    #         uDot=fgdvc.Rate(Name)
    #         vDot=model.calcRate(fgdvc,t,T,Name)
    #         return sum((uDot-vDot)**2)
    #     IC_ActEner=Parameter_Vector[2] #InitialCondition for Optimzation
    #     #Optimizing
    #     Optimized_E1=fmin(E_adjust,IC_ActEner,ftol=self.Pre_Tolerance,maxiter=self.PreMaxIter)  #shifts the rates curves by the Activation Energy together
    #     #reformate, because Vec has shape: [0.1, 5.0, array([ 26468.75])]
    #     Optimized_AE=Optimized_E1[0]
    #     return Optimized_AE
    #
    #
    # def improve_a(self,fgdvc,model,t,T,Parameter_Vector,Name):
    #     """Additional option: Only the preexponential factor in the Arrhenius Equation is optimized. Actual not necessary."""
    #     IC_a=Parameter_Vector
    #     def a_adjust(aPre):
    #         model.setParamVector([aPre,Parameter_Vector[1],Parameter_Vector[2]])
    #         print model.ParamVector()
    #         uDot=fgdvc.Rate(Name)
    #         vDot=model.calcRate(fgdvc,t,T,Name)
    #         return sum((uDot-vDot)**2)
    #     Optimized_a1=fmin(a_adjust,IC_a[0],ftol=self.Pre_Tolerance,maxiter=self.PreMaxIter)
    #     Optimized_a=Optimized_a1[0]
    #     return Optimized_a
    #
    # def maxLengthOfVectors(self,results):
    #     """Returns the minimum lenght of a all vectors from the several runs."""
    #     Len_tPointsL=[]
    #     for i in range(len(results)):
    #         Len_tPointsL.append(len(results[i].Time()))
    #     Len_tPoints=max(Len_tPointsL)
    #     return Len_tPoints

    def estimate(self, results, model, species):
        """ The main optimization method.

            Optimizes the pyrolysis parameters (i.e. t_start and k for constant rate)
            to match the preprocessor results.

            Inputs:
                results = a list of preprocessor result objects
                model = the pyrolysis model to optimise
                species = selected species to optimize

            How it works
            1. based on the selected optimsation procedure the input function
               is constructed
            2. for each run in the results the errors are evaluated and reduced
            3. when the optimiser converged the pyrolysis model with new parameters
               is returned
        """
        from PKP.src.Models import Model
        import scipy.optimize as scopt

        # NOTE here a block of functions is defined to evaluate the
        # the errors this should be moved to a separate class and
        # needs test methods
        def input_func(parameter, func, model, runs, species):
            """ The main function which returns the error, which serves as
                input for the optimiser

                this functions takes care on reducing the errors evaluated
                per run. If we get a list of scalars we just sum the error
                if get a list of lists containing errors per point we need
                to merge it
            """
            # rename it and make a class function
            model.updateParameter(parameter)
            # collect errors of individual runs
            ret = [errorPerRun(run, func, model, species) for run in runs]
            # if we have a simple scalar list just sum the errors
            # else we component wise sum the error
            return (sum(ret) if type(ret[0]) != list else map(np.add,ret))

        def errorPerRun(run, func, model, species):
            """ adapts final yield and  computes the modeled yield per run """
            times = run['time(ms)']*1e-3
            targetMass  = run[species]
            targetRate  = run[species] #FIXME
            model.final_yield = targetMass[-1] #FIXME does this make sense?
            modeledMass = model.calcMass(
                    init_mass = targetMass[0],
                    time = times,
                    temp = run['temp'],
                )
            dt = False # FIXME
            modeledRate = model.computeTimeDerivative(modeledMass, times = times)
            # normalisation factor
            def norm(weight, target):
                return weight/np.power(Model.yieldDelta(target), 2.0)
            normMass = norm(self.weightMass, targetMass)
            normRate = norm(self.weightRate, targetRate)
            return func(targetRate, modeledRate,
                        targetMass, modeledMass,
                        normRate, normMass, dt)

        def ls_input_func(tr, mr, tm, mm, nr, nm, dt):
            ErrorRate = Model.modelErrorSquared(tr, mr)*dt
            ErrorMass = Model.modelErrorSquared(tm, mm)*dt
            return (ErrorMass * nm + ErrorRate * nr) * self.scaleFactor * dt

        def min_input_func(tr, mr, tm, mm, nr, nm, dt):
            ErrorMass = Model.totModelErrorSquaredPerc(tm, mm)
            ErrorRate = Model.totModelErrorSquaredPerc(tr, mr)
            return (ErrorMass*nm + ErrorRate*nr)/len(tm)

        print 'start gradient based optimization, species: ' + species
        optimiser = getattr(scopt, self.optimizer)
        # select optimisation input function depending on the optimizer
        # this is needed since leastsq expects a list of errors where fmin
        # simply expects a global error
        error_func = (ls_input_func if self.optimizer == 'leastsq' else min_input_func)
        OptimizedVector = optimiser(
                func = input_func,
                x0   = model.parameter,
                args = (error_func, model, results, species),
                ftol = self.fitTolerance, # TODO GO what is the diff between gtol&ftol
                maxiter = self.maxIter
        ) # caLculates optimized vector
        # NOTE It seems unneccessary to reevaluate the model agian,
        #      but for now no better solution is in sight,
        #      so we will use it to store final yields and rates on the model
        self.deviation = input_func(OptimizedVector, error_func, model, results, species)
        return model


class GlobalOptimizer(object):
    """Makes runs over a defined range to look for global optimum. Local Optimizer is an LeastSquarsEstimator object, KineticModel is e.g. an constantRate model object, Fit_one_runObj is the List containing the Objects supporting the local fitting procedure with data."""
    def __init__(self,localOptimizer,KineticModel,Fit_one_runObj):
        self.LocOpt=localOptimizer
        self.KinModel=KineticModel
        self.FitInfo=Fit_one_runObj
        self.AllParameterList=[] #collects the final results of the local minima.

    def setParamList(self,ParameterList):
        """Sets the Parameter list."""
        self.__ParamList=ParameterList

    def ParamList(self):
        """Returns the Parameter list."""
        return self.__ParamList

    def GenerateOptima(self,Species,IndexListofParameterToOptimize,ArrayOfRanges,ListNrRuns):
        """This method makes a several number of runs and returns the Parameter having the lowest deviation of all local minima. The list contains 3 or 4 parameters. If e.g., the second and the third Parameter have to be optimized: IndexListofParameterToOptimize=[1,2]. If e.g. the range of the second is 10000 to 12000 and for the third 1 to 5,  ArrayOfRanges=[[10000,12000],[1,5]]. When ListNrRuns is e.g. [0,10,5,0] the the second Parameter will be optimizted eleven times between 10000 and 12000, the third six times between 1 and 5. Attention, the number of runs grows by NrRuns1*NrRuns2*NrRuns3*NrRuns4, which can lead to a very large number needing very much time!"""
        DevList=[]       #saves the deviation"""
        DevList=[]       #saves the deviation
        ParamArray=[]    #saves the corresponding Prameter lists
        Parameter=self.KinModel.ParamVector()
        for i in range(len(IndexListofParameterToOptimize)):
            Parameter[IndexListofParameterToOptimize[i]]=ArrayOfRanges[i][0]
        self.setParamList(Parameter)
        if len(self.__ParamList)==4:
            for i in range(ListNrRuns[0]+1):
                if ListNrRuns[0]!=0:
                    self.__ParamList[0] = ArrayOfRanges[IndexListofParameterToOptimize.index(0)][0] + (ArrayOfRanges[IndexListofParameterToOptimize.index(0)][1]-ArrayOfRanges[IndexListofParameterToOptimize.index(0)][0])*((float(i))/ListNrRuns[0])
                for j in range(ListNrRuns[1]+1):
                    if ListNrRuns[1]!=0:
                        self.__ParamList[1]=ArrayOfRanges[IndexListofParameterToOptimize.index(1)][0] + (ArrayOfRanges[IndexListofParameterToOptimize.index(1)][1]-ArrayOfRanges[IndexListofParameterToOptimize.index(1)][0])*((float(j))/ListNrRuns[1])
                    for k in range(ListNrRuns[2]+1):
                        if ListNrRuns[2]!=0:
                            self.__ParamList[2]=ArrayOfRanges[IndexListofParameterToOptimize.index(2)][0] + (ArrayOfRanges[IndexListofParameterToOptimize.index(2)][1]-ArrayOfRanges[IndexListofParameterToOptimize.index(2)][0])*((float(k))/ListNrRuns[2])
                        for l in range(ListNrRuns[3]+1):
                            if ListNrRuns[3]!=0:
                                self.__ParamList[3] = ArrayOfRanges[IndexListofParameterToOptimize.index(3)][0] + (ArrayOfRanges[IndexListofParameterToOptimize.index(3)][1]-ArrayOfRanges[IndexListofParameterToOptimize.index(3)][0])*((float(l))/ListNrRuns[3])
                            #
                            self.KinModel.setParamVector(self.LocOpt.estimate_T(self.FitInfo,self.KinModel,self.__ParamList,Species))
                            DevList.append(self.LocOpt.Deviation())
                            ParamArray.append(self.KinModel.ParamVector())
        ####
        if len(self.__ParamList)==3:
            for i in range(ListNrRuns[0]+1):
                if ListNrRuns[0]!=0:
                    self.__ParamList[0] = ArrayOfRanges[IndexListofParameterToOptimize.index(0)][0] + (ArrayOfRanges[IndexListofParameterToOptimize.index(0)][1]-ArrayOfRanges[IndexListofParameterToOptimize.index(0)][0])*((float(i))/ListNrRuns[0])
                for j in range(ListNrRuns[1]+1):
                    if ListNrRuns[1]!=0:
                            self.__ParamList[1]=ArrayOfRanges[IndexListofParameterToOptimize.index(1)][0] + (ArrayOfRanges[IndexListofParameterToOptimize.index(1)][1]-ArrayOfRanges[IndexListofParameterToOptimize.index(1)][0])*((float(j))/ListNrRuns[1])
                    for k in range(ListNrRuns[2]+1):
                        if ListNrRuns[2]!=0:
                            self.__ParamList[2]=ArrayOfRanges[IndexListofParameterToOptimize.index(2)][0] + (ArrayOfRanges[IndexListofParameterToOptimize.index(2)][1]-ArrayOfRanges[IndexListofParameterToOptimize.index(2)][0])*((float(k))/ListNrRuns[2])
                            #
                        self.KinModel.setParamVector(self.LocOpt.estimate_T(self.FitInfo,self.KinModel,self.__ParamList,Species))
                        DevList.append(self.LocOpt.Deviation())
                        ParamArray.append(self.KinModel.ParamVector())
        if len(self.__ParamList)==6:
            for i in range(ListNrRuns[0]+1):
                if ListNrRuns[0]!=0:
                    self.__ParamList[0] = ArrayOfRanges[IndexListofParameterToOptimize.index(0)][0] + (ArrayOfRanges[IndexListofParameterToOptimize.index(0)][1]-ArrayOfRanges[IndexListofParameterToOptimize.index(0)][0])*((float(i))/ListNrRuns[0])
                for j in range(ListNrRuns[1]+1):
                    if ListNrRuns[1]!=0:
                        self.__ParamList[1]=ArrayOfRanges[IndexListofParameterToOptimize.index(1)][0] + (ArrayOfRanges[IndexListofParameterToOptimize.index(1)][1]-ArrayOfRanges[IndexListofParameterToOptimize.index(1)][0])*((float(j))/ListNrRuns[1])
                    for k in range(ListNrRuns[2]+1):
                        if ListNrRuns[2]!=0:
                            self.__ParamList[2]=ArrayOfRanges[IndexListofParameterToOptimize.index(2)][0] + (ArrayOfRanges[IndexListofParameterToOptimize.index(2)][1]-ArrayOfRanges[IndexListofParameterToOptimize.index(2)][0])*((float(k))/ListNrRuns[2])
                        for l in range(ListNrRuns[3]+1):
                            if ListNrRuns[3]!=0:
                                self.__ParamList[3] = ArrayOfRanges[IndexListofParameterToOptimize.index(3)][0] + (ArrayOfRanges[IndexListofParameterToOptimize.index(3)][1]-ArrayOfRanges[IndexListofParameterToOptimize.index(3)][0])*((float(l))/ListNrRuns[3])
                            for m in range(ListNrRuns[4]+1):
                                if ListNrRuns[4]!=0:
                                    self.__ParamList[4]=ArrayOfRanges[IndexListofParameterToOptimize.index(4)][0] + (ArrayOfRanges[IndexListofParameterToOptimize.index(4)][1]-ArrayOfRanges[IndexListofParameterToOptimize.index(4)][0])*((float(m))/ListNrRuns[4])
                                for n in range(ListNrRuns[5]+1):
                                    if ListNrRuns[5]!=0:
                                        self.__ParamList[5] = ArrayOfRanges[IndexListofParameterToOptimize.index(5)][0] + (ArrayOfRanges[IndexListofParameterToOptimize.index(5)][1]-ArrayOfRanges[IndexListofParameterToOptimize.index(5)][0])*((float(n))/ListNrRuns[5])
                            #
                            self.KinModel.setParamVector(self.LocOpt.estimate_T(self.FitInfo,self.KinModel,self.__ParamList,Species))
                            DevList.append(self.LocOpt.Deviation())
                            ParamArray.append(self.KinModel.ParamVector())
        indexMinDeviation=DevList.index(np.min(DevList))
        return ParamArray[indexMinDeviation]





