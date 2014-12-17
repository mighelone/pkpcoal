import numpy as np
from scipy.optimize import fmin
from scipy.optimize import fmin_cg
from scipy.optimize import fmin_bfgs
from scipy.optimize import fmin_ncg
from scipy.optimize import leastsq
from scipy.optimize import fmin_slsqp

#not really precise, just for tests (Analytical solution)
def OptGradBased(inputs, model, results, finalYield, species):
    """ Starts a gradient Based Optimization and returns the final Fit.

    Parameter:
        inputs: dictionary containing the optimisation parameter
        model:  the pyrolsis model for which rate and yield defining
                parameters are optimised
        finalYield: the final yield of the considered species
        species: the name of the species to optimise

    Note:
        For Kobayashi Model set Final Yield to False (independent), for all 
        other set a value.  It will be excluded from the optimization. 
    """
    # Called by PyrolModelLauncher for every species
    ls = LeastSquaresEstimator(inputs['Optimisation'], finalYield)
    result = ls.estimate(results, model, species)
    print 'Final error= ' +  str(ls.deviation)
    return result

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

    def estimate(self, results, model, species, preLoopNumber=0):
        """ The main optimization method. 
            Optimizes the Fitting curve using the Least Squares for the weighted Yields 
            and the weighted Rates considering the temperature history.

            Requires at input: 
                The corresponding Fit_one_run object, the Model object, the kinetic parameter list, 
                a name (e.g. the species). 

            preLoopNumber is the number of running the improve_E and improve_a routines. 
            So the standard setting of preLoopNumber is equal zero. 
            It may be used if there is only a very bad convergence when optimize all three parameter.
            #TODO GO why estimateT?
        """
        def LeastSquaresFunction(parameter, model, run, species):
            """ The function which is to be optimised.
                
            """
            from PKP.src.Models import Model
            # TODO GO is it executed only for NRruns == 1
            # rename it and make a class function
            model.updateParameter(parameter)
            times = run['time(ms)']*1e-3
            modeled_mass = model.calcMass(
                    init_mass = run[species][0],
                    time = times,
                    temp = run['temp'], 
                )
            target_mass  = run[species] 
            target_rate  = run[species]
            modeled_rate = model.computeTimeDerivative(modeled_mass, times = times)
            massError = self.weightMass/np.power(Model.yieldDelta(target_mass), 2.0)
            if False: #self.selectedOptimizer == 'leastsq':
                w0 = self.scaleFactor * massError
                w1 = self.scaleFactor * self.weightRate/np.power(max(target_rate), 2.0)
                # Error = np.zeros(maxLen, dtype='d')
                # Dot2_ = w1*(
                #     ((uDot[runnedCaseNr]-vDot[runnedCaseNr])**2)
                #         *dt[runnedCaseNr]
                #     )
                # #the yield term
                # Dot1_ = w0*(
                #     ((u[runnedCaseNr]-v[runnedCaseNr])**2)
                #         *dt[runnedCaseNr]
                #     )
                # #makes an array, containing both, the rates and yields
                # Error[:len(Dot1_)] += Dot1_ + Dot2_
                # # print "deviation: ",np.sum(Error)
            else: # if self.inputs[] = 'fmin' ... 
                # NOTE: Removed loop since this part is only called if runs == 1 
                # TODO GO where does this come from?
                # TODO GO double check if its the rate or mass? propably its mass
                ErrorMass = (Model.totModelErrorSquaredPerc(target_mass, modeled_mass)
                            * massError)

                ErrorRate = (Model.totModelErrorSquaredPerc(target_rate, modeled_rate)
                            * self.weightRate/np.power(Model.yieldDelta(target_rate), 2.0))

                Error = (ErrorMass + ErrorRate)/len(target_mass)
                # print modeled_mass
                # print "ErrorMass " + str(ErrorMass)
            return Error

        # TODO GO Double check
        print 'start gradient based optimization, species: ' + species

        # add final yield to parameters
        # is this part of the optimized parameters too?
        # should it be passed as args to least squares
        #if model.final_yield == False: # TODO test if final_yield is needed
        model.final_yield = results[species][-1]
        #     # TODO where does the FinalY come from? 
        #     Parameter = list(Parameter)
        #     Parameter.append(self.FinalY)
        #     Parameter = np.array(Parameter)

        import scipy.optimize as scopt
        optimiser = getattr(scopt, self.optimizer)
        OptimizedVector = optimiser(
                func = LeastSquaresFunction,
                x0   = model.parameter,
                args = (model, results, species ), 
                ftol = self.fitTolerance, # TODO GO what is the diff between gtol&ftol
                maxiter = self.maxIter    
        ) # caLculates optimized vector
        self.deviation = LeastSquaresFunction(OptimizedVector, model, results, species)
        print self.deviation
        print OptimizedVector
        # appends final yield
        if self.FinalY != False:
            OptimizedVector = list(OptimizedVector)
            OptimizedVector.append(self.FinalY)
            OptimizedVector = np.array(OptimizedVector)
        return OptimizedVector


    def setPreTolerance(self, ToleranceForFminFunction):
        """ Sets the tolerance as a abort criterion for the prefitting procedure 
            (if preLoopNumber in estimate_T is not equal zero). """
        self.Pre_Tolerance=ToleranceForFminFunction

    def setPreMaxIter(self,MaxiumNumberOfIterationInPreProcedure):
        """ Sets the maximum number of iteration oin the optimizer as a abort 
            criterion for the prefitting procedure (if preLoopNumber in estimate_T 
            is not equal zero). """
        self.PreMaxIter=int(MaxiumNumberOfIterationInPreProcedure)


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





