import numpy as np
from scipy.optimize import fmin
from scipy.optimize import fmin_cg
from scipy.optimize import fmin_bfgs
from scipy.optimize import fmin_ncg
from scipy.optimize import leastsq
from scipy.optimize import fmin_slsqp
#
#not really precise, just for tests (Analytical solution)
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
class LeastSquarsEstimator(object):
    """Optimizes the Fitting curve using the Least Squares for the Yields and the Rates."""
    def __init__(self):
        print 'Least Square initialized'
        self.Pre_Tolerance=1.e-5   #Tolerance used for the self.improve_? functions
        self.Fit_Tolerance=1.e-10  #Tolerance for main run
        self.PreMaxIter=None       #Number of maximum iterations used for the self.improve_? functions
        self.MaxIter=None          #Number of maximum iterations for main run

    def improve_E(self,fgdvc,model,t,T,Parameter_Vector,Name,):
        """Additional option: Only the Activation Energy in the Arrhenius Equation is optimized. Actual not necessary."""
        def E_adjust(ActivationEnergie):
            model.setParamVector([Parameter_Vector[0],Parameter_Vector[1],ActivationEnergie])
            print model.ParamVector()
            uDot=fgdvc.Rate(Name)
            vDot=model.calcRate(fgdvc,t,T,Name)
            return sum((uDot-vDot)**2)
        IC_ActEner=Parameter_Vector[2] #InitialCondition for Optimzation
        #Optimizing
        Optimized_E1=fmin(E_adjust,IC_ActEner,ftol=self.Pre_Tolerance,maxiter=self.PreMaxIter)  #shifts the rates curves by the Activation Energy together
        #reformate, because Vec has shape: [0.1, 5.0, array([ 26468.75])]
        Optimized_AE=Optimized_E1[0]
        return Optimized_AE
       

    def improve_a(self,fgdvc,model,t,T,Parameter_Vector,Name):
        """Additional option: Only the preexponential factor in the Arrhenius Equation is optimized. Actual not necessary."""
        IC_a=Parameter_Vector
        def a_adjust(aPre):
            model.setParamVector([aPre,Parameter_Vector[1],Parameter_Vector[2]])
            print model.ParamVector()
            uDot=fgdvc.Rate(Name)
            vDot=model.calcRate(fgdvc,t,T,Name)
            return sum((uDot-vDot)**2)
        Optimized_a1=fmin(a_adjust,IC_a[0],ftol=self.Pre_Tolerance,maxiter=self.PreMaxIter)
        Optimized_a=Optimized_a1[0]
        return Optimized_a
    
    def minLengthOfVectors(self,fgdvc_list):
        """Returns the minimum lenght of a all vectors from the several runs."""
        Len_tPointsL=[]
        for i in range(len(fgdvc_list)):
            Len_tPointsL.append(len(fgdvc_list[i].Time()))
        Len_tPoints=min(Len_tPointsL)
        return Len_tPoints
        
        

    def estimate_T(self,fgdvc_list,model,Parameter_Vector,Name,preLoopNumber=0):
        """The main optimization method. Optimizes the Fitting curve using the Least Squares for the weighted Yields and the weighted Rates considering the temperatur history. Requires at input: The corresponding Fit_one_run object, the Model object, the kinetic parameter list, a name (e.g. the species). preLoopNumber is the number of running the  improve_E and improve_a routines. So the standard setting of preLoopNumber is equal zero. It may be used if there is only a very bad convergence when optimize all three parameter."""
        Len_tPoints=self.minLengthOfVectors(fgdvc_list)
        t=np.zeros([(Len_tPoints),len(fgdvc_list)]) #line index, time, column index: runned case
        dt=np.zeros([(Len_tPoints),len(fgdvc_list)]) #line index, time, column index: runned case
        T=list()
        u=np.zeros([(Len_tPoints),len(fgdvc_list)]) #line index, time, column index: runned case
        uDot=np.zeros([(Len_tPoints),len(fgdvc_list)]) #line index, time, column index: runned case
        v=np.zeros([(Len_tPoints),len(fgdvc_list)]) #line index, time, column index: runned case
        vDot=np.zeros([(Len_tPoints),len(fgdvc_list)]) #line index, time, column index: runned case
        for runnedCaseNr in range(len(fgdvc_list)):
            t[:,runnedCaseNr]=fgdvc_list[runnedCaseNr].Time()[:Len_tPoints]
            dt[:,runnedCaseNr]=fgdvc_list[runnedCaseNr].Dt()[:Len_tPoints]
            T.append(fgdvc_list[runnedCaseNr].Interpolate('Temp'))
            u[:,runnedCaseNr]=fgdvc_list[runnedCaseNr].Yield(Name)[:Len_tPoints]
        #updated Vector: CurrentVector
        #####improve all parameters
        w0=self.a0/(( max((fgdvc_list[0].Yield(Name))) -min((fgdvc_list[0].Yield(Name))) )**2)
        print fgdvc_list[0].SpeciesName(Name)
#        print 'a0: ', self.a0
#        print 'w0: ', w0
        w1=self.a1/(max( ((fgdvc_list[0].Rate(Name)))**2 ))
#        print 'a1: ', self.a1
#        print 'w1: ',w1
        print 'start gradient based optimization'
        #
        def LeastSquareFunction(Parameter):
            model.setParamVector(Parameter)
#            print model.ParamVector()
            for runnedCaseNr in range(len(fgdvc_list)):
                v[:,runnedCaseNr]=model.calcMass(fgdvc_list[runnedCaseNr],t[:,runnedCaseNr],T[runnedCaseNr],Name)[:Len_tPoints]
                uDot[:,runnedCaseNr]=fgdvc_list[runnedCaseNr].Rate(Name)[:Len_tPoints]
                vDot[:,runnedCaseNr]=model.deriveC(fgdvc_list[runnedCaseNr],v[:,runnedCaseNr],Len_tPoints)
            if self.selectedOptimizer=='leastsq':
                Dot1=np.zeros([(Len_tPoints),len(fgdvc_list)],dtype='d')
                Dot2=np.zeros([(Len_tPoints),len(fgdvc_list)],dtype='d')
                for runnedCaseNr in range(len(fgdvc_list)):
                    Dot2[:,runnedCaseNr]=w1*(((uDot[:,runnedCaseNr]-vDot[:,runnedCaseNr])**2)*dt[:,runnedCaseNr]) #the rate term
                    Dot1[:,runnedCaseNr]=w0*(((u[:,runnedCaseNr]-v[:,runnedCaseNr])**2)*dt[:,runnedCaseNr])        #the yield term
                Error=np.zeros((Len_tPoints),dtype='d')
                #makes a long array, containing both, the rates and yields
                Error[:]=np.sum(Dot1+Dot2,axis=1)
#                print np.sum(Error)
            else:
                sumYields_vec=(u-v)**2
                SumYields=np.sum(sumYields_vec*dt)
                SumRates_vec=(uDot-vDot)**2
                SumRates=np.sum(SumRates_vec*dt)
                Error= w0*SumYields+w1*SumRates
#                print Error
            return Error
        model.setParamVector(Parameter_Vector)
        if self.selectedOptimizer=='fmin':
            OptimizedVector=fmin(LeastSquareFunction,model.ParamVector(),ftol=self.Fit_Tolerance,maxiter=self.MaxIter)
            self.FinalDeviation=LeastSquareFunction(model.ParamVector())
            return OptimizedVector
        elif self.selectedOptimizer=='fmin_cg':
            OptimizedVector=fmin_cg(LeastSquareFunction,model.ParamVector(),gtol=self.Fit_Tolerance,maxiter=self.MaxIter)
            self.FinalDeviation=LeastSquareFunction(model.ParamVector())
            return OptimizedVector
        elif self.selectedOptimizer=='fmin_bfgs':
            OptimizedVector=fmin_bfgs(LeastSquareFunction,model.ParamVector(),gtol=self.Fit_Tolerance,maxiter=self.MaxIter)
            self.FinalDeviation=LeastSquareFunction(model.ParamVector())
            return OptimizedVector
        elif self.selectedOptimizer=='fmin_ncg':
            OptimizedVector=fmin_ncg(LeastSquareFunction,model.ParamVector(),avextol=self.Fit_Tolerance)#,maxiter=self.MaxIter)
            self.FinalDeviation=LeastSquareFunction(model.ParamVector())
            return OptimizedVector
        elif self.selectedOptimizer=='fmin_slsqp':
            OptimizedVector=fmin_slsqp(LeastSquareFunction,model.ParamVector(),acc=self.Fit_Tolerance)#,maxiter=self.MaxIter)
            self.FinalDeviation=LeastSquareFunction(model.ParamVector())
            return OptimizedVector
        elif self.selectedOptimizer=='leastsq':
            OptimizedVector=leastsq(LeastSquareFunction,model.ParamVector(),ftol=self.Fit_Tolerance,maxfev=self.MaxIter)
            self.FinalDeviation=LeastSquareFunction(model.ParamVector())
            return OptimizedVector[0]
        else:
            print "No Optimizer was selected. Please choose: 'fmin' or 'fmin_cg' or 'fmin_bfgs' or 'leastsq' or 'fmin_slsqp'\n"
        print  'Optm Vec:   ',OptimizedVector
        
    def Deviation(self):
        """Returns the Deviation after the optimization procedure."""
        return self.FinalDeviation
        

    def setWeights(self,WeightMass,WeightRates):
        """Sets the weights for the yields and the rates for the fitting procedure. See manual for equation."""
        self.a0=WeightMass
        self.a1=WeightRates
        
    def setTolerance(self, ToleranceForFminFunction):
        """Sets the tolerance as a abort criterion for the fitting procedure."""
        self.Fit_Tolerance=ToleranceForFminFunction

    def setPreTolerance(self, ToleranceForFminFunction):
        """Sets the tolerance as a abort criterion for the prefitting procedure (if preLoopNumber in estimate_T is not equal zero)."""
        self.Pre_Tolerance=ToleranceForFminFunction

    def setPreMaxIter(self,MaxiumNumberOfIterationInPreProcedure):
        """Sets the maximum number of iteration oin the optimizer as a abort criterion for the prefitting procedure (if preLoopNumber in estimate_T is not equal zero)."""
        self.PreMaxIter=int(MaxiumNumberOfIterationInPreProcedure)

    def setMaxIter(self,MaxiumNumberOfIterationInMainProcedure):
        """Sets the maximum number of iteration oin the optimizer as a abort criterion for the fitting procedure."""
        self.MaxIter=int(MaxiumNumberOfIterationInMainProcedure)

    def setOptimizer(self,ChosenOptimizer):
        """Select one optimizer of the scipy.optimizer library: 'fmin','fmin_cg','fmin_bfgs','fmin_ncg','fmin_slsqp' or 'leastsq'. According to experience 'fmin' (or also 'leastsq') generates at best the results."""
        self.selectedOptimizer=ChosenOptimizer

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





