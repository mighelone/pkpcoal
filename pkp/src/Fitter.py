import numpy as np
from scipy.optimize import fmin
from scipy.optimize import fmin_cg
from scipy.optimize import fmin_bfgs
from scipy.optimize import fmin_ncg
from scipy.optimize import leastsq
from scipy.optimize import fmin_slsqp

from pkp.src.Models import Model
from pkp.src.Models import ModelError

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
    optModel = ls.estimate(results, model, species)
    print 'Final error= ' +  str(ls.deviation)
    return optModel

def OptGenAlgBased(inputs, model, results, species):
    """ Starts a genetic algorithm and afterwards a gradient Based optimization.
        Sets the Final Fit result as the ParamVector in the Kinetic Model.
        Input are the Fit (Result Objects of the Detailed Models),
        the Parameter to initialize, the two Parameter vectors defining
        the range of the results and the Species index.
        TODO: Revise
    """
    from pkp.src import Evolve
    genAlg = Evolve.GenericOpt(inputs['Optimisation'])
    optModel = genAlg.estimate(results, model, species)
    # TODO we have to reset model.parameter to the values
    #      after genetic opt
    # TODO what is the reason for this test? 
    # afterwards grad based optimization
    # if GlobalOptParam.optimizGrad == True:
        #OptGradBased(Fit, ParameterVecInit, False, Species)
    return OptGradBased(inputs, optModel, results, species)


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
        import scipy.optimize as scopt

        print 'start gradient based optimization, species: ' + species
        print 'initial parameter: ' + str(model.parameter) 
        optimiser = getattr(scopt, self.optimizer)
        # select optimisation input function depending on the optimizer
        # this is needed since leastsq expects a list of errors where fmin
        # simply expects a global error
        error_func = (ModelError.ls_input_func if self.optimizer == 'leastsq' else ModelError.min_input_func)
        model_error = ModelError(results, model, species, error_func, self.weightMass, self.weightRate)
        OptimizedVector = optimiser(
                func = model_error.input_func,
                x0   = model.parameter,
                #args = (error_func, model, results, species),
                ftol = self.fitTolerance, # TODO GO what is the diff between gtol&ftol
                maxiter = self.maxIter
        ) # caLculates optimized vector
        # NOTE It seems unneccessary to reevaluate the model agian,
        #      but for now no better solution is in sight,
        #      so we will use it to store final yields and rates on the model
        self.deviation = model_error.error
        print 'final parameter: ' + str(model.parameter) 
        return model
