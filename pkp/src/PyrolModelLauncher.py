"""
    Launcher functions for the  different pyrolysis models



"""
import numpy as np
import pkp.src.Models as mdl

def constantRate(inputs, results):
    return genericRate(inputs, results, "constantRate")

def arrheniusRate(inputs, results):
    return genericRate(inputs, results, "arrheniusRate")

def C2SM(inputs, results):
    return genericRate(inputs, results, "C2SM")

def genericRate(inputs, results, pyrolModelName):
    """ Generates the results for constant Rate.

        Parameters:
            inputs: dictionary with input parameter
            results: a results  object containing data from
                     preprocessor runs
    """
    # NOTE: Following things have been removed and need to be reimplemented
    #       * plotting
    #       * write to data file
    #       * the species consider thing
    #       * the pccl dt thing
    #       * nrRuns > 1
    # init model
    #TODO implement corrector selector for models
    
    # uses for optimization gradient based (LS) optimizer if NrOfRuns == 1
    # otherwise an Evolutionary algorithm (GenAlg; global optimum) is used
    # Iterate species in Fit Object
    # NOTE it seems strange that we need a model object for every species
    # TODO Needs testing and cleaning of notes
    model = getattr(mdl, pyrolModelName)
    model_inputs = inputs['FIT']
    #FIXME very hacky
    species_names = results[results.keys()[0]].speciesNames
    presc_spec = inputs['FIT'].get('Species', [])
    spec_lst = (presc_spec if presc_spec else species_names)
    return FitResult({species_name : model(
                inputs  = model_inputs,
                runs    = results,
                species = species_name).fit()
        for species_name in species_names
        if species_name in spec_lst})

class FitResult(object):

    def __init__(self, res):
        self.res = res

    @property
    def _tsv(self):
        header = " ".join(self.res.keys())
        s = ""
        res = list(self.res.items())
        # res[n] is a tuple of species name and model of
        # the first species
        for i, _ in enumerate(res[0][1]):
        #     s += str(res[0][1]) + "\t"
             for name, mass in self.res.items():
                 s += str(mass[i]) + "\t"
             s += "\n"
        return "time\t" + header + "\n" + s


# def MakeResults_Arrh(self,PyrolProgram,File,Fit):
#     """Generates the results for Arrhenius Rate."""
#     #makes Species list which contains alls species to fit:
#     SpeciesList=[] # List containing the Species Index to fit
#     if self.ArrhSpec=='Total':
#         SpeciesList.append(Fit[0].SpeciesIndex('Total'))
#         if 'Total' not in self.SpeciesToConsider:
#             self.SpeciesToConsider.append('Total')
#     elif self.ArrhSpec=='MainSpecies':
#         SpeciesList.append(Fit[0].SpeciesIndex('Total'))
#         SpeciesList.append(Fit[0].SpeciesIndex('Tar'))
#         SpeciesList.append(Fit[0].SpeciesIndex('Gas'))
#         if 'Total' not in self.SpeciesToConsider:
#             self.SpeciesToConsider.append('Total')
#         if 'Tar' not in self.SpeciesToConsider:
#             self.SpeciesToConsider.append('Tar')
#         if 'Gas' not in self.SpeciesToConsider:
#             self.SpeciesToConsider.append('Gas')
#     elif self.ArrhSpec=='allSpecies':
#         for i in range(2,len(Fit[0].SpeciesNames()),1):
#             if Fit[0].SpeciesName(i) not in self.SpeciesToConsider:
#                 self.SpeciesToConsider.append(Fit[0].SpeciesName(i))
#             SpeciesList.append(i)
    # ##The single species:
    # for Species in SpeciesList:
    #     #
    #     m_final_prediction=Fit[0].Yield(Species)[-1]
    #     PredictionV0=[0.86e15,0.01,27700,m_final_prediction]  #for Standard Arrhenius
    #     #
    #     self.KinModel=Models.ArrheniusModel(PredictionV0)
    #     if PyrolProgram=='PCCL':
    #         self.KinModel.setDt4Intergrate(self.FG_dt)
    #     #
    #     print Fit[0].SpeciesName(Species)
    #     ParamInit = GlobalOptParam.EvAArrhInit
    #     if len(ParamInit) == 4:
    #         ParamInit.pop(-1)
    #     #
    #     if self.NrOfRuns == 1:
    #         self.OptGradBased(Fit,ParamInit,Fit[0].Yield(Species)[-1],Species)
    #     else:
    #         # init the Parameter for global optimization
    #         m_final_predictionAll=[] # final yields
    #         for i in range(len(Fit)):
    #             m_final_predictionAll.append(Fit[i].Yield(Species)[-1])
    #         ParamMin = GlobalOptParam.EvAArrhMin
    #         ParamMax = GlobalOptParam.EvAArrhMax
    #         if len(ParamMin) == 3:
    #             ParamMin.append(0.0)
    #         if len(ParamMax) == 3:
    #             ParamMax.append(0.0)
    #         if len(ParamInit) == 3:
    #             ParamInit.append(0.0)
    #         ParamInit[3] = (max(m_final_predictionAll)+min(m_final_predictionAll))/2.
    #         ParamMin[3] = (min(m_final_predictionAll))
    #         ParamMax[3] = (max(m_final_predictionAll))
    #         #
    #         self.OptGenAlgBased(Fit,ParamInit,ParamMin,ParamMax,Species)
#         #
#         self.KinModel.plot(Fit,Species)
#         self.Solution=self.KinModel.ParamVector()
#         #To avoid, a species with no yield is added to the parameter file
#         if np.sum(self.KinModel.ParamVector())!=np.sum(PredictionV0):
#             outfile.write(str(Fit[0].SpeciesName(Species))+'\t'+'%.6e  %6.4f  %11.4f  %7.4f  ' %(self.Solution[0],self.Solution[1],self.Solution[2],self.Solution[3])+'\n')
#     outfile.close()
#
#
# def MakeResults_ArrhNoB(self,PyrolProgram,File,Fit):
#     """Generates the results for Arrhenius Rate with no correction term T**b."""
#     outfile = open(PyrolProgram+'-Results_ArrheniusNoBRate.txt', 'w')
#     outfile.write("Species               A [1/s]                E_a [K]      FinalYield\n\n")
#     #makes Species list which contains alls species to fit:
#     SpeciesList=[]
#     if self.ArrhSpec=='Total':
#         SpeciesList.append(Fit[0].SpeciesIndex('Total'))
#         if 'Total' not in self.SpeciesToConsider:
#             self.SpeciesToConsider.append('Total')
#     elif self.ArrhSpec=='MainSpecies':
#         SpeciesList.append(Fit[0].SpeciesIndex('Total'))
#         SpeciesList.append(Fit[0].SpeciesIndex('Tar'))
#         SpeciesList.append(Fit[0].SpeciesIndex('Gas'))
#         if 'Total' not in self.SpeciesToConsider:
#             self.SpeciesToConsider.append('Total')
#         if 'Tar' not in self.SpeciesToConsider:
#             self.SpeciesToConsider.append('Tar')
#         if 'Gas' not in self.SpeciesToConsider:
#             self.SpeciesToConsider.append('Gas')
#     elif self.ArrhSpec=='allSpecies':
#         for i in range(2,len(Fit[0].SpeciesNames()),1):
#             if Fit[0].SpeciesName(i) not in self.SpeciesToConsider:
#                 self.SpeciesToConsider.append(Fit[0].SpeciesName(i))
#             SpeciesList.append(i)
#     ##The single species:
#     for Species in SpeciesList:
#         m_final_prediction=Fit[0].Yield(Species)[-1]
#         PredictionV0=[0.86e11,17700,m_final_prediction]  #for Standard Arrhenius
#         self.KinModel=Models.ArrheniusModelNoB(PredictionV0)
#         if PyrolProgram=='PCCL':
#             self.KinModel.setDt4Intergrate(self.FG_dt)
#         #
#         print Fit[0].SpeciesName(Species)
#         # gets final yields for all runs
#         m_final_predictionAll=[]
#         for i in range(len(Fit)):
#             m_final_predictionAll.append(Fit[i].Yield(Species)[-1])
#         ParamInit = GlobalOptParam.EvAArrhInit
#         if len(ParamInit) == 4:
#             ParamInit.pop(-1)
#         if len(ParamInit) == 3:
#             ParamInit.pop(1)
#         # optimization procedure
#         if self.NrOfRuns == 1:
#             self.OptGradBased(Fit,ParamInit,Fit[0].Yield(Species)[-1],Species)
#         else:
#             ParamMin = GlobalOptParam.EvAArrhMin
#             ParamMax = GlobalOptParam.EvAArrhMax
#             if len(ParamMin) == 3:
#                 ParamMin.pop(1)
#                 ParamMin.append(0.0)
#             if len(ParamMax) == 3:
#                 ParamMax.pop(1)
#                 ParamMax.append(0.0)
#             if len(ParamInit) == 3: # if final yield was appended
#                 ParamMax.pop(-1)
#             if len(ParamInit) == 2:
#                 ParamInit.append(0.0)
#             ParamInit[2] = (max(m_final_predictionAll)+min(m_final_predictionAll))/2.
#             ParamMin[2] = (min(m_final_predictionAll))
#             ParamMax[2] = (max(m_final_predictionAll))
#             #
#             self.OptGenAlgBased(Fit,ParamInit,ParamMin,ParamMax,Species)
#             #
#         self.Solution=self.KinModel.ParamVector()
#         self.KinModel.plot(Fit,Species)
#         if np.sum(self.KinModel.ParamVector())!=np.sum(PredictionV0): #To avoid, a species with no yield is added to the parameter file
#             outfile.write(str(Fit[0].SpeciesName(Species))+'\t'+'%.6e  %11.4f  %7.4f  ' %(self.Solution[0],self.Solution[1],self.Solution[2])+'\n')
#     outfile.close()
#     if oSystem=='Linux'  or oSystem == 'Darwin':
#         shutil.move(PyrolProgram+'-Results_ArrheniusNoBRate.txt','Result/'+PyrolProgram+'-Results_ArrheniusNoB.txt')
#     elif oSystem=='Windows':
#         shutil.move(PyrolProgram+'-Results_ArrheniusNoBRate.txt','Result\\'+PyrolProgram+'-Results_ArrheniusNoB.txt')
#     else:
#         print "The name of the operating system couldn't be found."
# #
# #
# def MakeResults_Kob(self,PyrolProgram,File,Fit):
#     """Generates the results for Kobayashi Rate."""
#     PredictionVKob0=[7e5,8e7/8314.33,2.3e8,1.6e8/8314.33,0.4,0.9]
#     outfile = open(PyrolProgram+'-Results_KobayashiRate.txt', 'w')
#     outfile.write("Species              A1 [1/s]              E_a1 [kJ/mol]             A2 [1/s]             E_a2 [kJ/mol]      alpha1     alpha2\n\n")
#     self.KinModel=Models.Kobayashi(GlobalOptParam.EvAKobInit) #(PredictionVKob0)
#     if PyrolProgram=='PCCL':
#             self.KinModel.setDt4Intergrate(self.FG_dt)
#     #######
#     ##The single species:
#     if 'Total' not in self.SpeciesToConsider:
#         self.SpeciesToConsider.append('Total')
#     for Species in [Fit[0].SpeciesIndex('Total')]:
#         # optimization procedure
#         print Fit[0].SpeciesName(Species)
#         if self.NrOfRuns == 1:
#             self.OptGradBased(Fit,self.KinModel.ParamVector(),False,Species)
#         else:
#             self.OptGenAlgBased(Fit,self.KinModel.ParamVector(),GlobalOptParam.EvAKobMin,GlobalOptParam.EvAKobMax,Species)
#         #
#         self.Solution=self.KinModel.ParamVector()
#         #
#         self.KinModel.plot(Fit,Species)
#         outfile.write(str(Fit[0].SpeciesName(Species))+'\t'+'%.6e  %11.4f  %.6e  %11.4f  %6.4f  %6.4f  ' %(self.Solution[0],self.Solution[1]*8314.33/1e6,self.Solution[2],self.Solution[3]*8314.33/1e6,self.Solution[4],self.Solution[5])+'\n')
#     outfile.close()
#     if oSystem=='Linux'  or oSystem == 'Darwin':
#         shutil.move(PyrolProgram+'-Results_KobayashiRate.txt','Result/'+PyrolProgram+'-Results_Kobayashi.txt')
#     elif oSystem=='Windows':
#         shutil.move(PyrolProgram+'-Results_KobayashiRate.txt','Result\\'+PyrolProgram+'-Results_Kobayashi.txt')
#     else:
#         print "The name of the operating system couldn't be found."
# #
# #
# def MakeResults_DEAM(self,PyrolProgram,File,Fit):
#     """Generates the results for DAEM model."""
# #    elif (CPD_FittingKineticParameter_Select=='DAEM' and PyrolProgram=='CPD') or (FG_FittingKineticParameter_Select=='DAEM' and PyrolProgram=='FGDVC'):
#     PredictionDAEM=[2e10,20e3,5e3,0.5]
#     outfile = open(PyrolProgram+'-Results_DAEM.txt', 'w')
#     outfile.write("Species                  A1 [1/s]               E_a1 [K]          sigma [K]   Final Yield\n\n")
#     self.KinModel=Models.DAEM(PredictionDAEM)
#     self.KinModel.setNrOfActivationEnergies(GlobalOptParam.NrOFActivtionEnergies)
#     if PyrolProgram=='PCCL':
#             self.KinModel.setDt4Intergrate(self.FG_dt)
#     #######
#     ParamInit = GlobalOptParam.EvADAEMInit
#     if len(ParamInit) == 4:
#             ParamInit.pop(-1)
#     ParamMin = GlobalOptParam.EvADAEMMin
#     if len(ParamMin) == 3:
#         ParamMin.append(0.0)
#     ParamMax = GlobalOptParam.EvADAEMMax
#     if len(ParamMax) == 3:
#         ParamMax.append(0.0)
#     ##The single species:
#     if 'Total' not in self.SpeciesToConsider:
#         self.SpeciesToConsider.append('Total')
#     for Species in [Fit[0].SpeciesIndex('Total')]:
#         print Fit[0].SpeciesName(Species)
#         m_final_predictionAll=[]
#         for i in range(len(Fit)):
#             m_final_predictionAll.append(Fit[i].Yield(Species)[-1])
#         # optimization procedure
#         if self.NrOfRuns == 1:
#             self.OptGradBased(Fit,ParamInit,Fit[0].Yield(Species)[-1],Species)
#         else:
#             if len(ParamInit) == 3:
#                 ParamInit.append(0.0)
#             ParamInit[3]= ((min(m_final_predictionAll)+max(m_final_predictionAll))/2.)
#             ParamMin[3]= (min(m_final_predictionAll))
#             ParamMax[3]= (max(m_final_predictionAll))
#             self.OptGenAlgBased(Fit,ParamInit,ParamMin,ParamMax,Species)
#         #
#         self.Solution=self.KinModel.ParamVector()
#         #
#         self.KinModel.plot(Fit,Species)
#         outfile.write(str(Fit[0].SpeciesName(Species))+'\t'+'%.6e  %11.4f  %11.4f  %6.4f  ' %(self.Solution[0],self.Solution[1],self.Solution[2],self.Solution[3])+'\n')
#     outfile.close()
#     if oSystem=='Linux' or oSystem == 'Darwin':
#         shutil.move(PyrolProgram+'-Results_DAEM.txt','Result/'+PyrolProgram+'-Results_DAEM.txt')
#     elif oSystem=='Windows':
#         shutil.move(PyrolProgram+'-Results_DAEM.txt','Result\\'+PyrolProgram+'-Results_DAEM.txt')
#     else:
#         print "The name of the operating system couldn't be found."
