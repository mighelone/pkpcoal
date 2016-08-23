from __future__ import division
import sys
import os
import pkp.CPD_SetAndLaunch as CPD_SetAndLaunch
import pkp.FGDVC_SetAndLaunch as FGDVC_SetAndLaunch
import pkp.PCCL_SetAndLaunch as PCCL_SetAndLaunch
import pkp.FGDVC_Result as FGDVC_Result
import pkp.CPD_Result as CPD_Result
import pkp.PCCL_Result as PCCL_Result
import pkp.Fitter as Fitter
import pkp.Models as Models
import pkp.FitInfo as FitInfo
import pkp.Compos_and_Energy as Compos_and_Energy
import pkp.InformationFiles as InformationFiles
import GlobalOptParam
import pkp.Evolve as Evolve
import numpy as np
import platform
import shutil
# import yaml
import ruamel_yaml as yaml
import logging
import argparse

# Which operating Sytem?
oSystem = platform.system()
# if oSystem == 'Darwin':
#    oSystem = 'Linux'
# Directories:
# gets the current directory:
workingDir = os.getcwd() + '/'
# FG-DVC Library coals folder name:
FG_LibCoalDir = 'input'
# FG-DVC coal generation (coalsd.exe) folder name:
FG_GenCoalDir = 'coals'
# FG-DVC fgdvcd.exe folder name:
FG_ExeCoalDir = 'FGDVC'
# Input file name for coalsd.exe:
FG_CoalGenFileName = 'CoalGen_FGDVC.txt'
# File name for generated coal file:
FG_CoalName = 'GenCoal'
#
#


class MainProcess(object):
    """Controls the whole process of generating input files, fitting etc."""

    def __init__(self):
        self.SpeciesToConsider = []  # for GUI
        self.ProgramModelDict = {}  # for GUI
        self.logger_root = 'main.' + self.__class__.__name__ + '.'
    #

    def ReadInputFiles(self, yml_file=None):
        """get parameters from input files"""
        #
        # Coal File:
        #
        logger = logging.getLogger(self.logger_root + 'ReadInputFiles')
        if not yml_file:
            yml_file = 'input.yml'
        logger.info('Reading %s', yml_file)
        with open(yml_file, 'r') as f:
            yml_input = yaml.load(f)

        # print 'Reading Coal.inp ...'
        logger.debug('Define coal properties')
        CoalInput = InformationFiles.ReadFile(workingDir + 'Coal.inp')
        coal_input = yml_input['Coal']
        proximate_analysis = coal_input['Proximate Analysis']
        logger.debug('PA read: %s', proximate_analysis)
        # normalize
        sum_pa = sum(proximate_analysis.itervalues())
        logger.debug('Sum PA %s', sum_pa)
        proximate_analysis = {el: (val / sum_pa * 100.)
                              for el, val in proximate_analysis.iteritems()}
        logger.debug('PA normalized: %s', proximate_analysis)

        self.PA_asrec = proximate_analysis

        self.PAFC_asrec = proximate_analysis['Fixed Carbon']
        self.PAVM_asrec = proximate_analysis['Volatile Matter']
        self.PAash = proximate_analysis['Ash']
        self.PAmoist = proximate_analysis['Moisture']
        logger.debug('PA FC: %s', self.PAFC_asrec)
        logger.debug('PA VM: %s', self.PAVM_asrec)
        #
        # gets daf values, as CPD needs daf as input:
        self.PAFC_daf, self.PAVM_daf = self.DAF(
            self.PAFC_asrec, self.PAVM_asrec)
        logger.debug('PA FC daf: %s', self.PAFC_daf)
        logger.debug('PA VM daf: %s', self.PAVM_daf)

        ultimate_analysis = coal_input['Ultimate Analysis']
        logger.debug('UA read %s', ultimate_analysis)
        sum_ua = sum(ultimate_analysis.itervalues())
        ultimate_analysis = {el: (val / sum_ua * 100)
                             for el, val in ultimate_analysis.iteritems()}
        logger.debug('UA normalized %s', ultimate_analysis)

        self.UAC, self.UAH, self.UAN, self.UAO, self.UAS = (
            ultimate_analysis[el] for el in ['C', 'H', 'N', 'O', 'S']
        )
        logger.debug('UAC %s', self.UAC)
        logger.debug('UAH %s', self.UAH)
        logger.debug('UAN %s', self.UAN)
        logger.debug('UAO %s', self.UAO)
        logger.debug('UAS %s', self.UAS)

        self.HHV = coal_input['HHV']
        logger.debug('HHV as rec: %s', self.HHV)
        self.densityDryCoal = coal_input['rhoDry']
        logger.debug('Rho dry: %s', self.densityDryCoal)

        self.WeightY = yml_input['Calibration']['weightYield']
        self.WeightR = yml_input['Calibration']['weightRate']
        logger.debug('Weight R %s Y %s', self.WeightR, self.WeightY)

        #
        # CPD Properties:
        #
        logger.info('Set CPD..')
        cpd_input = yml_input['CPD']
        self.MTar = cpd_input['MW_TAR']
        logger.debug('MW TAR (only for CPD) %s', self.MTar)
        self.CPDselect = cpd_input['active']
        logger.debug('CPDselect new %s', self.CPDselect)
        self.CPD_FittingKineticParameter_Select = cpd_input['fit']
        # self.CPD_FittingKineticParameter_Select = yml_input[
        #    'FIT']['Model']
        logger.debug('CPD fit new %s',
                     self.CPD_FittingKineticParameter_Select)

        # 0:initial dt, 1: print increment, 2: dt max
        self.CPDdt = [0, 1, 2]
        self.CPDdt[0:2] = (cpd_input['deltaT'], cpd_input['increment'])
        self.CPDdt[2] = 1e-5  # check the meaning of this option
        logger.debug('CPDdt %s', self.CPDdt)

        #
        #
        # FG-DVC Properties:
        #
        print 'Reading FGDVC.inp ...'
        FGDVCInput = InformationFiles.ReadFile(workingDir + 'FGDVC.inp')
        self.FG_select = FGDVCInput.UsePyrolProgr(
            InformationFiles.MF_sel)
        self.FG_FittingKineticParameter_Select = FGDVCInput.Fitting(
            InformationFiles.M_selFit)
        self.FG_CoalSelection = int(
            FGDVCInput.getValue(InformationFiles.MF_CoalSel))
        self.FG_MainDir = FGDVCInput.getText(InformationFiles.MF_dir[0])
        self.FG_DirOut = FGDVCInput.getText(InformationFiles.MF_dir[1])
        self.FG_TarCacking = FGDVCInput.getValue(
            InformationFiles.MF_TarCr)
        #
        # Polimi PMSKD model Properties
        # Predictive Multi-Step Kinetic Devolatilization
        #
        print 'Set PMSKD'
        pmskd_input = yml_input['PMSKD']
        self.PMSKD_select = pmskd_input['active']
        self.PMSKD_FittingKineticParameter_Select = pmskd_input['fit']
        self.PMSKD_npoint = pmskd_input['nStep']
        self.PMSKD_mechfile = pmskd_input['mechFile']
        self.PMSKD_ArrhSpec = 'Total'  # !!TODO check the meaning of this par
        logger.debug('PMSKD select %s', self.PMSKD_select)
        logger.debug('PMSKD Fitting %s',
                     self.PMSKD_FittingKineticParameter_Select)
        logger.debug('PMSKD ArrhSpecies %s', self.PMSKD_ArrhSpec)
        logger.debug('PMSKD npoints %s', self.PMSKD_npoint)
        logger.debug('PMSKD mechfile %s', self.PMSKD_mechfile)

        # BIO Polimi
        print 'Set bioPolimi ...'
        # BioInput = InformationFiles.ReadFile(
        #    workingDir + 'bioPolimi.inp')

        # self.bio_dict = {} # !! use consistent dictionary
        self.bio_dict = yml_input['BioPolimi']
        logger.debug('bio dict %s', self.bio_dict)

        # self.bio_dict['Cellulose'] = BioInput.getText('Cellulose:')
        # self.bio_dict['Hemicellulose'] = BioInput.getText('Hemicellulose:')
        # self.bio_dict['LigninC'] = BioInput.getText('LigninC:')
        # self.bio_dict['LigninH'] = BioInput.getText('LigninH:')
        # self.bio_dict['LigninO'] = BioInput.getText('LigninO:')

        #
        # PC Coal Lab Properties:
        #
        print 'Reading PCCL.inp ...'
        PCCLInput = InformationFiles.ReadFile(workingDir + 'PCCL.inp')
        self.PCCL_select = PCCLInput.UsePyrolProgr(
            InformationFiles.MPC_sel)
        self.PCCL_FittingKineticParameter_Select = PCCLInput.Fitting(
            InformationFiles.M_selFit)
        self.PCCL_Path = PCCLInput.getText(InformationFiles.MPC_dir)
        self.PCCL_Exe = PCCLInput.getText(InformationFiles.MPC_exe)
        try:
            self.PCCL_CoalCalFactor = float(
                PCCLInput.getText(InformationFiles.MPC_CoalCal))
        except ValueError:
            self.PCCL_CoalCalFactor = False
        self.PCCL_ParticleSize = PCCLInput.getValue(
            InformationFiles.MPC_partSize)
        #
        #
        # Operating Condition File:
        #
        logger.info('Set Operating condition')
        opconditions_input = yml_input['OperatingConditions']
        self.CPD_pressure = self.FG_pressure = opconditions_input[
            'pressure']  # atm
        self.NrOfRuns = opconditions_input['runs']
        for i in xrange(1, 6):
            time_temp = np.array(opconditions_input['run{}'.format(i)])
            setattr(self, 'TimeTemp{}'.format(i), time_temp)
            time_cpd = np.empty_like(time_temp)
            time_cpd[:, 0] = time_temp[:, 0] * 1e3
            time_cpd[:, 1] = time_temp[:, 1]
            setattr(self, 'CPD_TimeTemp{}'.format(i), time_cpd)
            logger.debug('Run %s: %s', i, getattr(
                self, 'TimeTemp{}'.format(i)))

        OpCondInp = InformationFiles.OperCondInput('OperCond.inp')
        # self.CPD_pressure = OpCondInp.getValue(
        #    InformationFiles.M_Pressure)
        # self.FG_pressure = OpCondInp.getValue(
        #    InformationFiles.M_Pressure)
        self.ArrhSpec = OpCondInp.getText(
            InformationFiles.M_selArrhSpec)
        logger.debug('ArrhSpec %s', self.ArrhSpec)

        # organize time temp for Polimi model
        self.timeHR = [self.TimeTemp1[:, 0],
                       self.TimeTemp2[:, 0],
                       self.TimeTemp3[:, 0],
                       self.TimeTemp4[:, 0],
                       self.TimeTemp5[:, 0]]
        self.temperatureHR = [self.TimeTemp1[:, 1],
                              self.TimeTemp2[:, 1],
                              self.TimeTemp3[:, 1],
                              self.TimeTemp4[:, 1],
                              self.TimeTemp5[:, 1]]

        # set when activate FG-DVC
        # self.FG_dt = OpCondInp.getValue(InformationFiles.M_dt)
        self.FG_T_t_History = self.FG_MainDir + 'tTHistory.txt'

        # self.CPD_TimeTemp1 = OpCondInp.getTimePoints(
        #    InformationFiles.M_TimePoints1[0], InformationFiles.M_TimePoints1[1])
        # logger.debug('CPD_timetemp %s', self.CPD_TimeTemp1)
        # logger.debug('TimeTemp1 %s', self.TimeTemp1)
        # sys.exit()
        # self.CPD_TimeTemp2 = OpCondInp.getTimePoints(
        #    InformationFiles.M_TimePoints2[0], InformationFiles.M_TimePoints2[1])
        # self.CPD_TimeTemp3 = OpCondInp.getTimePoints(
        #    InformationFiles.M_TimePoints3[0], InformationFiles.M_TimePoints3[1])
        # self.CPD_TimeTemp4 = OpCondInp.getTimePoints(
        #    InformationFiles.M_TimePoints4[0], InformationFiles.M_TimePoints4[1])
        # self.CPD_TimeTemp5 = OpCondInp.getTimePoints(
        #    InformationFiles.M_TimePoints5[0], InformationFiles.M_TimePoints5[1])
        # makes for CPD time in milliseconds:

        # tmax in s, not ms
        self.CPD_t_max1 = self.CPD_TimeTemp1[-1, 0] * 1.e-3
        # tmax in s, not ms
        self.CPD_t_max2 = self.CPD_TimeTemp2[-1, 0] * 1.e-3
        # tmax in s, not ms
        self.CPD_t_max3 = self.CPD_TimeTemp3[-1, 0] * 1.e-3
        # tmax in s, not ms
        self.CPD_t_max4 = self.CPD_TimeTemp4[-1, 0] * 1.e-3
        # tmax in s, not ms
        self.CPD_t_max5 = self.CPD_TimeTemp5[-1, 0] * 1.e-3

        #logger.debug('CPD_TimeTemp1 %s', self.CPD_TimeTemp1)
        #logger.debug('Max CPD t %s', self.CPD_t_max1)
        #
        #
        #

    def DAF(self, PAFC_asRecieved, PAVM_asRecieved):
        """calculates PAFC, PAVM  from the as recieved state to the daf state of coal"""
        fractionFC = PAFC_asRecieved / \
            (PAFC_asRecieved + PAVM_asRecieved)
        fractionVM = PAVM_asRecieved / \
            (PAFC_asRecieved + PAVM_asRecieved)
        return 100. * fractionFC, 100. * fractionVM
        #

    def CheckFGdt(self):
        """Aborts, if FG-DVC is selected and the timestep is lower than 1.e-3 (which is FG-DVC not able to read):"""
        if ((self.FG_select == True) and (self.FG_dt < 1e-4)):
            print "Please select for FG-DVC a time step greather equal 1e-4 in 'OperCond.inp'. FG-DVC would not be able to read the time history file for a dt<1e-4."
            sys.exit()

    def OptGradBased(self, Fit, ParameterVecInit, FinalYield, Species):
        """Starts a gradient Based Optimization. Sets the Final Fit result as the ParamVector in the Kinetic Model. Input are the Fit (Result Objects of the Detailed Models), the Parameter to initialize and the Final Yield (all dependent on the kinetic model). For Kobayashi Model set Final Yield to False (independent), for all other set a value. It will be excluded from the optimization. SPecies is the Species Index."""
        LS = Fitter.LeastSquarsEstimator()
        LS.setOptimizer(GlobalOptParam.selectedGradBasedOpt)
        LS.setTolerance(GlobalOptParam.Tolerance)
        LS.setMaxIter(GlobalOptParam.MaxIter)
        LS.setWeights(self.WeightY, self.WeightR)
        LS.setFinalYield(FinalYield)
        self.KinModel.setParamVector(LS.estimate_T(
            Fit, self.KinModel, ParameterVecInit, Species))
        print 'Final error=', LS.Deviation()

    def OptGenAlgBased(self, Fit, ParameterVecInit, ParameterVecMin, ParameterVecMax, Species):
        """Starts a genetic algorithm and afterwards a gradient Based optimization. Sets the Final Fit result as the ParamVector in the Kinetic Model. Input are the Fit (Result Objects of the Detailed Models), the Parameter to initialize, the two Parameter vectors defining the range of the results and the Species index."""
        GenAlg = Evolve.GenericOpt(self.KinModel, Fit, Species)
        GenAlg.setWeights(self.WeightY, self.WeightR)
        GenAlg.setParamRanges(
            ParameterVecInit, ParameterVecMin, ParameterVecMax)
        GenAlg.setNrPopulation(GlobalOptParam.NrOfPopulation)
        GenAlg.setNrGenerations(GlobalOptParam.NrOfGeneration)
        self.KinModel.setParamVector(GenAlg.mkResults())
        # afterwards grad based optimization
        if GlobalOptParam.optimizGrad == True:
            self.OptGradBased(Fit, ParameterVecInit, False, Species)

    def MakeResults_CR(self, PyrolProgram, File, Fit):
        """Generates the results for constant Rate."""
        outfile = open(PyrolProgram + '-Results_const_rate.txt', 'w')
        outfile.write(
            "Species     k [1/s]     t_start [s]   FinalYield\n\n")
        # init model
        ParamInit = GlobalOptParam.EvACRInit
        self.KinModel = Models.ConstantRateModel(ParamInit)
        if PyrolProgram == 'PCCL':
            self.KinModel.setDt4Intergrate(self.FG_dt)
        # uses for optimization gradient based (LS) optimizer if
        # NrOfRuns is one ; for more runs use Evolutionary algorithm
        # (GenAlg; global optimum)
        if self.NrOfRuns == 1:
            for Spec in range(2, len(Fit[0].SpeciesNames()), 1):
                #
                if Fit[0].SpeciesName(Spec) not in self.SpeciesToConsider:
                    self.SpeciesToConsider.append(
                        Fit[0].SpeciesName(Spec))
                #
                self.OptGradBased(Fit, ParamInit, Fit[
                                  0].Yield(Spec)[-1], Spec)
                self.KinModel.plot(Fit, Spec)
                self.Solution = self.KinModel.ParamVector()
                # if true nothing was optimized, no result to show
                if np.sum(self.Solution) != np.sum(ParamInit):
                    outfile.write(str(Fit[0].SpeciesName(Spec)) + '\t' + '%8.4f  %8.4f  %8.4f  ' % (
                        self.Solution[0], self.Solution[1], self.Solution[2]) + '\n')
        else:
            if len(ParamInit) == 2:
                ParamInit.append(0.0)
            ParamMin = GlobalOptParam.EvACRMin
            if len(ParamMin) == 2:
                ParamMin.append(0.0)
            ParamMax = GlobalOptParam.EvACRMax
            if len(ParamMax) == 2:
                ParamMax.append(0.0)
            for Spec in range(2, len(Fit[0].SpeciesNames()), 1):
                # max Yield, min Yield
                m_final_predictionAll = []
                for i in range(len(Fit)):
                    m_final_predictionAll.append(Fit[i].Yield(Spec)[-1])
                ParamInit[2] = (max(m_final_predictionAll) +
                                min(m_final_predictionAll)) / 2.
                ParamMin[2] = (min(m_final_predictionAll))
                ParamMax[2] = (max(m_final_predictionAll))
                #
                self.OptGenAlgBased(
                    Fit, ParamInit, ParamMin, ParamMax, Spec)
                self.Solution = self.KinModel.ParamVector()
                self.KinModel.plot(Fit, Spec)
                # if True nothing was optimized
                if np.sum(self.Solution) != np.sum(ParamInit):
                    outfile.write(str(Fit[0].SpeciesName(Spec)) + '\t' + '%8.4f  %8.4f  %8.4f  ' % (
                        self.Solution[0], self.Solution[1], self.Solution[2]) + '\n')
            #
            #
        outfile.close()
        if oSystem == 'Linux' or oSystem == 'Darwin':
            shutil.move(PyrolProgram + '-Results_const_rate.txt',
                        'Result/' + PyrolProgram + '-Results_constantRate.txt')
        elif oSystem == 'Windows':
            shutil.move(PyrolProgram + '-Results_const_rate.txt',
                        'Result\\' + PyrolProgram + '-Results_constantRate.txt')
        else:
            print "The name of the operating system couldn't be found."

    def MakeResults_Arrh(self, PyrolProgram, File, Fit):
        """Generates the results for Arrhenius Rate."""
        outfile = open(PyrolProgram + '-Results_ArrheniusRate.txt', 'w')
        outfile.write(
            "Species                         A [1/s]         b               E_a [K]    FinalYield\n\n")
        # makes Species list which contains alls species to fit:
        SpeciesList = []  # List containing the Species Index to fit
        if self.ArrhSpec == 'Total':
            SpeciesList.append(Fit[0].SpeciesIndex('Total'))
            if 'Total' not in self.SpeciesToConsider:
                self.SpeciesToConsider.append('Total')
        elif self.ArrhSpec == 'MainSpecies':
            SpeciesList.append(Fit[0].SpeciesIndex('Total'))
            SpeciesList.append(Fit[0].SpeciesIndex('Tar'))
            SpeciesList.append(Fit[0].SpeciesIndex('Gas'))
            if 'Total' not in self.SpeciesToConsider:
                self.SpeciesToConsider.append('Total')
            if 'Tar' not in self.SpeciesToConsider:
                self.SpeciesToConsider.append('Tar')
            if 'Gas' not in self.SpeciesToConsider:
                self.SpeciesToConsider.append('Gas')
        elif self.ArrhSpec == 'allSpecies':
            for i in range(2, len(Fit[0].SpeciesNames()), 1):
                if Fit[0].SpeciesName(i) not in self.SpeciesToConsider:
                    self.SpeciesToConsider.append(Fit[0].SpeciesName(i))
                SpeciesList.append(i)
        # The single species:
        for Species in SpeciesList:
            #
            m_final_prediction = Fit[0].Yield(Species)[-1]
            # for Standard Arrhenius
            PredictionV0 = [0.86e15, 0.01, 27700, m_final_prediction]
            #
            self.KinModel = Models.ArrheniusModel(PredictionV0)
            if PyrolProgram == 'PCCL':
                self.KinModel.setDt4Intergrate(self.FG_dt)
            #
            print Fit[0].SpeciesName(Species)
            ParamInit = GlobalOptParam.EvAArrhInit
            if len(ParamInit) == 4:
                ParamInit.pop(-1)
            #
            if self.NrOfRuns == 1:
                self.OptGradBased(Fit, ParamInit, Fit[
                                  0].Yield(Species)[-1], Species)
            else:
                # init the Parameter for global optimization
                m_final_predictionAll = []  # final yields
                for i in range(len(Fit)):
                    m_final_predictionAll.append(
                        Fit[i].Yield(Species)[-1])
                ParamMin = GlobalOptParam.EvAArrhMin
                ParamMax = GlobalOptParam.EvAArrhMax
                if len(ParamMin) == 3:
                    ParamMin.append(0.0)
                if len(ParamMax) == 3:
                    ParamMax.append(0.0)
                if len(ParamInit) == 3:
                    ParamInit.append(0.0)
                ParamInit[3] = (max(m_final_predictionAll) +
                                min(m_final_predictionAll)) / 2.
                ParamMin[3] = (min(m_final_predictionAll))
                ParamMax[3] = (max(m_final_predictionAll))
                #
                self.OptGenAlgBased(
                    Fit, ParamInit, ParamMin, ParamMax, Species)
            #
            self.KinModel.plot(Fit, Species)
            self.Solution = self.KinModel.ParamVector()
            # To avoid, a species with no yield is added to the
            # parameter file
            if np.sum(self.KinModel.ParamVector()) != np.sum(PredictionV0):
                outfile.write(str(Fit[0].SpeciesName(Species)) + '\t' + '%.6e  %6.4f  %11.4f  %7.4f  ' % (
                    self.Solution[0], self.Solution[1], self.Solution[2], self.Solution[3]) + '\n')
        outfile.close()
        if oSystem == 'Linux' or oSystem == 'Darwin':
            shutil.move(PyrolProgram + '-Results_ArrheniusRate.txt',
                        'Result/' + PyrolProgram + '-Results_Arrhenius.txt')
        elif oSystem == 'Windows':
            shutil.move(PyrolProgram + '-Results_ArrheniusRate.txt',
                        'Result\\' + PyrolProgram + '-Results_Arrhenius.txt')
        else:
            print "The name of the operating system couldn't be found."

    def MakeResults_ArrhNoB(self, PyrolProgram, File, Fit):
        """Generates the results for Arrhenius Rate with no correction term T**b."""
        outfile = open(
            PyrolProgram + '-Results_ArrheniusNoBRate.txt', 'w')
        outfile.write(
            "Species               A [1/s]                E_a [K]      FinalYield\n\n")
        # makes Species list which contains alls species to fit:
        SpeciesList = []
        if self.ArrhSpec == 'Total':
            SpeciesList.append(Fit[0].SpeciesIndex('Total'))
            if 'Total' not in self.SpeciesToConsider:
                self.SpeciesToConsider.append('Total')
        elif self.ArrhSpec == 'MainSpecies':
            SpeciesList.append(Fit[0].SpeciesIndex('Total'))
            SpeciesList.append(Fit[0].SpeciesIndex('Tar'))
            SpeciesList.append(Fit[0].SpeciesIndex('Gas'))
            if 'Total' not in self.SpeciesToConsider:
                self.SpeciesToConsider.append('Total')
            if 'Tar' not in self.SpeciesToConsider:
                self.SpeciesToConsider.append('Tar')
            if 'Gas' not in self.SpeciesToConsider:
                self.SpeciesToConsider.append('Gas')
        elif self.ArrhSpec == 'allSpecies':
            for i in range(2, len(Fit[0].SpeciesNames()), 1):
                if Fit[0].SpeciesName(i) not in self.SpeciesToConsider:
                    self.SpeciesToConsider.append(Fit[0].SpeciesName(i))
                SpeciesList.append(i)
        # The single species:
        for Species in SpeciesList:
            m_final_prediction = Fit[0].Yield(Species)[-1]
            # for Standard Arrhenius
            PredictionV0 = [0.86e11, 17700, m_final_prediction]
            self.KinModel = Models.ArrheniusModelNoB(PredictionV0)
            if PyrolProgram == 'PCCL':
                self.KinModel.setDt4Intergrate(self.FG_dt)
            #
            print Fit[0].SpeciesName(Species)
            # gets final yields for all runs
            m_final_predictionAll = []
            for i in range(len(Fit)):
                m_final_predictionAll.append(Fit[i].Yield(Species)[-1])
            ParamInit = GlobalOptParam.EvAArrhInit
            if len(ParamInit) == 4:
                ParamInit.pop(-1)
            if len(ParamInit) == 3:
                ParamInit.pop(1)
            # optimization procedure
            if self.NrOfRuns == 1:
                self.OptGradBased(Fit, ParamInit, Fit[
                                  0].Yield(Species)[-1], Species)
            else:
                ParamMin = GlobalOptParam.EvAArrhMin
                ParamMax = GlobalOptParam.EvAArrhMax
                if len(ParamMin) == 3:
                    ParamMin.pop(1)
                    ParamMin.append(0.0)
                if len(ParamMax) == 3:
                    ParamMax.pop(1)
                    ParamMax.append(0.0)
                if len(ParamInit) == 3:  # if final yield was appended
                    ParamMax.pop(-1)
                if len(ParamInit) == 2:
                    ParamInit.append(0.0)
                ParamInit[2] = (max(m_final_predictionAll) +
                                min(m_final_predictionAll)) / 2.
                ParamMin[2] = (min(m_final_predictionAll))
                ParamMax[2] = (max(m_final_predictionAll))
                #
                self.OptGenAlgBased(
                    Fit, ParamInit, ParamMin, ParamMax, Species)
                #
            self.Solution = self.KinModel.ParamVector()
            self.KinModel.plot(Fit, Species)
            # To avoid, a species with no yield is added to the
            # parameter file
            if np.sum(self.KinModel.ParamVector()) != np.sum(PredictionV0):
                outfile.write(str(Fit[0].SpeciesName(Species)) + '\t' + '%.6e  %11.4f  %7.4f  ' % (
                    self.Solution[0], self.Solution[1], self.Solution[2]) + '\n')
        outfile.close()
        if oSystem == 'Linux' or oSystem == 'Darwin':
            shutil.move(PyrolProgram + '-Results_ArrheniusNoBRate.txt',
                        'Result/' + PyrolProgram + '-Results_ArrheniusNoB.txt')
        elif oSystem == 'Windows':
            shutil.move(PyrolProgram + '-Results_ArrheniusNoBRate.txt',
                        'Result\\' + PyrolProgram + '-Results_ArrheniusNoB.txt')
        else:
            print "The name of the operating system couldn't be found."
    #
    #

    def MakeResults_Kob(self, PyrolProgram, File, Fit):
        """Generates the results for Kobayashi Rate."""
        PredictionVKob0 = [7e5, 8e7 / 8314.33,
                           2.3e8, 1.6e8 / 8314.33, 0.4, 0.9]
        outfile = open(PyrolProgram + '-Results_KobayashiRate.txt', 'w')
        outfile.write(
            "Species              A1 [1/s]              E_a1 [kJ/mol]             A2 [1/s]             E_a2 [kJ/mol]      alpha1     alpha2\n\n")
        self.KinModel = Models.Kobayashi(
            GlobalOptParam.EvAKobInit)  # (PredictionVKob0)
        if PyrolProgram == 'PCCL':
            self.KinModel.setDt4Intergrate(self.FG_dt)
        #######
        # The single species:
        if 'Total' not in self.SpeciesToConsider:
            self.SpeciesToConsider.append('Total')
        for Species in [Fit[0].SpeciesIndex('Total')]:
            # optimization procedure
            print Fit[0].SpeciesName(Species)
            if self.NrOfRuns == 1:
                self.OptGradBased(
                    Fit, self.KinModel.ParamVector(), False, Species)
            else:
                self.OptGenAlgBased(Fit, self.KinModel.ParamVector(
                ), GlobalOptParam.EvAKobMin, GlobalOptParam.EvAKobMax, Species)
            #
            self.Solution = self.KinModel.ParamVector()
            #
            self.KinModel.plot(Fit, Species)
            outfile.write(str(Fit[0].SpeciesName(Species)) + '\t' + '%.6e  %11.4f  %.6e  %11.4f  %6.4f  %6.4f  ' % (self.Solution[0], self.Solution[
                          1] * 8314.33 / 1e6, self.Solution[2], self.Solution[3] * 8314.33 / 1e6, self.Solution[4], self.Solution[5]) + '\n')
        outfile.close()
        if oSystem == 'Linux' or oSystem == 'Darwin':
            shutil.move(PyrolProgram + '-Results_KobayashiRate.txt',
                        'Result/' + PyrolProgram + '-Results_Kobayashi.txt')
        elif oSystem == 'Windows':
            shutil.move(PyrolProgram + '-Results_KobayashiRate.txt',
                        'Result\\' + PyrolProgram + '-Results_Kobayashi.txt')
        else:
            print "The name of the operating system couldn't be found."
#
#

    def MakeResults_DEAM(self, PyrolProgram, File, Fit):
        """Generates the results for DAEM model."""
# elif (CPD_FittingKineticParameter_Select=='DAEM' and
# PyrolProgram=='CPD') or (FG_FittingKineticParameter_Select=='DAEM' and
# PyrolProgram=='FGDVC'):
        PredictionDAEM = [2e10, 20e3, 5e3, 0.5]
        outfile = open(PyrolProgram + '-Results_DAEM.txt', 'w')
        outfile.write(
            "Species                  A1 [1/s]               E_a1 [kJ/kmol]          sigma [kJ/kmol]   Final Yield\n\n")
        self.KinModel = Models.DAEM(PredictionDAEM)
        # self.KinModel.setNrOfActivationEnergies(GlobalOptParam.NrOFActivtionEnergies)
        if PyrolProgram == 'PCCL':
            self.KinModel.setDt4Intergrate(self.FG_dt)
        #######
        ParamInit = GlobalOptParam.EvADAEMInit
        if len(ParamInit) == 4:
            ParamInit.pop(-1)
        ParamMin = GlobalOptParam.EvADAEMMin
        if len(ParamMin) == 3:
            ParamMin.append(0.0)
        ParamMax = GlobalOptParam.EvADAEMMax
        if len(ParamMax) == 3:
            ParamMax.append(0.0)
        # The single species:
        if 'Total' not in self.SpeciesToConsider:
            self.SpeciesToConsider.append('Total')
        for Species in [Fit[0].SpeciesIndex('Total')]:
            print Fit[0].SpeciesName(Species)
            m_final_predictionAll = []
            for i in range(len(Fit)):
                m_final_predictionAll.append(Fit[i].Yield(Species)[-1])
            # optimization procedure
            if self.NrOfRuns == 1:
                self.OptGradBased(Fit, ParamInit, Fit[
                                  0].Yield(Species)[-1], Species)
            else:
                if len(ParamInit) == 3:
                    ParamInit.append(0.0)
                ParamInit[3] = (
                    (min(m_final_predictionAll) + max(m_final_predictionAll)) / 2.)
                ParamMin[3] = (min(m_final_predictionAll))
                ParamMax[3] = (max(m_final_predictionAll))
                self.OptGenAlgBased(
                    Fit, ParamInit, ParamMin, ParamMax, Species)
            #
            self.Solution = self.KinModel.ParamVector()
            #
            self.KinModel.plot(Fit, Species)
            outfile.write(str(Fit[0].SpeciesName(Species)) + '\t' + '%.6e  %11.4f  %11.4f  %6.4f  ' % (
                self.Solution[0], self.Solution[1] / 1e6, self.Solution[2] / 1e6, self.Solution[3]) + '\n')
        outfile.close()
        if oSystem == 'Linux' or oSystem == 'Darwin':
            shutil.move(PyrolProgram + '-Results_DAEM.txt',
                        'Result/' + PyrolProgram + '-Results_DAEM.txt')
        elif oSystem == 'Windows':
            shutil.move(PyrolProgram + '-Results_DAEM.txt',
                        'Result\\' + PyrolProgram + '-Results_DAEM.txt')
        else:
            print "The name of the operating system couldn't be found."
#
#

    def SpeciesEnergy(self, PyrolProgram, File, FittingModel):
        """Carries out the species and Energy balance."""
        # SPECIES AND ENERGY BALANCE:
        for runNr in range(self.NrOfRuns):
            if PyrolProgram == 'CPD':
                print 'CPD energy and mass balance...'
                Compos_and_Energy.CPD_SpeciesBalance(File[runNr], self.UAC, self.UAH, self.UAN, self.UAO, self.UAS, self.PAVM_asrec,
                                                     self.PAFC_asrec, self.PAmoist, self.PAash, self.HHV, self.MTar, self.densityDryCoal, runNr)
            if PyrolProgram == 'FGDVC':
                print 'FG-DVC energy and mass balance...'
                Compos_and_Energy.FGPC_SpeciesBalance(File[runNr], self.UAC, self.UAH, self.UAN, self.UAO, self.UAS, self.PAVM_asrec,
                                                      self.PAFC_asrec, self.PAmoist, self.PAash, self.HHV, self.MTar, self.densityDryCoal, runNr, 'FGDVC')
                #    SpecCPD=Compos_and_Energy.CPD_SpeciesBalance(File[0],UAC,UAH,UAN,UAO,PAVM_asrec,PAFC_asrec,HHV,MTar,0)
#
            if PyrolProgram == 'PCCL':
                print 'PCCL-Flashchain energy and mass balance...TO FINISH'
                Compos_and_Energy.FGPC_SpeciesBalance(File[runNr], self.UAC, self.UAH, self.UAN, self.UAO, self.UAS, self.PAVM_asrec,
                                                      self.PAFC_asrec, self.PAmoist, self.PAash, self.HHV, self.MTar, self.densityDryCoal, runNr, 'PCCL')

        # new implementation of species energy for Kobayashi model
        # Michele Vascellari
        if FittingModel == 'Kobayashi':
            print FittingModel
            self.extrapolateYieldKoba(PyrolProgram, File)

    def extrapolateYieldKoba(self, PyrolProgram, File):
        '''
        extrapolate the results of Detailed model to the alpha1/alpha2 parameters
        in the Kobayashi model
        It is required for species/energy calculation
        TODO only CPD is working at the moment!!!
        '''
        yields = []
        lenFile = len(File)
        Yields2Cols = File[0].DictYields2Cols()
        iSol = Yields2Cols['Solid']
        for i in range(lenFile):
            yields.append(File[i].FinalYields())
        lenSpecies = len(yields[0])
        if PyrolProgram == 'CPD':
            volatileYield = []
            for i in range(lenFile):
                volatileYield.append(1. - yields[i][iSol])
            volatileYield = np.array(volatileYield)
            maxVol = volatileYield.max()
            iMax = volatileYield.argmax()
            minVol = volatileYield.min()
            iMin = volatileYield.argmin()
            alpha1 = self.Solution[4]
            alpha2 = self.Solution[5]
            yields1 = yields[
                iMin] + (yields[iMax] - yields[iMin]) / (maxVol - minVol) * (alpha1 - minVol)
            # check results
            iTar = Yields2Cols['Tar']
            if (yields1[iTar] < 0):
                ySolid = yields1[iSol]
                yTar = 0.5 * (1 - ySolid)
                index = [iSol, iTar, Yields2Cols['CO'], Yields2Cols['CO2'], Yields2Cols[
                    'H2O'], Yields2Cols['CH4'], Yields2Cols['Other']]
                sumy = yields1[index].sum()
                sumy = sumy - yields1[iTar] - yields1[iSol]
                factor = 1.0 + (yields1[iTar] - yTar) / sumy
                yields1[index] = yields1[index] * factor
                yields1[iTar] = yTar
                yields1[iSol] = ySolid

            yields2 = yields[
                iMin] + (yields[iMax] - yields[iMin]) / (maxVol - minVol) * (alpha2 - minVol)
            # print yields1
            # print yields2
            # print yields[iMin]
            # print yields[iMax]

            # print yields1[iTar],yields[iMin][iTar],yields[iMax][iTar],yields2[iTar]
            # print alpha1,minVol,maxVol,alpha2
            # print yields[iMin][iTar] +
            # (yields[iMax][iTar]-yields[iMin][iTar])/(maxVol-minVol)*(alpha1-minVol)
            fit1 = CPD_Result.CPD_ResultFake(yields1)
            Compos_and_Energy.CPD_SpeciesBalance(fit1, self.UAC, self.UAH, self.UAN, self.UAO, self.UAS, self.PAVM_asrec, self.PAFC_asrec, self.PAmoist, self.PAash,
                                                 self.HHV, self.MTar, self.densityDryCoal, 'Koba1')
            fit2 = CPD_Result.CPD_ResultFake(yields2)
            Compos_and_Energy.CPD_SpeciesBalance(fit2, self.UAC, self.UAH, self.UAN, self.UAO, self.UAS, self.PAVM_asrec, self.PAFC_asrec, self.PAmoist, self.PAash,
                                                 self.HHV, self.MTar, self.densityDryCoal, 'Koba2')


#
#
####CPD#####
    def MakeResults_CPD(self):
        """generates the result for CPD"""
        CPDFile = []
        CPDFit = []
        for runNr in range(self.NrOfRuns):
            # launches CPD
            CPD = CPD_SetAndLaunch.SetterAndLauncher()
            CPD.SetCoalParameter(
                self.UAC, self.UAH, self.UAN, self.UAO, self.PAVM_daf)
            CPD.CalcCoalParam()
            if runNr == 0:
                CPD.SetOperateCond(self.CPD_pressure,
                                   self.CPD_TimeTemp1)
                CPD.SetNumericalParam(self.CPDdt, self.CPD_t_max1)
            elif runNr == 1:
                CPD.SetOperateCond(self.CPD_pressure,
                                   self.CPD_TimeTemp2)
                CPD.SetNumericalParam(self.CPDdt, self.CPD_t_max2)
            elif runNr == 2:
                CPD.SetOperateCond(self.CPD_pressure,
                                   self.CPD_TimeTemp3)
                CPD.SetNumericalParam(self.CPDdt, self.CPD_t_max3)
            elif runNr == 3:
                CPD.SetOperateCond(self.CPD_pressure,
                                   self.CPD_TimeTemp4)
                CPD.SetNumericalParam(self.CPDdt, self.CPD_t_max4)
            elif runNr == 4:
                CPD.SetOperateCond(self.CPD_pressure,
                                   self.CPD_TimeTemp5)
                CPD.SetNumericalParam(self.CPDdt, self.CPD_t_max5)
            CPD.writeInstructFile(workingDir)
            print 'Running CPD ...', runNr
            if oSystem == 'Linux':
                # first Arg: CPD-executeable, second: Input data
                # containing CPD input file and the output files
                CPD.Run('./' + 'cpdnlg', 'IN.dat',
                        'CPD_' + str(runNr) + '_output.log')
            elif oSystem == 'Darwin':
                # first Arg: CPD-executeable, second: Input data
                # containing CPD input file and the output files
                CPD.Run('./' + 'cpdnlg.x', 'IN.dat',
                        'CPD_' + str(runNr) + '_output.log')
            elif oSystem == 'Windows':
                # first Arg: CPD-executeable, second: Input data
                # containing CPD input file and the output files
                CPD.Run('cpdnlg.exe', 'IN.dat', 'CPD_' +
                        str(runNr) + '_output.log')
            else:
                print "The name of the operating system couldn't be found."
            #
            # calibration of the kinetic parameter:
            # read result:
            CurrentCPDFile = CPD_Result.CPD_Result(workingDir)
            # creates object, required for fitting procedures
            CurrentCPDFit = FitInfo.Fit_one_run(CurrentCPDFile)
            CPDFile.append(CurrentCPDFile)
            CPDFit.append(CurrentCPDFit)
            #
            if oSystem == 'Linux' or oSystem == 'Darwin':
                shutil.move('CPD_Result1.dat', 'Result/' +
                            'CPD_' + str(runNr) + '_Result1.dat')
                shutil.move('CPD_Result2.dat', 'Result/' +
                            'CPD_' + str(runNr) + '_Result2.dat')
                shutil.move('CPD_Result3.dat', 'Result/' +
                            'CPD_' + str(runNr) + '_Result3.dat')
                shutil.move('CPD_Result4.dat', 'Result/' +
                            'CPD_' + str(runNr) + '_Result4.dat')
                shutil.move('CPD_' + str(runNr) + '_output.log',
                            'Result/' + 'CPD_' + str(runNr) + '_output.log')
            elif oSystem == 'Windows':
                shutil.move('CPD_Result1.dat', 'Result\\' +
                            'CPD_' + str(runNr) + '_Result1.dat')
                shutil.move('CPD_Result2.dat', 'Result\\' +
                            'CPD_' + str(runNr) + '_Result2.dat')
                shutil.move('CPD_Result3.dat', 'Result\\' +
                            'CPD_' + str(runNr) + '_Result3.dat')
                shutil.move('CPD_Result4.dat', 'Result\\' +
                            'CPD_' + str(runNr) + '_Result4.dat')
                shutil.move('CPD_' + str(runNr) + '_output.log',
                            'Result\\' + 'CPD_' + str(runNr) + '_output.log')
            else:
                print "The name of the operating system couldn't be found."
        #####
        M = Models.Model()
        for Species in CPDFit[0].SpeciesNames():
            M.mkSimpleResultFiles(CPDFit, Species)
            if (Species not in self.SpeciesToConsider) and (Species != 'Temp') and (Species != 'Time'):
                self.SpeciesToConsider.append(Species)
        if self.CPD_FittingKineticParameter_Select == 'constantRate':
            self.MakeResults_CR('CPD', CPDFile, CPDFit)
            currentDict = {'CPD': 'constantRate'}
        elif self.CPD_FittingKineticParameter_Select == 'Arrhenius':
            self.MakeResults_Arrh('CPD', CPDFile, CPDFit)
            currentDict = {'CPD': 'Arrhenius'}
        elif self.CPD_FittingKineticParameter_Select == 'ArrheniusNoB':
            self.MakeResults_ArrhNoB('CPD', CPDFile, CPDFit)
            currentDict = {'CPD': 'ArrheniusNoB'}
        elif self.CPD_FittingKineticParameter_Select == 'Kobayashi':
            self.MakeResults_Kob('CPD', CPDFile, CPDFit)
            currentDict = {'CPD': 'Kobayashi'}
        elif self.CPD_FittingKineticParameter_Select == 'DAEM':
            self.MakeResults_DEAM('CPD', CPDFile, CPDFit)
            currentDict = {'CPD': 'DAEM'}
        elif self.CPD_FittingKineticParameter_Select == None:
            currentDict = {'CPD': 'None'}
        else:
            print 'unspecified CPD_FittingKineticParameter_Select'
            currentDict = {}
        #
        self.ProgramModelDict.update(currentDict)
        #
        self.SpeciesEnergy(
            'CPD', CPDFile, self.CPD_FittingKineticParameter_Select)
        #
        #
    ####FG-DVC####

    def MakeResults_FG(self):
        """generates the result for FG-DVC"""
        # writes Time-Temperature file
        FG_TimeTemp1 = np.zeros(np.shape(self.CPD_TimeTemp1), order='F')
        FG_TimeTemp2 = np.zeros(np.shape(self.CPD_TimeTemp2), order='F')
        FG_TimeTemp3 = np.zeros(np.shape(self.CPD_TimeTemp3), order='F')
        FG_TimeTemp4 = np.zeros(np.shape(self.CPD_TimeTemp4), order='F')
        FG_TimeTemp5 = np.zeros(np.shape(self.CPD_TimeTemp5), order='F')
        FG_TimeTemp1[:, 0] = self.CPD_TimeTemp1[:, 0] * 1.e-3
        FG_TimeTemp2[:, 0] = self.CPD_TimeTemp2[:, 0] * 1.e-3
        FG_TimeTemp3[:, 0] = self.CPD_TimeTemp3[:, 0] * 1.e-3
        FG_TimeTemp4[:, 0] = self.CPD_TimeTemp4[:, 0] * 1.e-3
        FG_TimeTemp5[:, 0] = self.CPD_TimeTemp5[:, 0] * 1.e-3
        FG_TimeTemp1[:, 1] = self.CPD_TimeTemp1[:, 1]
        FG_TimeTemp2[:, 1] = self.CPD_TimeTemp2[:, 1]
        FG_TimeTemp3[:, 1] = self.CPD_TimeTemp3[:, 1]
        FG_TimeTemp4[:, 1] = self.CPD_TimeTemp4[:, 1]
        FG_TimeTemp5[:, 1] = self.CPD_TimeTemp5[:, 1]
        # initialize the launching object
        FGDVC = FGDVC_SetAndLaunch.SetterAndLauncher()
        # set and writes Coal Files:
        if self.FG_CoalSelection == 0:
            # deletes old generated file
            os.system('cd ' + self.FG_MainDir + FG_GenCoalDir + ' & del ' + FG_CoalName +
                      '_com.dat, ' + FG_CoalName + '_kin.dat, ' + FG_CoalName + '_pol.dat')
            # generates coalsd.exe input file
            MakeCoalGenFile = InformationFiles.WriteFGDVCCoalFile(
                FG_CoalGenFileName)
            MakeCoalGenFile.setCoalComp(
                self.UAC, self.UAH, self.UAO, self.UAN, self.UAS, 0)
            MakeCoalGenFile.write(
                self.FG_MainDir + FG_GenCoalDir + '\\', FG_CoalName, option=0)
            # makes new file
            try:
                os.system('cd ' + self.FG_MainDir + FG_GenCoalDir + ' & ' +
                          'coalsd.exe < ' + FG_CoalGenFileName + ' > coalsd_pkp.log')
            except OSError:
                print 'Problems with coalsd.exe'
            os.system('copy ' + self.FG_MainDir +
                      FG_GenCoalDir + '\coalsd_pkp.log . >> log.txt')
            # tests weather the coal file was genearated:
            if os.path.exists(self.FG_MainDir + '\\' + FG_GenCoalDir + '\\' + FG_CoalName + '_com.dat') == False:
                print 30 * '*', '\n', 'The coal is may outside the libraries coals. Select manually the closest library coal.', 30 * '*', '\n'
                MakeCoalGenFile.write(
                    self.FG_MainDir + FG_GenCoalDir + '\\', FG_CoalName, option=10)
                os.system('cd ' + self.FG_MainDir + FG_GenCoalDir + ' & ' +
                          'coalsd.exe < ' + FG_CoalGenFileName + ' > coalsd_pkp.log')
            # sets generated file for instruct.ini
            FGDVC.set1CoalLocation(
                self.FG_MainDir + FG_GenCoalDir + '\\' + FG_CoalName + '_com.dat')
            FGDVC.set2KinLocation(
                self.FG_MainDir + FG_GenCoalDir + '\\' + FG_CoalName + '_kin.dat')
            FGDVC.set3PolyLocation(
                self.FG_MainDir + FG_GenCoalDir + '\\' + FG_CoalName + '_pol.dat')
        elif self.FG_CoalSelection > 0 and self.FG_CoalSelection < 9:
            # sets library file for instruct.ini
            FGDVC.set1CoalLocation(
                self.FG_MainDir + FG_LibCoalDir + '\\coal.ar' + str(self.FG_CoalSelection))
            FGDVC.set2KinLocation(
                self.FG_MainDir + FG_LibCoalDir + '\\kin.ar' + str(self.FG_CoalSelection))
            FGDVC.set3PolyLocation(
                self.FG_MainDir + FG_LibCoalDir + '\\polymr.ar' + str(self.FG_CoalSelection))
        else:
            print "select Choose Coal: 0 interpolate between library coals and generate own coal. Set 1 to 8 for a library coal.' in FGDVC.inp equal a value between 0 and 8"
        # sets FG-DVC instruct.ini parameter
        FGDVC.set5Pressure(self.FG_pressure)
        if self.FG_TarCacking == 0.0:  # case: no tar cracking
            FGDVC.set6Theorie(13, 0.0)
        elif self.FG_TarCacking < 0.0:  # case: full tar cracking
            FGDVC.set6Theorie(15, 0.0)
        else:  # case: partial tar cracking
            FGDVC.set6Theorie(13, float(self.FG_TarCacking))
        #
        FGFile = []
        FGFit = []
        OpCondInp = InformationFiles.OperCondInput('OperCond.inp')
        for runNr in range(self.NrOfRuns):
            if runNr == 0:
                OpCondInp.writeFGDVCtTHist(
                    FG_TimeTemp1, self.FG_dt, self.FG_T_t_History)
            elif runNr == 1:
                OpCondInp.writeFGDVCtTHist(
                    FG_TimeTemp2, self.FG_dt, self.FG_T_t_History)
            elif runNr == 2:
                OpCondInp.writeFGDVCtTHist(
                    FG_TimeTemp3, self.FG_dt, self.FG_T_t_History)
            elif runNr == 3:
                OpCondInp.writeFGDVCtTHist(
                    FG_TimeTemp4, self.FG_dt, self.FG_T_t_History)
            elif runNr == 4:
                OpCondInp.writeFGDVCtTHist(
                    FG_TimeTemp5, self.FG_dt, self.FG_T_t_History)
            FGDVC.set7File(self.FG_T_t_History)
            FGDVC.set9AshMoisture(0.0, 0.0)
            # case: models temperature history with the file
            FGDVC.setTRamp_or_TFile('File')
            # writes the instruct.ini and launches FG-DVC (no graphical
            # user interface, only main file fgdvcd.exe)
            FGDVC.writeInstructFile(
                self.FG_MainDir + '\\' + FG_ExeCoalDir + '\\')
            FGDVC.Run('cd ' + self.FG_MainDir +
                      FG_ExeCoalDir + ' & ' + 'fgdvcd.exe')
            #
            # calibrate kinetic parameter:
            # read result:
            CurrentFGFile = FGDVC_Result.FGDVC_Result(self.FG_DirOut)
            # creates object, required for fitting procedures
            CurrentFGFit = FitInfo.Fit_one_run(CurrentFGFile)
            FGFile.append(CurrentFGFile)
            FGFit.append(CurrentFGFit)
            # copies file, keeping the name:
            if oSystem == 'Linux' or oSystem == 'Darwin':
                shutil.copyfile(self.FG_DirOut + 'gasyield.txt',
                                'Result/gasyield_' + str(runNr) + '.txt')
                shutil.copyfile(self.FG_DirOut + 'gasrate.txt',
                                'Result/gasrate_' + str(runNr) + '.txt')
            elif oSystem == 'Windows':
                shutil.copyfile(self.FG_DirOut + 'gasyield.txt',
                                'Result\\gasyield_' + str(runNr) + '.txt')
                shutil.copyfile(self.FG_DirOut + 'gasrate.txt',
                                'Result\\gasrate_' + str(runNr) + '.txt')
        #####
        M = Models.Model()
        for Species in FGFit[0].SpeciesNames():
            M.mkSimpleResultFiles(FGFit, Species)
            if (Species not in self.SpeciesToConsider) and (Species != 'Temp') and (Species != 'Time'):
                self.SpeciesToConsider.append(Species)
        if self.FG_FittingKineticParameter_Select == 'constantRate':
            self.MakeResults_CR('FGDVC', FGFile, FGFit)
            currentDict = {'FGDVC': 'constantRate'}
        elif self.FG_FittingKineticParameter_Select == 'Arrhenius':
            self.MakeResults_Arrh('FGDVC', FGFile, FGFit)
            currentDict = {'FGDVC': 'Arrhenius'}
        elif self.FG_FittingKineticParameter_Select == 'ArrheniusNoB':
            self.MakeResults_ArrhNoB('FGDVC', FGFile, FGFit)
            currentDict = {'FGDVC': 'ArrheniusNoB'}
        elif self.FG_FittingKineticParameter_Select == 'Kobayashi':
            self.MakeResults_Kob('FGDVC', FGFile, FGFit)
            currentDict = {'FGDVC': 'Kobayashi'}
        elif self.FG_FittingKineticParameter_Select == 'DAEM':
            self.MakeResults_DEAM('FGDVC', FGFile, FGFit)
            currentDict = {'FGDVC': 'DAEM'}
        elif self.FG_FittingKineticParameter_Select == None:
            currentDict = {'FGDVC': 'None'}
        else:
            print 'uspecified FG_FittingKineticParameter_Select'
            currentDict = {}
        #
        self.ProgramModelDict.update(currentDict)
        #
        self.SpeciesEnergy(
            'FGDVC', FGFile, self.FG_FittingKineticParameter_Select)
        #

    ####Pc Coal Lab####
    def MakeResults_PCCL(self):
        """generates the result for PC Coal Lab"""
        # writes Time-Temperature file
        PCCL_TimeTemp1 = np.zeros(
            np.shape(self.CPD_TimeTemp1), order='F')
        PCCL_TimeTemp2 = np.zeros(
            np.shape(self.CPD_TimeTemp2), order='F')
        PCCL_TimeTemp3 = np.zeros(
            np.shape(self.CPD_TimeTemp3), order='F')
        PCCL_TimeTemp4 = np.zeros(
            np.shape(self.CPD_TimeTemp4), order='F')
        PCCL_TimeTemp5 = np.zeros(
            np.shape(self.CPD_TimeTemp5), order='F')
        PCCL_TimeTemp1[:, 0] = self.CPD_TimeTemp1[:, 0] * 1.e-3
        PCCL_TimeTemp2[:, 0] = self.CPD_TimeTemp2[:, 0] * 1.e-3
        PCCL_TimeTemp3[:, 0] = self.CPD_TimeTemp3[:, 0] * 1.e-3
        PCCL_TimeTemp4[:, 0] = self.CPD_TimeTemp4[:, 0] * 1.e-3
        PCCL_TimeTemp5[:, 0] = self.CPD_TimeTemp5[:, 0] * 1.e-3
        PCCL_TimeTemp1[:, 1] = self.CPD_TimeTemp1[:, 1]
        PCCL_TimeTemp2[:, 1] = self.CPD_TimeTemp2[:, 1]
        PCCL_TimeTemp3[:, 1] = self.CPD_TimeTemp3[:, 1]
        PCCL_TimeTemp4[:, 1] = self.CPD_TimeTemp4[:, 1]
        PCCL_TimeTemp5[:, 1] = self.CPD_TimeTemp5[:, 1]
        # initialize the launching object
        PCCL = PCCL_SetAndLaunch.SetterAndLauncher()
        # set and writes Coal Files:
        PCCL.SetUACoalParameter(
            self.UAC, self.UAH, self.UAN, self.UAO, self.UAS)
        PCCL.SetPACoalParameter(
            self.PAVM_asrec, self.PAFC_asrec, self.PAmoist, self.PAash)
        if type(self.PCCL_CoalCalFactor) == float:
            PCCL.SetCoalCalibrationFactor(self.PCCL_CoalCalFactor)
        PCCL.SetPressure(self.FG_pressure)
        PCCL.SetParticleSize(self.PCCL_ParticleSize)
        #
        PCCL.writeCoalFiles(self.PCCL_Path)
        #
        PCCL.THist(PCCL_TimeTemp1[0, 1], PCCL_TimeTemp1[
                   1, 0], PCCL_TimeTemp1[-1, 1], PCCL_TimeTemp1[-1, 0])
        PCCL.writeInstructFiles(self.PCCL_Path, mkNewFile=True)
        PCCL.THist(PCCL_TimeTemp2[0, 1], PCCL_TimeTemp2[
                   1, 0], PCCL_TimeTemp2[-1, 1], PCCL_TimeTemp2[-1, 0])
        PCCL.writeInstructFiles(self.PCCL_Path, mkNewFile=False)
        PCCL.THist(PCCL_TimeTemp3[0, 1], PCCL_TimeTemp3[
                   1, 0], PCCL_TimeTemp3[-1, 1], PCCL_TimeTemp3[-1, 0])
        PCCL.writeInstructFiles(self.PCCL_Path, mkNewFile=False)
        PCCL.THist(PCCL_TimeTemp4[0, 1], PCCL_TimeTemp4[
                   1, 0], PCCL_TimeTemp4[-1, 1], PCCL_TimeTemp4[-1, 0])
        PCCL.writeInstructFiles(self.PCCL_Path, mkNewFile=False)
        PCCL.THist(PCCL_TimeTemp5[0, 1], PCCL_TimeTemp5[
                   1, 0], PCCL_TimeTemp5[-1, 1], PCCL_TimeTemp5[-1, 0])
        PCCL.writeInstructFiles(self.PCCL_Path, mkNewFile=False)
        PCCL.writeInstructFilesFinish()
        #
        PCCL.Run(self.PCCL_Path, self.PCCL_Exe)
        #
        PCCLFile = []
        PCCLFit = []
        for runNr in range(1, self.NrOfRuns + 1, 1):
            # read result:
            CurrentPCCLFile = PCCL_Result.PCCL_Result(
                self.PCCL_Path, runNr)
            # creates object, required for fitting procedures
            CurrentPCCLFit = FitInfo.Fit_one_run(CurrentPCCLFile)
            PCCLFile.append(CurrentPCCLFile)
            PCCLFit.append(CurrentPCCLFit)
            # copies file:
            if oSystem == 'Windows':
                shutil.copyfile(self.PCCL_Path + 'FDC1WT' + str(runNr) +
                                '.RPT', 'Result/PCCL_gasyield_wt_' + str(runNr) + '.txt')
                shutil.copyfile(self.PCCL_Path + 'FDC1NG' + str(runNr) +
                                '.RPT', 'Result/PCCL_gasyield_ng_' + str(runNr) + '.txt')
                shutil.copyfile(self.PCCL_Path + 'FDC1HC' + str(runNr) +
                                '.RPT', 'Result/PCCL_gasyield_hc_' + str(runNr) + '.txt')
        #####
        M = Models.Model()
        for Species in PCCLFit[0].SpeciesNames():
            M.mkSimpleResultFiles(PCCLFit, Species)
            if (Species not in self.SpeciesToConsider) and (Species != 'Temp') and (Species != 'Time'):
                self.SpeciesToConsider.append(Species)
        if self.PCCL_FittingKineticParameter_Select == 'constantRate':
            self.MakeResults_CR('PCCL', PCCLFile, PCCLFit)
            currentDict = {'PCCL': 'constantRate'}
        elif self.PCCL_FittingKineticParameter_Select == 'Arrhenius':
            self.MakeResults_Arrh('PCCL', PCCLFile, PCCLFit)
            currentDict = {'PCCL': 'Arrhenius'}
        elif self.PCCL_FittingKineticParameter_Select == 'ArrheniusNoB':
            self.MakeResults_ArrhNoB('PCCL', PCCLFile, PCCLFit)
            currentDict = {'PCCL': 'ArrheniusNoB'}
        elif self.PCCL_FittingKineticParameter_Select == 'Kobayashi':
            self.MakeResults_Kob('PCCL', PCCLFile, PCCLFit)
            currentDict = {'PCCL': 'Kobayashi'}
        elif self.PCCL_FittingKineticParameter_Select == 'DAEM':
            self.MakeResults_DEAM('PCCL', PCCLFile, PCCLFit)
            currentDict = {'PCCL': 'DAEM'}
        elif self.PCCL_FittingKineticParameter_Select == None:
            currentDict = {'PCCL': 'None'}
        else:
            print 'uspecified PCCL_FittingKineticParameter_Select'
            currentDict = {}
        #
        self.ProgramModelDict.update(currentDict)
        #
        self.SpeciesEnergy('PCCL', PCCLFile)
        #

    def RunPMSKD(self):
        '''
        run PMSKD
        '''
        # create object
        import pkp.coalPolimi
        try:
            coal = pkp.coalPolimi.coalPolimi(name='COAL', c=self.UAC, h=self.UAH, o=self.UAO, n=self.UAN, s=self.UAS,
                                             file=self.PMSKD_mechfile)
        except pkp.coalPolimi.compositionError:
            print 'Composition outside of triangle of definition'
            sys.exit()
        # organize TimeTemp
            PMSKDFile = []
        PMSKDFit = []
        for runNr in range(self.NrOfRuns):
            print 'Running PMSKD n. ' + str(runNr)
            # print self.timeHR[runNr]
            # print self.temperatureHR[runNr]
            # set heating rate
            coal.setHeatingRate(
                self.timeHR[runNr], self.temperatureHR[runNr])
            # coal.setTimeStep(self.PMSKD_npoint)
            coal.solvePyrolysis(
                filename='Result/PMSKD_{}.dat'.format(runNr))
            # plt.figure(runNr)
            # plt.plot(coal.getTemperature(),coal.getVolatile())
            # read result:
            # CurrentPMSKDFile=Coal
            # creates object, required for fitting procedures
            CurrentPMSKDFit = FitInfo.Fit_one_run(coal)
            # PMSKDFile.append(CurrentFGFile)
            PMSKDFit.append(CurrentPMSKDFit)
            # print coal.Yields_all()

            coal.reset()
            # print coal.timeHR
            # print coal.temperatureHR

        if self.PMSKD_FittingKineticParameter_Select == 'constantRate':
            self.MakeResults_CR('PMSKD', '', PMSKDFit)
            currentDict = {'PMSKD': 'constantRate'}
        elif self.PMSKD_FittingKineticParameter_Select == 'Arrhenius':
            self.MakeResults_Arrh('PMSKD', '', PMSKDFit)
            currentDict = {'PMSKD': 'Arrhenius'}
        elif self.PMSKD_FittingKineticParameter_Select == 'ArrheniusNoB':
            self.MakeResults_ArrhNoB('PMSKD', '', PMSKDFit)
            currentDict = {'PMSKD': 'ArrheniusNoB'}
        elif self.PMSKD_FittingKineticParameter_Select == 'Kobayashi':
            self.MakeResults_Kob('PMSKD', '', PMSKDFit)
            currentDict = {'PMSKD': 'Kobayashi'}
        elif self.PMSKD_FittingKineticParameter_Select == 'DAEM':
            self.MakeResults_DEAM('PMSKD', '', PMSKDFit)
            currentDict = {'PMSKD': 'DAEM'}
        elif self.PMSKD_FittingKineticParameter_Select == None:
            currentDict = {'PMSKD': 'None'}
            for Species in PMSKDFit[0].SpeciesNames():
                M = Models.Model()
                M.mkSimpleResultFiles(PMSKDFit, Species)
                if (Species not in self.SpeciesToConsider) and (Species != 'Temp') and (Species != 'Time'):
                    self.SpeciesToConsider.append(Species)
        else:
            print 'undefined PMSKD_FittingKineticParameter_Select'
            currentDict = {}
            #
        self.ProgramModelDict.update(currentDict)
        #
        # self.SpeciesEnergy('PMSKD',FGFile)

    def RunBioPolimi(self):
        '''
        run BioPolimi
        '''
        # create object
        print("Run BioPolimi...")
        import pkp.bioPolimi
        try:
            biomass = pkp.bioPolimi.bioPolimi(
                name='biomass', c=self.UAC, h=self.UAH, o=self.UAO, n=self.UAN, s=self.UAS, file=self.bio_dict['mechanism'])
        except pkp.bioPolimi.compositionError:
            print 'Composition outside of triangle of definition'
            sys.exit()

        # organize TimeTemp
        bioPolimiFile = []
        bioPolimiFit = []
        for runNr in range(self.NrOfRuns):
            print 'Running BioPolimi n. ' + str(runNr)
            biomass.setHeatingRate(
                self.timeHR[runNr], self.temperatureHR[runNr])
            biomass.solvePyrolysis(self.bio_dict['nPoints'])
            # print biomass.Yields_all()
            bioPolimiFit.append(FitInfo.Fit_one_run(biomass))
            biomass.reset()

        if self.bio_dict['FittingKineticParameter_Select'] == 'constantRate':
            self.MakeResults_CR('bioPolimi', '', bioPolimiFit)
            currentDict = {'bioPolimi': 'constantRate'}
        elif self.bio_dict['FittingKineticParameter_Select'] == 'Arrhenius':
            self.MakeResults_Arrh('bioPolimi', '', bioPolimiFit)
            currentDict = {'PMSKD': 'Arrhenius'}
        elif self.bio_dict['FittingKineticParameter_Select'] == 'ArrheniusNoB':
            self.MakeResults_ArrhNoB('bioPolimi', '', bioPolimiFit)
            currentDict = {'bioPolimi': 'ArrheniusNoB'}
        elif self.bio_dict['FittingKineticParameter_Select'] == 'Kobayashi':
            self.MakeResults_Kob('bioPolimi', '', bioPolimiFit)
            currentDict = {'bioPolimi': 'Kobayashi'}
        elif self.bio_dict['FittingKineticParameter_Select'] == 'DAEM':
            self.MakeResults_DEAM('bioPolimi', '', bioPolimiFit)
            currentDict = {'bioPolimi': 'DAEM'}
        elif self.bio_dict['FittingKineticParameter_Select'] == None:
            currentDict = {'bioPolimi': 'None'}
            for Species in bioPolimiFit[0].SpeciesNames():
                M = Models.Model()
                if Species != 'Temp' and Species != 'Time':
                    M.mkSimpleResultFiles(bioPolimiFit, Species)
                if (Species not in self.SpeciesToConsider) and (Species != 'Temp') and (Species != 'Time'):
                    self.SpeciesToConsider.append(Species)
        else:
            print 'undefined PMSKD_FittingKineticParameter_Select'
            currentDict = {}
            #
        self.ProgramModelDict.update(currentDict)
        #
        # self.SpeciesEnergy('PMSKD',FGFile)

# Main Part starting
if __name__ == "__main__":
    logger = logging.getLogger('main')
    parser = argparse.ArgumentParser(
        prog='PKP',
        description=('Python Pyrolysis Preprocessor PKP\n'
                     'Command line'))
    parser.add_argument('yml_file', action='store',
                        default=None, help='Input file')
    parser.add_argument('-d', action='store_true',
                        dest='debug', help='Pront debug messages')
    argument = parser.parse_args()
    if argument.debug:
        logger.setLevel(logging.DEBUG)
        logger.debug('Activate DEBUG messages')
    else:
        logger.setLevel(logging.INFO)
    yml_file = argument.yml_file
    logger.debug('Input file %s', yml_file)
    Case = MainProcess()
    Case.ReadInputFiles(yml_file=yml_file)
    if Case.CPDselect == True:
        Case.MakeResults_CPD()
    if Case.FG_select == True:
        Case.CheckFGdt()
        Case.MakeResults_FG()
    if Case.PMSKD_select == True:
        Case.RunPMSKD()
    if Case.PCCL_select == True:
        Case.MakeResults_PCCL()
    if Case.bio_dict['Use'] == True:
        Case.RunBioPolimi()
    print 'calculated Species: ', Case.SpeciesToConsider
