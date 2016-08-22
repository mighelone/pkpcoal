import shutil
import platform
import numpy
from numpy.linalg import solve
import numpy as np
from tabulate import tabulate
oSystem = platform.system()
# PARAMETER:
# Molecular weights [Formeln und Tabellen, Peatec, 8.Auflage] (g/mol):
MH = 2.01588 * 0.5
MC = 12.011
MN = 28.0134 * 0.5
MO = 31.9988 * 0.5
MS = 32.06
MRaw = 100.
# Latent Heat of Water in J/kg :
rH2O = 2263073
# Enthalpie of Formation in J/mol:
hfCO2 = -3.935428e8
hfSO2 = -2.968612e+08
hfH2O = -2.418428e8
hfCH4 = -7.4829298e7
hfCO = -1.105397e8
hfO2 = -847.6404
hfN2 = 1429.881
hfH2 = 2448.595
hfChar = -101.268
rhoH2O = 1000.0

# tablefmt = 'rst'
tablefmt = None
######################################


class SpeciesBalance(object):
    """
    This is the parent Class for the CPD and FG-DVC specific
    Species Balances, containing general methods like the Dulong
    formula.
    """

    def SpeciesIndex(self, species):
        """Returns the column number of the input species."""
        return self.Yields2Cols[species]

    def moistureVolumeFraction(self):
        """
        calculate volume fraction of moisture
        """
        volMoist = self.PAmoist / rhoH2O
        volDry = (1. - self.PAmoist) / self.densityDryCoal
        return volMoist / (volMoist + volDry)

        return self.densityDryCoal

    def Dulong(self):
        """
        Uses the Dulong formular to calculate the Higher heating value.
        The output is in J/kg.
        """
        # HHV = 32.79 MJ/kg fC + 150.4 (fH - fO/8) + 9.26 fS + 4.97 fO +
        # 2.42 fN
        HHV = 32.79 * self.UAC + 150.4 * \
            (self.UAH - self.UAO / 8) + 9.26 * \
            0.0 + 4.97 * self.UAO + 2.42 * self.UAN
        return HHV * 1.0e6

    def correctUA(self):
        """Scale Ultimate Analysis to have sum=1"""
        # UAC(new)=UAC+(1-UAH-UAN-UAO-UAC)
        # self.UAC=1.-self.UAH-self.UAN-self.UAO
        sumUA = self.UAC + self.UAH + self.UAO + self.UAN + self.UAS
        self.UAC = self.UAC / sumUA
        self.UAH = self.UAH / sumUA
        self.UAO = self.UAO / sumUA
        self.UAN = self.UAN / sumUA
        self.UAS = self.UAS / sumUA
        sumUA = self.UAC + self.UAH + self.UAO + self.UAN + self.UAS

# ##############CPD specific class#############################


class CPD_SpeciesBalance(SpeciesBalance):
    """
    This class calculates the Species and the Energy balance for CPD.
    See the manual for the formulas and more details.
    """

    def __init__(self, CPD_ResultObject, UAC, UAH, UAN, UAO, UAS, PAVM,
                 PAFC, PAmoist, PAash, HHV, MTar, densityDryCoal,
                 RunNr):
        # imports the dictionaries for the species defined in
        # 'CPD_Fit_one_run.py' -> CPD_Result()
        self.Yields2Cols = CPD_ResultObject.DictYields2Cols()
        self.Cols2Yields = CPD_ResultObject.DictCols2Yields()
        self.densityDryCoal = densityDryCoal
        # imports final yield composition defined in CPD_Fit_one_run.py'
        # -> CPD_Result()
        self.Yields = CPD_ResultObject.FinalYields()
        # 0:'Time', 1:'Temp', 2:'Tar', 3:'Gas', 4:'Solid', 5:'Total',
        # 6:'H2O', 7:'CO2', 8:'CH4', 9:'CO', 10:'Other'
        self.UAC = UAC / 100.
        self.UAH = UAH / 100.
        self.UAN = UAN / 100.
        self.UAO = UAO / 100.
        self.UAS = UAS / 100.
        self.PAVM = PAVM / 100.
        self.PAFC = PAFC / 100.
        self.PAash = PAash / 100.
        self.PAmoist = PAmoist / 100.
        if HHV == 0:
            HHV = self.Dulong()
        self.HHV = HHV
        self.MTar = MTar
        # corrects UA, if UA<1: Sulfur -> Carbon
        self.correctUA()
        self.CPDBalanceFile = open(
            'CPD-BalanceResults' + str(RunNr) + '.txt', 'w')
        self.CPDBalanceFile.write(
            'Coal Properties\n'
            '===============\n'
            'Ultimate Analysis\n'
            '-----------------\n\n')
        self.CPDBalanceFile.write(
            tabulate(
                [[sp, '{:4.3%}'.format(val)] for sp, val in zip(
                    ['C', 'H', 'O', 'N', 'S'],
                    [self.UAC, self.UAH, self.UAO, self.UAN,
                     self.UAS])],
                tablefmt=tablefmt)
        )
        self._addlines(2)
        self.CPDBalanceFile.write(
            'Proximate analysis\n'
            '------------------\n\n')
        PA = [v / 100. for v in [PAVM, PAFC, PAash, PAmoist]]
        PAdry = [PAi / (1. - PAmoist / 100.) for PAi in PA]
        PAdry[3] = 0
        PAdaf = [PAi / (1. - PAmoist / 100. - PAash / 100.)
                 for PAi in PA]
        PAdaf[2] = PAdaf[3] = 0
        self.CPDBalanceFile.write(
            tabulate(
                [
                    [sp, '{:4.3%}'.format(pa), '{:4.3%}'.format(
                        padry), '{:4.3%}'.format(padaf)]
                    for sp, pa, padry, padaf in
                    zip(['VM', 'FC', 'Ash', 'Moist.'], PA, PAdry, PAdaf)
                ],
                headers=['', 'AR%', 'dry%', 'daf%'],
                tablefmt=tablefmt
            )
        )
        self._addlines(2)
        self.CPDBalanceFile.write('\n')

        # if sum Yields != 1.0: scales every Yield up
        self.__correctYields()
        self.__Q_React()  # calculate heat of reaction for raw coal
        self.__CheckOxygen()
        self.__CheckOthers()
        self.__TarComp()
        self.__hfRaw()
        self.__hfTar()
        self.__QPyro()
        self.write_tar_as_benzene()
        self.__closeResultFile()
        if oSystem == 'Linux' or oSystem == 'Darwin':
            shutil.move('CPD-BalanceResults' + str(RunNr) + '.txt',
                        'Result/' + 'CPD-BalanceResults' + str(RunNr) + '.txt')
        elif oSystem == 'Windows':
            shutil.move('CPD-BalanceResults' + str(RunNr) + '.txt',
                        'Result\\' + 'CPD-BalanceResults' + str(RunNr) + '.txt')
        else:
            print "The name of the operating system couldn't be found."

    def __correctYields(self):
        """
        Correct the amount of the yields 'Other'.
        """
        SumYieldsWitoutN = (self.Yields[self.SpeciesIndex('Solid')] +
                            self.Yields[self.SpeciesIndex('Tar')] +
                            self.Yields[self.SpeciesIndex('CO')] +
                            self.Yields[self.SpeciesIndex('CO2')] +
                            self.Yields[self.SpeciesIndex('H2O')] +
                            self.Yields[self.SpeciesIndex('CH4')] +
                            self.Yields[self.SpeciesIndex('Other')])
        self.Yields[self.SpeciesIndex('Solid')] = self.Yields[
            self.SpeciesIndex('Solid')] / SumYieldsWitoutN
        self.Yields[self.SpeciesIndex('Tar')] = self.Yields[
            self.SpeciesIndex('Tar')] / SumYieldsWitoutN
        self.Yields[self.SpeciesIndex('CO')] = self.Yields[
            self.SpeciesIndex('CO')] / SumYieldsWitoutN
        self.Yields[self.SpeciesIndex('CO2')] = self.Yields[
            self.SpeciesIndex('CO2')] / SumYieldsWitoutN
        self.Yields[self.SpeciesIndex('H2O')] = self.Yields[
            self.SpeciesIndex('H2O')] / SumYieldsWitoutN
        self.Yields[self.SpeciesIndex('CH4')] = self.Yields[
            self.SpeciesIndex('CH4')] / SumYieldsWitoutN
        self.Yields[self.SpeciesIndex('Other')] = self.Yields[
            self.SpeciesIndex('Other')] / SumYieldsWitoutN

    def write_tar_as_benzene(self):
        composition = {sp: self.Yields[idx]
                       for sp, idx in self.Yields2Cols.iteritems()
                       if sp not in ('Time', 'Total', 'Temp', 'Gas')}
        ua = {el: val
              for el, val in zip(['C', 'H', 'O', 'N', 'S'],
                                 [self.UAC, self.UAH, self.UAO,
                                  self.UAN, self.UAS])}
        vol = VolatileComposition(
            composition=composition, ultimate_analysis=ua)
        vol.define_composition_empirical(CO=0.1)
        self.CPDBalanceFile.write(
            '\n'
            'Volatile Composition (C6H6 as TAR)\n'
            '==================================\n\n'
            'Use with Creck mechanism\n\n'
        )
        self.CPDBalanceFile.write(
            tabulate([[sp, val, vol.volatiles[sp]]
                      for sp, val in vol.composition.iteritems()],
                     headers=['Species', '% Daf', '%Vol'])
        )
        self._addlines(2)
        self.CPDBalanceFile.write(
            'Heat of pyrolysis: {:5.3f} MJ/kg'.format(
                vol.heat_of_pyrolysis(heat_coal=self.LHVdaf) / 1e6))

    def __CheckOxygen(self):
        """
        Checks weather the amount of Oxygen in the light gases
        is lower than in the Ultimate Analysis. If not, the amount
        of these species decreases. If yes, the tar contains oxygen.
        """
        # eq. (1) in Michele's report
        Sum1 = (1. * self.Yields[self.SpeciesIndex('H2O')]) /\
               (MO + 2 * MH) + (2. * self.Yields[self.SpeciesIndex('CO2')]) /\
               (MC + 2 * MO) + (1. * self.Yields[self.SpeciesIndex('CO')]) /\
               (MC + MO)
        Y0_cpd = MO * Sum1
        # gamma in eq. (2) in Michele's report
        self.gamma = (self.UAO / Y0_cpd)
        #
        if self.gamma < 1.:  # if gamma < 1 then is Y0_cpd > UAO
            # eq. (4) in Michele's report:
            Sum4 = self.Yields[self.SpeciesIndex(
                'H2O')] + self.Yields[self.SpeciesIndex('CO2')] + self.Yields[self.SpeciesIndex('CO')]
            self.Yields[self.SpeciesIndex('Other')] = self.Yields[
                self.SpeciesIndex('Other')] + (1 - self.gamma) * Sum4
            # eq. (3) in Michele's report:
            self.Yields[self.SpeciesIndex(
                'H2O')] = self.gamma * self.Yields[self.SpeciesIndex('H2O')]
            self.Yields[self.SpeciesIndex(
                'CO2')] = self.gamma * self.Yields[self.SpeciesIndex('CO2')]
            self.Yields[self.SpeciesIndex(
                'CO')] = self.gamma * self.Yields[self.SpeciesIndex('CO')]

    def __CheckOthers(self):
        """
        If the yield of nitrogen (is equal the UA of nitrogen) is lower 
        than the species 'Other' the difference is set equal Methane."""
        if self.UAN < self.Yields[self.SpeciesIndex('Other')]:
            # eq. (6) in Michele's report:
            self.Yields[self.SpeciesIndex('CH4')] = self.Yields[self.SpeciesIndex(
                'CH4')] + (self.Yields[self.SpeciesIndex('Other')] - self.UAN)
        self.__CheckHydrogen()
        # Result File:
        self.CPDBalanceFile.write(
            'Final yields\n'
            '------------\n\n')
        self.CPDBalanceFile.write(
            tabulate(
                [[sp, '{:7.5f}'.format(self.Yields[self.SpeciesIndex(sp_ind)])]
                 for sp, sp_ind in zip(
                    ['Char', 'Tar', 'H2O', 'CO2', 'CH4', 'CO'],
                    ['Solid', 'Tar', 'H2O', 'CO2', 'CH4', 'CO'])] +
                [['N2', '{:7.5f}'.format(self.UAN)]],
                tablefmt=tablefmt
            ))
        self._addlines(2)
        self.CPDBalanceFile.write(
            'Fluent\n'
            '======\n\n')
        ashdry = self.PAash / (1 - self.PAmoist)
        char = self.Yields[self.SpeciesIndex(
            'Solid')] * (1. - ashdry) * 100.
        vol = (
            1 - self.Yields[self.SpeciesIndex('Solid')]) * (1. - ashdry) * 100.
        ash = 100. - char - vol

        self.CPDBalanceFile.write(
            tabulate(
                [
                    ['Volatiles', '{:4.3f}'.format(vol)],
                    ['Char', '{:4.3f}'.format(char)],
                    ['Ash', '{:4.3f}'.format(ash)],
                    ['Moisture', '{:4.3f}'.format(
                        self.moistureVolumeFraction())]
                ],
                tablefmt=tablefmt))
        self._addlines(2)

    def __CheckHydrogen(self):
        '''
        check the Hydrogen balance... last step
        '''
        H = self.UAH / MH  # hydrogen content in the coal kmol H / kg DAF
        HCH4 = self.Yields[self.SpeciesIndex(
            'CH4')] * 4. / (MC + 4. * MH)  # H in CH4 kmol/kg DAF
        HH2O = self.Yields[self.SpeciesIndex(
            'H2O')] * 2. / (MO + 2. * MH)  # H in H2O kmol/kg DAF
        if H < HCH4 + HH2O:
            ratio = 1.  # molar ratio C/H
            A11 = 1. / self.MTar
            A12 = 1. / (4. * MH + MC)
            CCO = self.Yields[self.SpeciesIndex('CO')] / (MC + MO)
            CCO2 = self.Yields[
                self.SpeciesIndex('CO2')] / (MC + 2. * MO)
            Cchar = self.Yields[self.SpeciesIndex('Solid')] / MC
            B1 = self.UAC / MC - (Cchar + CCO + CCO2)
            A21 = ratio * A11
            A22 = 4. / (MC + 4. * MH)
            B2 = H - HH2O
            A = numpy.array([[A11, A12], [A21, A22]])
            B = numpy.array([B1, B2])
            res = solve(A, B)
            TAR = self.Yields[self.SpeciesIndex('Tar')]
            CH4 = self.Yields[self.SpeciesIndex('CH4')]
            TAR += CH4 - res[1]
            # set the new TAR yield
            self.Yields[self.SpeciesIndex('Tar')] += TAR
            self.Yields[self.SpeciesIndex('CH4')] = res[
                1]  # set the new CH4 yield

    def __TarComp(self):  # saves [m,n,o] from C_m H_n O_o
        """Calculates and returns the tar composition."""
        # eq. (7) in Michele's report:
        # check Sum UA, sums up carbon part: Sulfur->Carbon, newer,
        # higher amount of char (incresed by sulfur)
        # Carbon:
        sum = 0
        self.muCTar = (
            self.UAC / MC - self.Yields[self.SpeciesIndex('Solid')] /
            MC - self.Yields[self.SpeciesIndex('CO2')] / (MC + 2. * MO) -
            self.Yields[self.SpeciesIndex('CO')] / (MC + MO) -
            self.Yields[self.SpeciesIndex('CH4')] /
            (MC + 4. * MH)) * (self.MTar / self.Yields[self.SpeciesIndex('Tar')])
        # Hydrogen:
        self.muHTar = (self.UAH / MH - 2. * self.Yields[self.SpeciesIndex('H2O')] /
                       (2. * MH + MO) - 4. * self.Yields[self.SpeciesIndex('CH4')] /
                       (MC + 4. * MH)) * (self.MTar / self.Yields[self.SpeciesIndex('Tar')])
        # Oxygen:
        self.muOTar = (self.UAO / MO - self.Yields[self.SpeciesIndex('H2O')] /
                       (2. * MH + MO) - 2. * self.Yields[self.SpeciesIndex('CO2')] /
                       (MC + 2. * MO) - self.Yields[self.SpeciesIndex('CO')] /
                       (MC + MO)) * (self.MTar / self.Yields[self.SpeciesIndex('Tar')])
        if abs(self.muOTar) < 1.e-10:
            self.muOTar = 0.0
        self.muSTar = self.UAS / MS * \
            (self.MTar / self.Yields[self.SpeciesIndex('Tar')])
        self.muNTar = (self.UAN / MN - 2. * self.Yields[self.SpeciesIndex('Other')] / (
            2. * MN)) * (self.MTar / self.Yields[self.SpeciesIndex('Tar')])
        if abs(self.muNTar) < 1.e-10:
            self.muNTar = 0.0
        # print 'TarComposition: ', 'C=',self.muCTar,',H=',
        # self.muHTar,',O=',self.muOTar,',N=',self.muNTar,',S=',self.muSTar
        self.CPDBalanceFile.write(
            'Tar\n'
            '---\n'
            'Atomic composition\n'
            '~~~~~~~~~~~~~~~~~~\n\n')
        self.CPDBalanceFile.write(str('C_%5.3f ' % self.muCTar) +
                                  str('H_%5.3f ' % self.muHTar) +
                                  str('O_%5.3f ' % self.muOTar) +
                                  str('N_%5.3f ' % self.muNTar) +
                                  str('S_%5.3f\n' % self.muSTar))
        self._addlines(1)
        self.CPDBalanceFile.write(
            str('Molecular weigth = %5.3f kg/kmol\n\n' % self.MTar))
        self._addlines(1)
        # self.CPDBalanceFile.write('Carbon:   '+ str(self.muCTar)+'\n')
        # self.CPDBalanceFile.write('Hydrogen: '+ str(self.muHTar)+'\n')
        # self.CPDBalanceFile.write('Oxygen:   '+ str(self.muOTar)+'\n')
        # self.CPDBalanceFile.write('Nytrogen:   '+ str(self.muNTar)+'\n')
        # self.CPDBalanceFile.write('Sulphur:   '+
        # str(self.muSTar)+'\n\n')
        self.CPDBalanceFile.write(
            'TAR oxidation\n'
            '~~~~~~~~~~~~~\n\n')
        nuO2 = 0.5 * (2. * self.muCTar + self.muHTar /
                      2. + 2. * self.muSTar - self.muOTar)
        self.CPDBalanceFile.write(
            str('TAR + %5.3f O2 --> ' % nuO2) +
            str('%5.3f CO2 + ' % self.muCTar) +
            str('%5.3f H2O + ' % (self.muHTar / 2)) +
            str('%5.3f N2 + ' % (self.muNTar / 2)) +
            str('%5.3f SO2\n\n' % self.muSTar))
        self.CPDBalanceFile.write(
            'TAR partial oxidation\n'
            '~~~~~~~~~~~~~~~~~~~~~\n\n')
        nuO2 = 0.5 * (self.muCTar + 2. * self.muSTar - self.muOTar)
        self.CPDBalanceFile.write(
            str('TAR + %5.3f O2 -->' % nuO2) +
            str(' %5.3f CO + ' % self.muCTar) +
            str('%5.3f H2 + ' % (self.muHTar / 2)) +
            str('%5.3f N2 + ' % (self.muNTar / 2)) +
            str('%5.3f SO2\n\n' % self.muSTar))
        # check mass weight
        MTarCheck = self.muCTar * MC + self.muHTar * MH + \
            self.muOTar * MO + self.muNTar * MN + self.muSTar * MS
        return [self.muCTar, self.muHTar, self.muOTar, self.muSTar, self.muNTar]

    def __Q_React(self):
        """Calculates the Heat of Reaction of the coal cumbustion."""
        # eq. (8) in Michele's report
        HHVdaf = self.HHV / (self.PAFC + self.PAVM)
        # eq. (9) in Michele's report:
        LHVdaf = HHVdaf - ((MO + 2. * MH) / (2. * MH)) * self.UAH * rH2O
        self.LHVdaf = LHVdaf
        LHVar = LHVdaf * (self.PAFC + self.PAVM) - self.PAmoist * rH2O
        # Q_react in J/mol; eq. (10) in Michele's report:
        self.Q_react = LHVdaf * MRaw
        # self.CPDBalanceFile.write('== Heating value ==\n')
        self.CPDBalanceFile.write(
            '\n'
            'Heating Value\n'
            '-------------\n\n')
        self.CPDBalanceFile.write(
            tabulate(
                [['HHV', '{:4.3f}'.format(
                    self.HHV / 1e6), '{:4.3f}'.format(HHVdaf / 1e6)],
                 ['LHV', '{:4.3f}'.format(LHVar / 1e6), '{:4.3f}'.format(LHVdaf / 1e6)]],
                headers=['', 'ar MJ/kg', 'daf MJ/kg'],
                tablefmt=tablefmt
            ))
        self._addlines(2)
        # self.CPDBalanceFile.write('|   | ar MJ/kg|daf MJ/kg|\n')
        # self.CPDBalanceFile.write(
        #    '|HHV|' + str('%9.3f|' % (self.HHV / 1e6)) + str('%9.3f|\n' % (HHVdaf / 1e6)))
        # self.CPDBalanceFile.write(
        #    '|LHV|' + str('%9.3f|' % (LHVar / 1e6)) + str('%9.3f|\n\n' % (LHVdaf / 1e6)))

    def _addlines(self, n=1):
        for _ in xrange(n):
            self.CPDBalanceFile.write('\n')

    def __hfRaw(self):
        """Calculates the heat of formation of the coal molecule and writes it into the output file."""
        # eq. (12) in Michele's report:
        muCRaw = self.UAC * (MRaw / MC)
        muHRaw = self.UAH * (MRaw / MH)
        muORaw = self.UAO * (MRaw / MO)
        muNRaw = self.UAN * (MRaw / MN)
        muSRaw = self.UAS * (MRaw / MS)
# eq. (13) in Michele's report:
#        print 'Q_React =',self.Q_react/1e6,'MJ/kg'
#        print 'muCRaw =',muCRaw
#        print 'muHRaw =',muHRaw
#        print 'muORaw =',muORaw
#        print 'H2O =',muHRaw
#        print 'muORaw =',muORaw
        muCO2 = muCRaw
        muH2O = 0.5 * muHRaw
        muN2 = 0.5 * muNRaw
        muSO2 = muSRaw
        muO2 = 0.5 * (2 * muCO2 + muH2O + 2 * muSO2 - muORaw)
#        print 'muCO2=',muCO2
#        print 'muH2O=',muH2O
#        print 'muO2=',muO2
#        print 'muN2=',muN2
#        print 'mu*hfCO2=',muCO2*hfCO2
#        print 'muhfH2O=',muH2O*hfH2O
#        print 'muhfO2=',muO2*hfO2
#        print 'muhfN2=',muN2*hfN2
        self.hfraw = self.Q_react + muCO2 * hfCO2 + muH2O * \
            hfH2O + muN2 * hfN2 + muSO2 * hfSO2 - muO2 * hfO2
        # self.CPDBalanceFile.write('Heat of Formation in J/mol: \n')
        # self.CPDBalanceFile.write('Raw Coal Molecule: '+
        # str(self.hfraw)+'\n')

    def __nysEq15(self):
        """Calculates the stoichiometric coefficients of the devolatilization reaction."""
        # self.hfRaw()
        # eq. (15) in Michele's report:
        self.nyTar = self.Yields[
            self.SpeciesIndex('Tar')] * MRaw / self.MTar
        self.nyChar = self.Yields[
            self.SpeciesIndex('Solid')] * MRaw / MC
        self.nyH2O = self.Yields[
            self.SpeciesIndex('H2O')] * MRaw / (2. * MH + MO)
        self.nyCO2 = self.Yields[
            self.SpeciesIndex('CO2')] * MRaw / (MC + 2. * MO)
        self.nyCH4 = self.Yields[
            self.SpeciesIndex('CH4')] * MRaw / (MC + 4. * MH)
        self.nyCO = self.Yields[
            self.SpeciesIndex('CO')] * MRaw / (MC + MO)
        self.nyN2 = self.Yields[
            self.SpeciesIndex('Other')] * MRaw / (2. * MN)
#        print 'nyTar=',self.nyTar
#        print 'nyChar=',self.nyChar
#        print 'nyH2O=',self.nyH2O
#        print 'nyCO2=',self.nyCO2
#        print 'nyCH4=',self.nyCH4
#        print 'nyCO=',self.nyCO
#        print 'nyN2=',self.nyN2
      #
#        print 'Tar=',self.Yields[self.SpeciesIndex('Tar')]
#        print 'Char=',self.Yields[self.SpeciesIndex('Solid')]
#        print 'H2O=',self.Yields[self.SpeciesIndex('H2O') ]
#        print 'CO2=',self.Yields[self.SpeciesIndex('CO2')]
#        print 'CO=',self.Yields[self.SpeciesIndex('CO')]
#        print 'N2=',self.Yields[self.SpeciesIndex('Other')]

    def __hfTar(self):
        """Calculates the heat of formation of tar and writes it into the CPD Composition file."""
        self.__nysEq15()
        # eq. (16) in Michele's report:
        Sum16 = self.nyH2O * hfH2O + self.nyCO2 * hfCO2 + self.nyCH4 * \
            hfCH4 + self.nyCO * hfCO + self.nyChar * hfChar + self.nyN2 * hfN2
        hfTar = (self.hfraw - Sum16) / self.nyTar
        self.CPDBalanceFile.write(
            'Enthalpy of formation\n'
            '~~~~~~~~~~~~~~~~~~~~~\n\n'
            'hf = {:10.3f} J/kmol\n'.format(hfTar))

    def __QPyro(self):
        """Calculates the heat of the pyrolysis process, assuming heat of tar formation is equal zero."""
        self.__nysEq15()
        # eq. (17) in Michele's report:
        Sum17 = self.nyH2O * hfH2O + self.nyCO2 * \
            hfCO2 + self.nyCH4 * hfCH4 + self.nyCO * hfCO
        QPyro = -(self.hfraw - Sum17) / (MRaw)
        # eq. (18) in Michele's report:
        QPyroVM = QPyro / (1 - self.Yields[self.SpeciesIndex('Solid')])
        # self.CPDBalanceFile.write('Pyrolysis Heat in J/g: \n')
        # self.CPDBalanceFile.write('Q= '+ str(QPyro)+'\n')
        # self.CPDBalanceFile.write('refered to VM= '+
        # str(QPyroVM)+'\n')

    def __closeResultFile(self):
        """Closes the CPD Composition file."""
        self.CPDBalanceFile.close()


###############FG-DVC specific class###########################


class FGPC_SpeciesBalance(SpeciesBalance):
    """This class calculates the Species and the Energy balance for FG-DVC and PCCL. See the manual for the formulars and more details."""

    def __init__(self, FGDVC_ResultObject, UAC, UAH, UAN, UAO, UAS, PAVM, PAFC, PAmoist, PAash, HHV, MTar, densityDryCoal, RunNr, fileprefix):
        # imports the dictionaries for the species defined in
        # 'CPD_Fit_one_run.py' -> CPD_Result()
        self.densityDryCoal = densityDryCoal
        self.Yields2Cols = FGDVC_ResultObject.DictYields2Cols()
        self.Cols2Yields = FGDVC_ResultObject.DictCols2Yields()
        # imports final yield composition defined in CPD_Fit_one_run.py'
        # -> CPD_Result()
        self.Yields = FGDVC_ResultObject.FinalYields()
        # 0:'Time', 1:'Temp', 2:'Tar', 3:'Gas', 4:'Solid', 5:'Total',
        # 6:'H2O', 7:'CO2', 8:'CH4', 9:'CO', 10:'Other'
        self.UAC = UAC / 100.
        self.UAH = UAH / 100.
        self.UAN = UAN / 100.
        self.UAO = UAO / 100.
        self.UAS = UAS / 100.
        self.PAVM = PAVM / 100.
        self.PAFC = PAFC / 100.
        self.PAash = PAash / 100.
        self.PAmoist = PAmoist / 100.
        if HHV == 0:
            HHV = self.Dulong()
        self.HHV = HHV
        self.MTar = MTar
        # corrects UA, if UA<1: Sulfur -> Carbon
        self.correctUA()
        self.FGBalanceFile = open(
            fileprefix + '-BalanceResults' + str(RunNr) + '.txt', 'w')
        self.FGBalanceFile.write(
            '= Coal Properties =\n==Ultimate Analysis==\n')
        self.FGBalanceFile.write(
            '|C |' + str('%6.3f' % (self.UAC * 100.)) + '%|\n')
        self.FGBalanceFile.write(
            '|H |' + str('%6.3f' % (self.UAH * 100.)) + '%|\n')
        self.FGBalanceFile.write(
            '|N |' + str('%6.3f' % (self.UAN * 100.)) + '%|\n')
        self.FGBalanceFile.write(
            '|O |' + str('%6.3f' % (self.UAO * 100.)) + '%|\n')
        self.FGBalanceFile.write(
            '|S |' + str('%6.3f' % (self.UAS * 100.)) + '%|\n\n')
        self.FGBalanceFile.write('== Proximate analysis ==\n')
        self.FGBalanceFile.write('|   |  AR% | dry% | daf% |\n')
        self.FGBalanceFile.write('|VM |' + str('%6.3f|' % PAVM) + str('%6.3f|' % (100 * PAVM / (
            100. - PAmoist))) + str('%6.3f|' % (100 * PAVM / (100. - PAmoist - PAash))) + '\n')
        self.FGBalanceFile.write('|FC |' + str('%6.3f|' % PAFC) + str('%6.3f|' % (100 * PAFC / (
            100. - PAmoist))) + str('%6.3f|' % (100 * PAFC / (100. - PAmoist - PAash))) + '\n')
        self.FGBalanceFile.write('|ash|' + str('%6.3f|' % PAash) + str(
            '%6.3f|' % (100 * PAash / (100. - PAmoist))) + str('%6.3f|' % 0.0) + '\n')
        self.FGBalanceFile.write(
            '|H2O|' + str('%6.3f|' % PAmoist) + str('%6.3f|' % 0.0) + str('%6.3f|' % 0.0) + '\n\n')
        # self.FGBalanceFile.write('Moisture volume fraction ='+str(self.moistureVolumeFraction()))
        # considered yields: char, tar, CO, CO2, H2O, CH4, N2, H2, O2
        self.__correctYields()
        # the missing of the UA is included into carbon
        self.correctUA()
        # calculates the Tar composition:
        self.__TarComp()
        # writes Composition results into file:
        self.__writeSpeciesResults()
        # Energy Balance:
        self.__EnergyBalance()
        self.__writeEnergyResults()
        self.__closeFile()
        if oSystem == 'Linux':
            shutil.move(fileprefix + '-BalanceResults' + str(RunNr) + '.txt',
                        'Result/' + fileprefix + '-BalanceResults' + str(RunNr) + '.txt')
        elif oSystem == 'Windows':
            shutil.move(fileprefix + '-BalanceResults' + str(RunNr) + '.txt',
                        'Result\\' + fileprefix + '-BalanceResults' + str(RunNr) + '.txt')
        else:
            print "The name of the operating system couldn't be found."

    def __correctYields(self):
        """Further, only the yields of char, tar, CO, CO2, H2O, CH4, N2, H2, O2 are considered. Modifies the yields, merge species like Olefins parafins, HCN to tar."""
        self.Char = 1.0 - self.Yields[self.SpeciesIndex('Total')]
        self.CO = self.Yields[self.SpeciesIndex('CO')]
        self.CO2 = self.Yields[self.SpeciesIndex('CO2')]
        self.H2O = self.Yields[self.SpeciesIndex('H2O')]
        self.CH4 = self.Yields[self.SpeciesIndex('CH4')]
        # self.N2=self.UAN  ?? elemental Nitrogen or all in Tar
        self.H2 = self.Yields[self.SpeciesIndex('H2')]
        # self.Tar=self.Yields[self.SpeciesIndex('Tar')]
        Sum = self.Char + self.CO + self.CO2 + self.H2O + self.CH4 + self.H2
        # Tar is the difference between 1 and the rest of species
        self.Tar = 1. - Sum
        # print 'Sum of modified Species: ',
        # self.Char+self.CO+self.CO2+self.H2O+self.CH4+self.H2+self.Tar#+self.N2
        # #??

    def __TarComp(self):
        """Calculates the Tar composition using analyisis of the Ultimate Analysis."""
        # The species in tar (kg) per kg of coal
        # Carbon:
        SumC = self.Char + self.CO * MC / \
            (MC + MO) + self.CO2 * MC / (MC + 2. * MO) + \
            self.CH4 * MC / (MC + 4. * MH)
        TarC = self.UAC - SumC
        # Hydrogen:
        SumH = 2. * self.H2O * MH / \
            (2 * MH + MO) + 4. * self.CH4 * MH / (4 * MH + MC) + self.H2
        TarH = self.UAH - SumH
        # Oxygen:
        SumO = self.CO * MO / (MC + MO) + 2. * self.CO2 * \
            MO / (MC + 2. * MO) + self.H2O * MO / (2. * MH + MO)
        TarO = self.UAO - SumO
        # Nitrogen:
        TarN = self.UAN
        TarS = self.UAS
        #
        # now, tar coefficients (molar, not mass)
        self.TarC = (TarC * self.MTar) / (self.Tar * MC)
        self.TarH = (TarH * self.MTar) / (self.Tar * MH)
        self.TarO = (TarO * self.MTar) / (self.Tar * MO)
        self.TarN = (TarN * self.MTar) / (self.Tar * MN)
        self.TarS = (TarS * self.MTar) / (self.Tar * MS)
        # print "Mass Weight of Tar " ,  self.TarC*MC+self.TarH*MH
        # +self.TarO*MO+self.TarN

    def __writeSpeciesResults(self):
        """Writes the Species Balance results into the result file."""
        self.FGBalanceFile.write('== Final yields ==\n')
        self.FGBalanceFile.write('Char ' + str('%7.5f\n' % self.Char))
        self.FGBalanceFile.write('Tar  ' + str('%7.5f\n' % self.Tar))
        self.FGBalanceFile.write('H2O  ' + str('%7.5f\n' % self.H2O))
        self.FGBalanceFile.write('H2   ' + str('%7.5f\n' % self.H2))
        self.FGBalanceFile.write('CO2  ' + str('%7.5f\n' % self.CO2))
        self.FGBalanceFile.write('CH4  ' + str('%7.5f\n' % self.CH4))
        self.FGBalanceFile.write('CO   ' + str('%7.5f\n\n' % self.CO))
        self.FGBalanceFile.write('=== Fluent ===\n')
        char = self.Char * (1. - self.PAash) * 100.
        vol = (1 - self.Char) * (1. - self.PAash) * 100.
        ash = 100. - char - vol

        self.FGBalanceFile.write(
            'Volatile component fraction: ' + str('%7.3f\n' % vol))
        self.FGBalanceFile.write(
            'Combustibile fraction:       ' + str('%7.3f\n' % char))
        self.FGBalanceFile.write(
            'Ash fraction:                ' + str('%7.3f\n' % ash))
        self.FGBalanceFile.write(
            'Moisture volume fraction:    ' + str('%7.3f\n\n' % self.moistureVolumeFraction()))

        self.FGBalanceFile.write(
            '== Tar ==\n=== Atomic composition ===\n')
        self.FGBalanceFile.write(str('C_%5.3f ' % self.TarC) + str('H_%5.3f ' % self.TarH) + str(
            'O_%5.3f ' % self.TarO) + str('N_%5.3f ' % self.TarN) + str('S_%5.3f\n' % self.TarS))
        self.FGBalanceFile.write(
            str('Molecular weigth = %5.3f kg/kmol\n\n' % self.MTar))

        # self.FGBalanceFile.write('The tar molecule composition:\n')
        # self.FGBalanceFile.write('Carbon:   '+str(self.TarC)+'\n')
        # self.FGBalanceFile.write('Hydrogen: '+str(self.TarH)+'\n')
        # self.FGBalanceFile.write('Oxygen:   '+str(self.TarO)+'\n')
        # self.FGBalanceFile.write('Nitrogen: '+str(self.TarN)+'\n\n\n')

    def __EnergyBalance(self):
        """Calculates the heat of formation of tar."""
        # Reaction enthalpie (J/mol) for: (E..Energy, M...molar)
        CharEM = hfChar + hfO2 - hfCO2
        H2EM = hfH2 + 0.5 * hfO2 - hfH2O
        CH4EM = hfCH4 + 2. * hfO2 - hfCO2 - 2. * hfH2O
        COEM = hfCO + 0.5 * hfO2 - hfCO2
#        print 'CharEM ', CharEM
#        print 'H2EM', H2EM
#        print 'CH4EM', CH4EM
#        print 'COEM', COEM, '\n'
        # Reaction enthalpie (J/kg), abs. including yields:
#        print 'CharEM/MC',CharEM/MC
#        print 'H2EM/(2.*MH)',H2EM/(2.*MH)
#        print '(CH4EM)/(4.*MH+MC)',(CH4EM)/(4.*MH+MC)
#        print '(COEM)/(MC+MO)',(COEM)/(MC+MO), '\n'
        CharE = (CharEM * self.Char) / MC
        H2E = (H2EM * self.H2) / (2. * MH)
        CH4E = (CH4EM * self.CH4) / (4. * MH + MC)
        COE = (COEM * self.CO) / (MC + MO)
#        print 'CharE', CharE
#        print 'H2E', H2E
#        print 'CH4E', CH4E
#        print 'COE', COE, '\n'
        # Reaction enthalpie Tar (J/kg):
        print 'HHV as rec =', self.HHV / 1e6, 'MJ/kg'
        HHVdaf = self.HHV / (self.PAFC + self.PAVM)
        print 'HHV daf =', HHVdaf / 1e6, 'MJ/kg'
        LHVdaf = HHVdaf - ((MO + 2. * MH) / (2. * MH)) * self.UAH * rH2O
        self.LHVdaf = LHVdaf
        print 'LHV-daf =', LHVdaf / 1e6, 'MJ/kmol'
        TarE = LHVdaf - CharE - H2E - CH4E - COE
#        print 'TarE', TarE
        # Reaction enthalpie Tar in J/mol
        TarEM = (TarE * self.MTar) / self.Tar
#        print 'TarEM',TarEM
        # coeffs
        nuN2 = self.TarN * 0.5
        nuH2O = self.TarH * 0.5
        nuCO2 = self.TarC
        nuO2 = 0.5 * (2 * self.TarC + 0.5 * self.TarH - self.TarO)
        # heat of formation for tar:
#        print 'O2', nuO2*hfO2
#        print 'CO2', nuCO2*hfCO2
#        print 'H2O', nuH2O*hfH2O
#        print 'N2', nuN2*hfN2
        self.hfTar = TarEM - nuO2 * hfO2 + nuCO2 * \
            hfCO2 + nuH2O * hfH2O + nuN2 * hfN2

    def __writeEnergyResults(self):
        """Writes the Energy results into the result file."""
        # self.FGBalanceFile.write('The lower heating value of daf coal:\n')
        # self.FGBalanceFile.write('LHV:   %.4f MJ/mol\n' % (self.LHVdaf*1.e-6))
        # self.FGBalanceFile.write('The enthalpie of foramtion for the tar:\n')
        # self.FGBalanceFile.write('hfTar: %.4f MJ/mol\n' %
        # (self.hfTar*1.e-6))
        self.FGBalanceFile.write(
            str('=== Enthalpy of formation ===\nhf = %10.3f J/kmol\n' % self.hfTar))

    def __closeFile(self):
        """Closes the FG-DVC Composition file."""
        self.FGBalanceFile.close()


class VolatileComposition(object):
    '''
    Class to handle volatile composition
    '''

    def __init__(self, composition, ultimate_analysis):
        '''
        Parameters
        ----------
        Composition: dict

        '''
        import cantera
        self.gas = cantera.Solution('52.xml')
        sum_comp = sum(composition.values())
        self.composition = {sp: v / sum_comp for sp,
                            v in composition.iteritems()}
        self.composition['C6H6'] = self.composition.pop('Tar')
        sum_ua = sum(ultimate_analysis.values())
        self.ultimate_analysis = {
            sp: v / sum_ua for sp, v in ultimate_analysis.iteritems()}

    def calc_element_fraction(self, element, species):
        '''
        Calculate element fraction of the given element
        '''
        if species in ('Char', 'Solid'):
            if element == 'C':
                return 1.0
            else:
                return 0.0
        elif species == 'Other':
            return 0.0
        else:
            i = self.gas.element_index(element)
            return (self.gas.atomic_weights[i] *
                    self.gas.species(species).composition.get(element, 0) /
                    self.gas.molecular_weights[
                    self.gas.species_index(species)])

    def calc_tot_element_fraction(self, element):
        return np.sum([val * self.calc_element_fraction(element, sp)
                       for sp, val in self.composition.iteritems()])

    def assign_nitrogen(self):
        self.composition['N2'] = self.ultimate_analysis['N']
        self.composition['Other'] = self.composition[
            'Other'] - self.composition['N2']

    def fix_oxygen(self):
        '''
        Remove CO2 to balance O 
        '''
        sum_o = self.ultimate_analysis[
            'O'] - self.calc_tot_element_fraction('O')
        o_in_co2 = self.calc_element_fraction('O', 'CO2')
        remove_co2 = sum_o / o_in_co2
        self.composition['CO2'] = self.composition[
            'CO2'] + remove_co2
        self.composition['Other'] = self.composition[
            'Other'] - remove_co2
        # self.calc_others()

    def add_CO2_from_oxygen(self):
        sum_o = self.ultimate_analysis[
            'O'] - self.calc_tot_element_fraction('O')
        o_in_co2 = self.calc_element_fraction('O', 'CO2')
        add_co2 = sum_o / o_in_co2
        self.composition['CO2'] = self.composition[
            'CO2'] + add_co2
        self.composition['Other'] = self.composition[
            'Other'] - add_co2

    def define_composition(self):
        self.assign_nitrogen()
        print 'ua=', self.ultimate_analysis
        remaining = self.calc_remaining()
        print 'comp=', self.composition
        print 'Solid: ', self.calc_element_fraction('C', 'Solid')
        print 'remaining', remaining
        print 'sum remaining:', sum(remaining.itervalues())
        print 'others:', self.composition['Other']
        if remaining['O'] < 0:
            print 'fix oxygen'
            self.fix_oxygen()
            remaining = self.calc_remaining()
            print 'remaining=', remaining
            print 'sum remaining:', sum(remaining.itervalues())
            print 'others:', self.composition['Other']
            print 'sum composition: ', sum(self.composition.values())
        else:
            print 'Add CO2 from remaining O'
            self.add_CO2_from_oxygen()
            remaining = self.calc_remaining()
            print 'remaining=', remaining
            print 'sum remaining:', sum(remaining.itervalues())
            print 'others:', self.composition['Other']
            print 'sum composition: ', sum(self.composition.values())

        c_to_h_mass = remaining['C'] / remaining['H']
        c_to_h_molar = c_to_h_mass * \
            self.gas.atomic_weight('H') / self.gas.atomic_weight('C')
        print '\nC/H molar:', c_to_h_molar
        if 0 <= c_to_h_molar <= 0.5:
            # use C2H4 and H2
            print('\nAdd C2H4')
            c2h4 = remaining['C'] / \
                self.calc_element_fraction('C', 'C2H4')
            self.composition['C2H4'] = self.composition.get(
                'C2H4', 0) + c2h4
            self.composition['Other'] = self.composition[
                'Other'] - c2h4
            remaining = self.calc_remaining()
            print 'remaining=', remaining
            print 'sum remaining:', sum(remaining.itervalues())
            print 'others:', self.composition['Other']
            print 'sum composition: ', sum(self.composition.values())
        elif 0.5 < c_to_h_molar <= 1:
            # use C6H6 and H2
            pass
        else:
            raise ValueError(
                'C/H molar ratio: {} > 1 not possible'
                ' to assign a composition'.format(c_to_h_molar))
        print('\nAdd H2')
        self.composition['H2'] = remaining['H']
        self.composition['Other'] = self.composition[
            'other'] - remaining['H']
        remaining = self.calc_remaining()
        print 'remaining=', remaining
        print 'sum remaining:', sum(remaining.itervalues())
        print 'others:', self.composition['Other']
        print 'sum composition: ', sum(self.composition.values())
        self.calc_volatiles()

    def calc_volatiles(self):
        volatiles = self.composition.copy()
        volatiles['Solid'] = 0.0
        sum_volatiles = sum(volatiles.itervalues())

        self.volatiles = {
            sp: val / sum_volatiles for sp, val in volatiles.iteritems()}

    def calc_remaining(self):
        sum_elements = {
            el: self.calc_tot_element_fraction(el)
            for el in self.gas.element_names}
        return {el: self.ultimate_analysis[el] - se
                for el, se in sum_elements.iteritems()}

    def define_composition_empirical(self, CO=0.5, tar=None):
        '''
        http://www.sciencedirect.com/science/article/pii/S0255270104001916
        '''
        # print self.ultimate_analysis
        composition = {}
        composition['Solid'] = self.composition['Solid']
        # composition['C6H6'] = self.composition['C6H6']
        tar = self.composition['C6H6']
        composition['N2'] = self.ultimate_analysis['N']
        composition['CO'] = CO * self.ultimate_analysis['O'] / \
            self.calc_element_fraction('O', 'CO')
        composition['CO2'] = (1 - CO) * self.ultimate_analysis['O'] / \
            self.calc_element_fraction('O', 'CO2')
        # print composition

        self.composition = composition

        remaining = self.calc_remaining()
        # print 'Set CO/CO2', remaining
        # print 'tar', tar
        C_in_tar = tar * self.calc_element_fraction('C', 'C6H6')
        # print 'C in tar', C_in_tar

        if C_in_tar < remaining['C']:
            self.composition['C6H6'] = tar
        else:
            self.composition['C6H6'] = remaining['C'] / \
                self.calc_element_fraction('C', 'C6H6')

        # print self.composition
        remaining = self.calc_remaining()
        # print 'remaining', remaining

        c_to_h_mass = remaining['C'] / remaining['H']
        c_to_h_molar = c_to_h_mass * \
            self.gas.atomic_weight('H') / self.gas.atomic_weight('C')
        # print '\nC/H molar:', c_to_h_molar
        if 0 <= c_to_h_molar <= 0.5:
            # use C2H4 and H2
            #print('\nAdd C2H4')
            c2h4 = remaining['C'] / \
                self.calc_element_fraction('C', 'C2H4')
            self.composition['C2H4'] = self.composition.get(
                'C2H4', 0) + c2h4
            remaining = self.calc_remaining()
            # print 'remaining=', remaining
            # print 'sum remaining:', sum(remaining.itervalues())
            # print 'sum composition: ', sum(self.composition.values())
        elif 0.5 < c_to_h_molar <= 1:
            # use C6H6 and H2
            pass
        else:
            raise ValueError(
                'C/H molar ratio: {} > 1 not possible'
                ' to assign a composition'.format(c_to_h_molar))

        self.composition['H2'] = remaining['H']
        remaining = self.calc_remaining()
        # print 'remaining=', remaining
        # print 'sum remaining:', sum(remaining.itervalues())
        # print 'sum composition: ', sum(self.composition.values())
        # print 'Composition', self.composition
        self.calc_volatiles()

    def heat_of_reaction_species(self, sp):
        gas = self.gas
        T_ref = 273
        if sp in ('N2', 'Others'):
            return 0.0
        elif sp in ('Solid', 'Char'):
            hf = hfChar * 1e6
            n_o2 = 1
            n_co2 = 1
            n_h2o = 0
            mw = gas.atomic_weight('C')
        else:
            spc = self.gas.species(sp)
            hf = spc.thermo.h(T_ref)
            n_o2 = 0.5 * (2 * spc.composition.get('C', 0) + 0.5 *
                          spc.composition.get('H', 0) - spc.composition.get('O', 0))
            n_co2 = spc.composition.get('C', 0)
            # print 'n_co2', n_co2
            n_h2o = spc.composition.get('H', 0) * 0.5
            # print 'n_h2o', n_h2o
            mw = gas.molecular_weights[gas.species_index(sp)]
        h_o2 = gas.species('O2').thermo.h(T_ref)
        # print 'h_o2', h_o2 / 1e6
        h_co2 = gas.species('CO2').thermo.h(T_ref)
        # print 'h_co2', h_co2 / 1e6
        h_h2o = gas.species('H2O').thermo.h(T_ref)
        # print 'h_h2o', h_h2o / 1e6

        heat_molar = (hf + n_o2 * h_o2 - n_co2 *
                      h_co2 - n_h2o * h_h2o)
        return heat_molar / mw

    def heat_of_volatiles(self):
        return sum([self.heat_of_reaction_species(sp) * val
                    for sp, val in self.volatiles.iteritems()])

    def heat_of_pyrolysis(self, heat_coal):
        '''
        Parameters
        ----------
        heat_coal: float
            LHV coal dry
        '''
        heat_vol = self.heat_of_volatiles()
        heat_char = self.heat_of_reaction_species('Solid')
        print heat_vol, heat_char

        heat_vol_tot = (
            heat_coal - self.composition['Solid'] * heat_char) /\
            (1 - self.composition['Solid'])
        print heat_vol_tot

        return heat_vol - heat_vol_tot
