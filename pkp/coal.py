import numpy as np
from autologging import logged


@logged
class coal():

    def __init__(self, name='', c=0.845, h=0.05, o=0.1, n=0.05, s=0.):
        '''
        coal(name='',c=0.845,h=0.05,o=0.1,n=0.05,s=0.):
        create coal object, assigning the ultimate analysis
        '''
        self._Mc = 12.
        self._Mh = 1.
        self._Mo = 16.
        self.setUA(c, h, o, n, s)
        self.setName(name)
        '''
        self._cdaf = 0.845
        self._hdaf = 0.05
        self._odaf = 0.1
        self._ndaf = 0.05
        self._sdaf = 0.
        '''
        self._moist = 0.1
        self._ash = 0.1
        self._vm = 0.5
        self._char = 0.4
        self._rhoDry = 1400
        self._rhoH2O = 1000

    def setUA(self, c=0.845, h=0.05, o=0.1, n=0.05, s=0.):
        '''
        def setUA(self,c,h,o,n,s)
        set ultimate analysis of coal

        if sum is different from 1 values are sclaled
        '''
        sumUA = c + h + o + n + s
        self._c = c / sumUA
        self._h = h / sumUA
        self._o = o / sumUA
        self._n = n / sumUA
        self._s = s / sumUA

    def setPA(self, moisture=0.1, char=0.4, vm=0.5, ash=0.1):
        sumPA = moisture + char + vm + ash
        self._moist = moisture / sumPA
        self._ash = self._ash / sumPA
        self._vm = vm / sumPA
        self._char = char / sumPA
        self._daf = self._vm + self._char

    def setName(self,name):
        '''
        set name of the coal
        '''
        self._name = name

    def printName(self):
        '''
        print name of the coal
        '''
        print self._name

    def getMoist(self):
        return self._moist

    def getAsh(self):
        return self._ash

    def getChar(self):
        return self._char

    def getVM(self):
        return self._vm

    def __repr__(self):
        output = self._name+' coal\n\n'
        output += 'Ultimate analysis\n'
        output += 'C ='+str(self._c)+'\n'
        output += 'H ='+str(self._h)+'\n'
        output += 'O ='+str(self._o)+'\n'
        output += 'N ='+str(self._n)+'\n'
        output += 'S ='+str(self._s)+'\n\n'
        output += 'Proximate analysis\n'
        output += 'Char ='+str(self._char)+'\n'
        output += 'VM ='+str(self._vm)+'\n'
        output += 'Ash ='+str(self._ash)+'\n'
        output += 'Moist ='+str(self._moist)+'\n'

        return output

    def getMoistureVolFraction(self):
        vH2O = self._moist / self._rhoH2O
        vDry = (1-self._moist) / self._rhoDry
        return vH2O/(vH2O+vDry)

    def getVanKravelen(self):
        '''
        return coordinate O/C and H/C molar ratios
        in van Kravelen diagram
        '''
        oc = (self._o/self._Mo)/(self._c/self._Mc)
        hc = (self._h/self._Mh)/(self._c/self._Mc)
        return np.array([oc,hc])

    def getDensity(self):
        return self._rhoDry



