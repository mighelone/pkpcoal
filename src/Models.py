import numpy as np
import pylab as plt
import scipy as sp
import scipy.integrate
import scipy.interpolate
import platform
#
PngResolution=100
#
oSystem=platform.system()
#the general parent class
class Model(object):
    """Parent class of the children ConstantRateModel, the three Arrhenius Models (notations) and the Kobayashi models."""
    def pltYield(self,fgdvc_list,xValueToPlot,yValueToPlot):
        """Plots the yields (to select with yValueToPlot) over Time or Temperature (to slect with xValueToPlot)."""
        for runnedCaseNr in range(len(fgdvc_list)):
            plt.plot(fgdvc_list[runnedCaseNr].Yield(xValueToPlot),fgdvc_list[runnedCaseNr].Yield(yValueToPlot))
        if xValueToPlot=='Time':
            plt.xlabel('t in s')
        if xValueToPlot=='Temp':
            plt.xlabel('T in K')
        if type(yValueToPlot)==int:
            SpeciesForTitle=fgdvc_list[0].SpeciesName(yValueToPlot)
        if type(yValueToPlot)==str:
            SpeciesForTitle=yValueToPlot
        plt.title(SpeciesForTitle)
        plt.ylabel('yield in wt%')
        plt.legend()
        plt.grid()
        if oSystem=='Linux':
            plt.savefig('Result/'+'Yields_'+yValueToPlot+'VS'+xValueToPlot+'.pdf',format='pdf')
        elif oSystem=='Windows':
            plt.savefig('Result\\'+'Yields_'+yValueToPlot+'VS'+xValueToPlot+'.pdf',format='pdf')
        else:
            print 'Models: Operating Platform cannot be specified.'
        plt.clf(),plt.cla()

    def pltRate(self,fgdvc_list,xValueToPlot,yValueToPlot):
        """Plots the rates (to select with yValueToPlot) over Time or Temperature (to slect with xValueToPlot)."""
        for runnedCaseNr in range(len(fgdvc_list)):
            plt.plot(fgdvc_list[runnedCaseNr].Rate(fgdvc_list[runnedCaseNr].SpeciesIndex(xValueToPlot)),fgdvc_list[runnedCaseNr].Rate(fgdvc_list[runnedCaseNr].SpeciesIndex(yValueToPlot)),label=yValueToPlot)
        if xValueToPlot=='Time':
            plt.xlabel('t in s')
        if xValueToPlot=='Temp':
            plt.xlabel('T in K')
        plt.ylabel('rate in wt%/s')
        if type(yValueToPlot)==int:
            SpeciesForTitle=fgdvc_list[0].SpeciesName(yValueToPlot)
        if type(yValueToPlot)==str:
            SpeciesForTitle=yValueToPlot
        plt.title(SpeciesForTitle)
        plt.legend()
        plt.grid()
        if oSystem=='Linux':
            plt.savefig('Result/'+'Rates_'+yValueToPlot+'VS'+xValueToPlot+'.pdf',format='pdf')
        elif oSystem=='Windows':
            plt.savefig('Result\\'+'Rates_'+yValueToPlot+'VS'+xValueToPlot+'.pdf',format='pdf')
        else:
            print 'Models: Operating Platform cannot be specified.'
        plt.clf(),plt.cla()
        
    def minLengthOfVectors(self,fgdvc_list):
        """Returns the minimum lenght of a all vectors from the several runs."""
        Len_tPointsL=[]
        for i in range(len(fgdvc_list)):
            Len_tPointsL.append(len(fgdvc_list[i].Time()))
        Len_tPoints=min(Len_tPointsL)
        return Len_tPoints

    def plot(self,fgdvc_list,Species):
        """Plot the yield and the rates over time with two curves: one is the original data, the other the fitting curve. Also file 'PyrolysisProgramName-Species.out' (e.g. 'CPD-CO2.out') containing the time (s), yields (kg/kg), rates (kg/(kg s))."""
        #plots:
        colors=['r','b','g','black','purple']
        #Yields to compare
        Len_tPoints=self.minLengthOfVectors(fgdvc_list)
        u=np.zeros([(Len_tPoints),len(fgdvc_list)]) #line index, time, column index: runned case
        v=np.zeros([(Len_tPoints),len(fgdvc_list)]) #line index, time, column index: runned case
        for runnedCaseNr in range(len(fgdvc_list)):
            u[:,runnedCaseNr]=fgdvc_list[runnedCaseNr].Yield(Species)[:Len_tPoints]
            v[:,runnedCaseNr]=self.calcMass(fgdvc_list[runnedCaseNr],fgdvc_list[runnedCaseNr].Time(),fgdvc_list[runnedCaseNr].Interpolate('Temp'),Species)[:Len_tPoints]
        if type(Species)==int:
            SpeciesForTitle=fgdvc_list[0].SpeciesName(Species)
        if type(Species)==str:
            SpeciesForTitle=Species
        if SpeciesForTitle=='Solid':
            for runnedCaseNr in range(len(fgdvc_list)):
                plt.plot(fgdvc_list[runnedCaseNr].Time()[:len(u)],u[:,runnedCaseNr],'-',color=colors[runnedCaseNr],label=fgdvc_list[0].Name()+' '+str(runnedCaseNr))
                plt.plot(fgdvc_list[runnedCaseNr].Time()[:len(v)],v[:,runnedCaseNr],'--',color=colors[runnedCaseNr],label='fit')
                plt.plot(fgdvc_list[runnedCaseNr].Time()[:len(u)],(1.-u[:,runnedCaseNr]),'-',color=colors[runnedCaseNr],label='Sum yields'+' '+str(runnedCaseNr))
                plt.plot(fgdvc_list[runnedCaseNr].Time()[:len(v)],(1.-v[:,runnedCaseNr]),'--',color=colors[runnedCaseNr],label='fit')
        else:
            for runnedCaseNr in range(len(fgdvc_list)):
                plt.plot(fgdvc_list[runnedCaseNr].Time()[:len(u)],u[:,runnedCaseNr],'-',color=colors[runnedCaseNr],label=fgdvc_list[0].Name()+' '+str(runnedCaseNr))
                plt.plot(fgdvc_list[runnedCaseNr].Time()[:len(v)],v[:,runnedCaseNr],'--',color=colors[runnedCaseNr],label='fit')            
        plt.title(SpeciesForTitle)
        plt.xlabel('t in s')
        plt.ylabel('yield fraction in kg/kg_coal')
        plt.legend()
        plt.grid()
        if oSystem=='Linux':
            plt.savefig('Result/'+fgdvc_list[0].Name()+'-Fit_result_'+SpeciesForTitle+'_Y.pdf',format='pdf')
            plt.savefig('Result/'+fgdvc_list[0].Name()+'-Fit_result_'+SpeciesForTitle+'_Y.png',dpi=PngResolution,format='png')
        elif oSystem=='Windows':
            plt.savefig('Result\\'+fgdvc_list[0].Name()+'-Fit_result_'+SpeciesForTitle+'_Y.pdf',format='pdf')
            plt.savefig('Result\\'+fgdvc_list[0].Name()+'-Fit_result_'+SpeciesForTitle+'_Y.png',dpi=PngResolution,format='png')
        else:
            print 'Models: Operating Platform cannot be specified.'
        plt.clf(),plt.cla()
        #Rates to compare
        for runnedCaseNr in range(len(fgdvc_list)):
            ur=fgdvc_list[runnedCaseNr].Rate(Species)
            plt.plot(fgdvc_list[runnedCaseNr].Time(),ur,'-',color=colors[runnedCaseNr],label=fgdvc_list[runnedCaseNr].Name()+' '+str(runnedCaseNr))
            w=self.deriveC(fgdvc_list[runnedCaseNr],v[:,runnedCaseNr],Len_tPoints)
            plt.plot(fgdvc_list[runnedCaseNr].Time()[:Len_tPoints],w[:Len_tPoints],'--',color=colors[runnedCaseNr],label='fit')
        if type(Species)==int:
            SpeciesForTitle=fgdvc_list[0].SpeciesName(Species)
        if type(Species)==str:
            SpeciesForTitle=Species
        plt.title(SpeciesForTitle)
        plt.xlabel('t in s')
        plt.ylabel('rate in 1/s')#min')
        plt.legend()
        plt.grid()
        if oSystem=='Linux':
            plt.savefig('Result/'+fgdvc_list[0].Name()+'-Fit_result_'+SpeciesForTitle+'_R.pdf',format='pdf')
            plt.savefig('Result/'+fgdvc_list[0].Name()+'-Fit_result_'+SpeciesForTitle+'_R.png',dpi=PngResolution,format='png')
        elif oSystem=='Windows':
            plt.savefig('Result\\'+fgdvc_list[0].Name()+'-Fit_result_'+SpeciesForTitle+'_R.pdf',format='pdf')
            plt.savefig('Result\\'+fgdvc_list[0].Name()+'-Fit_result_'+SpeciesForTitle+'_R.png',dpi=PngResolution,format='png')
        else:
            print 'Models: Operating Platform cannot be specified.'
        plt.clf(),plt.cla()
        #writes result file
        for runnedCaseNr in range(len(fgdvc_list)):
            t=fgdvc_list[runnedCaseNr].Yield('Time')
            T=fgdvc_list[runnedCaseNr].Yield('Temp')
            w=self.deriveC(fgdvc_list[runnedCaseNr],v[:,runnedCaseNr],Len_tPoints)
            if oSystem=='Linux':
                resultFile=open('Result/'+fgdvc_list[runnedCaseNr].Name()+'-Fit_result_'+SpeciesForTitle+str(runnedCaseNr)+'.out','w')
            elif oSystem=='Windows':
                resultFile=open('Result\\'+fgdvc_list[runnedCaseNr].Name()+'-Fit_result_'+SpeciesForTitle+str(runnedCaseNr)+'.out','w')
            else:
                print 'Models: Operating Platform cannot be specified.'
            resultFile.write('    Time       Temperature    Yields       Rates    Yields(original) Rates(original) \n')
            for i in range(len(v)):
                resultFile.write('%7e  %11e %7e %8e %7e %8e \n' % (t[i], T[i], v[i,runnedCaseNr], w[i], u[i,runnedCaseNr], ur[i]))
            resultFile.close()

    def deriveC(self,fgdvc,yVector,maxVectorLenght=None):
        """Returns a CDS of the inputted yVector."""
        dt=fgdvc.Dt()
        if maxVectorLenght==None:
            yDot=np.zeros(fgdvc.NPoints())
        else:
            yDot=np.zeros(maxVectorLenght)
        yDot[0]=(yVector[1]-yVector[0])/dt[0]
#        print (yVector[2:]-yVector[:-2])
#        print (2*dt[1:-1])
        if maxVectorLenght!=None:
            dt=dt[:maxVectorLenght]
        yDot[1:-1]=(yVector[2:]-yVector[:-2])/(2*dt[1:-1])
        yDot[-1]=(yVector[-1]-yVector[-2])/dt[-1]
        return yDot
        
    def calcRate(self,fgdvc,time,T,Name):
        """Generates the Rates using the yields vector and a CDS."""
        m=self.calcMass(fgdvc,time,T,Name)
        mDot=self.deriveC(fgdvc,m)
        return mDot

    def setParamVector(self,ParameterList):
        """Sets the Vector containing the kinetic parameter of the Model (refering to the child model)."""
        self._ParamVector=ParameterList

    def ParamVector(self):
        """Returns the Vector containing the kinetic parameter of the Model (refering to the child model)."""
        return self._ParamVector

    def ErrorYield(self,fgdvc,Species):                                  #calculates avrg Error between curves
        """Returns the absolute deviation per point between the fitted and the original yield curve."""
        v=self.calcMass(fgdvc,fgdvc.Time(),fgdvc.Interpolate('Temp'),Species)
        u=fgdvc.Yield(Species)
        absE=(np.sum(u-v))/fgdvc.NPoints()
        return absE

    def ErrorRate(self,fgdvc,Species):                                  #calculates avrg Error between curves
        """Returns the absolute deviation per point between the fitted and the original rate curve."""
        v=self.calcMass(fgdvc,fgdvc.Time(),fgdvc.Interpolate('Temp'),Species)
        uDot=fgdvc.Rate(Species)
        vDot=self.deriveC(fgdvc,v)
        absE=(np.sum(uDot-vDot))/fgdvc.NPoints()
        return absE
    


################childrenclasses####################

class ConstantRateModel(Model):
    """The model calculating the mass with m(t)=m_s0+(m_s0-m_s,e)*e**(-k*(t-t_start)) from the ODE dm/dt = -k*(m-m_s,e). The Parameter to optimize are k and t_start."""
    def __init__(self,InitialParameterVector):
        print 'Constant rate initialized'
        self._ParamVector=InitialParameterVector

    def calcMass(self,fgdvc,t,T,SpeciesToCalc):
        ParamVector=self.ParamVector()
        u=fgdvc.Yield(SpeciesToCalc)
        u_0=u[0]; u_s0=ParamVector[2]
        if SpeciesToCalc=='Solid' or SpeciesToCalc==(fgdvc.SpeciesIndex('Solid')):
            v=u_s0+(u_0-u_s0)*np.exp(-ParamVector[0]*(t-ParamVector[1]))
        else:
            v=u_s0*(1.-np.exp(-ParamVector[0]*(t-ParamVector[1])))
        v=np.where(t>ParamVector[1],v,u_0)
        return v

                
class ArrheniusModel(Model):
    """The Arrhenius model in the standart notation: dm/dt=A*(T**b)*exp(-E/T)*(m_s-m) with the parameter a,b,E to optimize."""
    def __init__(self,InitialParameterVector):
        print 'Arrhenuis Model initialized'
        self._ParamVector=InitialParameterVector
        self.ODE_hmax=1.e-2
 
    def calcMass(self,fgdvc,time,T,Name):
        """Outputs the mass(t) using the model specific equation."""
        #numercal values:
        absoluteTolerance = 1.0e-8
        relativeTolerance = 1.0e-6
        ##################
        u=fgdvc.Yield(Name)
        ParamVec=self.ParamVector()
        m_s0=ParamVec[3]
        def dmdt(m,t):
            #A=parameter[0]
            #b=parameter[1]
            #E=parameter[2]
            if Name == 'Solid':# or Name == fgdvc.Yields2Cols['Solid']:
                dmdt_out=-ParamVec[0]*(T(t)**ParamVec[1])*np.exp(-ParamVec[2]/(T(t)))*(m-m_s0)
                dmdt_out=np.where(abs(dmdt_out)>1.e-300,dmdt_out,0.0) #sets values<0 =0.0, otherwise it will further cause problems (nan)
            else:
                dmdt_out= ParamVec[0]*(T(t)**ParamVec[1])*np.exp(-ParamVec[2]/(T(t)))*(m_s0-m)
                dmdt_out=np.where(abs(dmdt_out)>1.e-300,dmdt_out,0.0) #sets values<0 =0.0, otherwise it will further cause problems (nan)
                #print 'ParamVec[0]',ParamVec[0]
                #print '(T(t)**ParamVec[1])',(T(t)**ParamVec[1])
                #print 'np.exp(-ParamVec[2]/(T(t)))', np.exp(-ParamVec[2]/(T(t)))
                #print '(m_s0-m)',(m_s0-m)
            return dmdt_out
        InitialCondition=[u[0]]
        m_out=sp.integrate.odeint(dmdt,InitialCondition,time,atol=absoluteTolerance,rtol=relativeTolerance,hmax=self.ODE_hmax) 
        if (ParamVec[0]<0 or ParamVec[2]<0):
            m_out[:,0]=float('inf')
            return m_out[:,0]
        else:
            return m_out[:,0]


    def ConvertKinFactors(self,ParameterVector):
        """Dummy. Function actual has to convert the parameter into the standart Arrhenius notation."""
        #does nothing, just to have the same way of use for all notations
        return ParameterVector

        
        
class ArrheniusModelAlternativeNotation1(ArrheniusModel):
    """Arrhenius model with a notation having a better optimization behaviour: dm/dt=exp[k0-a*(T0/T(t)-1)]*(ms-m). See the documentation for the reference. The parameters to optimize are k0 and a."""
    def __init__(self,InitialParameterVector):
        print 'Arrhenuis Model class initialized'
        self._ParamVector=InitialParameterVector
        self.ODE_hmax=1.e-2
        
    def calcMass(self,fgdvc,time,T,Name):
        """Outputs the mass(t) using the model specific equation."""
        #numercal values:
        absoluteTolerance = 1.0e-8
        relativeTolerance = 1.0e-6
        ##################
        ParamVec=self.ParamVector()
        u=fgdvc.Yield(Name)
        m_s0=ParamVec[2]
        twhereTmax=time[fgdvc.LineNumberMaxRate(Name)]
        self.T0=T(twhereTmax)
        def dmdt(m,t):
            #k0=parameter[0]
            #n=0
            #a=parameter[1]
            if Name == 'Solid':# or Name == fgdvc.Yields2Cols['Solid']:
                dmdt_out=-np.exp( ParamVec[0]  - ParamVec[1]*( self.T0/T(t) - 1 ))*(m-m_s0) #-ParamVec[0]*( (T(t)/self.T0)**ParamVec[1] )*np.exp( -ParamVec[2]*(self.T0/(T(t))-1) )*(m)    #IC
                dmdt_out=np.where(abs(dmdt_out)>1.e-300,dmdt_out,0.0) #sets values<0 =0.0, otherwise it will further cuase problem(nan)
            else:
                dmdt_out= np.exp( ParamVec[0] - ParamVec[1]*( self.T0/T(t) - 1 ))*(m_s0-m)
            return dmdt_out
        InitialCondition=[u[0]]
        m_out=sp.integrate.odeint(dmdt,InitialCondition,time,atol=absoluteTolerance,rtol=relativeTolerance)
        return m_out[:,0]

    def ConvertKinFactors(self,ParameterVector):
        """Converts the own kinetic factors back to the standard Arrhenius kinetic factors."""
        #k=ParameterVector[0]
        #a=ParameterVector[1]
        A=np.exp( ParameterVector[0] + ParameterVector[1] )
        E=ParameterVector[1]*self.T0
        return [A,0.0,E,ParameterVector[2]]

    def ConvertKinFactorsToOwnNotation(self,fgdvc,ParameterVector,Species):
        """Converts the standard Arrhenius kinetic factors backk to the factors of the own notation."""
        time=fgdvc.Time()
        twhereTRmax=time[fgdvc.LineNumberMaxRate(Species)]
        T=fgdvc.Interpolate('Temp')
        self.T0=T(twhereTRmax)
        A=ParameterVector[0]
        if ParameterVector[1]!=0:
            print 'b is not considered in this Notation and set to 0'
        E=ParameterVector[2]
        a=E/self.T0
        k0=np.log(A)-a
        return [k0,a,ParameterVector[3]]

        
class ArrheniusModelAlternativeNotation2(ArrheniusModel):
    """Arrhenius model with a notation having a better optimization behaviour: dm/dt=exp[c*(b1*(1/T(t)-1/T_min)-b2*(1/T(t)-1/T_max))]*(ms-m); with c=(1/T_max-1/Tmin)**(-1). See the documentation for the reference. The parameters to optimize are b1 and b2."""
    def __init__(self,InitialParameterVector):
        print 'Arrhenius Model class initialized'
        self._ParamVector=InitialParameterVector
        self.T_min=300
        self.T_max=1500
        
    def setMinMaxTemp(self,Tmin,Tmax):
        """Sets the temperature constants, see the equation."""
        self.T_min=Tmin
        self.T_max=Tmax
        
    def calcMass(self,fgdvc,time,T,Name):
        """Outputs the mass(t) using the model specific equation."""
        self.c=1./(1./self.T_max-1./self.T_min)
        self.ODE_hmax=1.e-2
        self.fgdvc=fgdvc
        #numercal values:
        absoluteTolerance = 1.0e-8
        relativeTolerance = 1.0e-6
        ##################
        ParamVec=self.ParamVector()
        u=fgdvc.Yield(Name)
        #uDot=fgdvc.Rate(Name)
        m_s0=ParamVec[2]
        def dmdt(m,t):
            #b1=parameter[0]
            #b2=parameter[1]
            if Name == 'Solid':# or Name == fgdvc.Yields2Cols['Solid']:
                dmdt_out= -np.exp( self.c*( ParamVec[0]*(1./T(t)-1./self.T_min) - ParamVec[1]*(1./T(t)-1./self.T_max) ) )*(m-m_s0)
                dmdt_out=np.where(abs(dmdt_out)>1.e-300,dmdt_out,0.0) #sets values<0 =0.0, otherwise it will further cause problem(nan)
            else:
                dmdt_out= np.exp( self.c*( ParamVec[0]*(1./T(t)-1./self.T_min) - ParamVec[1]*(1./T(t)-1./self.T_max) ) )*(m_s0-m)
                dmdt_out=np.where(abs(dmdt_out)>1.e-300,dmdt_out,0.0) #sets values<0 =0.0, otherwise it will further cause problem(nan)
            return dmdt_out
        InitialCondition=[u[0]]
        m_out=sp.integrate.odeint(dmdt,InitialCondition,time,atol=absoluteTolerance,rtol=relativeTolerance,hmax=self.ODE_hmax)
        ArrStandartParam=self.ConvertKinFactors(ParamVec)
        if (ArrStandartParam[0]<0 or ArrStandartParam[2]<0):
            m_out[:,0]=float('inf')
            return m_out[:,0]
        else:
            return m_out[:,0]

    def ConvertKinFactors(self,ParameterVector):
        """Converts the own kinetic factors back to the standard Arrhenius kinetic factors."""
        #b1=ParameterVector[0]
        #b2=ParameterVector[1]
        self.c=1./(1./self.T_max-1./self.T_min)
        A=np.exp( -self.c*ParameterVector[0]/self.T_min + self.c*ParameterVector[1]/self.T_max )
        E=self.c*ParameterVector[1]-self.c*ParameterVector[0]
        m_s0=ParameterVector[2]
        return [A,0,E,m_s0]

    def ConvertKinFactorsToOwnNotation(self,ParameterVector):
        """Converts the standard Arrhenius kinetic factors back to the factors of the own notation."""
        self.c=1/(1/self.T_max-1/self.T_min)
        #A=ParameterVector[0]
        #b not used
        #E=ParameterVector[2]
        if ParameterVector[1]!=0:
            print 'b is not considered in this Notation and set to 0'
        b2=( (self.T_max/self.c)*np.log(ParameterVector[0])-(ParameterVector[2]*self.T_max)/(self.c*self.T_min) )*(1-self.T_max/self.T_min)**(-1)
        b1=b2-ParameterVector[2]/self.c
        return [b1,b2,ParameterVector[3]]

class Kobayashi(Model):
    """Calculates the devolatilization reaction using the Kobayashi model. The Arrhenius equation inside are in the standard notation."""
    def __init__(self,InitialParameterVector):
        print 'Kobayashi Model initialized'
        self._ParamVector=InitialParameterVector
        self.ODE_hmax=1.e-2
        
    def calcMass(self,fgdvc,time,T,Name):
        """Outputs the mass(t) using the model specific equation. The input Vector is [A1,E1,A2,E2,alpha1,alpha2]"""
        self.fgdvc=fgdvc
        time=np.delete(time,0)
        self.__Integral=0.0
        tList=[0.0]
        k1k2=[0.0]
        ParamVec=self.ParamVector()
        #
        def dmdt(m,t):
            k1=ParamVec[0]*np.exp(-ParamVec[1]/(T(t)))
            k2=ParamVec[2]*np.exp(-ParamVec[3]/(T(t)))
            tList.append(t)
            k1k2.append(k1+k2)
            self.__Integral+=0.5*(tList[-1]-tList[-2])*(k1k2[-1]+k1k2[-2])
            dmdt_out = ( (ParamVec[4]*k1+ParamVec[5]*k2)*np.exp(-self.__Integral) )
            dmdt_out=np.where(abs(dmdt_out)>1.e-300,dmdt_out,0.0) #sets values<0 =0.0, otherwise it will further cause problem(nan)
            return dmdt_out
        InitialCondition=[0]
        m_out=sp.integrate.odeint(dmdt,InitialCondition,time,atol=1.e-5,rtol=1.e-4,hmax=1.e-2)
        m_out=m_out[:,0]
        m_out=np.insert(m_out,0,0.0)
        if (ParamVec[0]<0 or ParamVec[1]<0 or ParamVec[2]<0 or ParamVec[3]<0 or ParamVec[4]<0  or ParamVec[5]>1  ):
            m_out[:]=float('inf')
            return m_out
        else:
            return m_out
        

#    def ConvertKinFactors(self,ParameterVector):
#        """Outputs the Arrhenius equation factors in the shape [A1,E1,A2,E2]. Here where the real Arrhenius model is in use only a dummy function."""
#        P=self.ParamVector()
#        return [P[0],P[1],P[2],P[3]]
#    
#    def setKobWeights(self,alpha1,alpha2):
#        """Sets the two Kobayashi weights alpha1 and alpha2."""
#        self.__alpha1=alpha1
#        self.__alpha2=alpha2
#        
#    def KobWeights(self):
#        """Returns the two Kobayashi weights alpha1 and alpha2."""
#        return self.__alpha1, self.__alpha2

class KobayashiPCCL(Model):
    """Calculates the devolatilization reaction using the Kobayashi model. The Arrhenius equation inside are in the standard notation. The fitting parameter are as in PCCL A1,A2,E1,alpha1."""
    def __init__(self,InitialParameterVector):
        print 'Kobayashi Model initialized'
        self._ParamVector=InitialParameterVector
        self.ODE_hmax=1.e-2
        
    def calcMass(self,fgdvc,time,T,Name):
        """Outputs the mass(t) using the model specific equation."""
        self.fgdvc=fgdvc
        time=np.delete(time,0)
        self.__Integral=0.0
        tList=[0.0]
        k1k2=[0.0]
        ParamVec=self.ParamVector()
        #
        def dmdt(m,t):
            k1=ParamVec[0]*np.exp(-ParamVec[2]/(T(t)))
            k2=ParamVec[1]*np.exp(-(ParamVec[2]+self.__E2diff)/(T(t)))
            tList.append(t)
            k1k2.append(k1+k2)
            self.__Integral+=0.5*(tList[-1]-tList[-2])*(k1k2[-1]+k1k2[-2])
            dmdt_out = ( (ParamVec[3]*k1+self.__alpha2*k2)*np.exp(-self.__Integral) )
            dmdt_out=np.where(abs(dmdt_out)>1.e-300,dmdt_out,0.0) #sets values<0 =0.0, otherwise it will further cause problem(nan)
            return dmdt_out
        InitialCondition=[0]
        m_out=sp.integrate.odeint(dmdt,InitialCondition,time,atol=1.e-5,rtol=1.e-4,hmax=1.e-2)
        m_out=m_out[:,0]
        m_out=np.insert(m_out,0,0.0)
        if (ParamVec[0]<0 or ParamVec[1]<0 or ParamVec[2]<0 or ParamVec[3]<0):
            m_out[:]=float('inf')
            return m_out
        else:
            return m_out
        

    def ConvertKinFactors(self,ParameterVector):
        """Outputs the Arrhenius equation factors in the shape [A1,E1,A2,E2]. Here where the real Arrhenius model is in use only a dummy function."""
        P=self.ParamVector()
        return [P[0],P[1],P[2],P[3]]
    
    def setKobWeights(self,alpha2):
        """Sets the two Kobayashi weights alpha2."""
        self.__alpha2=alpha2
        
    def KobWeights(self):
        """Returns the two Kobayashi weights alpha2."""
        return self.__alpha2
    
    def setE2Diff(self,DifferenceE1E2):
        """Sets the dE in E2=E1+dE."""
        self.__E2diff=DifferenceE1E2
        
    def E2Diff(self):
        """Returns the dE in E2=E1+dE."""
        return self.__E2diff

class KobayashiA2(Model):
    """Calculates the devolatilization reaction using the Kobayashi model. The Arrhenius equation inside are in the secend alternative notation (see class ArrheniusModelAlternativeNotation2)."""
    def __init__(self,InitialParameterVector):
        print 'Kobayashi Model initialized'
        self._ParamVector=InitialParameterVector
        self.ODE_hmax=1.e-2

    def calcMass(self,fgdvc,time,T,Name):
        """Outputs the mass(t) using the model specific equation."""
        self.fgdvc=fgdvc
        T_general=fgdvc.Rate('Temp')
        self.T_min=min(T_general)
        self.T_max=max(T_general)
        self.c=1./(1./self.T_max-1./self.T_min)
        #alpha has to be greater equal zero:
        absoluteTolerance = 1.0e-8
        relativeTolerance = 1.0e-6
        self.tList=[0.0]
        #deletes to solve trapezian rule:
        time=np.delete(time,0)
        self.Integral=0.0
        self.k1k2=[0.0]
        ParamVec=self.ParamVector()
        u=fgdvc.Yield(Name)
        #
        def dmdt(m,t):
            k1=np.exp( self.c*( ParamVec[0]*(1./T(t)-1./self.T_min) - ParamVec[1]*(1./T(t)-1./self.T_max) ) )
            k2=np.exp( self.c*( ParamVec[2]*(1./T(t)-1./self.T_min) - ParamVec[3]*(1./T(t)-1./self.T_max) ) )
            self.tList.append(t)
            self.k1k2.append(k1+k2)
            self.Integral+=0.5*(self.tList[-1]-self.tList[-2])*(self.k1k2[-1]+self.k1k2[-2])
            dmdt_out = ( (self.__alpha1*k1+self.__alpha2*k2)*np.exp(-self.Integral) )
            dmdt_out=np.where(abs(dmdt_out)>1.e-300,dmdt_out,0.0) #sets values<0 =0.0, otherwise it will further cause problem(nan)
            return dmdt_out
        InitialCondition=[u[0]]
        m_out=sp.integrate.odeint(dmdt,InitialCondition,time,atol=absoluteTolerance,rtol=relativeTolerance,hmax=self.ODE_hmax)
        m_out=m_out[:,0]
        m_out=np.insert(m_out,0,0.0)
        return m_out

    def ConvertKinFactors(self,ParameterVector):
        """Converts the alternative notaion Arrhenius factors into the satndard Arrhenius factors and return them in the shape  [A1,E1], [A2,E2]"""
        A1=np.exp( -self.c*ParameterVector[0]/self.T_min + self.c*ParameterVector[1]/self.T_max )
        E1=self.c*ParameterVector[1]-self.c*ParameterVector[0]
        A2=np.exp( -self.c*ParameterVector[2]/self.T_min + self.c*ParameterVector[3]/self.T_max )
        E2=self.c*ParameterVector[3]-self.c*ParameterVector[2]
        return [A1,E1,A2,E2]

    def setKobWeights(self,alpha1,alpha2):
        """Sets the two Kobayashi weights alpha1 and alpha2."""
        self.__alpha1=alpha1
        self.__alpha2=alpha2
        
    def KobWeights(self):
        """Returns the two Kobayashi weights alpha1 and alpha2."""
        return self.__alpha1, self.__alpha2
        

class DAEM(Model):
    """Calculates the devolatilization reaction using the Distributed Activation Energy Model."""
    def __init__(self,InitialParameterVector):
        print 'DAEM initialized'
        self._ParamVector=InitialParameterVector
        self.ODE_hmax=1.e-2
        self.NrOfActivationEnergies=50
    
    def setNrOfActivationEnergies(self,NrOfE):
        """Define for how many activation energies of the range of the whole distribution the integral shall be solved (using Simpson Rule)."""
        self.NrOfActivationEnergies=NrOfE
        
    def NrOfActivationEnergies(self):
        """Returns the number of activation enrgies the integral shall be solved for (using Simpson Rule)."""
        return self.NrOfActivationEnergies
        
    def calcMass(self,fgdvc,time,T,Name):
        """Outputs the mass(t) using the model specific equation."""
        self.E_List=np.arange(int(self._ParamVector[1]-3.*self._ParamVector[2]),int(self._ParamVector[1]+3.*self._ParamVector[2]),int((6.*self._ParamVector[2])/self.NrOfActivationEnergies)) #integration range E0 +- 3sigma, see [Cai 2008]
        #Inner Integral Funktion
        def II_dt(t,E_i):
            return np.exp( -E_i/T(t) )
        #outer Integral for one activation energy from t0 to tfinal
        #stores all values of the inner Integrals (time,ActivationEnergy) in a 2D-Array
        InnerInts=np.zeros([len(time),len(self.E_List)])
        CurrentInnerInt=np.zeros(len(time))
        for Ei in range(len(self.E_List)):
            CurrentInnerInt[:]=II_dt(time[:],self.E_List[Ei])
            InnerInts[1:,Ei] = sp.integrate.cumtrapz(CurrentInnerInt,time[:])
        #
        def OI_dE(EIndex,tIndex):
            m = np.exp(-self._ParamVector[0]*InnerInts[tIndex,EIndex])*(1./(self._ParamVector[2]*(2.*np.pi)**0.5))*np.exp(-(self.E_List[EIndex]-self._ParamVector[1])**2/(2.*self._ParamVector[2]**2)) 
#            print 'InnerInt',InnerInt,'mass',dm_dt
            return m
        m_out=np.zeros(np.shape(time))
        mE=np.zeros(np.shape(self.E_List))
        for ti in range(len(time)):
            for Ei in range(len(self.E_List)):
                mE[Ei]=OI_dE(Ei,ti)
            m_out[ti]=sp.integrate.simps(mE,self.E_List)
        #descaling
        m_out = self._ParamVector[3]*(1.-m_out)
        return m_out
