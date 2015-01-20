"""
plot functions from models

"""

def pltYield(self, fgdvc_list, xValueToPlot, yValueToPlot):
    """ Plots the yields (to select with yValueToPlot) over Time or 
        Temperature (to slect with xValueToPlot). """
    for runnedCaseNr in range(len(fgdvc_list)):
        plt.plot(
            fgdvc_list[runnedCaseNr].Yield(xValueToPlot),
            fgdvc_list[runnedCaseNr].Yield(yValueToPlot)
        )
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
    """ Plots the rates (to select with yValueToPlot) over Time or 
        Temperature (to slect with xValueToPlot).
    """
    for runnedCaseNr in range(len(fgdvc_list)):
        plt.plot(
            fgdvc_list[runnedCaseNr].Rate(fgdvc_list[runnedCaseNr].SpeciesIndex(xValueToPlot)),
            fgdvc_list[runnedCaseNr].Rate(fgdvc_list[runnedCaseNr].SpeciesIndex(yValueToPlot)),
            label=yValueToPlot
        )
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

def plot(self,fgdvc_list,Species):
    """ Plot the yield and the rates over time with two curves: 
        one is the original data, the other the fitting curve. 
        Also file 'PyrolysisProgramName-Species.out' (e.g. 'CPD-CO2.out') 
        containing the time (s), yields (kg/kg), rates (kg/(kg s)).
    """
    #plots:
    colors=['r','b','g','black','purple']
    #Yields to compare
    u=[] #line index, time, column index: runned case
    v=[] #line index, time, column index: runned case
    for runnedCaseNr in range(len(fgdvc_list)):
        u_=fgdvc_list[runnedCaseNr].Yield(Species)
        v_=self.calcMass(
            fgdvc_list[runnedCaseNr],
            fgdvc_list[runnedCaseNr].Time(),
            fgdvc_list[runnedCaseNr].Interpolate('Temp'),
            Species)
        u.append(u_)
        v.append(v_)
    if type(Species)==int:
        SpeciesForTitle=fgdvc_list[0].SpeciesName(Species)
    if type(Species)==str:
        SpeciesForTitle=Species
    if SpeciesForTitle=='Solid':
        for runnedCaseNr in range(len(fgdvc_list)):
            plt.plot(
                fgdvc_list[runnedCaseNr].Time()[:len(u[runnedCaseNr])],
                u[runnedCaseNr],
                '-',
                color=colors[runnedCaseNr],
                label=fgdvc_list[0].Name()+' '+str(runnedCaseNr))
            plt.plot(
                fgdvc_list[runnedCaseNr].Time()[:len(v[runnedCaseNr])],
                v[runnedCaseNr],
                '--',
                color=colors[runnedCaseNr],
                label='fit')
            plt.plot(
                fgdvc_list[runnedCaseNr].Time()[:len(u[runnedCaseNr])],
                (1.-u[runnedCaseNr]),
                '-',
                color=colors[runnedCaseNr],
                label='Sum yields' + ' ' +str(runnedCaseNr))
            plt.plot(
                fgdvc_list[runnedCaseNr].Time()[:len(v[runnedCaseNr])],
                (1.-v[runnedCaseNr]),
                '--',
                color=colors[runnedCaseNr],
                label='fit')
    else:
        for runnedCaseNr in range(len(fgdvc_list)):
            plt.plot(
                fgdvc_list[runnedCaseNr].Time()[:len(u[runnedCaseNr])],
                u[runnedCaseNr],
                '-',
                color=colors[runnedCaseNr],
                label=fgdvc_list[0].Name()+' '+str(runnedCaseNr)
            )
            plt.plot(
                fgdvc_list[runnedCaseNr].Time()[:len(v[runnedCaseNr])],
                v[runnedCaseNr],
                '--',
                color=colors[runnedCaseNr],
                label='fit'
            )
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
        w=self.deriveC(fgdvc_list[runnedCaseNr],v[runnedCaseNr])
        plt.plot(fgdvc_list[runnedCaseNr].Time(),w,'--',color=colors[runnedCaseNr],label='fit')
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
        w=self.deriveC(fgdvc_list[runnedCaseNr],v[runnedCaseNr])
        ur=fgdvc_list[runnedCaseNr].Rate(Species)
        if oSystem=='Linux':
            print fgdvc_list[runnedCaseNr]
            resultFile=open('Result/'+fgdvc_list[runnedCaseNr].Name()+'-Fit_result_'+SpeciesForTitle+'_'+str(runnedCaseNr)+'.out','w')
        elif oSystem=='Windows':
            print fgdvc_list[runnedCaseNr]
            resultFile=open('Result\\'+fgdvc_list[runnedCaseNr].Name()+'-Fit_result_'+SpeciesForTitle+'_'+str(runnedCaseNr)+'.out','w')
        else:
            print 'Models: Operating Platform cannot be specified.'
        resultFile.write('    Time       Temperature    Yields       Rates    Yields(original) Rates(original) \n')
        for i in range(len(t)):
            resultFile.write('%7e  %11e %7e %8e %7e %8e \n' % (t[i], T[i], v[runnedCaseNr][i], w[i], u[runnedCaseNr][i], ur[i]))
        resultFile.close()
