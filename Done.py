# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'Done.ui'

from PyQt4 import QtCore, QtGui
import platform
import os,sys
import numpy as np
import pylab as plt

from PyQt4.Qt import QWidget

OSys=platform.system()
sys.path.append('Result')


try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class Ui_Dialog(QWidget):
    def setupUi(self,SpeciesList,RunnedPyrolPrList,PyrolModelDict):
        QWidget.__init__(self)  #manually added
        self.SpecL=SpeciesList   #manually added
        self.QWidg= QWidget(self) #manually added
        self.PyrolPrL=RunnedPyrolPrList
        self.PyrModelsD=PyrolModelDict
        self.QWidg.setObjectName(_fromUtf8("Dialog"))
        self.QWidg.resize(286, 253)
        self.B_SpecAnal = QtGui.QPushButton(self.QWidg)
        self.B_SpecAnal.setGeometry(QtCore.QRect(60, 150, 171, 24))
        self.B_SpecAnal.setObjectName(_fromUtf8("B_SpecAnal"))
        self.B_KinParams = QtGui.QPushButton(self.QWidg)
        self.B_KinParams.setGeometry(QtCore.QRect(60, 180, 171, 24))
        self.B_KinParams.setObjectName(_fromUtf8("B_KinParams"))
        self.layoutWidget = QtGui.QWidget(self.QWidg)
        self.layoutWidget.setGeometry(QtCore.QRect(90, 50, 101, 73))
        self.layoutWidget.setObjectName(_fromUtf8("layoutWidget"))
        self.verticalLayout_2 = QtGui.QVBoxLayout(self.layoutWidget)
        self.verticalLayout_2.setMargin(0)
        self.verticalLayout_2.setObjectName(_fromUtf8("verticalLayout_2"))
        self.l_calcFinished_2 = QtGui.QLabel(self.layoutWidget)
        self.l_calcFinished_2.setObjectName(_fromUtf8("l_calcFinished_2"))
        self.verticalLayout_2.addWidget(self.l_calcFinished_2)
        self.cB_selSpec = QtGui.QComboBox(self.layoutWidget)
        self.cB_selSpec.setObjectName(_fromUtf8("cB_selSpec"))
        for specNr in range(len(self.SpecL)):       #manually added
            self.cB_selSpec.addItem(_fromUtf8(""))  #manually added
        self.verticalLayout_2.addWidget(self.cB_selSpec)
        self.B_ShowResults = QtGui.QPushButton(self.layoutWidget)
        self.B_ShowResults.setObjectName(_fromUtf8("B_ShowResults"))
        self.verticalLayout_2.addWidget(self.B_ShowResults)
        self.l_calcFinished = QtGui.QLabel(self.QWidg)
        self.l_calcFinished.setGeometry(QtCore.QRect(10, 10, 271, 31))
        self.l_calcFinished.setObjectName(_fromUtf8("l_calcFinished"))
#        self.B_Quit = QtGui.QPushButton(self.QWidg)
#        self.B_Quit.setGeometry(QtCore.QRect(110, 220, 51, 24))
#        self.B_Quit.setObjectName(_fromUtf8("B_Quit"))

        self.retranslateUi()
#        QtCore.QObject.connect(self.B_Quit, QtCore.SIGNAL(_fromUtf8("clicked()")), self.QWidg.close)
        QtCore.QObject.connect(self.B_SpecAnal, QtCore.SIGNAL(_fromUtf8("clicked()")), self.OpenSpecAnalysis)  #manually added
        QtCore.QObject.connect(self.B_KinParams, QtCore.SIGNAL(_fromUtf8("clicked()")), self.OpenKinParam)  #manually added
        QtCore.QObject.connect(self.B_ShowResults, QtCore.SIGNAL(_fromUtf8("clicked()")), self.PlotSpec)  #manually added
        QtCore.QMetaObject.connectSlotsByName(self.QWidg)

    def retranslateUi(self):
        self.QWidg.setWindowTitle(QtGui.QApplication.translate("Dialog", "Dialog", None, QtGui.QApplication.UnicodeUTF8))
        self.B_SpecAnal.setText(QtGui.QApplication.translate("Dialog", "Open Species Analysis Files", None, QtGui.QApplication.UnicodeUTF8))
        self.B_KinParams.setText(QtGui.QApplication.translate("Dialog", "Open kinetic Parameter File", None, QtGui.QApplication.UnicodeUTF8))
        self.l_calcFinished_2.setText(QtGui.QApplication.translate("Dialog", "<html><head/><body><p align=\"center\">Select Species</p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
#        self.cB_selSpec.setItemText(0, QtGui.QApplication.translate("Dialog", "Total", None, QtGui.QApplication.UnicodeUTF8))
        for specNr in range(len(self.SpecL)):#manually added
            self.cB_selSpec.setItemText(specNr, QtGui.QApplication.translate("Dialog", self.SpecL[specNr], None, QtGui.QApplication.UnicodeUTF8))#manually added
        self.B_ShowResults.setText(QtGui.QApplication.translate("Dialog", "Show Results", None, QtGui.QApplication.UnicodeUTF8))
        self.l_calcFinished.setText(QtGui.QApplication.translate("Dialog", "<html><head/><body><p align=\"center\"><span style=\" font-size:14pt;\">Calculation finished</span></p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
#        self.B_Quit.setText(QtGui.QApplication.translate("Dialog", "Quit", None, QtGui.QApplication.UnicodeUTF8))

    def OpenSpecAnalysis(self):
        """Opens the species analysis file with Kate/Windows Editor"""
        if OSys=='Linux':
            for PyrolPr in self.PyrModelsD:
                if os.path.exists('Result/'+PyrolPr+'-BalanceResults3.txt'):
                    os.system('kate Result/'+PyrolPr+'-BalanceResults0.txt Result/'+PyrolPr+'-BalanceResults1.txt Result/'+PyrolPr+'-BalanceResults2.txt Result/'+PyrolPr+'-BalanceResults3.txt')
                elif os.path.exists('Result/'+PyrolPr+'BalanceResults2.txt'):
                    os.system('kate Result/'+PyrolPr+'-BalanceResults0.txt Result/'+PyrolPr+'-BalanceResults1.txt Result/'+PyrolPr+'-BalanceResults2.txt')
                elif os.path.exists('Result/'+PyrolPr+'-BalanceResults1.txt'):
                    os.system('kate Result/'+PyrolPr+'-BalanceResults0.txt Result/'+PyrolPr+'-BalanceResults1.txt')
                elif os.path.exists('Result/'+PyrolPr+'-BalanceResults0.txt'):
                    os.system('kate Result/'+PyrolPr+'-BalanceResults0.txt')
        elif OSys=='Windows':
            for PyrolPr in self.PyrModelsD:
                if os.path.exists('Result\\'+PyrolPr+'-BalanceResults0.txt'):
                    os.system('Result\\'+PyrolPr+'-BalanceResults0.txt')
                if os.path.exists('Result\\'+PyrolPr+'-BalanceResults1.txt'):
                    os.system('Result\\'+PyrolPr+'-BalanceResults1.txt')
                if os.path.exists('Result'+PyrolPr+'-BalanceResults2.txt'):
                    os.system('Result\\'+PyrolPr+'-BalanceResults2.txt')
                if os.path.exists('Result\\'+PyrolPr+'-BalanceResults3.txt'):
                    os.system('Result\\'+PyrolPr+'-BalanceResults3.txt')

    
    def OpenKinParam(self):
        """Opens the kinetic analysis file with Kate/Windows Editor"""
        files=''
        if OSys=='Linux':
            for PyrolPr in self.PyrModelsD:
                if os.path.exists('Result/'+PyrolPr+'-Results_'+self.PyrModelsD[PyrolPr]+'.txt'):
                    files+='Result/'+PyrolPr+'-Results_'+self.PyrModelsD[PyrolPr]+'.txt '
            os.system('kate '+files)
        elif OSys=='Windows':
            for PyrolPr in self.PyrModelsD:
                if os.path.exists('Result\\'+PyrolPr+'-Results_'+self.PyrModelsD[PyrolPr]+'.txt'):
                    files='Result\\'+PyrolPr+'-Results_'+self.PyrModelsD[PyrolPr]+'.txt'
                    os.system(files)
        
    def PlotSpec(self):
        """Plots the results for all Programs and their fits."""
        SpecNr=self.cB_selSpec.currentIndex()
#        plt.clf(),plt.cla()
#        fig1 = plt.figure()
#        figSpec = fig1.add_subplot(111)
        plt.ylabel('Yields (mass fraction)')
        plt.xlabel('Time in s')
        plt.legend(loc='4')
        plt.grid()
        colors=['g','b','r','purple','c','m','y','b','purple','r','b','g','c','y','m']
        colorIndex=0
        Spec=self.SpecL[SpecNr]
        if OSys=='Linux':
            for PyrolPr in self.PyrModelsD:
                for i in range(4):
                    print 'Result/'+PyrolPr+'-Fit_result_'+Spec+str(i)+'.out'
                    if os.path.exists('Result/'+PyrolPr+'-Fit_result_'+Spec+str(i)+'.out'):
                        Y=np.genfromtxt('Result/'+PyrolPr+'-Fit_result_'+Spec+str(i)+'.out',skip_header=1)
                        if np.shape(Y)[1]==6:
                            plt.plot(Y[:,0],Y[:,2],'--',color=colors[colorIndex],label=PyrolPr+' fit '+str(i))
                            plt.plot(Y[:,0],Y[:,4],'-',color=colors[colorIndex],label=PyrolPr+' original '+str(i))
                            colorIndex+=1
                        elif np.shape(Y)[1]==4:
                            plt.plot(Y[:,0],Y[:,2],'-',color=colors[colorIndex],label=PyrolPr+' original '+str(i))
                            colorIndex+=1
        elif OSys=='Windows':
            for PyrolPr in self.PyrModelsD:
                for i in range(4):
                    if os.path.exists('Result\\'+PyrolPr+'-Fit_result_'+Spec+str(i)+'.out'):
                        Y=np.genfromtxt('Result\\'+PyrolPr+'-Fit_result_'+Spec+str(i)+'.out')
                        if np.shape(Y)[1]==6:
                            plt.plot(Y[:,0],Y[:,2],'--',color=colors[colorIndex],label=PyrolPr+' fit '+str(i))
                            plt.plot(Y[:,0],Y[:,4],'-',color=colors[colorIndex],label=PyrolPr+' original '+str(i))
                            colorIndex+=1
                        elif np.shape(Y)[1]==4:
                            plt.plot(Y[:,0],Y[:,2],'-',color=colors[colorIndex],label=PyrolPr+' original '+str(i))
                            colorIndex+=1
        plt.show()
        #
            
            

#SpecList=["Total","Tar","Gas"]
#RunnedPyrolPrL=["CPD"]#,"FGDVC"]
#PyrolModelD={'CPD':'ArrheniusRate'}#,'FGDVC':'DAEM'}
#if __name__ == "__main__":
#    import sys
#    app = QtGui.QApplication(sys.argv)
#    Dialog = QtGui.QDialog()
#    ui = Ui_Dialog()
#    ui.setupUi(Dialog,SpecList,RunnedPyrolPrL,PyrolModelD)
#    Dialog.show()
#    sys.exit(app.exec_())

