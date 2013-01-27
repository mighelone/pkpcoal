# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'Done2.ui'

import sys,os
sys.path.append('src')
from PyQt4 import QtCore, QtGui
from matplotlibwidgetFile import matplotlibWidget
import numpy as np
import platform
from PyQt4.Qt import QWidget

OSys=platform.system()
sys.path.append('Result')




try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class Ui_Dialog(QWidget):
    def setupUi(self, SpeciesList,RunnedPyrolPrList,PyrolModelDict, NrOfRuns):
	QWidget.__init__(self)  #manually added
	self.QWidg= QWidget(self) #manually added
	self.SpecL=SpeciesList   #manually added
	self.PyrolPrL=RunnedPyrolPrList
        self.PyrModelsD=PyrolModelDict
        self.NrRuns=NrOfRuns
        #Dialog.setObjectName(_fromUtf8("Dialog"))
        #Dialog.resize(820, 479)
        self.QWidg.setObjectName(_fromUtf8("Dialog"))
        self.QWidg.resize(820, 479)
        self.gridLayout = QtGui.QGridLayout(self.QWidg)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        #self.B_Quit = QtGui.QPushButton(self.QWidg)
        #self.B_Quit.setObjectName(_fromUtf8("B_Quit"))
        #self.gridLayout.addWidget(self.B_Quit, 2, 0, 1, 1)
        self.Tab_Main = QtGui.QTabWidget(self.QWidg)
        self.Tab_Main.setObjectName(_fromUtf8("Tab_Main"))
        self.tab_ResKin = QtGui.QWidget()
        self.tab_ResKin.setObjectName(_fromUtf8("tab_ResKin"))
        self.layoutWidget = QtGui.QWidget(self.tab_ResKin)
        self.layoutWidget.setGeometry(QtCore.QRect(11, 11, 771, 391))
        self.layoutWidget.setObjectName(_fromUtf8("layoutWidget"))
        self.gridLayout_2 = QtGui.QGridLayout(self.layoutWidget)
        self.gridLayout_2.setMargin(0)
        self.gridLayout_2.setObjectName(_fromUtf8("gridLayout_2"))
        self.cB_Kin = QtGui.QComboBox(self.layoutWidget)
        self.cB_Kin.setObjectName(_fromUtf8("cB_Kin"))
        for PyrolPr in range(len(self.PyrolPrL)):       #manually added
            self.cB_Kin.addItem(_fromUtf8(""))  #manually added
        #self.sB_Kin = QtGui.QSpinBox(self.layoutWidget)
        #self.sB_Kin.setObjectName(_fromUtf8("sB_Kin"))
        #self.gridLayout_2.addWidget(self.sB_Kin, 2, 1, 1, 1)
        self.gridLayout_2.addWidget(self.cB_Kin, 2, 1, 1, 1)
        self.l_NRKin = QtGui.QLabel(self.layoutWidget)
        self.l_NRKin.setObjectName(_fromUtf8("l_NRKin"))
        self.gridLayout_2.addWidget(self.l_NRKin, 1, 1, 1, 1)
        self.tE_Kin = QtGui.QTextEdit(self.layoutWidget)
        self.tE_Kin.setObjectName(_fromUtf8("tE_Kin"))
        self.gridLayout_2.addWidget(self.tE_Kin, 0, 0, 5, 1)
        spacerItem = QtGui.QSpacerItem(20, 158, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.gridLayout_2.addItem(spacerItem, 0, 1, 1, 1)
        spacerItem1 = QtGui.QSpacerItem(20, 158, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.gridLayout_2.addItem(spacerItem1, 4, 1, 1, 1)
        self.B_Kin = QtGui.QPushButton(self.layoutWidget)
        self.B_Kin.setObjectName(_fromUtf8("B_Kin"))
        self.gridLayout_2.addWidget(self.B_Kin, 3, 1, 1, 1)
        self.Tab_Main.addTab(self.tab_ResKin, _fromUtf8(""))
        self.tab_ResSpec = QtGui.QWidget()
        self.tab_ResSpec.setObjectName(_fromUtf8("tab_ResSpec"))
        self.layoutWidget_2 = QtGui.QWidget(self.tab_ResSpec)
        self.layoutWidget_2.setGeometry(QtCore.QRect(10, 10, 771, 391))
        self.layoutWidget_2.setObjectName(_fromUtf8("layoutWidget_2"))
        self.gridLayout_4 = QtGui.QGridLayout(self.layoutWidget_2)
        self.gridLayout_4.setMargin(0)
        self.gridLayout_4.setObjectName(_fromUtf8("gridLayout_4"))
        self.tE_Spec = QtGui.QTextEdit(self.layoutWidget_2)
        self.tE_Spec.setObjectName(_fromUtf8("tE_Spec"))
        self.gridLayout_4.addWidget(self.tE_Spec, 0, 0, 5, 1)
        spacerItem2 = QtGui.QSpacerItem(20, 158, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.gridLayout_4.addItem(spacerItem2, 0, 1, 1, 1)
        spacerItem3 = QtGui.QSpacerItem(20, 158, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.gridLayout_4.addItem(spacerItem3, 4, 1, 1, 1)
        self.l_NRSpec = QtGui.QLabel(self.layoutWidget_2)
        self.l_NRSpec.setObjectName(_fromUtf8("l_NRSpec"))
        self.gridLayout_4.addWidget(self.l_NRSpec, 1, 1, 1, 1)
        self.sB_Spec = QtGui.QSpinBox(self.layoutWidget_2)
        self.sB_Spec.setObjectName(_fromUtf8("sB_Spec"))
        self.gridLayout_4.addWidget(self.sB_Spec, 2, 1, 1, 1)
        self.cB_Spec = QtGui.QComboBox(self.layoutWidget)
        self.cB_Spec.setObjectName(_fromUtf8("cB_Spec"))
        for PyrolPr in range(len(self.PyrolPrL)):       #manually added
            self.cB_Spec.addItem(_fromUtf8(""))  #manually added
        self.gridLayout_4.addWidget(self.cB_Spec, 3, 1, 1, 1)
        self.B_Spec = QtGui.QPushButton(self.layoutWidget_2)
        self.B_Spec.setObjectName(_fromUtf8("B_Spec"))
        self.gridLayout_4.addWidget(self.B_Spec, 4, 1, 1, 1)
        self.Tab_Main.addTab(self.tab_ResSpec, _fromUtf8(""))
        self.tab_Plots = QtGui.QWidget()
        self.tab_Plots.setObjectName(_fromUtf8("tab_Plots"))
        self.widget = matplotlibWidget(self.tab_Plots)
        self.widget.setGeometry(QtCore.QRect(0, 0, 641, 381))
        self.widget.setObjectName(_fromUtf8("widget"))
        self.layoutWidget1 = QtGui.QWidget(self.tab_Plots)
        self.layoutWidget1.setGeometry(QtCore.QRect(650, 0, 103, 399))
        self.layoutWidget1.setObjectName(_fromUtf8("layoutWidget1"))
        self.verticalLayout = QtGui.QVBoxLayout(self.layoutWidget1)
        self.verticalLayout.setMargin(0)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        spacerItem4 = QtGui.QSpacerItem(20, 158, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem4)
        self.verticalLayout_2 = QtGui.QVBoxLayout()
        self.verticalLayout_2.setObjectName(_fromUtf8("verticalLayout_2"))
        self.l_Plot = QtGui.QLabel(self.layoutWidget1)
        self.l_Plot.setObjectName(_fromUtf8("l_Plot"))
        self.verticalLayout_2.addWidget(self.l_Plot)
        self.cB_Plot = QtGui.QComboBox(self.layoutWidget1)
        self.cB_Plot.setObjectName(_fromUtf8("cB_Plot"))
        for SpecNr in range(len(self.SpecL)):       #manually added
            self.cB_Plot.addItem(_fromUtf8(""))  #manually added
        self.verticalLayout_2.addWidget(self.cB_Plot)
        self.B_Plot = QtGui.QPushButton(self.layoutWidget1)
        self.B_Plot.setObjectName(_fromUtf8("B_Plot"))
        self.verticalLayout_2.addWidget(self.B_Plot)
        self.verticalLayout.addLayout(self.verticalLayout_2)
        spacerItem5 = QtGui.QSpacerItem(20, 158, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem5)
        self.Tab_Main.addTab(self.tab_Plots, _fromUtf8(""))
        self.gridLayout.addWidget(self.Tab_Main, 1, 0, 1, 1)
        #added
        self.sB_Spec.setMinimum(1)
        self.sB_Spec.setMaximum(self.NrRuns)
        #end added
        self.retranslateUi()
        self.Tab_Main.setCurrentIndex(1)
        #QtCore.QObject.connect(self.B_Quit, QtCore.SIGNAL(_fromUtf8("clicked()")), self.QWidg.close)
        QtCore.QObject.connect(self.B_Spec, QtCore.SIGNAL(_fromUtf8("clicked()")), self.OpenSpecAnalysis)  #manually added
        QtCore.QObject.connect(self.B_Kin, QtCore.SIGNAL(_fromUtf8("clicked()")), self.OpenKinParam)  #manually added
        QtCore.QObject.connect(self.B_Plot, QtCore.SIGNAL('clicked()'), self.PlotFunc)   #manually added
        QtCore.QMetaObject.connectSlotsByName(self.QWidg)


    def retranslateUi(self):
        self.QWidg.setWindowTitle(QtGui.QApplication.translate("Dialog", "Dialog", None, QtGui.QApplication.UnicodeUTF8))
        #self.B_Quit.setText(QtGui.QApplication.translate("Dialog", "Quit", None, QtGui.QApplication.UnicodeUTF8))
        self.l_NRKin.setText(QtGui.QApplication.translate("Dialog", "Program", None, QtGui.QApplication.UnicodeUTF8))
        self.B_Kin.setText(QtGui.QApplication.translate("Dialog", "Open", None, QtGui.QApplication.UnicodeUTF8))
        self.Tab_Main.setTabText(self.Tab_Main.indexOf(self.tab_ResKin), QtGui.QApplication.translate("Dialog", "Results - Kinetics", None, QtGui.QApplication.UnicodeUTF8))
        self.l_NRSpec.setText(QtGui.QApplication.translate("Dialog", "Run,Program", None, QtGui.QApplication.UnicodeUTF8))
        self.B_Spec.setText(QtGui.QApplication.translate("Dialog", "Open", None, QtGui.QApplication.UnicodeUTF8))
        self.Tab_Main.setTabText(self.Tab_Main.indexOf(self.tab_ResSpec), QtGui.QApplication.translate("Dialog", "Results - Species", None, QtGui.QApplication.UnicodeUTF8))
        self.l_Plot.setText(QtGui.QApplication.translate("Dialog", "<html><head/><body><p align=\"center\">Select Species</p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.cB_Plot.setItemText(0, QtGui.QApplication.translate("Dialog", "Total", None, QtGui.QApplication.UnicodeUTF8))
        self.B_Plot.setText(QtGui.QApplication.translate("Dialog", "Show Results", None, QtGui.QApplication.UnicodeUTF8))
        self.Tab_Main.setTabText(self.Tab_Main.indexOf(self.tab_Plots), QtGui.QApplication.translate("Dialog", "Plots", None, QtGui.QApplication.UnicodeUTF8))
        for specNr in range(len(self.SpecL)):#manually added
            self.cB_Plot.setItemText(specNr, QtGui.QApplication.translate("Dialog", self.SpecL[specNr], None, QtGui.QApplication.UnicodeUTF8))#manually added
        for PyrolPr in range(len(self.PyrolPrL)):       #manually added
            self.cB_Kin.setItemText(PyrolPr, QtGui.QApplication.translate("Dialog", self.PyrolPrL[PyrolPr], None, QtGui.QApplication.UnicodeUTF8))#manually added
            self.cB_Spec.setItemText(PyrolPr, QtGui.QApplication.translate("Dialog", self.PyrolPrL[PyrolPr], None, QtGui.QApplication.UnicodeUTF8))#manually added
        
        
    def OpenSpecAnalysis(self):
        """Opens the species analysis file in the TextEditor."""
        Nr=str(self.sB_Spec.value()-1) #number from SpinBox
        PyrolPr=self.PyrolPrL[self.cB_Spec.currentIndex()]
        if OSys=='Linux':
	  if os.path.exists('Result/'+PyrolPr+'-BalanceResults'+Nr+'.txt'):
	    SpecFileName='Result/'+PyrolPr+'-BalanceResults'+Nr+'.txt'
	  else:
	    print 'File cannot be found: ','Result/',PyrolPr,'-BalanceResults',Nr,'.txt'
        elif OSys=='Windows':
	  if os.path.exists('Result\\'+PyrolPr+'-BalanceResults'+Nr+'.txt'):
	    SpecFileName='Result\\'+PyrolPr+'-BalanceResults'+Nr+'.txt'
	  else:
	    print 'File cannot be found: ','Result\\',PyrolPr,'-BalanceResults',Nr,'.txt'
        SpecFile=open(SpecFileName,'r')
        data = SpecFile.read()
        self.tE_Spec.setText(data)
        SpecFile.close()
    
    def OpenKinParam(self):
        """Opens the species analysis file in the TextEditor."""
        Nr=(self.cB_Kin.currentIndex()) #number from ColumnsBar
        PyrolPr=self.PyrolPrL[Nr]
        if OSys=='Linux':
	  if os.path.exists('Result/'+PyrolPr+'-Results_'+self.PyrModelsD[PyrolPr]+'.txt'):
	    KinFileName='Result/'+PyrolPr+'-Results_'+self.PyrModelsD[PyrolPr]+'.txt'
	  else:
	    print 'File cannot be found: ','Result/'+PyrolPr+'-Results_'+self.PyrModelsD[PyrolPr]+'.txt'
        elif OSys=='Windows':
	  if os.path.exists('Result\\'+PyrolPr+'-Results_'+self.PyrModelsD[PyrolPr]+'.txt'):
	    KinFileName='Result\\'+PyrolPr+'-Results_'+self.PyrModelsD[PyrolPr]+'.txt'
	  else:
	    print 'File cannot be found: ','Result\\'+PyrolPr+'-Results_'+self.PyrModelsD[PyrolPr]+'.txt'
        KinFile=open(KinFileName,'r')
        data = KinFile.read()
        self.tE_Kin.setText(data)
        KinFile.close()     
        
    def PlotFunc(self):
	colors=['g','b','r','purple','c','m','y','b','purple','r','b','g','c','y','m']
        colorIndex=0
        #
        SpecNr=self.cB_Plot.currentIndex()
        Spec=self.SpecL[SpecNr]
        self.widget.canvas.ax.clear()
        self.widget.canvas.ax.grid()
        self.widget.canvas.ax.set_ylabel('Yields (mass fraction)')
        self.widget.canvas.ax.set_xlabel('Time in s')
        self.widget.canvas.ax.set_title(Spec)
        if OSys=='Linux':
            for PyrolPr in self.PyrModelsD:
                for i in range(self.NrRuns):
                    if os.path.exists('Result/'+PyrolPr+'-Fit_result_'+Spec+str(i)+'.out'):
                        Y=np.genfromtxt('Result/'+PyrolPr+'-Fit_result_'+Spec+str(i)+'.out',skip_header=1)
                        if np.shape(Y)[1]==6:
                            self.widget.canvas.ax.plot(Y[:,0],Y[:,2],'--',color=colors[colorIndex],label=PyrolPr+' fit '+str(i))
                            self.widget.canvas.ax.plot(Y[:,0],Y[:,4],'-',color=colors[colorIndex],label=PyrolPr+' original '+str(i))
                            colorIndex+=1
                        elif np.shape(Y)[1]==4:
                            self.widget.canvas.ax.plot(Y[:,0],Y[:,2],'-',color=colors[colorIndex],label=PyrolPr+' original '+str(i))
                            colorIndex+=1
        elif OSys=='Windows':
            for PyrolPr in self.PyrModelsD:
                for i in range(self.NrRuns):
                    if os.path.exists('Result\\'+PyrolPr+'-Fit_result_'+Spec+str(i)+'.out'):
                        Y=np.genfromtxt('Result\\'+PyrolPr+'-Fit_result_'+Spec+str(i)+'.out',skip_header=1)
                        if np.shape(Y)[1]==6:
                            self.widget.canvas.ax.plot(Y[:,0],Y[:,2],'--',color=colors[colorIndex],label=PyrolPr+' fit '+str(i))
                            self.widget.canvas.ax.plot(Y[:,0],Y[:,4],'-',color=colors[colorIndex],label=PyrolPr+' original '+str(i))
                            colorIndex+=1
                        elif np.shape(Y)[1]==4:
                            self.widget.canvas.ax.plot(Y[:,0],Y[:,2],'-',color=colors[colorIndex],label=PyrolPr+' original '+str(i))
                            colorIndex+=1
        self.widget.canvas.ax.legend(loc='4')
        self.widget.canvas.draw()
        
        

        
        

if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    myapp = GUIForm()
    myapp.show()
    sys.exit(app.exec_())

