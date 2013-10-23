import os
import sys
import platform
os.environ['QT_API'] = 'pyside'
import matplotlib
matplotlib.use('Qt4Agg')
matplotlib.rcParams['backend.qt4']='PySide'
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure

#from PySide import QtCore#, QtGui
#from PySide import *
#try:
#    _fromUtf8 = fromUtf8
#except AttributeError:
_fromUtf8 = lambda s: s


from PySide.QtCore import *
from PySide.QtGui import *

OSys=platform.system()



#class MainWindow(QMainWindow):
class Ui_Dialog(QWidget):#QMainWindow):
    #def __init__(self):
    def setupUi(self, SpeciesList,RunnedPyrolPrList,PyrolModelDict, NrOfRuns):
	# main stuff
        self.SpecL = SpeciesList   #manually added
        self.PyrolPrL = RunnedPyrolPrList
        self.PyrModelsD = PyrolModelDict
        self.NrRuns = NrOfRuns
        self.tab_Main = QTabWidget()
        self.tab_Main.setObjectName(_fromUtf8("tab_Main"))
        vbox = QVBoxLayout()
        vbox.addWidget(self.tab_Main)
        self.setLayout(vbox)
        self.setGeometry(900, 500, 250, 150)
        self.showMaximized()
        self.tab_ResKin = QWidget()
        self.tab_ResKin.setObjectName(_fromUtf8("tab_ResKin"))
        self.tab_Main.addTab(self.tab_ResKin, _fromUtf8(""))
        self.tab_ResSpec = QWidget()
        self.tab_ResSpec.setObjectName(_fromUtf8("tab_ResSpec"))
        self.tab_Main.addTab(self.tab_ResSpec, _fromUtf8(""))
        self.tab_Plots = QWidget()
        self.tab_Plots.setObjectName(_fromUtf8("tab_Plots"))
        self.tab_Main.addTab(self.tab_Plots, _fromUtf8(""))
        # kinetic tab
        self.gridLayout_2 = QGridLayout(self.tab_ResKin)
        self.gridLayout_2.setObjectName(_fromUtf8("gridLayout_2"))
        self.cB_Kin = QComboBox(self.tab_ResKin)
        self.cB_Kin.setObjectName(_fromUtf8("cB_Kin"))
        for PyrolPr in (self.PyrolPrL):       #manually added
            self.cB_Kin.addItem(_fromUtf8(""))  #manually added
        self.tE_Kin = QTextEdit(self.tab_ResKin)
        self.tE_Kin.setObjectName(_fromUtf8("tE_Kin"))
        self.l_NRKin = QLabel(self.tab_ResKin)
        self.l_NRKin.setObjectName(_fromUtf8("l_NRKin"))
        self.B_Kin = QPushButton(self.tab_ResKin)
        self.B_Kin.setObjectName(_fromUtf8("B_Kin"))
        self.gridLayout_2.addWidget(self.tE_Kin, 0, 0, 10, 6)
        self.gridLayout_2.addWidget(self.cB_Kin, 4, 10, 1, 1)
        self.gridLayout_2.addWidget(self.l_NRKin, 3, 10, 1, 1)
        self.gridLayout_2.addWidget(self.B_Kin, 6, 10, 1, 1)
        # species calc tab
        self.gridLayout_4 = QGridLayout(self.tab_ResSpec)
        self.gridLayout_4.setObjectName(_fromUtf8("gridLayout_4"))
        self.tE_Spec = QTextEdit(self.tab_ResSpec)
        self.tE_Spec.setObjectName(_fromUtf8("tE_Spec"))
        self.l_NRSpec = QLabel(self.tab_ResSpec)
        self.l_NRSpec.setObjectName(_fromUtf8("l_NRSpec"))
        self.sB_Spec = QSpinBox(self.tab_ResSpec)
        self.sB_Spec.setObjectName(_fromUtf8("sB_Spec"))
        self.sB_Spec.setMinimum(1)
        self.sB_Spec.setMaximum(self.NrRuns)
        self.cB_Spec = QComboBox(self.tab_ResSpec)
        self.cB_Spec.setObjectName(_fromUtf8("cB_Spec"))
        for PyrolPr in range(len(self.PyrolPrL)):       #manually added
            self.cB_Spec.addItem(_fromUtf8(""))  #manually added
        self.B_Spec = QPushButton(self.tab_ResSpec)
        self.B_Spec.setObjectName(_fromUtf8("B_Spec"))
        self.gridLayout_4.addWidget(self.tE_Spec, 0, 0, 10, 8)
        self.gridLayout_4.addWidget(self.l_NRSpec, 3, 10, 1, 1)
        self.gridLayout_4.addWidget(self.sB_Spec, 4, 10, 1, 1)
        self.gridLayout_4.addWidget(self.cB_Spec, 5, 10, 1, 1)
        self.gridLayout_4.addWidget(self.B_Spec, 7, 10, 1, 1)
        # Plot tab
        self.l_Plot = QLabel(self.tab_Plots)
        self.l_Plot.setObjectName(_fromUtf8("l_Plot"))
        self.cB_Plot = QComboBox(self.tab_Plots)
        self.cB_Plot.setObjectName(_fromUtf8("cB_Plot"))
        for SpecNr in range(len(self.SpecL)):       #manually added
            self.cB_Plot.addItem(_fromUtf8(""))  #manually added
        self.B_Plot = QPushButton(self.tab_Plots)
        self.B_Plot.setObjectName(_fromUtf8("B_Plot"))
        #plotting
        self.fig = Figure()
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.tab_Plots) 
        self.ax = self.fig.add_subplot(111)
        self.gridLayout_3 = QGridLayout(self.tab_Plots)
        self.gridLayout_3.setObjectName(_fromUtf8("gridLayout_3"))
        self.gridLayout_3.addWidget(self.canvas, 0, 0, 10, 10)
        self.gridLayout_3.addWidget(self.cB_Plot, 4, 10, 1, 1)
        self.gridLayout_3.addWidget(self.l_Plot, 3, 10, 1, 1)
        self.gridLayout_3.addWidget(self.B_Plot, 6, 10, 1, 1)
        self.retranslateUi()
        self.tab_Main.setCurrentIndex(1)
        QObject.connect(self.B_Spec, SIGNAL(_fromUtf8("clicked()")), self.OpenSpecAnalysis)  #manually added
        QObject.connect(self.B_Kin, SIGNAL(_fromUtf8("clicked()")), self.OpenKinParam)  #manually added
        QObject.connect(self.B_Plot, SIGNAL('clicked()'), self.PlotFunc)   #manually added
        QMetaObject.connectSlotsByName(self)
        # to close window
        exitAction = QAction('&Exit', self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.triggered.connect(self.close)


    def retranslateUi(self):
        self.setWindowTitle(QApplication.translate("Dialog", "Pyrolysis Kinetic Preprocessor - Results", None, QApplication.UnicodeUTF8))
        self.l_NRKin.setText(QApplication.translate("Dialog", "Program", None, QApplication.UnicodeUTF8))
        self.B_Kin.setText(QApplication.translate("Dialog", "Open", None, QApplication.UnicodeUTF8))
        self.tab_Main.setTabText(self.tab_Main.indexOf(self.tab_ResKin), QApplication.translate("Dialog", "Results - Kinetics", None, QApplication.UnicodeUTF8))
        self.l_NRSpec.setText(QApplication.translate("Dialog", "Run,Program", None, QApplication.UnicodeUTF8))
        self.B_Spec.setText(QApplication.translate("Dialog", "Open", None, QApplication.UnicodeUTF8))
        self.tab_Main.setTabText(self.tab_Main.indexOf(self.tab_ResSpec), QApplication.translate("Dialog", "Results - Species", None, QApplication.UnicodeUTF8))
        self.l_Plot.setText(QApplication.translate("Dialog", "<html><head/><body><p align=\"center\">Select Species</p></body></html>", None, QApplication.UnicodeUTF8))
        self.cB_Plot.setItemText(0, QApplication.translate("Dialog", "Total", None, QApplication.UnicodeUTF8))
        self.B_Plot.setText(QApplication.translate("Dialog", "Show Results", None, QApplication.UnicodeUTF8))
        self.tab_Main.setTabText(self.tab_Main.indexOf(self.tab_Plots), QApplication.translate("Dialog", "Plots", None, QApplication.UnicodeUTF8))
        for specNr in range(len(self.SpecL)):#manually added
            self.cB_Plot.setItemText(specNr, QApplication.translate("Dialog", self.SpecL[specNr], None, QApplication.UnicodeUTF8))#manually added
        for PyrolPr in range(len(self.PyrolPrL)):       #manually added
            self.cB_Kin.setItemText(PyrolPr, QApplication.translate("Dialog", self.PyrolPrL[PyrolPr], None, QApplication.UnicodeUTF8))#manually added
            self.cB_Spec.setItemText(PyrolPr, QApplication.translate("Dialog", self.PyrolPrL[PyrolPr], None, QApplication.UnicodeUTF8))#manually added


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
        print 'try open', KinFileName
        KinFile=open(KinFileName,'r')
        data = KinFile.read()
        self.tE_Kin.setText(data)
        KinFile.close()

    def PlotFunc(self):
        colors    = ['g','b','r','k']
        styles    = ['d','x','1','+','o']
        linewidths= [0.7, 1., 1., 1., 0.7 ] #size of the symbols in styles
        Index1=0
        Index2=0
        #
	#plotting
        SpecNr=self.cB_Plot.currentIndex()
        Spec=self.SpecL[SpecNr]
        self.ax.clear()
        self.ax.grid()
        self.ax.set_ylabel('Yields (mass fraction)')
        #self.ax.set_xlabel('Time in s')
        self.ax.set_xlabel('Temp in K')
        self.ax.set_title(Spec)
        if OSys=='Linux':
            for PyrolPr in self.PyrModelsD:
                Index2=0
                Index1+=1
                for i in range(self.NrRuns):
                    if os.path.exists('Result/'+PyrolPr+'-Fit_result_'+Spec+'_'+str(i)+'.out'):
                        Y=np.genfromtxt('Result/'+PyrolPr+'-Fit_result_'+Spec+'_'+str(i)+'.out',skip_header=1)
                        if np.shape(Y)[1]==6:
                            self.ax.plot(Y[:,0],Y[:,2],'-',color=colors[Index1],linewidth=1.5)
                            self.ax.plot(Y[:,0],Y[:,4],styles[Index2],color=colors[Index1],linewidth=linewidths[Index2],label=PyrolPr+' original '+str(i))
                            Index2+=1
                        elif np.shape(Y)[1]==4:
                            self.ax.plot(Y[:,0],Y[:,2],styles[Index2],color=colors[Index1],linewidth=linewidths[Index2],label=PyrolPr+' original '+str(i))
                            Index2+=1
        elif OSys=='Windows':
            for PyrolPr in self.PyrModelsD:
                Index2=0
                Index1+=1
                for i in range(self.NrRuns):
                    if os.path.exists('Result\\'+PyrolPr+'-Fit_result_'+Spec+'_'+str(i)+'.out'):
                        Y=np.genfromtxt('Result\\'+PyrolPr+'-Fit_result_'+Spec+'_'+str(i)+'.out',skip_header=1)
                        if np.shape(Y)[1]==6:
                            self.ax.plot(Y[:,0],Y[:,2],'-',color=colors[Index1],linewidth=1.5)
                            self.ax.plot(Y[:,0],Y[:,4],styles[Index2],color=colors[Index1],linewidth=linewidths[Index2],label=PyrolPr+' original '+str(i))
                            Index2+=1
                        elif np.shape(Y)[1]==4:
                            self.ax.plot(Y[:,0],Y[:,2],styles[Index2],color=colors[Index1],linewidth=linewidths[Index2],label=PyrolPr+' original '+str(i))
                            Index2+=1
        self.ax.legend(loc='lower right')#'4')
        self.canvas.draw()




