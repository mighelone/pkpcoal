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

from PySide import QtCore, QtGui


try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

OSys=platform.system()



#class MainWindow(QMainWindow):
class Ui_Dialog(QtGui.QWidget):#QMainWindow):
    #def __init__(self):
    def setupUi(self, SpeciesList,RunnedPyrolPrList,PyrolModelDict, NrOfRuns):
	# main stuff
        self.SpecL = SpeciesList   #manually added
        self.PyrolPrL = RunnedPyrolPrList
        self.PyrModelsD = PyrolModelDict
        self.NrRuns = NrOfRuns
        self.tab_Main = QtGui.QTabWidget()
        self.tab_Main.setObjectName(_fromUtf8("tab_Main"))
	vbox = QtGui.QVBoxLayout()
        vbox.addWidget(self.tab_Main)
        self.setLayout(vbox)
	self.setGeometry(900, 500, 250, 150)
	self.showMaximized()
        self.tab_ResKin = QtGui.QWidget()
        self.tab_ResKin.setObjectName(_fromUtf8("tab_ResKin"))
        self.tab_Main.addTab(self.tab_ResKin, _fromUtf8(""))
        self.tab_ResSpec = QtGui.QWidget()
        self.tab_ResSpec.setObjectName(_fromUtf8("tab_ResSpec"))
        self.tab_Main.addTab(self.tab_ResSpec, _fromUtf8(""))
        self.tab_Plots = QtGui.QWidget()
        self.tab_Plots.setObjectName(_fromUtf8("tab_Plots"))
        self.tab_Main.addTab(self.tab_Plots, _fromUtf8(""))
	# kinetic tab
        self.gridLayout_2 = QtGui.QGridLayout(self.tab_ResKin)
        self.gridLayout_2.setObjectName(_fromUtf8("gridLayout_2"))
        self.cB_Kin = QtGui.QComboBox(self.tab_ResKin)
        self.cB_Kin.setObjectName(_fromUtf8("cB_Kin"))
        for PyrolPr in range(len(self.PyrolPrL)):       #manually added
            self.cB_Kin.addItem(_fromUtf8(""))  #manually added
        self.tE_Kin = QtGui.QTextEdit(self.tab_ResKin)
        self.tE_Kin.setObjectName(_fromUtf8("tE_Kin"))
        self.l_NRKin = QtGui.QLabel(self.tab_ResKin)
        self.l_NRKin.setObjectName(_fromUtf8("l_NRKin"))
        self.B_Kin = QtGui.QPushButton(self.tab_ResKin)
        self.B_Kin.setObjectName(_fromUtf8("B_Kin"))
        self.gridLayout_2.addWidget(self.tE_Kin, 0, 0, 10, 6)
        self.gridLayout_2.addWidget(self.cB_Kin, 4, 10, 1, 1)
        self.gridLayout_2.addWidget(self.l_NRKin, 3, 10, 1, 1)
        self.gridLayout_2.addWidget(self.B_Kin, 6, 10, 1, 1)
	# species calc tab
        self.gridLayout_4 = QtGui.QGridLayout(self.tab_ResSpec)
        self.gridLayout_4.setObjectName(_fromUtf8("gridLayout_4"))
        self.tE_Spec = QtGui.QTextEdit(self.tab_ResSpec)
        self.tE_Spec.setObjectName(_fromUtf8("tE_Spec"))
        self.l_NRSpec = QtGui.QLabel(self.tab_ResSpec)
        self.l_NRSpec.setObjectName(_fromUtf8("l_NRSpec"))
        self.sB_Spec = QtGui.QSpinBox(self.tab_ResSpec)
        self.sB_Spec.setObjectName(_fromUtf8("sB_Spec"))
	self.sB_Spec.setMinimum(1)
        self.sB_Spec.setMaximum(self.NrRuns)
        self.cB_Spec = QtGui.QComboBox(self.tab_ResSpec)
        self.cB_Spec.setObjectName(_fromUtf8("cB_Spec"))
        for PyrolPr in range(len(self.PyrolPrL)):       #manually added
            self.cB_Spec.addItem(_fromUtf8(""))  #manually added
        self.B_Spec = QtGui.QPushButton(self.tab_ResSpec)
        self.B_Spec.setObjectName(_fromUtf8("B_Spec"))
        self.gridLayout_4.addWidget(self.tE_Spec, 0, 0, 10, 8)
        self.gridLayout_4.addWidget(self.l_NRSpec, 3, 10, 1, 1)
        self.gridLayout_4.addWidget(self.sB_Spec, 4, 10, 1, 1)
        self.gridLayout_4.addWidget(self.cB_Spec, 5, 10, 1, 1)
        self.gridLayout_4.addWidget(self.B_Spec, 7, 10, 1, 1)
	# Plot tab
        self.l_Plot = QtGui.QLabel(self.tab_Plots)
        self.l_Plot.setObjectName(_fromUtf8("l_Plot"))
        self.cB_Plot = QtGui.QComboBox(self.tab_Plots)
        self.cB_Plot.setObjectName(_fromUtf8("cB_Plot"))
        for SpecNr in range(len(self.SpecL)):       #manually added
            self.cB_Plot.addItem(_fromUtf8(""))  #manually added
        self.B_Plot = QtGui.QPushButton(self.tab_Plots)
        self.B_Plot.setObjectName(_fromUtf8("B_Plot"))
	#plotting
        self.fig = Figure()
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self)
        self.ax = self.fig.add_subplot(111)
	self.gridLayout_3 = QtGui.QGridLayout(self.tab_Plots)
        self.gridLayout_3.setObjectName(_fromUtf8("gridLayout_3"))
	self.gridLayout_3.addWidget(self.canvas, 0, 0, 10, 10)
	self.gridLayout_3.addWidget(self.cB_Plot, 4, 10, 1, 1)
	self.gridLayout_3.addWidget(self.l_Plot, 3, 10, 1, 1)
	self.gridLayout_3.addWidget(self.B_Plot, 6, 10, 1, 1)
        self.retranslateUi()
        self.tab_Main.setCurrentIndex(1)
        QtCore.QObject.connect(self.B_Spec, QtCore.SIGNAL(_fromUtf8("clicked()")), self.OpenSpecAnalysis)  #manually added
        QtCore.QObject.connect(self.B_Kin, QtCore.SIGNAL(_fromUtf8("clicked()")), self.OpenKinParam)  #manually added
        QtCore.QObject.connect(self.B_Plot, QtCore.SIGNAL('clicked()'), self.PlotFunc)   #manually added
        QtCore.QMetaObject.connectSlotsByName(self)
	# to close window
	exitAction = QtGui.QAction('&Exit', self)
        exitAction.setShortcut('Ctrl+Q')
	exitAction.triggered.connect(self.close)


    def retranslateUi(self):
        self.setWindowTitle(QtGui.QApplication.translate("Dialog", "Dialog", None, QtGui.QApplication.UnicodeUTF8))
        self.l_NRKin.setText(QtGui.QApplication.translate("Dialog", "Program", None, QtGui.QApplication.UnicodeUTF8))
        self.B_Kin.setText(QtGui.QApplication.translate("Dialog", "Open", None, QtGui.QApplication.UnicodeUTF8))
        self.tab_Main.setTabText(self.tab_Main.indexOf(self.tab_ResKin), QtGui.QApplication.translate("Dialog", "Results - Kinetics", None, QtGui.QApplication.UnicodeUTF8))
        self.l_NRSpec.setText(QtGui.QApplication.translate("Dialog", "Run,Program", None, QtGui.QApplication.UnicodeUTF8))
        self.B_Spec.setText(QtGui.QApplication.translate("Dialog", "Open", None, QtGui.QApplication.UnicodeUTF8))
        self.tab_Main.setTabText(self.tab_Main.indexOf(self.tab_ResSpec), QtGui.QApplication.translate("Dialog", "Results - Species", None, QtGui.QApplication.UnicodeUTF8))
        self.l_Plot.setText(QtGui.QApplication.translate("Dialog", "<html><head/><body><p align=\"center\">Select Species</p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.cB_Plot.setItemText(0, QtGui.QApplication.translate("Dialog", "Total", None, QtGui.QApplication.UnicodeUTF8))
        self.B_Plot.setText(QtGui.QApplication.translate("Dialog", "Show Results", None, QtGui.QApplication.UnicodeUTF8))
        self.tab_Main.setTabText(self.tab_Main.indexOf(self.tab_Plots), QtGui.QApplication.translate("Dialog", "Plots", None, QtGui.QApplication.UnicodeUTF8))
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
	#plotting
        SpecNr=self.cB_Plot.currentIndex()
        Spec=self.SpecL[SpecNr]
        self.ax.clear()
        self.ax.grid()
	self.ax.set_ylabel('Yields (mass fraction)')
        self.ax.set_xlabel('Time in s')
        self.ax.set_title(Spec)
        if OSys=='Linux':
            for PyrolPr in self.PyrModelsD:
                for i in range(self.NrRuns):
                    if os.path.exists('Result/'+PyrolPr+'-Fit_result_'+Spec+str(i)+'.out'):
                        Y=np.genfromtxt('Result/'+PyrolPr+'-Fit_result_'+Spec+str(i)+'.out',skip_header=1)
                        if np.shape(Y)[1]==6:
                            self.ax.plot(Y[:,0],Y[:,2],'--',color=colors[colorIndex],label=PyrolPr+' fit '+str(i))
                            self.ax.plot(Y[:,0],Y[:,4],'-',color=colors[colorIndex],label=PyrolPr+' original '+str(i))
                            colorIndex+=1
                        elif np.shape(Y)[1]==4:
                            self.ax.plot(Y[:,0],Y[:,2],'-',color=colors[colorIndex],label=PyrolPr+' original '+str(i))
                            colorIndex+=1
        elif OSys=='Windows':
            for PyrolPr in self.PyrModelsD:
                for i in range(self.NrRuns):
                    if os.path.exists('Result\\'+PyrolPr+'-Fit_result_'+Spec+str(i)+'.out'):
                        Y=np.genfromtxt('Result\\'+PyrolPr+'-Fit_result_'+Spec+str(i)+'.out',skip_header=1)
                        if np.shape(Y)[1]==6:
                            self.ax.plot(Y[:,0],Y[:,2],'--',color=colors[colorIndex],label=PyrolPr+' fit '+str(i))
                            self.ax.plot(Y[:,0],Y[:,4],'-',color=colors[colorIndex],label=PyrolPr+' original '+str(i))
                            colorIndex+=1
                        elif np.shape(Y)[1]==4:
                            self.ax.plot(Y[:,0],Y[:,2],'-',color=colors[colorIndex],label=PyrolPr+' original '+str(i))
                            colorIndex+=1
        self.ax.legend(loc='lower right')#'4')
        self.canvas.draw()





if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    myapp = GUIForm()
    myapp.show()
    sys.exit(app.exec_())


