import os
os.environ['QT_API'] = 'pyside'
import matplotlib
matplotlib.use('Qt4Agg')
matplotlib.rcParams['backend.qt4']='PySide'
import pylab as plt
import numpy as np
import writeInfoFiles
import sys
import platform
import shutil

OSys=platform.system()
sys.path.append('src')
import PKP
import InformationFiles

import Done   #second GUI

from PySide.QtCore import * #QWidget, QMainWindow
from PySide.QtGui import *


#try:
#    _fromUtf8 = QtCore.QString.fromUtf8
#except AttributeError:
_fromUtf8 = lambda s: s

#class Ui_PKP(QMainWindow):
class Ui_PKP(QMainWindow):
    def setupUi(self, PKP, SaveInfoObj):
        #QMainWindow.__init__(self) #manually added
        self.QWidgMain= QWidget(self) #manually added
        self.setCentralWidget(self.QWidgMain) #manually added
        self.TheCalculationsAreDone=False #manually added
        self.svInfo=SaveInfoObj
        PKP.setObjectName(_fromUtf8("PKP"))
        PKP.resize(1033, 753)
        self.centralwidget = QWidget(PKP)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.Header1 = QLabel(self.centralwidget)
        self.Header1.setGeometry(QRect(-10, 50, 721, 31))
        self.Header1.setObjectName(_fromUtf8("Header1"))
        self.Header2 = QLabel(self.centralwidget)
        self.Header2.setGeometry(QRect(30, 230, 721, 31))
        self.Header2.setObjectName(_fromUtf8("Header2"))
        self.L_HHV = QLabel(self.centralwidget)
        self.L_HHV.setGeometry(QRect(570, 350, 191, 51))
        self.L_HHV.setObjectName(_fromUtf8("L_HHV"))
        self.lE_HHV = QLineEdit(self.centralwidget)
        self.lE_HHV.setGeometry(QRect(570, 400, 113, 23))
        self.lE_HHV.setObjectName(_fromUtf8("lE_HHV"))
        self.Header3 = QLabel(self.centralwidget)
        self.Header3.setGeometry(QRect(20, 520, 721, 31))
        self.Header3.setObjectName(_fromUtf8("Header3"))
        self.L_pressure = QLabel(self.centralwidget)
        self.L_pressure.setGeometry(QRect(30, 560, 141, 24))
        self.L_pressure.setObjectName(_fromUtf8("L_pressure"))
        self.lE_pressure = QLineEdit(self.centralwidget)
        self.lE_pressure.setGeometry(QRect(220, 560, 113, 23))
        self.lE_pressure.setObjectName(_fromUtf8("lE_pressure"))
        self.L_THist = QLabel(self.centralwidget)
        self.L_THist.setGeometry(QRect(30, 590, 191, 71))
        self.L_THist.setObjectName(_fromUtf8("L_THist"))
        self.tE_THist_1 = QTextEdit(self.centralwidget)
        self.tE_THist_1.setGeometry(QRect(220, 590, 111, 76))
        self.tE_THist_1.setObjectName(_fromUtf8("tE_THist_1"))
        self.layoutWidget = QWidget(self.centralwidget)
        self.layoutWidget.setGeometry(QRect(40, 100, 241, 91))
        self.layoutWidget.setObjectName(_fromUtf8("layoutWidget"))
        self.formLayout = QFormLayout(self.layoutWidget)
        self.formLayout.setFieldGrowthPolicy(QFormLayout.ExpandingFieldsGrow)
        #self.formLayout.setContentsMargins (0)#setMargin(0)
        self.formLayout.setObjectName(_fromUtf8("formLayout"))
        self.L_CPD = QLabel(self.layoutWidget)
        self.L_CPD.setObjectName(_fromUtf8("L_CPD"))
        self.formLayout.setWidget(0, QFormLayout.LabelRole, self.L_CPD)
        self.cB_CPD = QComboBox(self.layoutWidget)
        self.cB_CPD.setObjectName(_fromUtf8("cB_CPD"))
        self.cB_CPD.addItem(_fromUtf8(""))
        self.cB_CPD.addItem(_fromUtf8(""))
        self.cB_CPD.addItem(_fromUtf8(""))
        self.cB_CPD.addItem(_fromUtf8(""))
        self.cB_CPD.addItem(_fromUtf8(""))
        self.cB_CPD.addItem(_fromUtf8(""))
        self.cB_CPD.addItem(_fromUtf8(""))
        self.formLayout.setWidget(0, QFormLayout.FieldRole, self.cB_CPD)
        self.L_FGDVC = QLabel(self.layoutWidget)
        self.L_FGDVC.setObjectName(_fromUtf8("L_FGDVC"))
        self.formLayout.setWidget(1, QFormLayout.LabelRole, self.L_FGDVC)
        self.cB_FGDVC = QComboBox(self.layoutWidget)
        self.cB_FGDVC.setObjectName(_fromUtf8("cB_FGDVC"))
        self.cB_FGDVC.addItem(_fromUtf8(""))
        self.cB_FGDVC.addItem(_fromUtf8(""))
        self.cB_FGDVC.addItem(_fromUtf8(""))
        self.cB_FGDVC.addItem(_fromUtf8(""))
        self.cB_FGDVC.addItem(_fromUtf8(""))
        self.cB_FGDVC.addItem(_fromUtf8(""))
        self.cB_FGDVC.addItem(_fromUtf8(""))
        self.formLayout.setWidget(1, QFormLayout.FieldRole, self.cB_FGDVC)
        self.L_PCCL = QLabel(self.layoutWidget)
        self.L_PCCL.setObjectName(_fromUtf8("L_PCCL"))
        self.formLayout.setWidget(2, QFormLayout.LabelRole, self.L_PCCL)
        self.cB_PCCL = QComboBox(self.layoutWidget)
        self.cB_PCCL.setObjectName(_fromUtf8("cB_PCCL"))
        self.cB_PCCL.addItem(_fromUtf8(""))
        self.cB_PCCL.addItem(_fromUtf8(""))
        self.cB_PCCL.addItem(_fromUtf8(""))
        self.cB_PCCL.addItem(_fromUtf8(""))
        self.cB_PCCL.addItem(_fromUtf8(""))
        self.cB_PCCL.addItem(_fromUtf8(""))
        self.cB_PCCL.addItem(_fromUtf8(""))
        self.formLayout.setWidget(2, QFormLayout.FieldRole, self.cB_PCCL)
        self.layoutWidget1 = QWidget(self.centralwidget)
        self.layoutWidget1.setGeometry(QRect(540, 100, 201, 81))
        self.layoutWidget1.setObjectName(_fromUtf8("layoutWidget1"))
        self.formLayout_2 = QFormLayout(self.layoutWidget1)
        self.formLayout_2.setFieldGrowthPolicy(QFormLayout.ExpandingFieldsGrow)
        #self.formLayout_2.setMargin(0)
        self.formLayout_2.setObjectName(_fromUtf8("formLayout_2"))
        self.L_WeightParam = QLabel(self.layoutWidget1)
        self.L_WeightParam.setObjectName(_fromUtf8("L_WeightParam"))
        self.formLayout_2.setWidget(0, QFormLayout.SpanningRole, self.L_WeightParam)
        self.L_Yweight = QLabel(self.layoutWidget1)
        self.L_Yweight.setObjectName(_fromUtf8("L_Yweight"))
        self.formLayout_2.setWidget(1, QFormLayout.LabelRole, self.L_Yweight)
        self.lE_Yweight = QLineEdit(self.layoutWidget1)
        self.lE_Yweight.setObjectName(_fromUtf8("lE_Yweight"))
        self.formLayout_2.setWidget(1, QFormLayout.FieldRole, self.lE_Yweight)
        self.L_Rweight = QLabel(self.layoutWidget1)
        self.L_Rweight.setObjectName(_fromUtf8("L_Rweight"))
        self.formLayout_2.setWidget(2, QFormLayout.LabelRole, self.L_Rweight)
        self.lE_Rweight = QLineEdit(self.layoutWidget1)
        self.lE_Rweight.setObjectName(_fromUtf8("lE_Rweight"))
        self.formLayout_2.setWidget(2, QFormLayout.FieldRole, self.lE_Rweight)
        self.layoutWidget2 = QWidget(self.centralwidget)
        self.layoutWidget2.setGeometry(QRect(570, 280, 154, 51))
        self.layoutWidget2.setObjectName(_fromUtf8("layoutWidget2"))
        self.gridLayout_4 = QGridLayout(self.layoutWidget2)
        #self.gridLayout_4.setMargin(0)
        self.gridLayout_4.setObjectName(_fromUtf8("gridLayout_4"))
        self.L_MW = QLabel(self.layoutWidget2)
        self.L_MW.setObjectName(_fromUtf8("L_MW"))
        self.gridLayout_4.addWidget(self.L_MW, 0, 0, 1, 2)
        self.L_MWTar = QLabel(self.layoutWidget2)
        self.L_MWTar.setObjectName(_fromUtf8("L_MWTar"))
        self.gridLayout_4.addWidget(self.L_MWTar, 1, 0, 1, 1)
        self.lE_MWTar = QLineEdit(self.layoutWidget2)
        self.lE_MWTar.setObjectName(_fromUtf8("lE_MWTar"))
        self.gridLayout_4.addWidget(self.lE_MWTar, 1, 1, 1, 1)
        self.layoutWidget3 = QWidget(self.centralwidget)
        self.layoutWidget3.setGeometry(QRect(160, 450, 341, 53))
        self.layoutWidget3.setObjectName(_fromUtf8("layoutWidget3"))
        self.gridLayout_5 = QGridLayout(self.layoutWidget3)
        #self.gridLayout_5.setMargin(0)
        self.gridLayout_5.setObjectName(_fromUtf8("gridLayout_5"))
        self.L_FGDVCcoal = QLabel(self.layoutWidget3)
        self.L_FGDVCcoal.setObjectName(_fromUtf8("L_FGDVCcoal"))
        self.gridLayout_5.addWidget(self.L_FGDVCcoal, 0, 0, 1, 1)
        self.cB_FGDVCcoal = QComboBox(self.layoutWidget3)
        self.cB_FGDVCcoal.setObjectName(_fromUtf8("cB_FGDVCcoal"))
        self.cB_FGDVCcoal.addItem(_fromUtf8(""))
        self.cB_FGDVCcoal.addItem(_fromUtf8(""))
        self.cB_FGDVCcoal.addItem(_fromUtf8(""))
        self.cB_FGDVCcoal.addItem(_fromUtf8(""))
        self.cB_FGDVCcoal.addItem(_fromUtf8(""))
        self.cB_FGDVCcoal.addItem(_fromUtf8(""))
        self.cB_FGDVCcoal.addItem(_fromUtf8(""))
        self.cB_FGDVCcoal.addItem(_fromUtf8(""))
        self.cB_FGDVCcoal.addItem(_fromUtf8(""))
        self.gridLayout_5.addWidget(self.cB_FGDVCcoal, 0, 1, 1, 1)
        self.L_FGDVCtarCr = QLabel(self.layoutWidget3)
        self.L_FGDVCtarCr.setObjectName(_fromUtf8("L_FGDVCtarCr"))
        self.gridLayout_5.addWidget(self.L_FGDVCtarCr, 1, 0, 1, 1)
        self.lE_FGDVCtarCr = QLineEdit(self.layoutWidget3)
        self.lE_FGDVCtarCr.setObjectName(_fromUtf8("lE_FGDVCtarCr"))
        self.gridLayout_5.addWidget(self.lE_FGDVCtarCr, 1, 1, 1, 1)
        self.lE_numTimeStep = QLineEdit(self.centralwidget)
        self.lE_numTimeStep.setGeometry(QRect(580, 560, 113, 23))
        self.lE_numTimeStep.setObjectName(_fromUtf8("lE_numTimeStep"))
        self.L_numTimeStep = QLabel(self.centralwidget)
        self.L_numTimeStep.setGeometry(QRect(390, 560, 181, 24))
        self.L_numTimeStep.setObjectName(_fromUtf8("L_numTimeStep"))
        self.label = QLabel(self.centralwidget)
        self.label.setGeometry(QRect(780, 200, 53, 15))
        self.label.setText(_fromUtf8(""))
        self.label.setObjectName(_fromUtf8("label"))
        self.label_2 = QLabel(self.centralwidget)
        self.label_2.setGeometry(QRect(750, 380, 231, 171))
        self.label_2.setText(_fromUtf8(""))
        self.label_2.setPixmap(QPixmap(_fromUtf8("Logos/virtuhcon_logo.png")))
        self.label_2.setScaledContents(True)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.label_3 = QLabel(self.centralwidget)
        self.label_3.setGeometry(QRect(750, 80, 231, 171))
        self.label_3.setText(_fromUtf8(""))
        self.label_3.setPixmap(QPixmap(_fromUtf8("Logos/ntfd_rgb.png")))
        self.label_3.setScaledContents(True)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.tE_THist_2 = QTextEdit(self.centralwidget)
        self.tE_THist_2.setGeometry(QRect(340, 590, 111, 76))
        self.tE_THist_2.setObjectName(_fromUtf8("tE_THist_2"))
        self.tE_THist_3 = QTextEdit(self.centralwidget)
        self.tE_THist_3.setGeometry(QRect(460, 590, 111, 76))
        self.tE_THist_3.setObjectName(_fromUtf8("tE_THist_3"))
        self.tE_THist_4 = QTextEdit(self.centralwidget)
        self.tE_THist_4.setGeometry(QRect(580, 590, 111, 76))
        self.tE_THist_4.setObjectName(_fromUtf8("tE_THist_4"))
        self.tE_THist_5 = QTextEdit(self.centralwidget)
        self.tE_THist_5.setGeometry(QRect(700, 590, 111, 76))
        self.tE_THist_5.setObjectName(_fromUtf8("tE_THist_5"))
        self.sB_Nr_THist = QSpinBox(self.centralwidget)
        self.sB_Nr_THist.setGeometry(QRect(140, 630, 51, 23))
        self.sB_Nr_THist.setObjectName(_fromUtf8("sB_Nr_THist"))
        #manually added:
        self.sB_Nr_THist.setMinimum(1)
        self.sB_Nr_THist.setMaximum(5)
        #end manually added:
        self.B_Launch = QPushButton(self.centralwidget)
        self.B_Launch.setGeometry(QRect(880, 613, 101, 31))
        self.B_Launch.setObjectName(_fromUtf8("B_Launch"))
        self.B_Plot1 = QPushButton(self.centralwidget)
        self.B_Plot1.setGeometry(QRect(280, 670, 51, 24))
        self.B_Plot1.setObjectName(_fromUtf8("B_Plot1"))
        self.B_Plot2 = QPushButton(self.centralwidget)
        self.B_Plot2.setGeometry(QRect(400, 670, 51, 24))
        self.B_Plot2.setObjectName(_fromUtf8("B_Plot2"))
        self.B_Plot3 = QPushButton(self.centralwidget)
        self.B_Plot3.setGeometry(QRect(520, 670, 51, 24))
        self.B_Plot3.setObjectName(_fromUtf8("B_Plot3"))
        self.B_Plot4 = QPushButton(self.centralwidget)
        self.B_Plot4.setGeometry(QRect(640, 670, 51, 24))
        self.B_Plot4.setObjectName(_fromUtf8("B_Plot4"))
        self.B_Plot5 = QPushButton(self.centralwidget)
        self.B_Plot5.setGeometry(QRect(760, 670, 51, 24))
        self.B_Plot5.setObjectName(_fromUtf8("B_Plot5"))
        self.B_Open4 = QPushButton(self.centralwidget)
        self.B_Open4.setGeometry(QRect(580, 670, 51, 24))
        self.B_Open4.setObjectName(_fromUtf8("B_Open4"))
        self.B_Open2 = QPushButton(self.centralwidget)
        self.B_Open2.setGeometry(QRect(340, 670, 51, 24))
        self.B_Open2.setObjectName(_fromUtf8("B_Open2"))
        self.B_Open5 = QPushButton(self.centralwidget)
        self.B_Open5.setGeometry(QRect(700, 670, 51, 24))
        self.B_Open5.setObjectName(_fromUtf8("B_Open5"))
        self.B_Open3 = QPushButton(self.centralwidget)
        self.B_Open3.setGeometry(QRect(460, 670, 51, 24))
        self.B_Open3.setObjectName(_fromUtf8("B_Open3"))
        self.B_Open1 = QPushButton(self.centralwidget)
        self.B_Open1.setGeometry(QRect(220, 670, 51, 24))
        self.B_Open1.setObjectName(_fromUtf8("B_Open1"))
        self.Header1_2 = QLabel(self.centralwidget)
        self.Header1_2.setGeometry(QRect(170, 0, 721, 31))
        self.Header1_2.setObjectName(_fromUtf8("Header1_2"))
        self.layoutWidget4 = QWidget(self.centralwidget)
        self.layoutWidget4.setGeometry(QRect(31, 281, 202, 159))
        self.layoutWidget4.setObjectName(_fromUtf8("layoutWidget4"))
        self.formLayout_3 = QFormLayout(self.layoutWidget4)
        #self.formLayout_3.setMargin(0)
        self.formLayout_3.setObjectName(_fromUtf8("formLayout_3"))
        self.L_UA = QLabel(self.layoutWidget4)
        self.L_UA.setObjectName(_fromUtf8("L_UA"))
        self.formLayout_3.setWidget(0, QFormLayout.SpanningRole, self.L_UA)
        self.L_UAC = QLabel(self.layoutWidget4)
        self.L_UAC.setObjectName(_fromUtf8("L_UAC"))
        self.formLayout_3.setWidget(1, QFormLayout.LabelRole, self.L_UAC)
        self.lE_UAC = QLineEdit(self.layoutWidget4)
        self.lE_UAC.setObjectName(_fromUtf8("lE_UAC"))
        self.formLayout_3.setWidget(1, QFormLayout.FieldRole, self.lE_UAC)
        self.L_UAH = QLabel(self.layoutWidget4)
        self.L_UAH.setObjectName(_fromUtf8("L_UAH"))
        self.formLayout_3.setWidget(2, QFormLayout.LabelRole, self.L_UAH)
        self.lE_UAH = QLineEdit(self.layoutWidget4)
        self.lE_UAH.setObjectName(_fromUtf8("lE_UAH"))
        self.formLayout_3.setWidget(2, QFormLayout.FieldRole, self.lE_UAH)
        self.L_UAN = QLabel(self.layoutWidget4)
        self.L_UAN.setObjectName(_fromUtf8("L_UAN"))
        self.formLayout_3.setWidget(3, QFormLayout.LabelRole, self.L_UAN)
        self.lE_UAN = QLineEdit(self.layoutWidget4)
        self.lE_UAN.setObjectName(_fromUtf8("lE_UAN"))
        self.formLayout_3.setWidget(3, QFormLayout.FieldRole, self.lE_UAN)
        self.L_UAO = QLabel(self.layoutWidget4)
        self.L_UAO.setObjectName(_fromUtf8("L_UAO"))
        self.formLayout_3.setWidget(4, QFormLayout.LabelRole, self.L_UAO)
        self.lE_UAO = QLineEdit(self.layoutWidget4)
        self.lE_UAO.setObjectName(_fromUtf8("lE_UAO"))
        self.formLayout_3.setWidget(4, QFormLayout.FieldRole, self.lE_UAO)
        self.L_UAS = QLabel(self.layoutWidget4)
        self.L_UAS.setObjectName(_fromUtf8("L_UAS"))
        self.formLayout_3.setWidget(5, QFormLayout.LabelRole, self.L_UAS)
        self.lE_UAS = QLineEdit(self.layoutWidget4)
        self.lE_UAS.setObjectName(_fromUtf8("lE_UAS"))
        self.formLayout_3.setWidget(5, QFormLayout.FieldRole, self.lE_UAS)
        self.layoutWidget5 = QWidget(self.centralwidget)
        self.layoutWidget5.setGeometry(QRect(281, 281, 253, 134))
        self.layoutWidget5.setObjectName(_fromUtf8("layoutWidget5"))
        self.verticalLayout = QVBoxLayout(self.layoutWidget5)
        #self.verticalLayout.setMargin(0)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.L_PA = QLabel(self.layoutWidget5)
        self.L_PA.setObjectName(_fromUtf8("L_PA"))
        self.verticalLayout.addWidget(self.L_PA)
        self.gridLayout = QGridLayout()
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.L_PAFC = QLabel(self.layoutWidget5)
        self.L_PAFC.setObjectName(_fromUtf8("L_PAFC"))
        self.gridLayout.addWidget(self.L_PAFC, 0, 0, 1, 1)
        self.lE_PAFC = QLineEdit(self.layoutWidget5)
        self.lE_PAFC.setObjectName(_fromUtf8("lE_PAFC"))
        self.gridLayout.addWidget(self.lE_PAFC, 0, 1, 1, 2)
        self.L_PAVM = QLabel(self.layoutWidget5)
        self.L_PAVM.setObjectName(_fromUtf8("L_PAVM"))
        self.gridLayout.addWidget(self.L_PAVM, 1, 0, 1, 1)
        self.lE_PAVM = QLineEdit(self.layoutWidget5)
        self.lE_PAVM.setObjectName(_fromUtf8("lE_PAVM"))
        self.gridLayout.addWidget(self.lE_PAVM, 1, 1, 1, 2)
        self.L_PAMoi = QLabel(self.layoutWidget5)
        self.L_PAMoi.setObjectName(_fromUtf8("L_PAMoi"))
        self.gridLayout.addWidget(self.L_PAMoi, 2, 0, 1, 1)
        self.lE_PAMoi = QLineEdit(self.layoutWidget5)
        self.lE_PAMoi.setObjectName(_fromUtf8("lE_PAMoi"))
        self.gridLayout.addWidget(self.lE_PAMoi, 2, 2, 1, 1)
        self.L_PAAsh = QLabel(self.layoutWidget5)
        self.L_PAAsh.setObjectName(_fromUtf8("L_PAAsh"))
        self.gridLayout.addWidget(self.L_PAAsh, 3, 0, 1, 2)
        self.lE_PAAsh = QLineEdit(self.layoutWidget5)
        self.lE_PAAsh.setObjectName(_fromUtf8("lE_PAAsh"))
        self.gridLayout.addWidget(self.lE_PAAsh, 3, 2, 1, 1)
        self.verticalLayout.addLayout(self.gridLayout)
	self.cB_ArrhSpec = QComboBox(self.centralwidget)                         
        self.cB_ArrhSpec.setGeometry(QRect(330, 160, 118, 24))                  
        self.cB_ArrhSpec.setObjectName(_fromUtf8("cB_ArrhSpec"))                       
        self.cB_ArrhSpec.addItem(_fromUtf8(""))                                        
        self.cB_ArrhSpec.addItem(_fromUtf8(""))                                        
        self.cB_ArrhSpec.addItem(_fromUtf8(""))                                        
        self.L_ArrhSpec = QLabel(self.centralwidget)                             
        self.L_ArrhSpec.setGeometry(QRect(290, 110, 199, 41))                   
        self.L_ArrhSpec.setObjectName(_fromUtf8("L_ArrhSpec"))
        PKP.setCentralWidget(self.centralwidget)
        self.menubar = QMenuBar(PKP)
        self.menubar.setGeometry(QRect(0, 0, 1033, 21))
        self.menubar.setObjectName(_fromUtf8("menubar"))
        self.menuFile = QMenu(self.menubar)
        self.menuFile.setObjectName(_fromUtf8("menuFile"))
        self.menuHelp = QMenu(self.menubar)
        self.menuHelp.setObjectName(_fromUtf8("menuHelp"))
        PKP.setMenuBar(self.menubar)
        self.statusbar = QStatusBar(PKP)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        PKP.setStatusBar(self.statusbar)
        self.actionWrite_into_File = QAction(PKP)
        self.actionWrite_into_File.setObjectName(_fromUtf8("actionWrite_into_File"))
        self.actionWrite_and_Run = QAction(PKP)
        self.actionWrite_and_Run.setObjectName(_fromUtf8("actionWrite_and_Run"))
        self.actionExit = QAction(PKP)
        self.actionExit.setObjectName(_fromUtf8("actionExit"))
        self.actionShow_generated_Results = QAction(PKP)
        self.actionShow_generated_Results.setObjectName(_fromUtf8("actionShow_generated_Results"))
        self.actionShow_saved_state = QAction(PKP)
        self.actionShow_saved_state.setObjectName(_fromUtf8("actionShow_saved_state"))
        self.actionManual = QAction(PKP)
        self.actionManual.setObjectName(_fromUtf8("actionManual"))
        self.menuFile.addAction(self.actionWrite_into_File)
        self.menuFile.addAction(self.actionWrite_and_Run)
        self.menuFile.addAction(self.actionShow_saved_state)
        self.menuFile.addAction(self.actionShow_generated_Results)
        self.menuFile.addSeparator()
        self.menuFile.addAction(self.actionExit)
        self.menuHelp.addAction(self.actionManual)
        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuHelp.menuAction())

        self.retranslateUi(PKP)
        QObject.connect(self.actionExit, SIGNAL(_fromUtf8("triggered()")), PKP.close)
        QObject.connect(self.B_Launch, SIGNAL(_fromUtf8("clicked()")), self.WriteRun)
        QObject.connect(self.actionWrite_and_Run, SIGNAL(_fromUtf8("triggered()")), self.WriteRun)
        #manually added
        QObject.connect(self.actionWrite_into_File, SIGNAL(_fromUtf8("triggered()")), self.SaveInfos)
        QObject.connect(self.actionShow_saved_state, SIGNAL(_fromUtf8("triggered()")), self.LoadState)
        QObject.connect(self.actionShow_generated_Results, SIGNAL(_fromUtf8("triggered()")), self.ReOpenResultWindow)
        QObject.connect(self.actionManual, SIGNAL(_fromUtf8("triggered()")), self.OpenManual)
        #
        QObject.connect(self.B_Open1, SIGNAL(_fromUtf8("clicked()")), self.LoadTtFile1)
        QObject.connect(self.B_Open2, SIGNAL(_fromUtf8("clicked()")), self.LoadTtFile2)
        QObject.connect(self.B_Open3, SIGNAL(_fromUtf8("clicked()")), self.LoadTtFile3)
        QObject.connect(self.B_Open4, SIGNAL(_fromUtf8("clicked()")), self.LoadTtFile4)
        QObject.connect(self.B_Open5, SIGNAL(_fromUtf8("clicked()")), self.LoadTtFile5)
        #
        QObject.connect(self.B_Plot1, SIGNAL(_fromUtf8("clicked()")), self.Plot1)
        QObject.connect(self.B_Plot2, SIGNAL(_fromUtf8("clicked()")), self.Plot2)
        QObject.connect(self.B_Plot3, SIGNAL(_fromUtf8("clicked()")), self.Plot3)
        QObject.connect(self.B_Plot4, SIGNAL(_fromUtf8("clicked()")), self.Plot4)
        QObject.connect(self.B_Plot5, SIGNAL(_fromUtf8("clicked()")), self.Plot5)        
        #end manually added
        QMetaObject.connectSlotsByName(PKP)
        PKP.setTabOrder(self.cB_CPD, self.cB_FGDVC)
        PKP.setTabOrder(self.cB_FGDVC, self.cB_PCCL)
        PKP.setTabOrder(self.cB_PCCL, self.lE_Yweight)
        PKP.setTabOrder(self.lE_Yweight, self.lE_Rweight)
        PKP.setTabOrder(self.lE_Rweight, self.lE_UAC)
        PKP.setTabOrder(self.lE_UAC, self.lE_UAH)
        PKP.setTabOrder(self.lE_UAH, self.lE_UAN)
        PKP.setTabOrder(self.lE_UAN, self.lE_UAO)
        PKP.setTabOrder(self.lE_UAO, self.lE_UAS)
        PKP.setTabOrder(self.lE_UAS, self.lE_PAFC)
        PKP.setTabOrder(self.lE_PAFC, self.lE_PAVM)
        PKP.setTabOrder(self.lE_PAVM, self.lE_PAMoi)
        PKP.setTabOrder(self.lE_PAMoi, self.lE_PAAsh)
        PKP.setTabOrder(self.lE_PAAsh, self.lE_MWTar)
        PKP.setTabOrder(self.lE_MWTar, self.lE_HHV)
        PKP.setTabOrder(self.lE_HHV, self.cB_FGDVCcoal)
        PKP.setTabOrder(self.cB_FGDVCcoal, self.lE_FGDVCtarCr)
        PKP.setTabOrder(self.lE_FGDVCtarCr, self.lE_pressure)
        PKP.setTabOrder(self.lE_pressure, self.lE_numTimeStep)
        PKP.setTabOrder(self.lE_numTimeStep, self.tE_THist_1)
        PKP.setTabOrder(self.tE_THist_1, self.B_Open1)
        PKP.setTabOrder(self.B_Open1, self.B_Plot1)
        PKP.setTabOrder(self.B_Plot1, self.tE_THist_2)
        PKP.setTabOrder(self.tE_THist_2, self.B_Open2)
        PKP.setTabOrder(self.B_Open2, self.B_Plot2)
        PKP.setTabOrder(self.B_Plot2, self.tE_THist_3)
        PKP.setTabOrder(self.tE_THist_3, self.B_Open3)
        PKP.setTabOrder(self.B_Open3, self.B_Plot3)
        PKP.setTabOrder(self.B_Plot3, self.tE_THist_4)
        PKP.setTabOrder(self.tE_THist_4, self.B_Open4)
        PKP.setTabOrder(self.B_Open4, self.B_Plot4)
        PKP.setTabOrder(self.B_Plot4, self.tE_THist_5)
        PKP.setTabOrder(self.tE_THist_5, self.B_Open5)
        PKP.setTabOrder(self.B_Open5, self.B_Plot5)
        PKP.setTabOrder(self.B_Plot5, self.sB_Nr_THist)
        PKP.setTabOrder(self.sB_Nr_THist, self.B_Launch)

    def retranslateUi(self, PKP):
        PKP.setWindowTitle(QApplication.translate("PKP", "MainWindow", None, QApplication.UnicodeUTF8))
        self.Header1.setText(QApplication.translate("PKP", "<html><head/><body><p align=\"center\"><span style=\" font-size:18pt; font-weight:600;\">Pyrolysis Programs and Fitting</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.Header2.setText(QApplication.translate("PKP", "<html><head/><body><p align=\"center\"><span style=\" font-size:18pt; font-weight:600;\">Coal Properties (as received)</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.L_HHV.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:14pt;\">Higher Heating Value<br/>in MJ/kg</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.Header3.setText(QApplication.translate("PKP", "<html><head/><body><p align=\"center\"><span style=\" font-size:18pt; font-weight:600;\">Operating Conditions</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.L_pressure.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:14pt;\">Pressure in atm</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.lE_pressure.setText(QApplication.translate("PKP", "1", None, QApplication.UnicodeUTF8))
        self.L_THist.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:14pt;\">Temperature History</span></p><p><span style=\" font-size:10pt;\">t in s    T in K</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.L_CPD.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:14pt;\">CPD</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.cB_CPD.setItemText(0, QApplication.translate("PKP", "None", None, QApplication.UnicodeUTF8))
        self.cB_CPD.setItemText(1, QApplication.translate("PKP", "Run", None, QApplication.UnicodeUTF8))
        self.cB_CPD.setItemText(2, QApplication.translate("PKP", "constant Rate", None, QApplication.UnicodeUTF8))
        self.cB_CPD.setItemText(3, QApplication.translate("PKP", "Arrhenius", None, QApplication.UnicodeUTF8))
        self.cB_CPD.setItemText(4, QApplication.translate("PKP", "Arrhenius no B", None, QApplication.UnicodeUTF8))
        self.cB_CPD.setItemText(5, QApplication.translate("PKP", "Kobayashi", None, QApplication.UnicodeUTF8))
        self.cB_CPD.setItemText(6, QApplication.translate("PKP", "DAEM", None, QApplication.UnicodeUTF8))
        self.L_FGDVC.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:14pt;\">FG-DVC</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.cB_FGDVC.setItemText(0, QApplication.translate("PKP", "None", None, QApplication.UnicodeUTF8))
        self.cB_FGDVC.setItemText(1, QApplication.translate("PKP", "Run", None, QApplication.UnicodeUTF8))
        self.cB_FGDVC.setItemText(2, QApplication.translate("PKP", "constant Rate", None, QApplication.UnicodeUTF8))
        self.cB_FGDVC.setItemText(3, QApplication.translate("PKP", "Arrhenius", None, QApplication.UnicodeUTF8))
        self.cB_FGDVC.setItemText(4, QApplication.translate("PKP", "Arrhenius No B", None, QApplication.UnicodeUTF8))
        self.cB_FGDVC.setItemText(5, QApplication.translate("PKP", "Kobayashi", None, QApplication.UnicodeUTF8))
        self.cB_FGDVC.setItemText(6, QApplication.translate("PKP", "DAEM", None, QApplication.UnicodeUTF8))
        self.L_PCCL.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:14pt;\">PC Coal Lab</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.cB_PCCL.setItemText(0, QApplication.translate("PKP", "None", None, QApplication.UnicodeUTF8))
        self.cB_PCCL.setItemText(1, QApplication.translate("PKP", "Run", None, QApplication.UnicodeUTF8))
        self.cB_PCCL.setItemText(2, QApplication.translate("PKP", "constant Rate", None, QApplication.UnicodeUTF8))
        self.cB_PCCL.setItemText(3, QApplication.translate("PKP", "Arrhenius", None, QApplication.UnicodeUTF8))
        self.cB_PCCL.setItemText(4, QApplication.translate("PKP", "Arrhenius No B", None, QApplication.UnicodeUTF8))
        self.cB_PCCL.setItemText(5, QApplication.translate("PKP", "Kobayashi", None, QApplication.UnicodeUTF8))
        self.cB_PCCL.setItemText(6, QApplication.translate("PKP", "DAEM", None, QApplication.UnicodeUTF8))
        self.L_WeightParam.setText(QApplication.translate("PKP", "<html><head/><body><p align=\"center\"><span style=\" font-size:14pt;\">Weight Parameter</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.L_Yweight.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:14pt;\">Yields</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.lE_Yweight.setText(QApplication.translate("PKP", "1", None, QApplication.UnicodeUTF8))
        self.L_Rweight.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:14pt;\">Rates</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.lE_Rweight.setText(QApplication.translate("PKP", "1", None, QApplication.UnicodeUTF8))
        self.L_MW.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:14pt;\">Molecule Weight</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.L_MWTar.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:14pt;\">Tar</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.L_FGDVCcoal.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:14pt;\">FG-DVC Coal #</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.cB_FGDVCcoal.setItemText(0, QApplication.translate("PKP", "0 - Interpolate", None, QApplication.UnicodeUTF8))
        self.cB_FGDVCcoal.setItemText(1, QApplication.translate("PKP", "1 - Beulah-Zap", None, QApplication.UnicodeUTF8))
        self.cB_FGDVCcoal.setItemText(2, QApplication.translate("PKP", "2 - Woydak-Anderson", None, QApplication.UnicodeUTF8))
        self.cB_FGDVCcoal.setItemText(3, QApplication.translate("PKP", "3 - Illinois # 6", None, QApplication.UnicodeUTF8))
        self.cB_FGDVCcoal.setItemText(4, QApplication.translate("PKP", "4 - Bind Canyon, UT", None, QApplication.UnicodeUTF8))
        self.cB_FGDVCcoal.setItemText(5, QApplication.translate("PKP", "5 - Lewis-Stockton, WV", None, QApplication.UnicodeUTF8))
        self.cB_FGDVCcoal.setItemText(6, QApplication.translate("PKP", "6 - Pittsburgh # 8", None, QApplication.UnicodeUTF8))
        self.cB_FGDVCcoal.setItemText(7, QApplication.translate("PKP", "7 - Upper Freeport, PA", None, QApplication.UnicodeUTF8))
        self.cB_FGDVCcoal.setItemText(8, QApplication.translate("PKP", "8 - Pocahontas, VA", None, QApplication.UnicodeUTF8))
        self.L_FGDVCtarCr.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:14pt;\">FG-DVC Tar Cracking</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.lE_FGDVCtarCr.setText(QApplication.translate("PKP", "0", None, QApplication.UnicodeUTF8))
        self.lE_numTimeStep.setText(QApplication.translate("PKP", "1.e-4", None, QApplication.UnicodeUTF8))
        self.L_numTimeStep.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:14pt;\">numerical time step</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.B_Launch.setText(QApplication.translate("PKP", "Launch", None, QApplication.UnicodeUTF8))
        self.B_Plot1.setText(QApplication.translate("PKP", "Plot", None, QApplication.UnicodeUTF8))
        self.B_Plot2.setText(QApplication.translate("PKP", "Plot", None, QApplication.UnicodeUTF8))
        self.B_Plot3.setText(QApplication.translate("PKP", "Plot", None, QApplication.UnicodeUTF8))
        self.B_Plot4.setText(QApplication.translate("PKP", "Plot", None, QApplication.UnicodeUTF8))
        self.B_Plot5.setText(QApplication.translate("PKP", "Plot", None, QApplication.UnicodeUTF8))
        self.B_Open4.setText(QApplication.translate("PKP", "Open", None, QApplication.UnicodeUTF8))
        self.B_Open2.setText(QApplication.translate("PKP", "Open", None, QApplication.UnicodeUTF8))
        self.B_Open5.setText(QApplication.translate("PKP", "Open", None, QApplication.UnicodeUTF8))
        self.B_Open3.setText(QApplication.translate("PKP", "Open", None, QApplication.UnicodeUTF8))
        self.B_Open1.setText(QApplication.translate("PKP", "Open", None, QApplication.UnicodeUTF8))
        self.Header1_2.setText(QApplication.translate("PKP", "<html><head/><body><p align=\"center\"><span style=\" font-size:22pt; font-weight:600;\">PKP</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.L_UA.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:14pt;\">Ultimate Analysis in %</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.L_UAC.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:14pt;\">Carbon</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.L_UAH.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:14pt;\">Hydrogen</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.L_UAN.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:14pt;\">Nitrogen</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.L_UAO.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:14pt;\">Oxygen</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.L_UAS.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:14pt;\">Sulphur</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.L_PA.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:14pt;\">Proximate Analysis in %</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.L_PAFC.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:14pt;\">Fixed Carbon</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.L_PAVM.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:14pt;\">Volatile Matter</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.L_PAMoi.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:14pt;\">Moisture</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.L_PAAsh.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:14pt;\">Ash</span></p></body></html>", None, QApplication.UnicodeUTF8))
 	self.cB_ArrhSpec.setItemText(0, QApplication.translate("PKP", "Total", None, QApplication.UnicodeUTF8))
        self.cB_ArrhSpec.setItemText(1, QApplication.translate("PKP", "Main Species", None, QApplication.UnicodeUTF8))
        self.cB_ArrhSpec.setItemText(2, QApplication.translate("PKP", "all Species", None, QApplication.UnicodeUTF8))
        self.L_ArrhSpec.setText(QApplication.translate("PKP", "<html><head/><body><p align=\"center\"><span style=\" font-size:14pt;\">selected Fit Species<br/>(Arrhenius)</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.menuFile.setTitle(QApplication.translate("PKP", "File", None, QApplication.UnicodeUTF8))
        self.menuHelp.setTitle(QApplication.translate("PKP", "Help", None, QApplication.UnicodeUTF8))
        self.actionWrite_into_File.setText(QApplication.translate("PKP", "Write into File", None, QApplication.UnicodeUTF8))
        self.actionWrite_into_File.setShortcut(QApplication.translate("PKP", "Ctrl+S", None, QApplication.UnicodeUTF8))
        self.actionWrite_and_Run.setText(QApplication.translate("PKP", "Write and Run", None, QApplication.UnicodeUTF8))
        self.actionWrite_and_Run.setShortcut(QApplication.translate("PKP", "Ctrl+R", None, QApplication.UnicodeUTF8))
        self.actionExit.setText(QApplication.translate("PKP", "Exit", None, QApplication.UnicodeUTF8))
        self.actionExit.setShortcut(QApplication.translate("PKP", "Ctrl+Q", None, QApplication.UnicodeUTF8))
        self.actionShow_generated_Results.setText(QApplication.translate("PKP", "Show generated Results", None, QApplication.UnicodeUTF8))
        self.actionShow_generated_Results.setShortcut(QApplication.translate("PKP", "Ctrl+T", None, QApplication.UnicodeUTF8))
        self.actionShow_saved_state.setText(QApplication.translate("PKP", "Load saved state", None, QApplication.UnicodeUTF8))
        self.actionShow_saved_state.setShortcut(QApplication.translate("PKP", "Ctrl+O", None, QApplication.UnicodeUTF8))
	self.actionManual.setText(QApplication.translate("PKP", "Manual", None, QApplication.UnicodeUTF8))

    #manually added
    def OpenManual(self):
        """Opens the manual."""
        if OSys=='Linux':
            os.system('okular Documentation/Manual/PKPDocumentation.pdf')
        elif OSys=='Windows':
            os.system('Documentation\\Manual\\PKPDocumentation.pdf')
    
    def SaveInfos(self):
        """Saves the Information when the save or the run -option is used."""
        #Pyrolysis Program selection
        CPDsel = self.cB_CPD.currentIndex()
        FGsel  = self.cB_FGDVC.currentIndex()
        PCCLsel= self.cB_PCCL.currentIndex()
        self.svInfo.SetRunPyrolProg(CPDsel,FGsel,PCCLsel)
        #saves the Species for Arrhenius to fit
        ArrSpec=self.cB_ArrhSpec.currentIndex()
        self.svInfo.SetArrhSpec(ArrSpec)
        #saves the WeightParameter
        weightY = str(self.lE_Yweight.text())
        weightR = str(self.lE_Rweight.text())
        self.svInfo.setWeightYR(weightY,weightR)
        #Coal Properties
        UAC = str(self.lE_UAC.text()) #information UA
        UAH = str(self.lE_UAH.text())
        UAN = str(self.lE_UAN.text())
        UAO = str(self.lE_UAO.text())
        UAS = str(self.lE_UAS.text())
        self.svInfo.setUA(UAC,UAH,UAN,UAO,UAS)
        PAFC = str(self.lE_PAFC.text())
        PAVM = str(self.lE_PAVM.text())
        PAMoi = str(self.lE_PAMoi.text())
        PAAsh = str(self.lE_PAAsh.text())
        self.svInfo.setPA(PAFC,PAVM,PAMoi,PAAsh)
        HHV = str(self.lE_HHV.text()) #	information HHV
        MWTar = str(self.lE_MWTar.text())
        self.svInfo.setMwsHHV(MWTar,HHV)
        FGCoal = str(self.cB_FGDVCcoal.currentIndex())  # information FG-DVC coal interpolation
        FGTar  = self.lE_FGDVCtarCr.text()         # information FG-DVC tar cracking
        self.svInfo.setFGCoalProp(FGCoal,FGTar)
        #operating condition
        pressure = str(self.lE_pressure.text()) # operating pressure
        dt = str(self.lE_numTimeStep.text()) # numerical time step
        self.svInfo.setOperCond(pressure,dt)
        # setTime Histories
        NrOfT = str(self.sB_Nr_THist.value())
#        T1 = str(self.tE_THist_1.toPlainText())
#        T2 = str(self.tE_THist_2.toPlainText())
#        T3 = str(self.tE_THist_3.toPlainText())
#        T4 = str(self.tE_THist_4.toPlainText())
#        T5 = str(self.tE_THist_5.toPlainText())
        file1=open('TempHist1.dat','w')
        file1.write(self.tE_THist_1.toPlainText())
        file1.close()
        file2=open('TempHist2.dat','w')
        file2.write(self.tE_THist_2.toPlainText())
        file2.close()
        file3=open('TempHist3.dat','w')
        file3.write(self.tE_THist_3.toPlainText())
        file3.close()
        file4=open('TempHist4.dat','w')
        file4.write(self.tE_THist_4.toPlainText())
        file4.close()
        file5=open('TempHist5.dat','w')
        file5.write(self.tE_THist_5.toPlainText())
        file5.close()
        self.svInfo.setTimeHistories(NrOfT)
        #
        #write Files
        writeInfoFiles.WriteCoalFile(self.svInfo)
        writeInfoFiles.WriteCPDFile(self.svInfo)
        writeInfoFiles.WriteFGFile(self.svInfo)
        writeInfoFiles.WriteOCFile(self.svInfo)
        
    def LoadTtFile1(self):
        """Loads the temperature history nr 1 file via file browser"""
        filename = QFileDialog.getOpenFileName()
        file1=open(filename[0],'r')
        data = file1.read()
        self.tE_THist_1.setText(data)
        file1.close()
    def LoadTtFile2(self):
        """Loads the temperature history nr 2 file via file browser"""
        filename = QFileDialog.getOpenFileName()
        file2=open(filename[0],'r')
        data = file2.read()
        self.tE_THist_2.setText(data)
        file2.close()
    def LoadTtFile3(self):
        """Loads the temperature history nr 3 file via file browser"""
        filename = QFileDialog.getOpenFileName()
        file3=open(filename[0],'r')
        data = file3.read()
        self.tE_THist_3.setText(data)
        file3.close()
    def LoadTtFile4(self):
        """Loads the temperature history nr 4 file via file browser"""
        filename = QFileDialog.getOpenFileName()
        file4=open(filename[0],'r')
        data = file4.read()
        self.tE_THist_4.setText(data)
        file4.close()
    def LoadTtFile5(self):
        """Loads the temperature history nr 5 file via file browser"""
        filename = QFileDialog.getOpenFileName()
        file5=open(filename[0],'r')
        data = file5.read()
        self.tE_THist_5.setText(data)
        file5.close()
        
    def Plot1(self):
        """Plots the temperature over time history (temperature history nr 1) and saves temperatuer history in "TempHist1.dat"."""
        file1=open('TempHist1.dat','w')
        file1.write(self.tE_THist_1.toPlainText())
        file1.close()
#        plt.clf(),plt.cla()
        Tt=np.genfromtxt('TempHist1.dat')
#        figT1=plt.figure()
#        figTt1=figT1.add_subplot(111)
        plt.ylabel('Temperature in K')
        plt.xlabel('Time in s')
        plt.legend()
        plt.grid()
        if len(Tt[:,0])<50:
            plt.plot(Tt[:,0],Tt[:,1],'-x')
        else:
            plt.plot(Tt[:,0],Tt[:,1],'-')
        plt.show()
    def Plot2(self):
        """Plots the temperature over time history (temperature history nr 2) and saves temperatuer history in "TempHist2.dat"."""
        file2=open('TempHist2.dat','w')
        file2.write(self.tE_THist_2.toPlainText())
        file2.close()
        Tt=np.genfromtxt('TempHist2.dat')
#        plt.clf(),plt.cla()
#        figT2=plt.figure()
#        figTt2=figT2.add_subplot(111)
        plt.ylabel('Temperature in K')
        plt.xlabel('Time in s')
        plt.legend()
        plt.grid()
        if len(Tt[:,0])<50:
            plt.plot(Tt[:,0],Tt[:,1],'-x')
        else:
            plt.plot(Tt[:,0],Tt[:,1],'-')
        plt.show()
    def Plot3(self):
        """Plots the temperature over time history (temperature history nr 3) and saves temperatuer history in "TempHist3.dat"."""
        file3=open('TempHist3.dat','w')
        file3.write(self.tE_THist_3.toPlainText())
        file3.close()
        Tt=np.genfromtxt('TempHist3.dat')
#        plt.clf(),plt.cla()
#        figT3=plt.figure()
#        figTt3=figT3.add_subplot(111)
        plt.ylabel('Temperature in K')
        plt.xlabel('Time in s')
        plt.legend()
        plt.grid()
        if len(Tt[:,0])<50:
            plt.plot(Tt[:,0],Tt[:,1],'-x')
        else:
            plt.plot(Tt[:,0],Tt[:,1],'-')
        plt.show()
    def Plot4(self):
        """Plots the temperature over time history (temperature history nr 4) and saves temperatuer history in "TempHist4.dat"."""
        file4=open('TempHist4.dat','w')
        file4.write(self.tE_THist_4.toPlainText())
        file4.close()
        Tt=np.genfromtxt('TempHist4.dat')
#        plt.clf(),plt.cla()
#        figT4=plt.figure()
#        figTt4=figT4.add_subplot(111)
        plt.ylabel('Temperature in K')
        plt.xlabel('Time in s')
        plt.legend()
        plt.grid()
        if len(Tt[:,0])<50:
            plt.plot(Tt[:,0],Tt[:,1],'-x')
        else:
            plt.plot(Tt[:,0],Tt[:,1],'-')
        plt.show()
    def Plot5(self):
        """Plots the temperature over time history (temperature history nr 5) and saves temperatuer history in "TempHist5.dat"."""
        file5=open('TempHist5.dat','w')
        file5.write(self.tE_THist_5.toPlainText())
        file5.close()
        Tt=np.genfromtxt('TempHist5.dat')
#        plt.clf(),plt.cla()
#        figT5=plt.figure()
#        figTt5=figT5.add_subplot(111)
        plt.ylabel('Temperature in K')
        plt.xlabel('Time in s')
        plt.legend()
        plt.grid()
        if len(Tt[:,0])<50:
            plt.plot(Tt[:,0],Tt[:,1],'-x')
        else:
            plt.plot(Tt[:,0],Tt[:,1],'-')
        plt.show()

    def WriteRun(self):
        """Writes the *.inp files and launch the process."""
        self.SaveInfos()
        #os.system('python Pyrolysis.py')
        #print 'START executable'
        #
        #removes date from Result directory
        for filename in os.listdir('Result'):
            filepath = os.path.join('Result', filename)
            try:
                shutil.rmtree(filepath)
            except OSError:
                    os.remove(filepath)
        #
        Case=PKP.MainProcess()
        Case.ReadInputFiles()
        ProgramL=[]
        if Case.CPDselect==True:
            Case.MakeResults_CPD()
            ProgramL.append('CPD')
        if Case.FG_select==True:
            Case.CheckFGdt()
            Case.MakeResults_FG()
            ProgramL.append('FGDVC')
        SpeciesL=Case.SpeciesToConsider #e.g. ["Total","Tar","Gas"]
        ProgrFitD=Case.ProgramModelDict #e.g. {'CPD':'ArrheniusRate'}
        #
        self.Done=Done.Ui_Dialog()
        self.Done.setupUi(SpeciesL,ProgramL,ProgrFitD,(self.sB_Nr_THist.value()))
        self.Done.show()
        self.TheCalculationsAreDone=True
        #
#        os.system('python Done.py')
#        self.actionExit
        
    def ReOpenResultWindow(self):
        """If the calculation were done once, the window showing the results can be opened again."""
        if self.TheCalculationsAreDone==True:
            self.Done.show()
        else:
            print 'Please launch the process first to show the results.'
        
        
    def LoadState(self):
        CoalInput=InformationFiles.ReadFile('Coal.inp')
        PAFC_asrec=CoalInput.getText(InformationFiles.M_PA[0])
        PAVM_asrec=CoalInput.getText(InformationFiles.M_PA[1])
        PAmoist = CoalInput.getText(InformationFiles.M_PA[2])
        PAash = CoalInput.getText(InformationFiles.M_PA[3])
        #
        #gets daf values, as CPD needs daf as input:
        UAC=CoalInput.getText(InformationFiles.M_UA[0])
        UAH=CoalInput.getText(InformationFiles.M_UA[1])
        UAN=CoalInput.getText(InformationFiles.M_UA[2])
        UAO=CoalInput.getText(InformationFiles.M_UA[3])
        UAS=CoalInput.getText(InformationFiles.M_UA[4])
        # scale ultimate analysis
        HHV=str(CoalInput.getValue(InformationFiles.M_HHV)*1e-6)
        MTar=CoalInput.getText(InformationFiles.M_MTar)
        DryDens=CoalInput.getText(InformationFiles.M_density)
        WeightY=CoalInput.getText(InformationFiles.M_Weight[0])
        WeightR=CoalInput.getText(InformationFiles.M_Weight[1])
        #
        #CPD Properties:
        CPDInput=InformationFiles.ReadFile('CPD.inp')
#        CPDselect=CPDInput.UsePyrolProgr(InformationFiles.MC_sel)
        CPD_FittingKineticParameter_Select=CPDInput.Fitting(InformationFiles.M_selFit)
        CPD_ArrhSpec=CPDInput.getText(InformationFiles.M_selArrhSpec)
        CPDdt=[0,1,2] #0:initila dt, 1: print increment, 2: dt max
        CPDdt[0]=(CPDInput.getText(InformationFiles.MC_dt[0]))
        CPDdt[1]=(CPDInput.getText(InformationFiles.MC_dt[1]))
        #
        #
        #FG-DVC Properties:
        FGDVCInput=InformationFiles.ReadFile('FGDVC.inp')
#        FG_select=FGDVCInput.UsePyrolProgr(InformationFiles.MF_sel)
        FG_FittingKineticParameter_Select=FGDVCInput.Fitting(InformationFiles.M_selFit)
        FG_CoalSelection=int(FGDVCInput.getValue(InformationFiles.MF_CoalSel))
        FG_TarCacking=FGDVCInput.getText(InformationFiles.MF_TarCr)
        #
        #
        #Operating Condition File:
        OpCondInp=InformationFiles.OperCondInput('OperCond.inp')
        CPD_pressure=OpCondInp.getText(InformationFiles.M_Pressure)
        #Number of FG-DVC/CPD/PCCL runs:
        NrOfRuns=int(OpCondInp.getValue(InformationFiles.M_NrRuns))
        CPDdt[2]=OpCondInp.getText(InformationFiles.M_dt)
        FG_dt=OpCondInp.getText(InformationFiles.M_dt)
        #
        #
        #
        #set it into the GUI:
        if CPD_FittingKineticParameter_Select == None:
            CPD_FittingKineticParameter_Select='None'
        if FG_FittingKineticParameter_Select == None:
            FG_FittingKineticParameter_Select='None'
        self.cB_CPD.setCurrentIndex(self.svInfo.RunPyrolProgReverse(CPD_FittingKineticParameter_Select))
        self.cB_FGDVC.setCurrentIndex(self.svInfo.RunPyrolProgReverse(FG_FittingKineticParameter_Select))
#        self.cB_PCCL.setCurrentIndex(self.svInfo.RunPyrolProgReverse(PCCL_FittingKineticParameter_Select))
        #
        self.cB_ArrhSpec.setCurrentIndex(self.svInfo.ArrhSpecReverse(CPD_ArrhSpec))
        #saves the WeightParameter
        self.lE_Yweight.setText(WeightY)
        self.lE_Rweight.setText(WeightR)
        #Coal Properties
        self.lE_UAC.setText(UAC) #information UA
        self.lE_UAH.setText(UAH)
        self.lE_UAN.setText(UAN)
        self.lE_UAO.setText(UAO)
        self.lE_UAS.setText(UAS)
        self.lE_PAFC.setText(PAFC_asrec)
        self.lE_PAVM.setText(PAVM_asrec)
        self.lE_PAMoi.setText(PAmoist)
        self.lE_PAAsh.setText(PAash)
        #s
        self.lE_HHV.setText(HHV) #	information HHV
        self.lE_MWTar.setText(MTar)
        self.cB_FGDVCcoal.setCurrentIndex(FG_CoalSelection)  # information FG-DVC coal interpolation
        self.lE_FGDVCtarCr.setText(FG_TarCacking)         # information FG-DVC tar cracking
        #operating condition
        self.lE_pressure.setText(CPD_pressure) # operating pressure
        self.lE_numTimeStep.setText(FG_dt) # numerical time step
        # setTime Histories
        self.sB_Nr_THist.setValue(NrOfRuns)
        file1=open('TempHist1.dat','r')
        self.tE_THist_1.setText(file1.read())
        file1.close()
        file2=open('TempHist2.dat','r')
        self.tE_THist_2.setText(file2.read())
        file2.close()
        file3=open('TempHist3.dat','r')
        self.tE_THist_3.setText(file3.read())
        file3.close()
        file4=open('TempHist4.dat','r')
        self.tE_THist_4.setText(file4.read())
        file4.close()
        file5=open('TempHist5.dat','r')
        self.tE_THist_5.setText(file5.read())
        file5.close()
        #



############################################################################


class InfosFromGUI(object):
    """Saves the information from the GUI."""
    
    def SetRunPyrolProg(self,CPDIndex,FGDVCIndex,PCCLIndex):
        """Saves which options of the three Pyrolysis programs are used."""
        FitDict={0:'None',1:'Run',2:'constantRate',3:'Arrhenius',4:'ArrheniusNoB',5:'Kobayashi',6:'DAEM'}
        self.__CPDsel = FitDict[CPDIndex]
        self.__FGsel  = FitDict[FGDVCIndex]
        self.__PCCLsel= FitDict[PCCLIndex]
        
    def RunPyrolProgReverse(self,ModelName):
        """Returns the GUI column bar index of the models selected."""
        FitDict={'None':0,'Run':1,'constantRate':2,'Arrhenius':3,'ArrheniusNoB':4,'Kobayashi':5,'DAEM':6}
        return FitDict[ModelName]

    def RunPyrolProg(self):
        """Returns which options of the three Pyrolysis programs are used."""
        return self.__CPDsel, self.__FGsel, self.__PCCLsel

    def SetArrhSpec(self,SpeciesIndex):
        """Sets which species shall be fitted for Arrhenius."""
        SpecIndexDict={0:'Total', 1:'MainSpecies', 2:'allSpecies'}
        Species=SpecIndexDict[SpeciesIndex]
        self.__ArrSpec=Species
        
    def ArrhSpecReverse(self,SpeciesName):
        """returns the UI columns bar index of species fitted for Arrhenius."""
        SpecIndexDict={'Total':0, 'MainSpecies':1, 'allSpecies':2}
        return SpecIndexDict[SpeciesName]
        
    def ArrhSpec(self):
        """Returns which species shall be fitted for Arrhenius."""
        return self.__ArrSpec

    def setWeightYR(self,Y,R):
        """Sets the weghts for Yields and Rates"""
        self.__WeightY=Y
        self.__WeightR=R
        
    def WeightYR(self):
        """Returns the weghts for Yields and Rates"""
        return self.__WeightY, self.__WeightR

    def setUA(self,UAC,UAH,UAN,UAO,UAS):
        """Saves the coal UA properties."""
        self.__UAC=UAC
        self.__UAH=UAH
        self.__UAN=UAN
        self.__UAO=UAO
        self.__UAS=UAS
    
    def UA(self):
        """Returns the coal UA properties."""
        return self.__UAC,self.__UAH,self.__UAN,self.__UAO,self.__UAS

    def setPA(self,FixedCarbon,VolatileMatter,Moisture,Ash):
        """Saves the coal PA properties."""
        self.__PAFC = FixedCarbon
        self.__PAVM = VolatileMatter
        self.__PAMoi = Moisture
        self.__PAASh = Ash
    
    def PA(self):
        """Returns the coal PA properties."""
        return self.__PAFC,self.__PAVM,self.__PAMoi,self.__PAASh
    
    def setMwsHHV(self,MoleWweightTar,HHV):
        """Sets the Molar Weight of Tar and sets the Higher Heating Value."""
        self.__MwTar = MoleWweightTar
        self.__HHV   = str(float(HHV)*1e6)

    def MwsHHV(self):
        """Retruns the Molar Weight of Tar and sets the Higher Heating Value."""
        return self.__MwTar,self.__HHV

    def setFGCoalProp(self,FGCoalFit,FGTarModeling):
        """Defines the way of the FG-DVC Coal Fitting and the Tar Modeling."""
        self.__FGCoalFit = FGCoalFit
        self.__FGTarModeling = FGTarModeling

    def FGCoalProp(self):
        """Defines the way of the FG-DVC Coal Fitting and the Tar Modeling."""
        return self.__FGCoalFit, self.__FGTarModeling
        
    def setOperCond(self,pressure,timestep):
        """Sets the pressure and the time step."""
        self.__p = pressure
        self.__dt = timestep

    def OperCond(self):
        """Sets the pressure and the time step."""
        return self.__p, self.__dt
        
    def setTimeHistories(self,NrOfTs):
        """Saves the Number time history"""
        self.__nrT=NrOfTs

    def TimeHistories(self):
        """Returns the number of time history."""
        return self.__nrT
    
    def setCoalDens(self,CoalDensity):
        """Saves the Coal Density"""
        self.__CoalDens = CoalDensity
        
    def CoalDens(self):
        """Returns the Coal Density"""
        return self.__CoalDens

if __name__ == "__main__":
    import sys
    Infosaver=InfosFromGUI()
    app = QApplication(sys.argv)
    PKPWindow = QMainWindow()
    ui = Ui_PKP()
    ui.setupUi(PKPWindow,Infosaver)
    PKPWindow.show()
    sys.exit(app.exec_())

