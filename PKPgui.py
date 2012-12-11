# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'MainWindow.ui'

from PyQt4 import QtCore, QtGui
import pylab as plt
import numpy as np
import os
import writeInfoFiles
import sys
import platform

from PyQt4.Qt import QWidget, QMainWindow

OSys=platform.system()
sys.path.append('src')
import PKP
import InformationFiles

import Done   #second GUI

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class Ui_PKP(QMainWindow):
    def setupUi(self, PKP, SaveInfoObj):
        QMainWindow.__init__(self) #manually added
        self.QWidgMain= QWidget(self) #manually added
        self.setCentralWidget(self.QWidgMain) #manually added
        self.TheCalculationsAreDone=False #manually added
        self.svInfo=SaveInfoObj
        PKP.setObjectName(_fromUtf8("PKP"))
        PKP.resize(1033, 753)
        self.centralwidget = QtGui.QWidget(PKP)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.Header1 = QtGui.QLabel(self.centralwidget)
        self.Header1.setGeometry(QtCore.QRect(-10, 50, 721, 31))
        self.Header1.setObjectName(_fromUtf8("Header1"))
        self.Header2 = QtGui.QLabel(self.centralwidget)
        self.Header2.setGeometry(QtCore.QRect(30, 230, 721, 31))
        self.Header2.setObjectName(_fromUtf8("Header2"))
        self.L_HHV = QtGui.QLabel(self.centralwidget)
        self.L_HHV.setGeometry(QtCore.QRect(570, 350, 191, 51))
        self.L_HHV.setObjectName(_fromUtf8("L_HHV"))
        self.lE_HHV = QtGui.QLineEdit(self.centralwidget)
        self.lE_HHV.setGeometry(QtCore.QRect(570, 400, 113, 23))
        self.lE_HHV.setObjectName(_fromUtf8("lE_HHV"))
        self.Header3 = QtGui.QLabel(self.centralwidget)
        self.Header3.setGeometry(QtCore.QRect(20, 520, 721, 31))
        self.Header3.setObjectName(_fromUtf8("Header3"))
        self.L_pressure = QtGui.QLabel(self.centralwidget)
        self.L_pressure.setGeometry(QtCore.QRect(30, 560, 141, 24))
        self.L_pressure.setObjectName(_fromUtf8("L_pressure"))
        self.lE_pressure = QtGui.QLineEdit(self.centralwidget)
        self.lE_pressure.setGeometry(QtCore.QRect(220, 560, 113, 23))
        self.lE_pressure.setObjectName(_fromUtf8("lE_pressure"))
        self.L_THist = QtGui.QLabel(self.centralwidget)
        self.L_THist.setGeometry(QtCore.QRect(30, 590, 191, 71))
        self.L_THist.setObjectName(_fromUtf8("L_THist"))
        self.tE_THist_1 = QtGui.QTextEdit(self.centralwidget)
        self.tE_THist_1.setGeometry(QtCore.QRect(220, 590, 111, 76))
        self.tE_THist_1.setObjectName(_fromUtf8("tE_THist_1"))
        self.layoutWidget = QtGui.QWidget(self.centralwidget)
        self.layoutWidget.setGeometry(QtCore.QRect(40, 100, 241, 91))
        self.layoutWidget.setObjectName(_fromUtf8("layoutWidget"))
        self.formLayout = QtGui.QFormLayout(self.layoutWidget)
        self.formLayout.setFieldGrowthPolicy(QtGui.QFormLayout.ExpandingFieldsGrow)
        self.formLayout.setMargin(0)
        self.formLayout.setObjectName(_fromUtf8("formLayout"))
        self.L_CPD = QtGui.QLabel(self.layoutWidget)
        self.L_CPD.setObjectName(_fromUtf8("L_CPD"))
        self.formLayout.setWidget(0, QtGui.QFormLayout.LabelRole, self.L_CPD)
        self.cB_CPD = QtGui.QComboBox(self.layoutWidget)
        self.cB_CPD.setObjectName(_fromUtf8("cB_CPD"))
        self.cB_CPD.addItem(_fromUtf8(""))
        self.cB_CPD.addItem(_fromUtf8(""))
        self.cB_CPD.addItem(_fromUtf8(""))
        self.cB_CPD.addItem(_fromUtf8(""))
        self.cB_CPD.addItem(_fromUtf8(""))
        self.cB_CPD.addItem(_fromUtf8(""))
        self.cB_CPD.addItem(_fromUtf8(""))
        self.formLayout.setWidget(0, QtGui.QFormLayout.FieldRole, self.cB_CPD)
        self.L_FGDVC = QtGui.QLabel(self.layoutWidget)
        self.L_FGDVC.setObjectName(_fromUtf8("L_FGDVC"))
        self.formLayout.setWidget(1, QtGui.QFormLayout.LabelRole, self.L_FGDVC)
        self.cB_FGDVC = QtGui.QComboBox(self.layoutWidget)
        self.cB_FGDVC.setObjectName(_fromUtf8("cB_FGDVC"))
        self.cB_FGDVC.addItem(_fromUtf8(""))
        self.cB_FGDVC.addItem(_fromUtf8(""))
        self.cB_FGDVC.addItem(_fromUtf8(""))
        self.cB_FGDVC.addItem(_fromUtf8(""))
        self.cB_FGDVC.addItem(_fromUtf8(""))
        self.cB_FGDVC.addItem(_fromUtf8(""))
        self.cB_FGDVC.addItem(_fromUtf8(""))
        self.formLayout.setWidget(1, QtGui.QFormLayout.FieldRole, self.cB_FGDVC)
        self.L_PCCL = QtGui.QLabel(self.layoutWidget)
        self.L_PCCL.setObjectName(_fromUtf8("L_PCCL"))
        self.formLayout.setWidget(2, QtGui.QFormLayout.LabelRole, self.L_PCCL)
        self.cB_PCCL = QtGui.QComboBox(self.layoutWidget)
        self.cB_PCCL.setObjectName(_fromUtf8("cB_PCCL"))
        self.cB_PCCL.addItem(_fromUtf8(""))
        self.cB_PCCL.addItem(_fromUtf8(""))
        self.cB_PCCL.addItem(_fromUtf8(""))
        self.cB_PCCL.addItem(_fromUtf8(""))
        self.cB_PCCL.addItem(_fromUtf8(""))
        self.cB_PCCL.addItem(_fromUtf8(""))
        self.cB_PCCL.addItem(_fromUtf8(""))
        self.formLayout.setWidget(2, QtGui.QFormLayout.FieldRole, self.cB_PCCL)
        self.layoutWidget1 = QtGui.QWidget(self.centralwidget)
        self.layoutWidget1.setGeometry(QtCore.QRect(540, 100, 201, 81))
        self.layoutWidget1.setObjectName(_fromUtf8("layoutWidget1"))
        self.formLayout_2 = QtGui.QFormLayout(self.layoutWidget1)
        self.formLayout_2.setFieldGrowthPolicy(QtGui.QFormLayout.ExpandingFieldsGrow)
        self.formLayout_2.setMargin(0)
        self.formLayout_2.setObjectName(_fromUtf8("formLayout_2"))
        self.L_WeightParam = QtGui.QLabel(self.layoutWidget1)
        self.L_WeightParam.setObjectName(_fromUtf8("L_WeightParam"))
        self.formLayout_2.setWidget(0, QtGui.QFormLayout.SpanningRole, self.L_WeightParam)
        self.L_Yweight = QtGui.QLabel(self.layoutWidget1)
        self.L_Yweight.setObjectName(_fromUtf8("L_Yweight"))
        self.formLayout_2.setWidget(1, QtGui.QFormLayout.LabelRole, self.L_Yweight)
        self.lE_Yweight = QtGui.QLineEdit(self.layoutWidget1)
        self.lE_Yweight.setObjectName(_fromUtf8("lE_Yweight"))
        self.formLayout_2.setWidget(1, QtGui.QFormLayout.FieldRole, self.lE_Yweight)
        self.L_Rweight = QtGui.QLabel(self.layoutWidget1)
        self.L_Rweight.setObjectName(_fromUtf8("L_Rweight"))
        self.formLayout_2.setWidget(2, QtGui.QFormLayout.LabelRole, self.L_Rweight)
        self.lE_Rweight = QtGui.QLineEdit(self.layoutWidget1)
        self.lE_Rweight.setObjectName(_fromUtf8("lE_Rweight"))
        self.formLayout_2.setWidget(2, QtGui.QFormLayout.FieldRole, self.lE_Rweight)
        self.layoutWidget2 = QtGui.QWidget(self.centralwidget)
        self.layoutWidget2.setGeometry(QtCore.QRect(570, 280, 154, 51))
        self.layoutWidget2.setObjectName(_fromUtf8("layoutWidget2"))
        self.gridLayout_4 = QtGui.QGridLayout(self.layoutWidget2)
        self.gridLayout_4.setMargin(0)
        self.gridLayout_4.setObjectName(_fromUtf8("gridLayout_4"))
        self.L_MW = QtGui.QLabel(self.layoutWidget2)
        self.L_MW.setObjectName(_fromUtf8("L_MW"))
        self.gridLayout_4.addWidget(self.L_MW, 0, 0, 1, 2)
        self.L_MWTar = QtGui.QLabel(self.layoutWidget2)
        self.L_MWTar.setObjectName(_fromUtf8("L_MWTar"))
        self.gridLayout_4.addWidget(self.L_MWTar, 1, 0, 1, 1)
        self.lE_MWTar = QtGui.QLineEdit(self.layoutWidget2)
        self.lE_MWTar.setObjectName(_fromUtf8("lE_MWTar"))
        self.gridLayout_4.addWidget(self.lE_MWTar, 1, 1, 1, 1)
        self.layoutWidget3 = QtGui.QWidget(self.centralwidget)
        self.layoutWidget3.setGeometry(QtCore.QRect(160, 450, 341, 53))
        self.layoutWidget3.setObjectName(_fromUtf8("layoutWidget3"))
        self.gridLayout_5 = QtGui.QGridLayout(self.layoutWidget3)
        self.gridLayout_5.setMargin(0)
        self.gridLayout_5.setObjectName(_fromUtf8("gridLayout_5"))
        self.L_FGDVCcoal = QtGui.QLabel(self.layoutWidget3)
        self.L_FGDVCcoal.setObjectName(_fromUtf8("L_FGDVCcoal"))
        self.gridLayout_5.addWidget(self.L_FGDVCcoal, 0, 0, 1, 1)
        self.cB_FGDVCcoal = QtGui.QComboBox(self.layoutWidget3)
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
        self.L_FGDVCtarCr = QtGui.QLabel(self.layoutWidget3)
        self.L_FGDVCtarCr.setObjectName(_fromUtf8("L_FGDVCtarCr"))
        self.gridLayout_5.addWidget(self.L_FGDVCtarCr, 1, 0, 1, 1)
        self.lE_FGDVCtarCr = QtGui.QLineEdit(self.layoutWidget3)
        self.lE_FGDVCtarCr.setObjectName(_fromUtf8("lE_FGDVCtarCr"))
        self.gridLayout_5.addWidget(self.lE_FGDVCtarCr, 1, 1, 1, 1)
        self.lE_numTimeStep = QtGui.QLineEdit(self.centralwidget)
        self.lE_numTimeStep.setGeometry(QtCore.QRect(580, 560, 113, 23))
        self.lE_numTimeStep.setObjectName(_fromUtf8("lE_numTimeStep"))
        self.L_numTimeStep = QtGui.QLabel(self.centralwidget)
        self.L_numTimeStep.setGeometry(QtCore.QRect(390, 560, 181, 24))
        self.L_numTimeStep.setObjectName(_fromUtf8("L_numTimeStep"))
        self.label = QtGui.QLabel(self.centralwidget)
        self.label.setGeometry(QtCore.QRect(780, 200, 53, 15))
        self.label.setText(_fromUtf8(""))
        self.label.setObjectName(_fromUtf8("label"))
        self.label_2 = QtGui.QLabel(self.centralwidget)
        self.label_2.setGeometry(QtCore.QRect(750, 380, 231, 171))
        self.label_2.setText(_fromUtf8(""))
        self.label_2.setPixmap(QtGui.QPixmap(_fromUtf8("Logos/virtuhcon_logo.png")))
        self.label_2.setScaledContents(True)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.label_3 = QtGui.QLabel(self.centralwidget)
        self.label_3.setGeometry(QtCore.QRect(750, 80, 231, 171))
        self.label_3.setText(_fromUtf8(""))
        self.label_3.setPixmap(QtGui.QPixmap(_fromUtf8("Logos/ntfd_rgb.png")))
        self.label_3.setScaledContents(True)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.tE_THist_2 = QtGui.QTextEdit(self.centralwidget)
        self.tE_THist_2.setGeometry(QtCore.QRect(340, 590, 111, 76))
        self.tE_THist_2.setObjectName(_fromUtf8("tE_THist_2"))
        self.tE_THist_3 = QtGui.QTextEdit(self.centralwidget)
        self.tE_THist_3.setGeometry(QtCore.QRect(460, 590, 111, 76))
        self.tE_THist_3.setObjectName(_fromUtf8("tE_THist_3"))
        self.tE_THist_4 = QtGui.QTextEdit(self.centralwidget)
        self.tE_THist_4.setGeometry(QtCore.QRect(580, 590, 111, 76))
        self.tE_THist_4.setObjectName(_fromUtf8("tE_THist_4"))
        self.tE_THist_5 = QtGui.QTextEdit(self.centralwidget)
        self.tE_THist_5.setGeometry(QtCore.QRect(700, 590, 111, 76))
        self.tE_THist_5.setObjectName(_fromUtf8("tE_THist_5"))
        self.sB_Nr_THist = QtGui.QSpinBox(self.centralwidget)
        self.sB_Nr_THist.setGeometry(QtCore.QRect(140, 630, 51, 23))
        self.sB_Nr_THist.setObjectName(_fromUtf8("sB_Nr_THist"))
        #manually added:
        self.sB_Nr_THist.setMinimum(1)
        self.sB_Nr_THist.setMaximum(5)
        #end manually added:
        self.B_Launch = QtGui.QPushButton(self.centralwidget)
        self.B_Launch.setGeometry(QtCore.QRect(880, 613, 101, 31))
        self.B_Launch.setObjectName(_fromUtf8("B_Launch"))
        self.B_Plot1 = QtGui.QPushButton(self.centralwidget)
        self.B_Plot1.setGeometry(QtCore.QRect(280, 670, 51, 24))
        self.B_Plot1.setObjectName(_fromUtf8("B_Plot1"))
        self.B_Plot2 = QtGui.QPushButton(self.centralwidget)
        self.B_Plot2.setGeometry(QtCore.QRect(400, 670, 51, 24))
        self.B_Plot2.setObjectName(_fromUtf8("B_Plot2"))
        self.B_Plot3 = QtGui.QPushButton(self.centralwidget)
        self.B_Plot3.setGeometry(QtCore.QRect(520, 670, 51, 24))
        self.B_Plot3.setObjectName(_fromUtf8("B_Plot3"))
        self.B_Plot4 = QtGui.QPushButton(self.centralwidget)
        self.B_Plot4.setGeometry(QtCore.QRect(640, 670, 51, 24))
        self.B_Plot4.setObjectName(_fromUtf8("B_Plot4"))
        self.B_Plot5 = QtGui.QPushButton(self.centralwidget)
        self.B_Plot5.setGeometry(QtCore.QRect(760, 670, 51, 24))
        self.B_Plot5.setObjectName(_fromUtf8("B_Plot5"))
        self.B_Open4 = QtGui.QPushButton(self.centralwidget)
        self.B_Open4.setGeometry(QtCore.QRect(580, 670, 51, 24))
        self.B_Open4.setObjectName(_fromUtf8("B_Open4"))
        self.B_Open2 = QtGui.QPushButton(self.centralwidget)
        self.B_Open2.setGeometry(QtCore.QRect(340, 670, 51, 24))
        self.B_Open2.setObjectName(_fromUtf8("B_Open2"))
        self.B_Open5 = QtGui.QPushButton(self.centralwidget)
        self.B_Open5.setGeometry(QtCore.QRect(700, 670, 51, 24))
        self.B_Open5.setObjectName(_fromUtf8("B_Open5"))
        self.B_Open3 = QtGui.QPushButton(self.centralwidget)
        self.B_Open3.setGeometry(QtCore.QRect(460, 670, 51, 24))
        self.B_Open3.setObjectName(_fromUtf8("B_Open3"))
        self.B_Open1 = QtGui.QPushButton(self.centralwidget)
        self.B_Open1.setGeometry(QtCore.QRect(220, 670, 51, 24))
        self.B_Open1.setObjectName(_fromUtf8("B_Open1"))
        self.Header1_2 = QtGui.QLabel(self.centralwidget)
        self.Header1_2.setGeometry(QtCore.QRect(170, 0, 721, 31))
        self.Header1_2.setObjectName(_fromUtf8("Header1_2"))
        self.layoutWidget4 = QtGui.QWidget(self.centralwidget)
        self.layoutWidget4.setGeometry(QtCore.QRect(31, 281, 202, 159))
        self.layoutWidget4.setObjectName(_fromUtf8("layoutWidget4"))
        self.formLayout_3 = QtGui.QFormLayout(self.layoutWidget4)
        self.formLayout_3.setMargin(0)
        self.formLayout_3.setObjectName(_fromUtf8("formLayout_3"))
        self.L_UA = QtGui.QLabel(self.layoutWidget4)
        self.L_UA.setObjectName(_fromUtf8("L_UA"))
        self.formLayout_3.setWidget(0, QtGui.QFormLayout.SpanningRole, self.L_UA)
        self.L_UAC = QtGui.QLabel(self.layoutWidget4)
        self.L_UAC.setObjectName(_fromUtf8("L_UAC"))
        self.formLayout_3.setWidget(1, QtGui.QFormLayout.LabelRole, self.L_UAC)
        self.lE_UAC = QtGui.QLineEdit(self.layoutWidget4)
        self.lE_UAC.setObjectName(_fromUtf8("lE_UAC"))
        self.formLayout_3.setWidget(1, QtGui.QFormLayout.FieldRole, self.lE_UAC)
        self.L_UAH = QtGui.QLabel(self.layoutWidget4)
        self.L_UAH.setObjectName(_fromUtf8("L_UAH"))
        self.formLayout_3.setWidget(2, QtGui.QFormLayout.LabelRole, self.L_UAH)
        self.lE_UAH = QtGui.QLineEdit(self.layoutWidget4)
        self.lE_UAH.setObjectName(_fromUtf8("lE_UAH"))
        self.formLayout_3.setWidget(2, QtGui.QFormLayout.FieldRole, self.lE_UAH)
        self.L_UAN = QtGui.QLabel(self.layoutWidget4)
        self.L_UAN.setObjectName(_fromUtf8("L_UAN"))
        self.formLayout_3.setWidget(3, QtGui.QFormLayout.LabelRole, self.L_UAN)
        self.lE_UAN = QtGui.QLineEdit(self.layoutWidget4)
        self.lE_UAN.setObjectName(_fromUtf8("lE_UAN"))
        self.formLayout_3.setWidget(3, QtGui.QFormLayout.FieldRole, self.lE_UAN)
        self.L_UAO = QtGui.QLabel(self.layoutWidget4)
        self.L_UAO.setObjectName(_fromUtf8("L_UAO"))
        self.formLayout_3.setWidget(4, QtGui.QFormLayout.LabelRole, self.L_UAO)
        self.lE_UAO = QtGui.QLineEdit(self.layoutWidget4)
        self.lE_UAO.setObjectName(_fromUtf8("lE_UAO"))
        self.formLayout_3.setWidget(4, QtGui.QFormLayout.FieldRole, self.lE_UAO)
        self.L_UAS = QtGui.QLabel(self.layoutWidget4)
        self.L_UAS.setObjectName(_fromUtf8("L_UAS"))
        self.formLayout_3.setWidget(5, QtGui.QFormLayout.LabelRole, self.L_UAS)
        self.lE_UAS = QtGui.QLineEdit(self.layoutWidget4)
        self.lE_UAS.setObjectName(_fromUtf8("lE_UAS"))
        self.formLayout_3.setWidget(5, QtGui.QFormLayout.FieldRole, self.lE_UAS)
        self.layoutWidget5 = QtGui.QWidget(self.centralwidget)
        self.layoutWidget5.setGeometry(QtCore.QRect(281, 281, 253, 134))
        self.layoutWidget5.setObjectName(_fromUtf8("layoutWidget5"))
        self.verticalLayout = QtGui.QVBoxLayout(self.layoutWidget5)
        self.verticalLayout.setMargin(0)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.L_PA = QtGui.QLabel(self.layoutWidget5)
        self.L_PA.setObjectName(_fromUtf8("L_PA"))
        self.verticalLayout.addWidget(self.L_PA)
        self.gridLayout = QtGui.QGridLayout()
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.L_PAFC = QtGui.QLabel(self.layoutWidget5)
        self.L_PAFC.setObjectName(_fromUtf8("L_PAFC"))
        self.gridLayout.addWidget(self.L_PAFC, 0, 0, 1, 1)
        self.lE_PAFC = QtGui.QLineEdit(self.layoutWidget5)
        self.lE_PAFC.setObjectName(_fromUtf8("lE_PAFC"))
        self.gridLayout.addWidget(self.lE_PAFC, 0, 1, 1, 2)
        self.L_PAVM = QtGui.QLabel(self.layoutWidget5)
        self.L_PAVM.setObjectName(_fromUtf8("L_PAVM"))
        self.gridLayout.addWidget(self.L_PAVM, 1, 0, 1, 1)
        self.lE_PAVM = QtGui.QLineEdit(self.layoutWidget5)
        self.lE_PAVM.setObjectName(_fromUtf8("lE_PAVM"))
        self.gridLayout.addWidget(self.lE_PAVM, 1, 1, 1, 2)
        self.L_PAMoi = QtGui.QLabel(self.layoutWidget5)
        self.L_PAMoi.setObjectName(_fromUtf8("L_PAMoi"))
        self.gridLayout.addWidget(self.L_PAMoi, 2, 0, 1, 1)
        self.lE_PAMoi = QtGui.QLineEdit(self.layoutWidget5)
        self.lE_PAMoi.setObjectName(_fromUtf8("lE_PAMoi"))
        self.gridLayout.addWidget(self.lE_PAMoi, 2, 2, 1, 1)
        self.L_PAAsh = QtGui.QLabel(self.layoutWidget5)
        self.L_PAAsh.setObjectName(_fromUtf8("L_PAAsh"))
        self.gridLayout.addWidget(self.L_PAAsh, 3, 0, 1, 2)
        self.lE_PAAsh = QtGui.QLineEdit(self.layoutWidget5)
        self.lE_PAAsh.setObjectName(_fromUtf8("lE_PAAsh"))
        self.gridLayout.addWidget(self.lE_PAAsh, 3, 2, 1, 1)
        self.verticalLayout.addLayout(self.gridLayout)
	self.cB_ArrhSpec = QtGui.QComboBox(self.centralwidget)                         
        self.cB_ArrhSpec.setGeometry(QtCore.QRect(330, 160, 118, 24))                  
        self.cB_ArrhSpec.setObjectName(_fromUtf8("cB_ArrhSpec"))                       
        self.cB_ArrhSpec.addItem(_fromUtf8(""))                                        
        self.cB_ArrhSpec.addItem(_fromUtf8(""))                                        
        self.cB_ArrhSpec.addItem(_fromUtf8(""))                                        
        self.L_ArrhSpec = QtGui.QLabel(self.centralwidget)                             
        self.L_ArrhSpec.setGeometry(QtCore.QRect(290, 110, 199, 41))                   
        self.L_ArrhSpec.setObjectName(_fromUtf8("L_ArrhSpec"))
        PKP.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(PKP)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1033, 21))
        self.menubar.setObjectName(_fromUtf8("menubar"))
        self.menuFile = QtGui.QMenu(self.menubar)
        self.menuFile.setObjectName(_fromUtf8("menuFile"))
        self.menuHelp = QtGui.QMenu(self.menubar)
        self.menuHelp.setObjectName(_fromUtf8("menuHelp"))
        PKP.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(PKP)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        PKP.setStatusBar(self.statusbar)
        self.actionWrite_into_File = QtGui.QAction(PKP)
        self.actionWrite_into_File.setObjectName(_fromUtf8("actionWrite_into_File"))
        self.actionWrite_and_Run = QtGui.QAction(PKP)
        self.actionWrite_and_Run.setObjectName(_fromUtf8("actionWrite_and_Run"))
        self.actionExit = QtGui.QAction(PKP)
        self.actionExit.setObjectName(_fromUtf8("actionExit"))
        self.actionShow_generated_Results = QtGui.QAction(PKP)
        self.actionShow_generated_Results.setObjectName(_fromUtf8("actionShow_generated_Results"))
        self.actionShow_saved_state = QtGui.QAction(PKP)
        self.actionShow_saved_state.setObjectName(_fromUtf8("actionShow_saved_state"))
        self.menuFile.addAction(self.actionWrite_into_File)
        self.menuFile.addAction(self.actionWrite_and_Run)
        self.menuFile.addAction(self.actionShow_saved_state)
        self.menuFile.addAction(self.actionShow_generated_Results)
        self.menuFile.addSeparator()
        self.menuFile.addAction(self.actionExit)
        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuHelp.menuAction())

        self.retranslateUi(PKP)
        QtCore.QObject.connect(self.actionExit, QtCore.SIGNAL(_fromUtf8("activated()")), PKP.close)
        QtCore.QObject.connect(self.B_Launch, QtCore.SIGNAL(_fromUtf8("clicked()")), self.WriteRun)
        QtCore.QObject.connect(self.actionWrite_and_Run, QtCore.SIGNAL(_fromUtf8("activated()")), self.WriteRun)
        #manually added
        QtCore.QObject.connect(self.actionWrite_into_File, QtCore.SIGNAL(_fromUtf8("activated()")), self.SaveInfos)
        QtCore.QObject.connect(self.actionShow_saved_state, QtCore.SIGNAL(_fromUtf8("activated()")), self.LoadState)
        QtCore.QObject.connect(self.actionShow_generated_Results, QtCore.SIGNAL(_fromUtf8("activated()")), self.ReOpenResultWindow)
        #
        QtCore.QObject.connect(self.B_Open1, QtCore.SIGNAL(_fromUtf8("clicked()")), self.LoadTtFile1)
        QtCore.QObject.connect(self.B_Open2, QtCore.SIGNAL(_fromUtf8("clicked()")), self.LoadTtFile2)
        QtCore.QObject.connect(self.B_Open3, QtCore.SIGNAL(_fromUtf8("clicked()")), self.LoadTtFile3)
        QtCore.QObject.connect(self.B_Open4, QtCore.SIGNAL(_fromUtf8("clicked()")), self.LoadTtFile4)
        QtCore.QObject.connect(self.B_Open5, QtCore.SIGNAL(_fromUtf8("clicked()")), self.LoadTtFile5)
        #
        QtCore.QObject.connect(self.B_Plot1, QtCore.SIGNAL(_fromUtf8("clicked()")), self.Plot1)
        QtCore.QObject.connect(self.B_Plot2, QtCore.SIGNAL(_fromUtf8("clicked()")), self.Plot2)
        QtCore.QObject.connect(self.B_Plot3, QtCore.SIGNAL(_fromUtf8("clicked()")), self.Plot3)
        QtCore.QObject.connect(self.B_Plot4, QtCore.SIGNAL(_fromUtf8("clicked()")), self.Plot4)
        QtCore.QObject.connect(self.B_Plot5, QtCore.SIGNAL(_fromUtf8("clicked()")), self.Plot5)        
        #end manually added
        QtCore.QMetaObject.connectSlotsByName(PKP)
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
        PKP.setWindowTitle(QtGui.QApplication.translate("PKP", "MainWindow", None, QtGui.QApplication.UnicodeUTF8))
        self.Header1.setText(QtGui.QApplication.translate("PKP", "<html><head/><body><p align=\"center\"><span style=\" font-size:18pt; font-weight:600;\">Pyrolysis Programs and Fitting</span></p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.Header2.setText(QtGui.QApplication.translate("PKP", "<html><head/><body><p align=\"center\"><span style=\" font-size:18pt; font-weight:600;\">Coal Properties (as received)</span></p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.L_HHV.setText(QtGui.QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:14pt;\">Higher Heating Value<br/>in MJ/kg</span></p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.Header3.setText(QtGui.QApplication.translate("PKP", "<html><head/><body><p align=\"center\"><span style=\" font-size:18pt; font-weight:600;\">Operating Conditions</span></p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.L_pressure.setText(QtGui.QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:14pt;\">Pressure in atm</span></p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.lE_pressure.setText(QtGui.QApplication.translate("PKP", "1", None, QtGui.QApplication.UnicodeUTF8))
        self.L_THist.setText(QtGui.QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:14pt;\">Temperature History</span></p><p><span style=\" font-size:10pt;\">t in s    T in K</span></p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.L_CPD.setText(QtGui.QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:14pt;\">CPD</span></p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.cB_CPD.setItemText(0, QtGui.QApplication.translate("PKP", "None", None, QtGui.QApplication.UnicodeUTF8))
        self.cB_CPD.setItemText(1, QtGui.QApplication.translate("PKP", "Run", None, QtGui.QApplication.UnicodeUTF8))
        self.cB_CPD.setItemText(2, QtGui.QApplication.translate("PKP", "constant Rate", None, QtGui.QApplication.UnicodeUTF8))
        self.cB_CPD.setItemText(3, QtGui.QApplication.translate("PKP", "Arrhenius", None, QtGui.QApplication.UnicodeUTF8))
        self.cB_CPD.setItemText(4, QtGui.QApplication.translate("PKP", "Arrhenius no B", None, QtGui.QApplication.UnicodeUTF8))
        self.cB_CPD.setItemText(5, QtGui.QApplication.translate("PKP", "Kobayashi", None, QtGui.QApplication.UnicodeUTF8))
        self.cB_CPD.setItemText(6, QtGui.QApplication.translate("PKP", "DAEM", None, QtGui.QApplication.UnicodeUTF8))
        self.L_FGDVC.setText(QtGui.QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:14pt;\">FG-DVC</span></p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.cB_FGDVC.setItemText(0, QtGui.QApplication.translate("PKP", "None", None, QtGui.QApplication.UnicodeUTF8))
        self.cB_FGDVC.setItemText(1, QtGui.QApplication.translate("PKP", "Run", None, QtGui.QApplication.UnicodeUTF8))
        self.cB_FGDVC.setItemText(2, QtGui.QApplication.translate("PKP", "constant Rate", None, QtGui.QApplication.UnicodeUTF8))
        self.cB_FGDVC.setItemText(3, QtGui.QApplication.translate("PKP", "Arrhenius", None, QtGui.QApplication.UnicodeUTF8))
        self.cB_FGDVC.setItemText(4, QtGui.QApplication.translate("PKP", "Arrhenius No B", None, QtGui.QApplication.UnicodeUTF8))
        self.cB_FGDVC.setItemText(5, QtGui.QApplication.translate("PKP", "Kobayashi", None, QtGui.QApplication.UnicodeUTF8))
        self.cB_FGDVC.setItemText(6, QtGui.QApplication.translate("PKP", "DAEM", None, QtGui.QApplication.UnicodeUTF8))
        self.L_PCCL.setText(QtGui.QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:14pt;\">PC Coal Lab</span></p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.cB_PCCL.setItemText(0, QtGui.QApplication.translate("PKP", "None", None, QtGui.QApplication.UnicodeUTF8))
        self.cB_PCCL.setItemText(1, QtGui.QApplication.translate("PKP", "Run", None, QtGui.QApplication.UnicodeUTF8))
        self.cB_PCCL.setItemText(2, QtGui.QApplication.translate("PKP", "constant Rate", None, QtGui.QApplication.UnicodeUTF8))
        self.cB_PCCL.setItemText(3, QtGui.QApplication.translate("PKP", "Arrhenius", None, QtGui.QApplication.UnicodeUTF8))
        self.cB_PCCL.setItemText(4, QtGui.QApplication.translate("PKP", "Arrhenius No B", None, QtGui.QApplication.UnicodeUTF8))
        self.cB_PCCL.setItemText(5, QtGui.QApplication.translate("PKP", "Kobayashi", None, QtGui.QApplication.UnicodeUTF8))
        self.cB_PCCL.setItemText(6, QtGui.QApplication.translate("PKP", "DAEM", None, QtGui.QApplication.UnicodeUTF8))
        self.L_WeightParam.setText(QtGui.QApplication.translate("PKP", "<html><head/><body><p align=\"center\"><span style=\" font-size:14pt;\">Weight Parameter</span></p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.L_Yweight.setText(QtGui.QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:14pt;\">Yields</span></p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.lE_Yweight.setText(QtGui.QApplication.translate("PKP", "1", None, QtGui.QApplication.UnicodeUTF8))
        self.L_Rweight.setText(QtGui.QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:14pt;\">Rates</span></p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.lE_Rweight.setText(QtGui.QApplication.translate("PKP", "1", None, QtGui.QApplication.UnicodeUTF8))
        self.L_MW.setText(QtGui.QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:14pt;\">Molecule Weight</span></p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.L_MWTar.setText(QtGui.QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:14pt;\">Tar</span></p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.L_FGDVCcoal.setText(QtGui.QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:14pt;\">FG-DVC Coal #</span></p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.cB_FGDVCcoal.setItemText(0, QtGui.QApplication.translate("PKP", "0 - Interpolate", None, QtGui.QApplication.UnicodeUTF8))
        self.cB_FGDVCcoal.setItemText(1, QtGui.QApplication.translate("PKP", "1 - Beulah-Zap", None, QtGui.QApplication.UnicodeUTF8))
        self.cB_FGDVCcoal.setItemText(2, QtGui.QApplication.translate("PKP", "2 - Woydak-Anderson", None, QtGui.QApplication.UnicodeUTF8))
        self.cB_FGDVCcoal.setItemText(3, QtGui.QApplication.translate("PKP", "3 - Illinois # 6", None, QtGui.QApplication.UnicodeUTF8))
        self.cB_FGDVCcoal.setItemText(4, QtGui.QApplication.translate("PKP", "4 - Bind Canyon, UT", None, QtGui.QApplication.UnicodeUTF8))
        self.cB_FGDVCcoal.setItemText(5, QtGui.QApplication.translate("PKP", "5 - Lewis-Stockton, WV", None, QtGui.QApplication.UnicodeUTF8))
        self.cB_FGDVCcoal.setItemText(6, QtGui.QApplication.translate("PKP", "6 - Pittsburgh # 8", None, QtGui.QApplication.UnicodeUTF8))
        self.cB_FGDVCcoal.setItemText(7, QtGui.QApplication.translate("PKP", "7 - Upper Freeport, PA", None, QtGui.QApplication.UnicodeUTF8))
        self.cB_FGDVCcoal.setItemText(8, QtGui.QApplication.translate("PKP", "8 - Pocahontas, VA", None, QtGui.QApplication.UnicodeUTF8))
        self.L_FGDVCtarCr.setText(QtGui.QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:14pt;\">FG-DVC Tar Cracking</span></p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.lE_FGDVCtarCr.setText(QtGui.QApplication.translate("PKP", "0", None, QtGui.QApplication.UnicodeUTF8))
        self.lE_numTimeStep.setText(QtGui.QApplication.translate("PKP", "1.e-4", None, QtGui.QApplication.UnicodeUTF8))
        self.L_numTimeStep.setText(QtGui.QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:14pt;\">numerical time step</span></p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.B_Launch.setText(QtGui.QApplication.translate("PKP", "Launch", None, QtGui.QApplication.UnicodeUTF8))
        self.B_Plot1.setText(QtGui.QApplication.translate("PKP", "Plot", None, QtGui.QApplication.UnicodeUTF8))
        self.B_Plot2.setText(QtGui.QApplication.translate("PKP", "Plot", None, QtGui.QApplication.UnicodeUTF8))
        self.B_Plot3.setText(QtGui.QApplication.translate("PKP", "Plot", None, QtGui.QApplication.UnicodeUTF8))
        self.B_Plot4.setText(QtGui.QApplication.translate("PKP", "Plot", None, QtGui.QApplication.UnicodeUTF8))
        self.B_Plot5.setText(QtGui.QApplication.translate("PKP", "Plot", None, QtGui.QApplication.UnicodeUTF8))
        self.B_Open4.setText(QtGui.QApplication.translate("PKP", "Open", None, QtGui.QApplication.UnicodeUTF8))
        self.B_Open2.setText(QtGui.QApplication.translate("PKP", "Open", None, QtGui.QApplication.UnicodeUTF8))
        self.B_Open5.setText(QtGui.QApplication.translate("PKP", "Open", None, QtGui.QApplication.UnicodeUTF8))
        self.B_Open3.setText(QtGui.QApplication.translate("PKP", "Open", None, QtGui.QApplication.UnicodeUTF8))
        self.B_Open1.setText(QtGui.QApplication.translate("PKP", "Open", None, QtGui.QApplication.UnicodeUTF8))
        self.Header1_2.setText(QtGui.QApplication.translate("PKP", "<html><head/><body><p align=\"center\"><span style=\" font-size:22pt; font-weight:600;\">PKP</span></p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.L_UA.setText(QtGui.QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:14pt;\">Ultimate Analysis in %</span></p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.L_UAC.setText(QtGui.QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:14pt;\">Carbon</span></p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.L_UAH.setText(QtGui.QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:14pt;\">Hydrogen</span></p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.L_UAN.setText(QtGui.QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:14pt;\">Nitrogen</span></p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.L_UAO.setText(QtGui.QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:14pt;\">Oxygen</span></p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.L_UAS.setText(QtGui.QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:14pt;\">Sulphur</span></p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.L_PA.setText(QtGui.QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:14pt;\">Proximate Analysis in %</span></p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.L_PAFC.setText(QtGui.QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:14pt;\">Fixed Carbon</span></p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.L_PAVM.setText(QtGui.QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:14pt;\">Volatile Matter</span></p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.L_PAMoi.setText(QtGui.QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:14pt;\">Moisture</span></p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.L_PAAsh.setText(QtGui.QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:14pt;\">Ash</span></p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
 	self.cB_ArrhSpec.setItemText(0, QtGui.QApplication.translate("PKP", "Total", None, QtGui.QApplication.UnicodeUTF8))
        self.cB_ArrhSpec.setItemText(1, QtGui.QApplication.translate("PKP", "Main Species", None, QtGui.QApplication.UnicodeUTF8))
        self.cB_ArrhSpec.setItemText(2, QtGui.QApplication.translate("PKP", "all Species", None, QtGui.QApplication.UnicodeUTF8))
        self.L_ArrhSpec.setText(QtGui.QApplication.translate("PKP", "<html><head/><body><p align=\"center\"><span style=\" font-size:14pt;\">selected Fit Species<br/>(Arrhenius)</span></p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.menuFile.setTitle(QtGui.QApplication.translate("PKP", "File", None, QtGui.QApplication.UnicodeUTF8))
        self.menuHelp.setTitle(QtGui.QApplication.translate("PKP", "Help", None, QtGui.QApplication.UnicodeUTF8))
        self.actionWrite_into_File.setText(QtGui.QApplication.translate("PKP", "Write into File", None, QtGui.QApplication.UnicodeUTF8))
        self.actionWrite_into_File.setShortcut(QtGui.QApplication.translate("PKP", "Ctrl+S", None, QtGui.QApplication.UnicodeUTF8))
        self.actionWrite_and_Run.setText(QtGui.QApplication.translate("PKP", "Write and Run", None, QtGui.QApplication.UnicodeUTF8))
        self.actionWrite_and_Run.setShortcut(QtGui.QApplication.translate("PKP", "Ctrl+R", None, QtGui.QApplication.UnicodeUTF8))
        self.actionExit.setText(QtGui.QApplication.translate("PKP", "Exit", None, QtGui.QApplication.UnicodeUTF8))
        self.actionExit.setShortcut(QtGui.QApplication.translate("PKP", "Ctrl+Q", None, QtGui.QApplication.UnicodeUTF8))
        self.actionShow_generated_Results.setText(QtGui.QApplication.translate("PKP", "Show generated Results", None, QtGui.QApplication.UnicodeUTF8))
        self.actionShow_generated_Results.setShortcut(QtGui.QApplication.translate("PKP", "Ctrl+T", None, QtGui.QApplication.UnicodeUTF8))
        self.actionShow_saved_state.setText(QtGui.QApplication.translate("PKP", "Load saved state", None, QtGui.QApplication.UnicodeUTF8))
        self.actionShow_saved_state.setShortcut(QtGui.QApplication.translate("PKP", "Ctrl+O", None, QtGui.QApplication.UnicodeUTF8))

    #manually added
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
        filename = QtGui.QFileDialog.getOpenFileName()
        file1=open(filename)
        data = file1.read()
        self.tE_THist_1.setText(data)
        file1.close()
    def LoadTtFile2(self):
        """Loads the temperature history nr 2 file via file browser"""
        filename = QtGui.QFileDialog.getOpenFileName()
        file2=open(filename)
        data = file2.read()
        self.tE_THist_2.setText(data)
        file2.close()
    def LoadTtFile3(self):
        """Loads the temperature history nr 3 file via file browser"""
        filename = QtGui.QFileDialog.getOpenFileName()
        file3=open(filename)
        data = file3.read()
        self.tE_THist_3.setText(data)
        file3.close()
    def LoadTtFile4(self):
        """Loads the temperature history nr 4 file via file browser"""
        filename = QtGui.QFileDialog.getOpenFileName()
        file4=open(filename)
        data = file4.read()
        self.tE_THist_4.setText(data)
        file4.close()
    def LoadTtFile5(self):
        """Loads the temperature history nr 5 file via file browser"""
        filename = QtGui.QFileDialog.getOpenFileName()
        file5=open(filename)
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
        self.Done.setupUi(SpeciesL,ProgramL,ProgrFitD)
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


if __name__ == "__main__":
    import sys
    Infosaver=InfosFromGUI()
    app = QtGui.QApplication(sys.argv)
    PKP = QtGui.QMainWindow()
    ui = Ui_PKP()
    ui.setupUi(PKP,Infosaver)
    PKP.show()
    sys.exit(app.exec_())

