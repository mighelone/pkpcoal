import os
os.environ['QT_API'] = 'pyside'
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Qt4Agg')
matplotlib.rcParams['backend.qt4']='PySide'
import numpy as np
import writeInfoFiles
import sys
import platform
import shutil

OSys=platform.system()
sys.path.append('src')
import PKP
import InformationFiles
import GUI_ErrorPrompt

import Done   #second GUI

from PySide.QtCore import * #QWidget, QMainWindow
from PySide.QtGui import *


#try:
#    _fromUtf8 = QtCore.QString.fromUtf8
#except AttributeError:
_fromUtf8 = lambda s: s

#class Ui_PKP(QMainWindow):
class Ui_PKP(QMainWindow):
    """ The general notation:
        - L_ is a label
        - cB_ is a column bar
        - lE_ is a line Edit
        - gL_ is a grid Layout
        - tE_ is a text Edit
        - sB_ is a scroll Bar
        - B_ is a bottom
        - lW_ is a (layout) QWidget
        """
    def setupUi(self, PKP, SaveInfoObj):
        #QMainWindow.__init__(self) #manually added
        self.QWidgMain= QWidget(self) # Main Window
        self.setCentralWidget(self.QWidgMain) # Central Widget
        # information
        self.TheCalculationsAreDone=False # Checking if new Window showing the results should be added
        self.svInfo=SaveInfoObj
        # General
        PKP.setObjectName(_fromUtf8("PKP"))
#        PKP.resize(1033, 753)
        PKP.showMaximized()
        self.centralwidget = QWidget(PKP)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        # Headlines
        self.Header1_2 = QLabel(self.centralwidget) # the central 'PKP' header
        self.Header1_2.setGeometry(QRect(170, 0, 721, 31))
        self.Header1_2.setObjectName(_fromUtf8("Header1_2"))
        self.Header1 = QLabel(self.centralwidget)   # the subheader 'Pyrolysis Programs and Fitting'
        self.Header1.setGeometry(QRect(-10, 50, 721, 31))
        self.Header1.setObjectName(_fromUtf8("Header1"))
        self.Header2 = QLabel(self.centralwidget)   # the subheader 'Coal Properties (as received)'
        self.Header2.setGeometry(QRect(30, 230, 721, 31))
        self.Header2.setObjectName(_fromUtf8("Header2"))
        self.Header3 = QLabel(self.centralwidget)    # the subheader 'Operating Conditions'
        self.Header3.setGeometry(QRect(20, 520, 721, 31))
        self.Header3.setObjectName(_fromUtf8("Header3"))
        # all Labels (standard text in the GUI)
        self.lW_PyrolPrFit = QWidget(self.centralwidget)
        self.lW_PyrolPrFit.setObjectName(_fromUtf8("lW_PyrolPrFit"))
        self.lW_PyrolPrFit1 = QWidget(self.lW_PyrolPrFit)
        self.lW_PyrolPrFit1.setObjectName(_fromUtf8("lW_PyrolPrFit1"))
        self.lW_PyrolPrFit2 = QWidget(self.lW_PyrolPrFit)
        self.lW_PyrolPrFit2.setObjectName(_fromUtf8("lW_PyrolPrFit2"))
        self.lW_PyrolPrFit3 = QWidget(self.lW_PyrolPrFit)
        self.lW_PyrolPrFit3.setObjectName(_fromUtf8("lW_PyrolPrFit3"))
        self.L_FGDVC = QLabel(self.lW_PyrolPrFit)
        self.L_FGDVC.setObjectName(_fromUtf8("L_FGDVC"))
        self.L_CPD = QLabel(self.lW_PyrolPrFit)
        self.L_CPD.setObjectName(_fromUtf8("L_CPD"))
        self.L_PCCL = QLabel(self.lW_PyrolPrFit)
        self.L_PCCL.setObjectName(_fromUtf8("L_PCCL"))
        self.L_Polimi = QLabel(self.lW_PyrolPrFit)
        self.L_Polimi.setObjectName(_fromUtf8("L_Polimi"))
        self.L_WeightParam = QLabel(self.lW_PyrolPrFit)
        self.L_WeightParam.setObjectName(_fromUtf8("L_WeightParam"))
        self.L_Yweight = QLabel(self.lW_PyrolPrFit)
        self.L_Yweight.setObjectName(_fromUtf8("L_Yweight"))
        self.L_Rweight = QLabel(self.lW_PyrolPrFit)
        self.L_Rweight.setObjectName(_fromUtf8("L_Rweight"))
        self.L_ArrhSpec1 = QLabel(self.lW_PyrolPrFit)                             
        self.L_ArrhSpec1.setObjectName(_fromUtf8("L_ArrhSpec1"))
        self.L_ArrhSpec2 = QLabel(self.lW_PyrolPrFit)                             
        self.L_ArrhSpec2.setObjectName(_fromUtf8("L_ArrhSpec2"))
        #
        # The input objects, field 1 'Pyrolysis Programs and Fitting'
        self.cB_CPD = QComboBox(self.lW_PyrolPrFit)
        self.cB_CPD.setObjectName(_fromUtf8("cB_CPD"))
        self.cB_CPD.addItem(_fromUtf8(""))
        self.cB_CPD.addItem(_fromUtf8(""))
        self.cB_CPD.addItem(_fromUtf8(""))
        self.cB_CPD.addItem(_fromUtf8(""))
        self.cB_CPD.addItem(_fromUtf8(""))
        self.cB_CPD.addItem(_fromUtf8(""))
        self.cB_CPD.addItem(_fromUtf8(""))
        self.cB_FGDVC = QComboBox(self.lW_PyrolPrFit)
        self.cB_FGDVC.setObjectName(_fromUtf8("cB_FGDVC"))
        self.cB_FGDVC.addItem(_fromUtf8(""))
        self.cB_FGDVC.addItem(_fromUtf8(""))
        self.cB_FGDVC.addItem(_fromUtf8(""))
        self.cB_FGDVC.addItem(_fromUtf8(""))
        self.cB_FGDVC.addItem(_fromUtf8(""))
        self.cB_FGDVC.addItem(_fromUtf8(""))
        self.cB_FGDVC.addItem(_fromUtf8(""))
        self.cB_PCCL = QComboBox(self.lW_PyrolPrFit)
        self.cB_PCCL.setObjectName(_fromUtf8("cB_PCCL"))
        self.cB_PCCL.addItem(_fromUtf8(""))
        self.cB_PCCL.addItem(_fromUtf8(""))
        self.cB_PCCL.addItem(_fromUtf8(""))
        self.cB_PCCL.addItem(_fromUtf8(""))
        self.cB_PCCL.addItem(_fromUtf8(""))
        self.cB_PCCL.addItem(_fromUtf8(""))
        self.cB_PCCL.addItem(_fromUtf8(""))
        self.cB_Polimi = QComboBox(self.lW_PyrolPrFit)
        self.cB_Polimi.setObjectName(_fromUtf8("cB_Polimi"))
        self.cB_Polimi.addItem(_fromUtf8(""))
        self.cB_Polimi.addItem(_fromUtf8(""))
        self.cB_Polimi.addItem(_fromUtf8(""))
        self.cB_Polimi.addItem(_fromUtf8(""))
        self.cB_Polimi.addItem(_fromUtf8(""))
        self.cB_Polimi.addItem(_fromUtf8(""))
        self.cB_Polimi.addItem(_fromUtf8(""))
        self.cB_ArrhSpec = QComboBox(self.lW_PyrolPrFit)                         
        self.cB_ArrhSpec.setObjectName(_fromUtf8("cB_ArrhSpec"))                       
        self.cB_ArrhSpec.addItem(_fromUtf8(""))                                        
        self.cB_ArrhSpec.addItem(_fromUtf8(""))                                        
        self.cB_ArrhSpec.addItem(_fromUtf8(""))
        self.lE_Yweight = QLineEdit(self.lW_PyrolPrFit)
        self.lE_Yweight.setObjectName(_fromUtf8("lE_Yweight"))
        self.lE_Rweight = QLineEdit(self.lW_PyrolPrFit)
        self.lE_Rweight.setObjectName(_fromUtf8("lE_Rweight"))
        # grid for Part 'Pyrolyiss Program and Fitting'
        self.gL_PyrolPrFit = QGridLayout(self.lW_PyrolPrFit)
        self.gL_PyrolPrFit1 = QGridLayout(self.lW_PyrolPrFit1)
        self.gL_PyrolPrFit2 = QGridLayout(self.lW_PyrolPrFit2)
        self.gL_PyrolPrFit3 = QGridLayout(self.lW_PyrolPrFit3)
        self.gL_PyrolPrFit1.addWidget(self.L_CPD,1,1,1,1)
        self.gL_PyrolPrFit1.addWidget(self.cB_CPD,1,2,1,1)
        self.gL_PyrolPrFit1.addWidget(self.L_FGDVC,2,1,1,1)
        self.gL_PyrolPrFit1.addWidget(self.cB_FGDVC,2,2,1,1)
        self.gL_PyrolPrFit1.addWidget(self.L_PCCL,3,1,1,1)
        self.gL_PyrolPrFit1.addWidget(self.cB_PCCL,3,2,1,1)
        self.gL_PyrolPrFit1.addWidget(self.L_Polimi,4,1,1,1)
        self.gL_PyrolPrFit1.addWidget(self.cB_Polimi,4,2,1,1)
        self.gL_PyrolPrFit2.addWidget(self.L_ArrhSpec1,1,4,1,2)
        self.gL_PyrolPrFit2.addWidget(self.cB_ArrhSpec,2,4,1,1)
        self.gL_PyrolPrFit2.addWidget(self.L_ArrhSpec2,3,4,1,2)
        self.gL_PyrolPrFit3.addWidget(self.L_WeightParam,1,7,1,2)
        self.gL_PyrolPrFit3.addWidget(self.L_Yweight,2,7,1,1)
        self.gL_PyrolPrFit3.addWidget(self.L_Rweight,3,7,1,1)
        self.gL_PyrolPrFit3.addWidget(self.lE_Yweight,2,8,1,1)
        self.gL_PyrolPrFit3.addWidget(self.lE_Rweight,3,8,1,1)
        self.gL_PyrolPrFit.addWidget(self.lW_PyrolPrFit1,1,1,1,2)
        self.gL_PyrolPrFit.setColumnMinimumWidth(3,50)
        self.gL_PyrolPrFit.addWidget(self.lW_PyrolPrFit2,1,4,1,1)
        self.gL_PyrolPrFit.setColumnMinimumWidth(5,50)
        self.gL_PyrolPrFit.addWidget(self.lW_PyrolPrFit3,1,6,1,1)
        #
        # The input objects, field 2 'Coal Properties'
        self.lW_Coal = QWidget(self.centralwidget)
        self.lW_Coal.setObjectName(_fromUtf8("lW_Coal"))
        self.L_UA = QLabel(self.lW_Coal)
        self.L_UA.setObjectName(_fromUtf8("L_UA"))
        self.L_UAC = QLabel(self.lW_Coal)
        self.L_UAC.setObjectName(_fromUtf8("L_UAC"))
        self.lE_UAC = QLineEdit(self.lW_Coal)
        self.lE_UAC.setObjectName(_fromUtf8("lE_UAC"))
        self.L_UAH = QLabel(self.lW_Coal)
        self.L_UAH.setObjectName(_fromUtf8("L_UAH"))
        self.lE_UAH = QLineEdit(self.lW_Coal)
        self.lE_UAH.setObjectName(_fromUtf8("lE_UAH"))
        self.L_UAN = QLabel(self.lW_Coal)
        self.L_UAN.setObjectName(_fromUtf8("L_UAN"))
        self.lE_UAN = QLineEdit(self.lW_Coal)
        self.lE_UAN.setObjectName(_fromUtf8("lE_UAN"))
        self.L_UAO = QLabel(self.lW_Coal)
        self.L_UAO.setObjectName(_fromUtf8("L_UAO"))
        self.lE_UAO = QLineEdit(self.lW_Coal)
        self.lE_UAO.setObjectName(_fromUtf8("lE_UAO"))
        self.L_UAS = QLabel(self.lW_Coal)
        self.L_UAS.setObjectName(_fromUtf8("L_UAS"))
        self.lE_UAS = QLineEdit(self.lW_Coal)
        self.lE_UAS.setObjectName(_fromUtf8("lE_UAS"))
        self.L_PA = QLabel(self.lW_Coal)
        self.L_PA.setObjectName(_fromUtf8("L_PA"))
        self.L_PAFC = QLabel(self.lW_Coal)
        self.L_PAFC.setObjectName(_fromUtf8("L_PAFC"))
        self.lE_PAFC = QLineEdit(self.lW_Coal)
        self.lE_PAFC.setObjectName(_fromUtf8("lE_PAFC"))
        self.L_PAVM = QLabel(self.lW_Coal)
        self.L_PAVM.setObjectName(_fromUtf8("L_PAVM"))
        self.lE_PAVM = QLineEdit(self.lW_Coal)
        self.lE_PAVM.setObjectName(_fromUtf8("lE_PAVM"))
        self.L_PAMoi = QLabel(self.lW_Coal)
        self.L_PAMoi.setObjectName(_fromUtf8("L_PAMoi"))
        self.lE_PAMoi = QLineEdit(self.lW_Coal)
        self.lE_PAMoi.setObjectName(_fromUtf8("lE_PAMoi"))
        self.L_PAAsh = QLabel(self.lW_Coal)
        self.L_PAAsh.setObjectName(_fromUtf8("L_PAAsh"))
        self.lE_PAAsh = QLineEdit(self.lW_Coal)
        self.lE_PAAsh.setObjectName(_fromUtf8("lE_PAAsh"))
        self.L_MW = QLabel(self.lW_Coal)
        self.L_MW.setObjectName(_fromUtf8("L_MW"))
        self.L_MWTar = QLabel(self.lW_Coal)
        self.L_MWTar.setObjectName(_fromUtf8("L_MWTar"))
        self.lE_MWTar = QLineEdit(self.lW_Coal)
        self.lE_MWTar.setObjectName(_fromUtf8("lE_MWTar"))
        self.lE_HHV = QLineEdit(self.lW_Coal)
        self.L_HHV = QLabel(self.lW_Coal)
        self.L_HHV.setObjectName(_fromUtf8("L_HHV"))
        self.lE_HHV.setObjectName(_fromUtf8("lE_HHV"))
        self.lE_FGDVCtarCr = QLineEdit(self.lW_Coal)
        self.lE_FGDVCtarCr.setObjectName(_fromUtf8("lE_FGDVCtarCr"))
        self.L_FGDVCcoal = QLabel(self.lW_Coal)
        self.L_FGDVCcoal.setObjectName(_fromUtf8("L_FGDVCcoal"))
        self.L_FGDVCtarCr = QLabel(self.lW_Coal)
        self.L_FGDVCtarCr.setObjectName(_fromUtf8("L_FGDVCtarCr"))
        self.cB_FGDVCcoal = QComboBox(self.lW_Coal)
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
        self.L_PartDia1 = QLabel(self.lW_Coal)
        self.L_PartDia1.setObjectName(_fromUtf8("L_PartDia1"))
        self.L_PartDia2 = QLabel(self.lW_Coal)
        self.L_PartDia2.setObjectName(_fromUtf8("L_PartDia2"))
        self.lE_PartDia = QLineEdit(self.lW_Coal)
        self.lE_PartDia.setObjectName(_fromUtf8("lE_PartDia"))
        # grid for Part 'Coal Properties'
        self.lW_Coal1 = QWidget(self.lW_Coal)
        self.lW_Coal1.setObjectName(_fromUtf8("lW_Coal1"))
        self.lW_Coal2 = QWidget(self.lW_Coal)
        self.lW_Coal2.setObjectName(_fromUtf8("lW_Coal2"))
        self.lW_Coal3 = QWidget(self.lW_Coal)
        self.lW_Coal3.setObjectName(_fromUtf8("lW_Coal3"))
        self.lW_Coal4 = QWidget(self.lW_Coal)
        self.lW_Coal4.setObjectName(_fromUtf8("lW_Coal4"))
        self.gL_Coal = QGridLayout(self.lW_Coal)
        self.gL_Coal1 = QGridLayout(self.lW_Coal1)
#        self.gL_Coal2 = QGridLayout(self.lW_Coal2)
        self.gL_Coal3 = QGridLayout(self.lW_Coal3)
        self.gL_Coal4 = QGridLayout(self.lW_Coal4)
        self.gL_Coal1.addWidget(self.L_UA,1,1,1,2)
        self.gL_Coal1.addWidget(self.L_UAC,2,1,1,1)
        self.gL_Coal1.addWidget(self.lE_UAC,2,2,1,1)
        self.gL_Coal1.addWidget(self.L_UAH,3,1,1,1)
        self.gL_Coal1.addWidget(self.lE_UAH,3,2,1,1)
        self.gL_Coal1.addWidget(self.L_UAO,4,1,1,1)
        self.gL_Coal1.addWidget(self.lE_UAO,4,2,1,1)
        self.gL_Coal1.addWidget(self.L_UAN,5,1,1,1)
        self.gL_Coal1.addWidget(self.lE_UAN,5,2,1,1)
        self.gL_Coal1.addWidget(self.L_UAS,6,1,1,1)
        self.gL_Coal1.addWidget(self.lE_UAS,6,2,1,1)
        self.gL_Coal1.setColumnMinimumWidth(3,50)
        self.gL_Coal1.addWidget(self.L_PA,1,4,1,2)
        self.gL_Coal1.addWidget(self.L_PAFC,2,4,1,1)
        self.gL_Coal1.addWidget(self.lE_PAFC,2,5,1,1)
        self.gL_Coal1.addWidget(self.L_PAVM,3,4,1,1)
        self.gL_Coal1.addWidget(self.lE_PAVM,3,5,1,1)
        self.gL_Coal1.addWidget(self.L_PAMoi,4,4,1,1)
        self.gL_Coal1.addWidget(self.lE_PAMoi,4,5,1,1)
        self.gL_Coal1.addWidget(self.L_PAAsh,5,4,1,1)
        self.gL_Coal1.addWidget(self.lE_PAAsh,5,5,1,1)
        self.gL_Coal3.addWidget(self.L_MW,1,1,1,2)
        self.gL_Coal3.addWidget(self.L_MWTar,2,1,1,1)
        self.gL_Coal3.addWidget(self.lE_MWTar,2,2,1,1)
        self.gL_Coal3.addWidget(self.L_HHV,4,1,1,2)
        self.gL_Coal3.addWidget(self.lE_HHV,5,2,1,1)
        self.gL_Coal4.addWidget(self.L_FGDVCcoal,1,1,1,2)
        self.gL_Coal4.addWidget(self.cB_FGDVCcoal,1,3,1,2)
        self.gL_Coal4.addWidget(self.L_FGDVCtarCr,2,1,1,2)
        self.gL_Coal4.addWidget(self.lE_FGDVCtarCr,2,3,1,2)
        self.gL_Coal4.setColumnMinimumWidth(5,50)
        self.gL_Coal4.addWidget(self.L_PartDia1,1,6,1,2)
        self.gL_Coal4.addWidget(self.L_PartDia2,2,6,1,1)
        self.gL_Coal4.addWidget(self.lE_PartDia,2,7,1,1)
        self.gL_Coal.addWidget(self.lW_Coal1,1,1,1,2)
#        self.gL_Coal.addWidget(self.lW_Coal2,1,3,1,2)
        self.gL_Coal.addWidget(self.lW_Coal3,1,5,1,1)
        self.gL_Coal.addWidget(self.lW_Coal4,2,1,1,3)
        self.gL_Coal4.setColumnMinimumWidth(6,20)
        #
        # The input objects, field 3 'Operating Conditions'
        self.lW_OperCond = QWidget(self.centralwidget)
        self.lW_OperCond.setObjectName(_fromUtf8("lW_OperCond"))
        #
        self.L_pressure = QLabel(self.lW_OperCond)
        self.L_pressure.setObjectName(_fromUtf8("L_pressure"))
        self.lE_pressure = QLineEdit(self.lW_OperCond)
        self.lE_pressure.setObjectName(_fromUtf8("lE_pressure"))
        self.L_THist1 = QLabel(self.lW_OperCond)
        self.L_THist1.setObjectName(_fromUtf8("L_THist1"))
        self.L_THist2 = QLabel(self.lW_OperCond)
        self.L_THist2.setObjectName(_fromUtf8("L_THist2"))
        self.tE_THist_1 = QTextEdit(self.lW_OperCond)
        self.tE_THist_1.setObjectName(_fromUtf8("tE_THist_1"))
        self.tE_THist_2 = QTextEdit(self.lW_OperCond)
        self.tE_THist_2.setObjectName(_fromUtf8("tE_THist_2"))
        self.tE_THist_3 = QTextEdit(self.lW_OperCond)
        self.tE_THist_3.setObjectName(_fromUtf8("tE_THist_3"))
        self.tE_THist_4 = QTextEdit(self.lW_OperCond)
        self.tE_THist_4.setObjectName(_fromUtf8("tE_THist_4"))
        self.tE_THist_5 = QTextEdit(self.lW_OperCond)
        self.tE_THist_5.setObjectName(_fromUtf8("tE_THist_5"))
        self.sB_Nr_THist = QSpinBox(self.lW_OperCond)
        self.sB_Nr_THist.setObjectName(_fromUtf8("sB_Nr_THist"))
        self.sB_Nr_THist.setMinimum(1)
        self.sB_Nr_THist.setMaximum(5)
        self.L_numTimeStep = QLabel(self.lW_PyrolPrFit) 
        self.L_numTimeStep.setObjectName(_fromUtf8("L_numTimeStep"))
        self.lE_numTimeStep = QLineEdit(self.lW_OperCond)
        self.lE_numTimeStep.setObjectName(_fromUtf8("lE_numTimeStep"))
        self.B_Plot1 = QPushButton(self.lW_OperCond)
        self.B_Plot1.setObjectName(_fromUtf8("B_Plot1"))
        self.B_Plot2 = QPushButton(self.lW_OperCond)
        self.B_Plot2.setObjectName(_fromUtf8("B_Plot2"))
        self.B_Plot3 = QPushButton(self.lW_OperCond)
        self.B_Plot3.setObjectName(_fromUtf8("B_Plot3"))
        self.B_Plot4 = QPushButton(self.lW_OperCond)
        self.B_Plot4.setObjectName(_fromUtf8("B_Plot4"))
        self.B_Plot5 = QPushButton(self.lW_OperCond)
        self.B_Plot5.setObjectName(_fromUtf8("B_Plot5"))
        self.B_Open4 = QPushButton(self.lW_OperCond)
        self.B_Open4.setObjectName(_fromUtf8("B_Open4"))
        self.B_Open2 = QPushButton(self.lW_OperCond)
        self.B_Open2.setObjectName(_fromUtf8("B_Open2"))
        self.B_Open5 = QPushButton(self.lW_OperCond)
        self.B_Open5.setObjectName(_fromUtf8("B_Open5"))
        self.B_Open3 = QPushButton(self.lW_OperCond)
        self.B_Open3.setObjectName(_fromUtf8("B_Open3"))
        self.B_Open1 = QPushButton(self.lW_OperCond)
        self.B_Open1.setObjectName(_fromUtf8("B_Open1"))
        self.B_Launch = QPushButton(self.lW_OperCond)
        self.B_Launch.setObjectName(_fromUtf8("B_Launch"))
        # grid for Part 'Operating Conditions'
        self.lW_OperCond1 = QWidget(self.lW_OperCond)
        self.lW_OperCond2 = QWidget(self.lW_OperCond)
        self.gL_OperCond = QGridLayout(self.lW_OperCond)
        self.gL_OperCond1 = QGridLayout(self.lW_OperCond1)
        self.gL_OperCond2 = QGridLayout(self.lW_OperCond2)
        self.gL_OperCond1.addWidget(self.L_pressure,1,1,1,1)
        self.gL_OperCond1.addWidget(self.lE_pressure,2,1,1,1)
#        self.gL_OperCond1.setColumnMinimumWidth(3,260)
        self.gL_OperCond1.setRowMinimumHeight(3,10)
        self.gL_OperCond1.addWidget(self.L_numTimeStep,4,1,1,1)
        self.gL_OperCond1.addWidget(self.lE_numTimeStep,5,1,1,1)
#        self.gL_OperCond1.setColumnMinimumWidth(6,260)
        self.gL_OperCond2.addWidget(self.L_THist1,2,1,1,2)
        self.gL_OperCond2.addWidget(self.L_THist2,3,1,1,1)
        self.gL_OperCond2.addWidget(self.sB_Nr_THist,3,2,1,1)
        self.gL_OperCond2.addWidget(self.tE_THist_1,1,3,3,2)
        self.gL_OperCond2.addWidget(self.B_Open1,4,3,1,1)
        self.gL_OperCond2.addWidget(self.B_Plot1,4,4,1,1)
        self.gL_OperCond2.addWidget(self.tE_THist_2,1,5,3,2)
        self.gL_OperCond2.addWidget(self.B_Open2,4,5,1,1)
        self.gL_OperCond2.addWidget(self.B_Plot2,4,6,1,1)
        self.gL_OperCond2.addWidget(self.tE_THist_3,1,7,3,2)
        self.gL_OperCond2.addWidget(self.B_Open3,4,7,1,1)
        self.gL_OperCond2.addWidget(self.B_Plot3,4,8,1,1)
        self.gL_OperCond2.addWidget(self.tE_THist_4,1,9,3,2)
        self.gL_OperCond2.addWidget(self.B_Open4,4,9,1,1)
        self.gL_OperCond2.addWidget(self.B_Plot4,4,10,1,1)
        self.gL_OperCond2.addWidget(self.tE_THist_5,1,11,3,2)
        self.gL_OperCond2.addWidget(self.B_Open5,4,11,1,1)
        self.gL_OperCond2.addWidget(self.B_Plot5,4,12,1,1)
        self.gL_OperCond2.addWidget(self.B_Launch,2,13,2,2)
        self.gL_OperCond.addWidget(self.lW_OperCond1,1,1,1,1)
        self.gL_OperCond.addWidget(self.lW_OperCond2,1,2,3,6)
        #
        #
        # Pictures
        self.lW_Figs = QWidget(self.centralwidget)
        self.lW_Figs.setObjectName(_fromUtf8("lW_Figs"))
        #
        self.PictureVirtuh = QLabel(self.centralwidget)
#        self.PictureVirtuh.setGeometry(QRect(750, 380, 231, 171))
        self.PictureVirtuh.setText(_fromUtf8(""))
        self.PictureVirtuh.setPixmap(QPixmap(_fromUtf8("Logos/virtuhcon_logo.png")))
        self.PictureVirtuh.setScaledContents(True)
        self.PictureVirtuh.setObjectName(_fromUtf8("PictureVirtuh"))
        self.PictureNTFD = QLabel(self.centralwidget)
#        self.PictureNTFD.setGeometry(QRect(750, 80, 231, 171))
        self.PictureNTFD.setText(_fromUtf8(""))
        self.PictureNTFD.setPixmap(QPixmap(_fromUtf8("Logos/ntfd_rgb.png")))
        self.PictureNTFD.setScaledContents(True)
        self.PictureNTFD.setObjectName(_fromUtf8("PictureNTFD"))
        self.gL_Figs = QGridLayout(self.lW_Figs)
        self.gL_Figs.addWidget(self.Header1_2,1,1,2,6)
        self.gL_Figs.addWidget(self.PictureNTFD,2,2,4,4)
        self.gL_Figs.addWidget(self.PictureVirtuh,9,2,4,4)
        #
        #
        # Menubar
        self.menubar = QMenuBar(PKP)
#        self.menubar.setGeometry(QRect(0, 0, 1033, 21))
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
        #
        # Grid Layout (for whole GUI)        
        self.MainGrid = QGridLayout(self.centralwidget)
        self.MainGrid.setObjectName(_fromUtf8("MainGrid"))
#        self.MainGrid.addWidget(self.Header1_2,1,1,1,4)
        self.MainGrid.addWidget(self.Header1,1,1,1,2)
        self.MainGrid.addWidget(self.lW_PyrolPrFit,2,1,3,2)
        self.MainGrid.addWidget(self.lW_Figs,1,4,10,1)
        self.MainGrid.addWidget(self.Header2,5,1,1,2)
        self.MainGrid.addWidget(self.lW_Coal,6,1,5,3)
        self.MainGrid.addWidget(self.Header3,11,1,1,3)
        self.MainGrid.addWidget(self.lW_OperCond,12,1,3,4)
        #
        PKP.setCentralWidget(self.centralwidget)
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
        PKP.setTabOrder(self.cB_PCCL, self.cB_Polimi)
        PKP.setTabOrder(self.cB_Polimi, self.cB_ArrhSpec)
        PKP.setTabOrder(self.cB_ArrhSpec, self.lE_Yweight)
        PKP.setTabOrder(self.lE_Yweight, self.lE_Rweight)
        PKP.setTabOrder(self.lE_Rweight, self.lE_UAC)
        PKP.setTabOrder(self.lE_UAC, self.lE_UAH)
        PKP.setTabOrder(self.lE_UAH, self.lE_UAO)
        PKP.setTabOrder(self.lE_UAO, self.lE_UAN)
        PKP.setTabOrder(self.lE_UAN, self.lE_UAS)
        PKP.setTabOrder(self.lE_UAS, self.lE_PAFC)
        PKP.setTabOrder(self.lE_PAFC, self.lE_PAVM)
        PKP.setTabOrder(self.lE_PAVM, self.lE_PAMoi)
        PKP.setTabOrder(self.lE_PAMoi, self.lE_PAAsh)
        PKP.setTabOrder(self.lE_PAAsh, self.lE_MWTar)
        PKP.setTabOrder(self.lE_MWTar, self.lE_HHV)
        PKP.setTabOrder(self.lE_HHV, self.cB_FGDVCcoal)
        PKP.setTabOrder(self.cB_FGDVCcoal, self.lE_FGDVCtarCr)
        PKP.setTabOrder(self.lE_FGDVCtarCr, self.lE_PartDia)
        PKP.setTabOrder(self.lE_PartDia, self.lE_pressure)
        PKP.setTabOrder(self.lE_pressure, self.lE_numTimeStep)
        PKP.setTabOrder(self.lE_numTimeStep, self.sB_Nr_THist)
        PKP.setTabOrder(self.sB_Nr_THist, self.B_Open1)
        PKP.setTabOrder(self.B_Open1, self.B_Plot1)
        PKP.setTabOrder(self.B_Plot1, self.B_Open2)
        PKP.setTabOrder(self.B_Open2, self.B_Plot2)
        PKP.setTabOrder(self.B_Plot2, self.B_Open3)
        PKP.setTabOrder(self.B_Open3, self.B_Plot3)
        PKP.setTabOrder(self.B_Plot3, self.B_Open4)
        PKP.setTabOrder(self.B_Open4, self.B_Plot4)
        PKP.setTabOrder(self.B_Plot4, self.B_Open5)
        PKP.setTabOrder(self.B_Open5, self.B_Plot5)
        PKP.setTabOrder(self.B_Plot5, self.B_Launch)


    def retranslateUi(self, PKP):
        """Sets all the text style and connects the bottons with their individual functions"""
        PKP.setWindowTitle(QApplication.translate("PKP", "Pyrolysis Kinetic Preprocessor - Central Input Window", None, QApplication.UnicodeUTF8))
        self.Header1_2.setText(QApplication.translate("PKP", "<html><head/><body><p align=\"center\"><span style=\" font-size:40pt; font-weight:600;\">PKP</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.Header1.setText(QApplication.translate("PKP", "<html><head/><body><p align=\"center\"><span style=\" font-size:14pt; font-weight:600;\">Pyrolysis Programs and Fitting</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.Header2.setText(QApplication.translate("PKP", "<html><head/><body><p align=\"center\"><span style=\" font-size:14pt; font-weight:600;\">Coal Properties (as received)</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.L_HHV.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:13pt;\">Higher Heating Value<br/>in MJ/kg</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.Header3.setText(QApplication.translate("PKP", "<html><head/><body><p align=\"center\"><span style=\" font-size:14pt; font-weight:600;\">Operating Conditions</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.L_pressure.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:13pt;\">Pressure in atm</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.lE_pressure.setText(QApplication.translate("PKP", "1", None, QApplication.UnicodeUTF8))
        self.L_THist1.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:13pt;\">Temperature History</span></p><p></span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.L_THist2.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:10pt;\">t in s ,   T in K</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.L_CPD.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:13pt;\">CPD</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.cB_CPD.setItemText(0, QApplication.translate("PKP", "None", None, QApplication.UnicodeUTF8))
        self.cB_CPD.setItemText(1, QApplication.translate("PKP", "Run", None, QApplication.UnicodeUTF8))
        self.cB_CPD.setItemText(2, QApplication.translate("PKP", "constant Rate", None, QApplication.UnicodeUTF8))
        self.cB_CPD.setItemText(3, QApplication.translate("PKP", "Arrhenius", None, QApplication.UnicodeUTF8))
        self.cB_CPD.setItemText(4, QApplication.translate("PKP", "Arrhenius no B", None, QApplication.UnicodeUTF8))
        self.cB_CPD.setItemText(5, QApplication.translate("PKP", "Kobayashi", None, QApplication.UnicodeUTF8))
        self.cB_CPD.setItemText(6, QApplication.translate("PKP", "DAEM", None, QApplication.UnicodeUTF8))
        self.L_FGDVC.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:13pt;\">FG-DVC</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.cB_FGDVC.setItemText(0, QApplication.translate("PKP", "None", None, QApplication.UnicodeUTF8))
        self.cB_FGDVC.setItemText(1, QApplication.translate("PKP", "Run", None, QApplication.UnicodeUTF8))
        self.cB_FGDVC.setItemText(2, QApplication.translate("PKP", "constant Rate", None, QApplication.UnicodeUTF8))
        self.cB_FGDVC.setItemText(3, QApplication.translate("PKP", "Arrhenius", None, QApplication.UnicodeUTF8))
        self.cB_FGDVC.setItemText(4, QApplication.translate("PKP", "Arrhenius No B", None, QApplication.UnicodeUTF8))
        self.cB_FGDVC.setItemText(5, QApplication.translate("PKP", "Kobayashi", None, QApplication.UnicodeUTF8))
        self.cB_FGDVC.setItemText(6, QApplication.translate("PKP", "DAEM", None, QApplication.UnicodeUTF8))
        self.L_PCCL.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:13pt;\">PC Coal Lab</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.cB_PCCL.setItemText(0, QApplication.translate("PKP", "None", None, QApplication.UnicodeUTF8))
        self.cB_PCCL.setItemText(1, QApplication.translate("PKP", "Run", None, QApplication.UnicodeUTF8))
        self.cB_PCCL.setItemText(2, QApplication.translate("PKP", "constant Rate", None, QApplication.UnicodeUTF8))
        self.cB_PCCL.setItemText(3, QApplication.translate("PKP", "Arrhenius", None, QApplication.UnicodeUTF8))
        self.cB_PCCL.setItemText(4, QApplication.translate("PKP", "Arrhenius No B", None, QApplication.UnicodeUTF8))
        self.cB_PCCL.setItemText(5, QApplication.translate("PKP", "Kobayashi", None, QApplication.UnicodeUTF8))
        self.cB_PCCL.setItemText(6, QApplication.translate("PKP", "DAEM", None, QApplication.UnicodeUTF8))
        self.L_Polimi.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:13pt;\">Polimi</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.cB_Polimi.setItemText(0, QApplication.translate("PKP", "None", None, QApplication.UnicodeUTF8))
        self.cB_Polimi.setItemText(1, QApplication.translate("PKP", "Run", None, QApplication.UnicodeUTF8))
        self.cB_Polimi.setItemText(2, QApplication.translate("PKP", "constant Rate", None, QApplication.UnicodeUTF8))
        self.cB_Polimi.setItemText(3, QApplication.translate("PKP", "Arrhenius", None, QApplication.UnicodeUTF8))
        self.cB_Polimi.setItemText(4, QApplication.translate("PKP", "Arrhenius No B", None, QApplication.UnicodeUTF8))
        self.cB_Polimi.setItemText(5, QApplication.translate("PKP", "Kobayashi", None, QApplication.UnicodeUTF8))
        self.cB_Polimi.setItemText(6, QApplication.translate("PKP", "DAEM", None, QApplication.UnicodeUTF8))        
        self.L_WeightParam.setText(QApplication.translate("PKP", "<html><head/><body><p align=\"center\"><span style=\" font-size:13pt;\">Weight Parameter</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.L_Yweight.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:13pt;\">Yields</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.lE_Yweight.setText(QApplication.translate("PKP", "1", None, QApplication.UnicodeUTF8))
        self.L_Rweight.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:13pt;\">Rates</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.lE_Rweight.setText(QApplication.translate("PKP", "1", None, QApplication.UnicodeUTF8))
        self.lE_PartDia.setText(QApplication.translate("PKP", "100.", None, QApplication.UnicodeUTF8))
        self.L_MW.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:13pt;\">Molecule Weight</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.L_MWTar.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:13pt;\">Tar</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.L_PartDia1.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:13pt;\">Particle Diameter (PCCL)</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.L_PartDia2.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:10pt;\">in micrometer</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.L_FGDVCcoal.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:13pt;\">FG-DVC Coal #</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.cB_FGDVCcoal.setItemText(0, QApplication.translate("PKP", "0 - Interpolate", None, QApplication.UnicodeUTF8))
        self.cB_FGDVCcoal.setItemText(1, QApplication.translate("PKP", "1 - Beulah-Zap", None, QApplication.UnicodeUTF8))
        self.cB_FGDVCcoal.setItemText(2, QApplication.translate("PKP", "2 - Woydak-Anderson", None, QApplication.UnicodeUTF8))
        self.cB_FGDVCcoal.setItemText(3, QApplication.translate("PKP", "3 - Illinois # 6", None, QApplication.UnicodeUTF8))
        self.cB_FGDVCcoal.setItemText(4, QApplication.translate("PKP", "4 - Bind Canyon, UT", None, QApplication.UnicodeUTF8))
        self.cB_FGDVCcoal.setItemText(5, QApplication.translate("PKP", "5 - Lewis-Stockton, WV", None, QApplication.UnicodeUTF8))
        self.cB_FGDVCcoal.setItemText(6, QApplication.translate("PKP", "6 - Pittsburgh # 8", None, QApplication.UnicodeUTF8))
        self.cB_FGDVCcoal.setItemText(7, QApplication.translate("PKP", "7 - Upper Freeport, PA", None, QApplication.UnicodeUTF8))
        self.cB_FGDVCcoal.setItemText(8, QApplication.translate("PKP", "8 - Pocahontas, VA", None, QApplication.UnicodeUTF8))
        self.L_FGDVCtarCr.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:13pt;\">FG-DVC Tar Cracking</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.lE_FGDVCtarCr.setText(QApplication.translate("PKP", "0", None, QApplication.UnicodeUTF8))
        self.lE_numTimeStep.setText(QApplication.translate("PKP", "1.e-4", None, QApplication.UnicodeUTF8))
        self.L_numTimeStep.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:13pt;\">numerical time step</span></p></body></html>", None, QApplication.UnicodeUTF8))
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
        self.L_UA.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:13pt;\">Ultimate Analysis in %</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.L_UAC.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:13pt;\">Carbon</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.L_UAH.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:13pt;\">Hydrogen</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.L_UAN.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:13pt;\">Nitrogen</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.L_UAO.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:13pt;\">Oxygen</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.L_UAS.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:13pt;\">Sulphur</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.L_PA.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:13pt;\">Proximate Analysis in %</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.L_PAFC.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:13pt;\">Fixed Carbon</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.L_PAVM.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:13pt;\">Volatile Matter</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.L_PAMoi.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:13pt;\">Moisture</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.L_PAAsh.setText(QApplication.translate("PKP", "<html><head/><body><p><span style=\" font-size:13pt;\">Ash</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.cB_ArrhSpec.setItemText(0, QApplication.translate("PKP", "Total", None, QApplication.UnicodeUTF8))
        self.cB_ArrhSpec.setItemText(1, QApplication.translate("PKP", "Main Species", None, QApplication.UnicodeUTF8))
        self.cB_ArrhSpec.setItemText(2, QApplication.translate("PKP", "all Species", None, QApplication.UnicodeUTF8))
        self.L_ArrhSpec1.setText(QApplication.translate("PKP", "<html><head/><body><p align=\"center\"><span style=\" font-size:13pt;\">selected Species</span></p></body></html>", None, QApplication.UnicodeUTF8))
        self.L_ArrhSpec2.setText(QApplication.translate("PKP", "<html><head/><body><p align=\"center\"><span style=\" font-size:11pt;\">Arrhenius Fits only</span></p></body></html>", None, QApplication.UnicodeUTF8))
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
        Polimisel= self.cB_Polimi.currentIndex()
        self.svInfo.SetRunPyrolProg(CPDsel,FGsel,PCCLsel,Polimisel)
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
        PCCLParticleDiam = str(self.lE_PartDia.text())
        self.svInfo.setPCCLParticleSize(PCCLParticleDiam)
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
        writeInfoFiles.WritePCCLFile(self.svInfo)
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
        Tt=np.genfromtxt('TempHist1.dat')
        fig = plt.figure()
        ax = fig.add_subplot(111)
        if len(Tt[:,0])<50:
            ax.plot(Tt[:,0],Tt[:,1],'-x')
        else:
            ax.plot(Tt[:,0],Tt[:,1],'-')
        ax.set_title('Temperature History')
        ax.set_xlabel('Time in s')
        ax.set_ylabel('Temperature in K')
        ax.grid()
        plt.show()
    def Plot2(self):
        """Plots the temperature over time history (temperature history nr 2) and saves temperatuer history in "TempHist2.dat"."""
        file2=open('TempHist2.dat','w')
        file2.write(self.tE_THist_2.toPlainText())
        file2.close()
        Tt=np.genfromtxt('TempHist2.dat')
        fig = plt.figure()
        ax = fig.add_subplot(111)
        if len(Tt[:,0])<50:
            ax.plot(Tt[:,0],Tt[:,1],'-x')
        else:
            ax.plot(Tt[:,0],Tt[:,1],'-')
        ax.set_title('Temperature History')
        ax.set_xlabel('Time in s')
        ax.set_ylabel('Temperature in K')
        ax.grid()
        plt.show()
    def Plot3(self):
        """Plots the temperature over time history (temperature history nr 3) and saves temperatuer history in "TempHist3.dat"."""
        file3=open('TempHist3.dat','w')
        file3.write(self.tE_THist_3.toPlainText())
        file3.close()
        Tt=np.genfromtxt('TempHist3.dat')
        fig = plt.figure()
        ax = fig.add_subplot(111)
        if len(Tt[:,0])<50:
            ax.plot(Tt[:,0],Tt[:,1],'-x')
        else:
            ax.plot(Tt[:,0],Tt[:,1],'-')
        ax.set_title('Temperature History')
        ax.set_xlabel('Time in s')
        ax.set_ylabel('Temperature in K')
        ax.grid()
        plt.show()
    def Plot4(self):
        """Plots the temperature over time history (temperature history nr 4) and saves temperatuer history in "TempHist4.dat"."""
        file4=open('TempHist4.dat','w')
        file4.write(self.tE_THist_4.toPlainText())
        file4.close()
        Tt=np.genfromtxt('TempHist4.dat')
        fig = plt.figure()
        ax = fig.add_subplot(111)
        if len(Tt[:,0])<50:
            ax.plot(Tt[:,0],Tt[:,1],'-x')
        else:
            ax.plot(Tt[:,0],Tt[:,1],'-')
        ax.set_title('Temperature History')
        ax.set_xlabel('Time in s')
        ax.set_ylabel('Temperature in K')
        ax.grid()
        plt.show()
    def Plot5(self):
        """Plots the temperature over time history (temperature history nr 5) and saves temperatuer history in "TempHist5.dat"."""
        file5=open('TempHist5.dat','w')
        file5.write(self.tE_THist_5.toPlainText())
        file5.close()
        Tt=np.genfromtxt('TempHist5.dat')
        fig = plt.figure()
        ax = fig.add_subplot(111)
        if len(Tt[:,0])<50:
            ax.plot(Tt[:,0],Tt[:,1],'-x')
        else:
            ax.plot(Tt[:,0],Tt[:,1],'-')
        ax.set_title('Temperature History')
        ax.set_xlabel('Time in s')
        ax.set_ylabel('Temperature in K')
        ax.grid()
        plt.show()

    def WriteRun(self):
        """Writes the *.inp files and launch the process."""
        # checks whether the input is consistent
        InpIsOk = self.__checkInput()
        PCCLIsOk = self.__PCCLInputIsOk()
        if InpIsOk != 'True':
            self.RaiseError(InpIsOk)
        # checks whether PCCL input is ok, raise otherwise error
        elif (self.cB_PCCL.currentIndex()!=0) and (PCCLIsOk != 'True'):
                self.RaiseError('Please insert for PC Coal Lab\n-a linear heating ramp\n-with a const. hold temperature.\nError in Temperature input '+PCCLIsOk)
        else:
            # tells user if pressure for PCCL is too low
            if self.cB_PCCL.currentIndex()!=0:
                if float(str(self.lE_pressure.text())) < 0.101325:
                    self.RaiseError('The selected pressure for PC Coal Lab is too low. Minimum input pressure is 0.01MPa (0.101325atm). Continue with this value.')
            print 'Checked Input in GUI. Was ok'
            self.SaveInfos()
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
            if Case.PCCL_select==True:
                Case.MakeResults_PCCL()
                ProgramL.append('PCCL')
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
        
    def RaiseError(self,ErrorMassage):
        """Raise an Window with the Error massage inputted"""
        GUI_ErrorPrompt.PromptError(ErrorMassage)
    
    def __checkInput(self):
        """Checks wheather the general input is ok. If not returns the wrong T-t input field info."""
        isOk = 'True'
        # checks input fields:
        try:
            float(str(self.lE_Yweight.text()))
        except ValueError:
            return 'No valid input for Weight-Parameter (Yields).'
        try:
            float(str(self.lE_Rweight.text()))
        except ValueError:
            return 'No valid input for Weight-Parameter (Rates).'
        try:
            float(str(self.lE_UAC.text()))
        except ValueError:
            return 'No valid input for Ultimate Analysis (Carbon).'
        try:
            float(str(self.lE_UAH.text()))
        except ValueError:
            return 'No valid input for Ultimate Analysis (Hydrogen).'
        try:
            float(str(self.lE_UAN.text()))
        except ValueError:
            return 'No valid input for Ultimate Analysis (Nitrogen).'
        try:
            float(str(self.lE_UAO.text()))
        except ValueError:
            return 'No valid input for Ultimate Analysis (Oxygen).'
        try:
            float(str(self.lE_UAS.text()))
        except ValueError:
            return 'No valid input for Ultimate Analysis (Sulfur).'
        try:
            float(str(self.lE_PAFC.text()))
        except ValueError:
            return 'No valid input for Proximate Analysis (Fixed Carbon).'
        try:
            float(str(self.lE_PAVM.text()))
        except ValueError:
            return 'No valid input for Proximate Analysis (Volatile Matter).'
        try:
            float(str(self.lE_PAMoi.text()))
        except ValueError:
            return 'No valid input for Proximate Analysis (Moisture).'
        try:
            float(str(self.lE_PAAsh.text()))
        except ValueError:
            return 'No valid input for Proximate Analysis (Ash).'
        try:
            float(str(self.lE_MWTar.text()))
        except ValueError:
            return 'No valid input for Tar Molecule Weight.'
        try:
            float(str(self.lE_HHV.text()))
        except ValueError:
            return 'No valid input for Higher Heating Value.'
        try:
            float(str(self.lE_FGDVCtarCr.text()))
        except ValueError:
            return 'No valid input for FG-DVC Tar cracking option. Enter -1 for Full, 0 for no tar cracking and a value larger zero for the tar residence time.'
        try:
            float(str(self.lE_PartDia.text()))
        except ValueError:
            return 'No valid input for Particle Diameter.'
        try:
            float(str(self.lE_pressure.text()))
        except ValueError:
            return 'No valid input for pressure.'
        try:
            float(str(self.lE_numTimeStep.text()))
        except ValueError:
            return 'No valid input for numerical time step.'
        # checks temperature histories
        i = 0 ; j = 0
        L = [self.tE_THist_1,self.tE_THist_2,self.tE_THist_3,self.tE_THist_4,self.tE_THist_5]
        for i in range(self.sB_Nr_THist.value()):
            tT = str(L[i].toPlainText()).split() #e.g. ['0.0', '300', '0.01', '1500', '0.08', '1500']
            if len(tT)%2 != 0:
                isOk = 'Temperture history #'+str(i+1)+':Odd number of input values.Please Input t[s] T[K] linewise separated by space.'
                return isOk
            elif True in map(lambda j: ',' in tT[j] , range(len(tT))): # checks wheather comma is in input
                isOk = 'Temperture history #'+str(i+1)+':Please use no comma, only a space character to separate time and temperature.'
            elif float(tT[0])!=0.0: # t_{i} != 0
                isOk = 'Temperture history #'+str(i+1)+': first time must be zero'
            elif False in (map(lambda j: float(tT[j])>=float(tT[j-2]) , range(2,len(tT),2))): # t_{i} < t_{i-2} -> False
                isOk = 'Temperture history #'+str(i+1)+': time in line '+str((j+1)/2+2)+' > '+str((j+3)/2+2)
                return isOk
        return isOk
          
    def __PCCLInputIsOk(self):
        """Checks wheather the input is ok for PC Coal Lab. If not returns the wrong input field info."""
        isOk = 'True'
        i = 0
        L = [self.tE_THist_1,self.tE_THist_2,self.tE_THist_3,self.tE_THist_4,self.tE_THist_5]
        for i in range(self.sB_Nr_THist.value()):
            tT = str(L[i].toPlainText()).split() #e.g. ['0.0', '300', '0.01', '1500', '0.08', '1500']
            if (len(tT)>6):
                isOk = str(i+1)+':\n Too many data points.'
                return isOk
            elif  len(tT)<4:
                isOk = str(i+1)+':\n Too few data points.'
                return isOk
            elif len(tT)>3:
                if tT[3] != tT[5]:
                    isOk = str(i+1)+':\n No constant hold temperature.'
                    return isOk
        return isOk

        
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
        # PCCL Properties
        PCCLInput=InformationFiles.ReadFile('PCCL.inp')
        PCCL_FittingKineticParameter_Select=PCCLInput.Fitting(InformationFiles.M_selFit)
        # !!! Update
        #
        #
        #Operating Condition File:
        OpCondInp=InformationFiles.OperCondInput('OperCond.inp')
        ArrhSpec=OpCondInp.getText(InformationFiles.M_selArrhSpec)
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
        if PCCL_FittingKineticParameter_Select == None:
            PCCL_FittingKineticParameter_Select='None'
        self.cB_CPD.setCurrentIndex(self.svInfo.RunPyrolProgReverse(CPD_FittingKineticParameter_Select))
        self.cB_FGDVC.setCurrentIndex(self.svInfo.RunPyrolProgReverse(FG_FittingKineticParameter_Select))
        self.cB_PCCL.setCurrentIndex(self.svInfo.RunPyrolProgReverse(PCCL_FittingKineticParameter_Select))
#        self.cB_PCCL.setCurrentIndex(self.svInfo.RunPyrolProgReverse(PCCL_FittingKineticParameter_Select))
        #
        self.cB_ArrhSpec.setCurrentIndex(self.svInfo.ArrhSpecReverse(ArrhSpec))
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
    
    def SetRunPyrolProg(self,CPDIndex,FGDVCIndex,PCCLIndex,PolimiIndex):
        """Saves which options of the three Pyrolysis programs are used."""
        FitDict={0:'None',1:'Run',2:'constantRate',3:'Arrhenius',4:'ArrheniusNoB',5:'Kobayashi',6:'DAEM'}
        self.__CPDsel    = FitDict[CPDIndex]
        self.__FGsel     = FitDict[FGDVCIndex]
        self.__PCCLsel   = FitDict[PCCLIndex]
        self.__Polimisel = FitDict[PolimiIndex]
        
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
        
#    def setPCCLCoalCalibrFactor(self,CoalCalibrFactor):
#        """Sets the PC Coal Lab coal Calibration Factor. If input is float or string of float then it is saved as str(float). Otherwise (e.g. input None or False) is saved as None."""
#        try:
#            self.__PCCLCoalCalibr = str(float(CoalCalibrFactor))
#        except ValueError:
#            self.__PCCLCoalCalibr = 'None'
#        
#    def PCCLCoalCalibrFactor(self):
#        """Returns the PC Coal Lab coal Calibration Factor."""
#        return self.__PCCLCoalCalibr
    
    def setPCCLParticleSize(self,ParticleDiamInMicrometer):
        """Defines the Pc Coal Lab Particle Size. Input is in microMeter"""
        self.__PCCLParticleDiam = ParticleDiamInMicrometer
    
    def PCCLParticleSize(self):
        """Returns the Pc Coal Lab Particle Size. Output is in microMeter"""
        return self.__PCCLParticleDiam
        
        
        


if __name__ == "__main__":
    import sys
    Infosaver=InfosFromGUI()
    app = QApplication(sys.argv)
    PKPWindow = QMainWindow()
    ui = Ui_PKP()
    ui.setupUi(PKPWindow,Infosaver)
    PKPWindow.show()
    sys.exit(app.exec_())

