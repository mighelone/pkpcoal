import sys
from PySide.QtGui import *
from PySide.QtCore import *
 

def PromptError(ErrorMassage=''):
  #"""Generates Window to raise an Error corresponding to the input Massage. Just use this and apply .show() to plot result."""
  """Raise an Error corresponding to the input Massage"""
  ErrorMassage=str(ErrorMassage)
  msgBox = QMessageBox()
  msgBox.setText('Input Error')
  msgBox.setInformativeText(ErrorMassage)
  msgBox.exec_()