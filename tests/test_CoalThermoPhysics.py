import numpy as np
import pytest

import pkp.src.CPD as cpdsl
from mock_data import *
from test_main import (linear_preProc, 
    cpd_linear, cpd_run_single, 
    cpd_run_multi, cpd_run_ifrf,
    almost)


def test_PostulateSubstanceMethane():
    from pkp.src.CoalThermoPhysics import Coal, R
    from pkp.src.CoalThermoPhysics import PostulateSubstance
    coal = Coal(methane_as_coal['Coal'])
    methane = PostulateSubstance(coal=coal,qFactor=1.0)
    massComp = methane.VolatileCompositionMass()
    molComp =  methane.VolatileCompositionMol()
    assert almost(massComp['Carbon'], 0.75)
    assert almost(massComp['Hydrogen'], 0.25)
    assert massComp['Oxygen'] == 0.0
    assert almost(molComp['Carbon'], 1.0)
    assert almost(molComp['Hydrogen'], 4.0)
    assert molComp['Oxygen'] == 0.0
    molCompProd = methane.ProductCompositionMol(partialOx=False)
    assert almost(molCompProd['H2O'], 2.0)
    assert almost(molCompProd['CO2'], 1.0)
    molCompProdPart = methane.ProductCompositionMol(partialOx=True)
    assert almost(molCompProdPart['H2O'], 2.0)
    assert almost(molCompProdPart['CO'], 1.0)
    print ""
    print "methane mechanism partial: ", methane.mechanism()
    print "methane mechanism full: ", methane.mechanism(partialOx=False)
    h0f_CH4 = methane.EnthalpyOfFormation()
    # check if less then 5% deviation from book value
    assert abs((-74831.0 - h0f_CH4)/h0f_CH4) < 0.05
    print "methane a6 for OpenFOAM: ", h0f_CH4 / R

def test_PostulateSubstancePropane():
    from pkp.src.CoalThermoPhysics import Coal, R
    from pkp.src.CoalThermoPhysics import PostulateSubstance
    coal = Coal(propane_as_coal['Coal'])
    methane = PostulateSubstance(coal=coal,qFactor=1.0)
    massComp = methane.VolatileCompositionMass()
    molComp =  methane.VolatileCompositionMol()
    assert almost(massComp['Carbon'], 0.8181)
    assert almost(massComp['Hydrogen'], 0.1818)
    assert massComp['Oxygen'] == 0.0
    assert almost(molComp['Carbon'], 3.0)
    assert almost(molComp['Hydrogen'], 8.0)
    assert molComp['Oxygen'] == 0.0
    molCompProd = methane.ProductCompositionMol(partialOx=False)
    molCompProdPart = methane.ProductCompositionMol(partialOx=True)
    print ""
    print "Propane mechanism partial: ", methane.mechanism()
    print "Propane mechanism full: ", methane.mechanism(partialOx=False)
    h0f_CH4 = methane.EnthalpyOfFormation()
    # check if less then 5% deviation from book value
    print "propane a6 for OpenFOAM: ", h0f_CH4 / R
