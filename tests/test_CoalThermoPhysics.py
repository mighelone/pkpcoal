import numpy as np
import pytest

import pkp.src.CPD as cpdsl
from mock_data import *
from test_main import (linear_preProc, 
    cpd_linear, cpd_run_single, 
    cpd_run_multi, cpd_run_ifrf,
    almost)

def test_PostulateSubstanceQfactor():
    """ test if changing the qfactor conserves the mass """
    from pkp.src.CoalThermoPhysics import (Coal, R,
                PostulateSubstance)
    from pkp.src.PreProc import ManualQfactor 
    coal = ManualQfactor(
            Coal(mock_ps['Coal']),
            qFactor=1.0)
    post_subs = PostulateSubstance(coal)
    massComp = post_subs.VolatileCompositionMass
    assert massComp['Carbon'] < 1e-64
    for q in [1.1, 1.5, 2.0]:
        coal = ManualQfactor(
                Coal(mock_ps['Coal']),
                qFactor=q)
        post_subs = PostulateSubstance(coal)
        massComp = post_subs.VolatileCompositionMass
        #TODO FIX THIS
        # assert q == 100.0/(massComp['Hydrogen'] 
        #                   + massComp['Oxygen'])
    

def test_PostulateSubstanceCO():
    from pkp.src.CoalThermoPhysics import (Coal, R,
             PostulateSubstance, enthalpy_balance,
             enthalpy_balance_mol)
    from pkp.src.PreProc import ManualQfactor 
    coal = ManualQfactor(
                Coal(co_as_coal['Coal']),
                qFactor=1.0)
    methane = PostulateSubstance(coal)
    massComp = methane.VolatileCompositionMass
    molComp =  methane.VolatileCompositionMol
    assert almost(massComp['Carbon'], 12.0/28.0*100.0)
    assert massComp['Hydrogen'] == 0.0
    assert almost(massComp['Oxygen'], 16.0/28.0*100.0)
    # FULL OX: CH4 + 2O2 -> CO2 + 2H2O
    # mass basis
    massCompProd = methane.ProductCompositionMass(partialOx=False)
    assert massCompProd['N2'] == 0.0
    assert massCompProd['H2O'] == 0.0 # 36/80
    assert almost(massCompProd['CO2'], 100.0) # 44/80
    # 1kg CH4 produces 5kg Products
    FTP = methane.Fuel_to_Product_Mass
    assert almost(FTP, 28.0/44.0) 
    assert almost(methane.Oxygen_from_yield, 0.5)
    FOx = methane.FOx
    assert almost(FOx, 28.0/16.0)
    hp = methane.enthalpy_balance_products
    
    hp_CH4 = (-393546.0/44.0)*FTP
    hp_CH4_mol = -393546.0
    assert almost(enthalpy_balance_mol(
         {'CO2': 1.0}), hp_CH4_mol)
    
    assert almost(enthalpy_balance(
         {'CO2': 1.0})*FTP, hp_CH4)
    assert almost(hp, hp_CH4)

    # 
    #print ""
    #print "methane mechanism partial: ", methane.mechanism()
    #print "methane mechanism full: ", methane.mechanism(partialOx=False)
    h0f_CH4 = methane.EnthalpyOfFormation()
    # check if less then 5% deviation from book value
    #assert abs((-74831.0 - h0f_CH4)/h0f_CH4) < 0.05
    #print "methane a6 for OpenFOAM: ", h0f_CH4 / R

# def test_PostulateSubstancePropane():
#     from pkp.src.CoalThermoPhysics import Coal, R
#     from pkp.src.CoalThermoPhysics import PostulateSubstance
#     coal = Coal(propane_as_coal['Coal'])
#     methane = PostulateSubstance(coal=coal, qFactor=1.0)
#     massComp = methane.VolatileCompositionMass
#     molComp =  methane.VolatileCompositionMol
#     assert almost(massComp['Carbon'], 0.8181)
#     assert almost(massComp['Hydrogen'], 0.1818)
#     assert massComp['Oxygen'] == 0.0
#     assert almost(molComp['Carbon'], 3.0)
#     assert almost(molComp['Hydrogen'], 8.0)
#     assert molComp['Oxygen'] == 0.0
#     molCompProd = methane.ProductCompositionMol(partialOx=False)
#     molCompProdPart = methane.ProductCompositionMol(partialOx=True)
#     print ""
#     print "Propane mechanism partial: ", methane.mechanism()
#     print "Propane mechanism full: ", methane.mechanism(partialOx=False)
#     h0f_CH4 = methane.EnthalpyOfFormation()
#     # check if less then 5% deviation from book value
#     print "propane a6 for OpenFOAM: ", h0f_CH4 / R
