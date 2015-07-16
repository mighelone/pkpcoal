import numpy as np
import pytest

import pkp.src.CPD as cpdsl
from mock_data import *
from test_main import (linear_preProc,
    cpd_linear, cpd_run_single,
    cpd_run_multi, cpd_run_ifrf,
    almost)

def test_PostulateSubstance_conservativeMassComposition():
    """ test if changing the qfactor conserves the mass
        and composition """
    from pkp.src.CoalThermoPhysics import YieldFromCoalInp
    postSubs = YieldFromCoalInp(HydrogenOnlyVM['Coal'], 1.0, [''])
    # All the coal carbon is fixed carbon for a q of 1.0
    # hence the carbon in the volatiles should be zero
    assert postSubs.VolatileCompositionMass['Carbon'] < 1e-64
    assert postSubs.VolatileCompositionMass['Hydrogen'] == 100.0

    # increasing the q-factor will introduce carbon mass
    # to the yield, initially carbon in VM is 0.0.
    for q in [1.05, 1.1, 1.2, 1.5]:
        ps = YieldFromCoalInp(HydrogenOnlyVM['Coal'], q, [''])
        massComp = ps.VolatileCompositionMass
        # changing the q factor should not produce any oxygen
        assert massComp['Oxygen'] == 0.0
        # carbon and hydrogen should add up to 100%
        assert almost(massComp['Carbon'] + massComp['Hydrogen'], 100.0)
        # the q-Factor is the ratio between old and new yield
        # mass. q=m_yield_new/m_yield_old

        sc = 100.0/massComp['Hydrogen']
        assert almost(sc, q)

def test_PostulateSubstance_conservativeMassCompositionII():
    """ test if changing the qfactor conserves the mass
        and composition """
    from pkp.src.CoalThermoPhysics import YieldFromCoalInp
    HydrogenOnlyVM['Coal']['Ultimate Analysis']['Hydrogen'] = 1.0/3.0
    HydrogenOnlyVM['Coal']['Ultimate Analysis']['Oxygen'] = 1.0/3.0
    HydrogenOnlyVM['Coal']['Ultimate Analysis']['Carbon'] = 1.0/3.0
    HydrogenOnlyVM['Coal']['Proximate Analysis']['Fixed Carbon'] = 1.0/3.0
    HydrogenOnlyVM['Coal']['Proximate Analysis']['Volatile Matter'] = 2.0/3.0
    postSubs = YieldFromCoalInp(HydrogenOnlyVM['Coal'], 1.0, [''])
    # All the coal carbon is fixed carbon for a q of 1.0
    # hence the carbon in the volatiles should be zero
    assert postSubs.VolatileCompositionMass['Carbon'] < 1e-64
    assert postSubs.VolatileCompositionMass['Hydrogen'] == 50.0
    assert postSubs.VolatileCompositionMass['Oxygen'] == 50.0

    # increasing the q-factor will introduce carbon mass
    # to the yield, initially carbon in VM is 0.0.
    for q in [1.05, 1.1, 1.2, 1.5]:
        ps = YieldFromCoalInp(HydrogenOnlyVM['Coal'], q, [''])
        massComp = ps.VolatileCompositionMass
        # changing the q factor should not produce any oxygen
        # carbon and hydrogen should add up to 100%
        assert almost(
                  massComp['Oxygen']
                + massComp['Carbon']
                + massComp['Hydrogen'], 100.0)
        # the q-Factor is the ratio between old and new yield
        # mass. q=m_yield_new/m_yield_old

        sc = 100.0/(massComp['Hydrogen']+massComp['Oxygen'])
        assert almost(sc, q)

@pytest.mark.xfail
def test_GeneratePostulatSubstance():
    from pkp.src.CoalThermoPhysics import (EnthOfForm,
        generatePostulateSubstanceSpecies, PreProcFromCoalInp)
    coPreProc = PreProcFromCoalInp(co_as_coal['Coal'])
    CxHyOz = generatePostulateSubstanceSpecies(coPreProc)
    print CxHyOz.components
    # test if molecular has been read from coal dict
    assert CxHyOz.molecular_weight == co_as_coal['Coal']['MW_PS']
    # the mock co data has 100 percents volatiles
    # hence the postulate substance should consist
    # only of 12/28 C and 16/28 O

    assert CxHyOz.components['Carbon'] == 12.0/28.0*100.0
    assert CxHyOz.components['Oxygen'] == 16.0/28.0*100.0
    assert CxHyOz.components['Hydrogen'] == 0.0

    # Failing needs fixing
    assert CxHyOz.hf_perMol == EnthOfForm['CO']

def test_YieldCOComposition():
    """ Test standard methods of the postulate substance class
        for a known species """
    from pkp.src.CoalThermoPhysics import (
        enthalpy_balance, enthalpy_balance_mol, YieldFromCoalInp, CO)
    co = YieldFromCoalInp(co_as_coal['Coal'], 1.0, [CO])
    massComp = co.VolatileCompositionMass
    molComp = co.VolatileCompositionMol
    # Preliminary test if compostion is identical to expected
    # co composition from
    # test_PostulateSubstance_conservativeMassComposition
    # it should be clear that postulate substance conserves the
    # composition, but it shoould also be clear that co is matched
    assert almost(massComp['Carbon'], 12.0/28.0*100.0)
    assert massComp['Hydrogen'] == 0.0
    assert almost(massComp['Oxygen'], 16.0/28.0*100.0)

    co.composition = [1.0] # here we mock a composition fit
    assert almost(co.MW, 28.0)

    # Test if the predicted product composition is correct
    # FULL OX: CO + 0.5O2 -> CO2
    # mass basis
    massCompProd = co.ProductCompositionMass(partialOx=False)
    assert massCompProd['N2'] == 0.0
    assert massCompProd['H2O'] == 0.0 # 36/80
    assert almost(massCompProd['CO2'], 100.0) # 4_        4/80
    # TODO test also molar composition

    # 1 kg CO forms 44.0/28.0 kg of products (CO2)
    # FTP = [kg_yield/kg_Products]
    assert almost(co.Fuel_to_Product_Mass, 28.0/44.0)

    assert almost(co.ProductsMass(False)['CO2'], 44.0/28.0)
    assert co.ProductsMass(False)['H2O'] == 0.0

    # 1 mol CO forms 1 mol of CO2 and 0 H2O
    assert almost(co.ProductMols(False)['CO2'], 1.0)
    assert co.ProductMols(False)['H2O'] == 0.0

    # 50 % of CO2 oxygen is from the fuel
    assert almost(co.Oxygen_from_yield, 0.5)

    # 1 kg CO needs 0.5 * 32.0/28.0 O2
    assert almost(co.FOx, 28.0/16.0)

    # 1 mol CO needs 0.5 mol O2
    assert almost(co.FOx_molar, 0.5)

    # print "CO mechanism partial: ", co.mechanism()
    print "CO mechanism full: ", co.mechanism(partialOx=False)


def test_YieldCOEnthalpies():
    """ Test standard methods of the postulate substance class
        for a known species """
    from pkp.src.CoalThermoPhysics import (
        enthalpy_balance, enthalpy_balance_mol, YieldFromCoalInp, CO)
    co = YieldFromCoalInp(co_as_coal['Coal'], 1.0, [CO])
    co.composition = [1.0] # here we mock a composition fit
    # compute the enthalpy balance of the products
    # and compare with known enthalpy of CO
    # per kg fuel
    h_products = co.EnthalpyBalanceProds
    # per mol fuel
    h_products_mol = co.EnthalpyBalanceProdsMol
    FTP = co.Fuel_to_Product_Mass
    # The enthalpy balance of CO products is
    # h_fCO2 * 1 mol
    h_products_co_mol = -393546.0
    assert almost(enthalpy_balance_mol({'CO2': 1.0}), h_products_co_mol)
    # to get get enthalpy the enthalpy balance
    # on the per kg_yield basis we have h_fCO/MW_CO/FTP
    #  kJ/kg_yield
    h_products_co = h_products_co_mol/44.0/FTP

    assert almost(enthalpy_balance({'CO2': 1.0})/FTP, h_products_co)

    assert almost(h_products, h_products_co)
    assert almost(h_products_mol, h_products_co_mol)

    h_fCO = -110541.0 # kJ/kmol
    assert almost(co.EnthalpyOfFormation, h_fCO/28.0)
    assert almost(co.EnthalpyOfFormationMol, h_fCO)
    # test if alias delivers same result
    assert co.pre.coal.hhv_daf == co.hhv_daf
    # test if CO has expected HHV in kJ/kg
    assert co.hhv_daf == 10159.968
    # CO produces no H2O -> HHV and LHV are identical
    assert co.hhv_daf == co.lhv_daf
    # Since the mock species has no char LHV_daf and
    # LHV_yield should be identical
    assert co.lhv_daf == co.lhv_yield
    # and the heat released through char conversion
    # should be zero
    assert co.pre.ftot == 1.0
    assert co.q_release_char_yield == 0.0
    # Finally the enthalpy of formation should close to zero
    # allow less then 0.01 percent deviation between
    # LHV_yield and LHV_yield_species
    print "co.EnthalpyOfDevol" , co.EnthalpyOfDevol
    assert abs(co.EnthalpyOfDevol)/co.lhv_yield < 0.01

def test_YieldCharCH4Enthalpies():
    """ Test standard methods of the postulate substance class
        for a known species """
    from pkp.src.CoalThermoPhysics import (
        enthalpy_balance, enthalpy_balance_mol, YieldFromCoalInp, CH4)
    methane = YieldFromCoalInp(charCH4_as_coal['Coal'], 1.0, [CH4])
    methane.composition = [1.0] # here we mock a composition fit
    massComp = methane.VolatileCompositionMass
    molComp = methane.VolatileCompositionMol
    massCompProd = methane.ProductCompositionMass(partialOx=False)
    assert almost(massComp['Carbon'], 12.0/16.0*100.0)
    assert massComp['Oxygen'] == 0.0
    assert almost(massComp['Hydrogen'], 4.0/16.0*100.0)

    methane.composition = [1.0] # here we mock a composition fit
    assert almost(methane.MW, 16.0)

    # Test if the predicted product composition is correct
    # FULL OX: CH4 + 2O2 -> CO2 + 2H2O
    # mass basis
    massCompProd = methane.ProductCompositionMass(partialOx=False)
    assert massCompProd['N2'] == 0.0
    assert almost(massCompProd['H2O'], 36.0/80.0*100.0) # 4_        4/80
    assert almost(massCompProd['CO2'], 44.0/80.0*100.0) # 4_        4/80

    # 1 kg CH4 forms 80.0/16.0 kg of products (CO2)
    # FTP = [kg_yield/kg_Products]
    assert almost(methane.Fuel_to_Product_Mass, 16.0/80.0)

    assert almost(methane.ProductsMass(False)['CO2'], 44.0/16.0)
    assert almost(methane.ProductsMass(False)['H2O'], 36.0/16.0)

    # 1 mol CO forms 1 mol of CO2 and 0 H2O
    assert almost(methane.ProductMols(False)['CO2'], 1.0)
    assert almost(methane.ProductMols(False)['H2O'], 2.0)

    # 0 % of CO2 oxygen is from the fuel
    assert methane.Oxygen_from_yield == 0.0

    # 1 kg CH4 needs 2.0 * 32.0/16.0 O2
    assert almost(methane.FOx, 16.0/64.0)

    # 1 mol CH4 needs 2 mol O2
    assert almost(methane.FOx_molar, 2.0)

    # compute the enthalpy balance of the products
    # and compare with known enthalpy of CO
    # per kg fuel
    h_products = methane.EnthalpyBalanceProds
    # per mol fuel
    h_products_mol = methane.EnthalpyBalanceProdsMol
    FTP = methane.Fuel_to_Product_Mass
    # The enthalpy balance of CO products is
    # h_fCO2 * 1 mol
    h_products_methane_mol = -393546.0 + 2.0 * -241845.0
    assert almost(enthalpy_balance_mol({'CO2': 1.0, 'H2O':2.0}),
                  h_products_methane_mol)
    # to get get enthalpy the enthalpy balance
    # on the per kg_yield basis we have h_fCO/MW_CO/FTP
    #  kJ/kg_yield
    #h_products_methane = h_products_co_mol/44.0/FTP

    # assert almost(enthalpy_balance({'CO2': 1.0})/FTP, h_products_co)
    #
    # assert almost(h_products, h_products_co)
    assert almost(h_products_mol, h_products_methane_mol)

    h_fCH4 = -74831 # kJ/kmol
    assert almost(methane.EnthalpyOfFormation, h_fCH4/16.0)
    assert almost(methane.EnthalpyOfFormationMol, h_fCH4)

    print "CH4 mechanism partial: ", methane.mechanism()
    print "CH4 mechanism full: ", methane.mechanism(partialOx=False)
    # here we have the mixture 0.5 char and 0.5 methane
    assert methane.hhv_daf == 0.5*55528.0 + 0.5*32666.6
    # allow less then 100.0 kJ/kg deviation between
    # LHV_yield and LHV_yield_species
    assert almost(methane.lhv_yield, 50016.0)
    assert almost(methane.lhv_yield_species, 50016.0)
    assert abs(methane.EnthalpyOfDevol)/55528.0 < 0.01

def test_YieldCH4():
    """ Test standard methods of the postulate substance class
        for a known species """
    from pkp.src.CoalThermoPhysics import (
        enthalpy_balance, enthalpy_balance_mol, YieldFromCoalInp, CH4)
    methane = YieldFromCoalInp(methane_as_coal['Coal'], 1.0, [CH4])
    massComp = methane.VolatileCompositionMass
    molComp = methane.VolatileCompositionMol
    # Preliminary test if compostion is identical to expected
    # co composition from
    # test_PostulateSubstance_conservativeMassComposition
    # it should be clear that postulate substance conserves the
    # composition, but it shoould also be clear that co is matched
    assert almost(massComp['Carbon'], 12.0/16.0*100.0)
    assert massComp['Oxygen'] == 0.0
    assert almost(massComp['Hydrogen'], 4.0/16.0*100.0)

    methane.composition = [1.0] # here we mock a composition fit
    assert almost(methane.MW, 16.0)

    # Test if the predicted product composition is correct
    # FULL OX: CH4 + 2O2 -> CO2 + 2H2O
    # mass basis
    massCompProd = methane.ProductCompositionMass(partialOx=False)
    assert massCompProd['N2'] == 0.0
    assert almost(massCompProd['H2O'], 36.0/80.0*100.0) # 4_        4/80
    assert almost(massCompProd['CO2'], 44.0/80.0*100.0) # 4_        4/80

    # 1 kg CH4 forms 80.0/16.0 kg of products (CO2)
    # FTP = [kg_yield/kg_Products]
    assert almost(methane.Fuel_to_Product_Mass, 16.0/80.0)

    assert almost(methane.ProductsMass(False)['CO2'], 44.0/16.0)
    assert almost(methane.ProductsMass(False)['H2O'], 36.0/16.0)

    # 1 mol CO forms 1 mol of CO2 and 0 H2O
    assert almost(methane.ProductMols(False)['CO2'], 1.0)
    assert almost(methane.ProductMols(False)['H2O'], 2.0)

    # 0 % of CO2 oxygen is from the fuel
    assert methane.Oxygen_from_yield == 0.0

    # 1 kg CH4 needs 2.0 * 32.0/16.0 O2
    assert almost(methane.FOx, 16.0/64.0)

    # 1 mol CH4 needs 2 mol O2
    assert almost(methane.FOx_molar, 2.0)

    # compute the enthalpy balance of the products
    # and compare with known enthalpy of CO
    # per kg fuel
    h_products = methane.EnthalpyBalanceProds
    # per mol fuel
    h_products_mol = methane.EnthalpyBalanceProdsMol
    FTP = methane.Fuel_to_Product_Mass
    # The enthalpy balance of CO products is
    # h_fCO2 * 1 mol
    h_products_methane_mol = -393546.0 + 2.0 * -241845.0
    assert almost(enthalpy_balance_mol({'CO2': 1.0, 'H2O':2.0}),
                  h_products_methane_mol)
    # to get get enthalpy the enthalpy balance
    # on the per kg_yield basis we have h_fCO/MW_CO/FTP
    #  kJ/kg_yield
    #h_products_methane = h_products_co_mol/44.0/FTP

    # assert almost(enthalpy_balance({'CO2': 1.0})/FTP, h_products_co)
    #
    # assert almost(h_products, h_products_co)
    assert almost(h_products_mol, h_products_methane_mol)

    h_fCH4 = -74831 # kJ/kmol
    assert almost(methane.EnthalpyOfFormation, h_fCH4/16.0)
    assert almost(methane.EnthalpyOfFormationMol, h_fCH4)

    print "CH4 mechanism partial: ", methane.mechanism()
    print "CH4 mechanism full: ", methane.mechanism(partialOx=False)
    assert methane.hhv_daf == 55528.0
    assert almost(methane.lhv_daf, 50016.0)
    assert almost(methane.lhv_yield, 50016.0)
    assert almost(methane.hhv_yield, 55528.0)
    # allow less then 100.0 kJ/kg deviation between
    # LHV_yield and LHV_yield_species
    assert abs(methane.EnthalpyOfDevol)/55528.0 < 0.01

def test_YieldEnthalpyOfFormation():
    """ Test standard methods of the postulate substance class
        for a known species """
    from pkp.src.CoalThermoPhysics import (
        enthalpy_balance, enthalpy_balance_mol, YieldFromCoalInp, CH4)
    from copy import deepcopy
    mod_methane = deepcopy(methane_as_coal['Coal'])
    hhv_orig = methane_as_coal['Coal']['hhv']
    h_devol_orig = 0.0
    for delta in [0, 1, 10, 1000, 10000]:
        # if we lower the target heating value
        # the enthalpy of devolatilisation should increase
        # since ..
        mod_methane['hhv'] = hhv_orig - delta
        methane = YieldFromCoalInp(mod_methane, 1.0, [CH4])
        methane.composition = [1.0] # here we mock a composition fit
        # The HHV change shouldn't affect the enthalpy of formation
        h_fCH4 = -74831 # kJ/kmol
        assert almost(methane.EnthalpyOfFormation, h_fCH4/16.0)
        assert almost(methane.EnthalpyOfFormationMol, h_fCH4)

        assert almost(methane.hhv_daf, hhv_orig-delta)
        if delta == 0:
            h_devol_orig  = methane.EnthalpyOfDevol
        print methane.EnthalpyOfDevol
        assert h_devol_orig - methane.EnthalpyOfDevol - delta < 0.01

def test_YieldCO2H2O():
    """ Test standard methods of the postulate substance class
        for a known species """
    from pkp.src.CoalThermoPhysics import (
        enthalpy_balance, enthalpy_balance_mol, YieldFromCoalInp, H2O, CO2)
    coal = YieldFromCoalInp(inert_coal['Coal'], 1.0, [CO2, H2O])

    coal.composition = [0.5, 0.5] # here we mock a composition fit
    print "coal mechanism partial: ", coal.mechanism(False)
    print coal.hhv_daf
    print coal.lhv_daf
    print coal.lhv_yield
    print coal.lhv_yield_species
    # h0f_CH4 = co.EnthalpyOfFormation()
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
