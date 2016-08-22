from pkp.Compos_and_Energy import VolatileComposition


composition = {
    'Tar': 0.3069,
    'H2O': 0.0732,
    'CO2': 0.0239,
    'CH4': 0.0386,
    'CO':  0.0471,
    'other': 0.0450
}

composition['Solid'] = 1 - sum(composition.values())
print 'sum of comp', sum(composition.values())

ultimate_analysis = {'C': 75.23, 'H': 5.16,
                     'O': 9.83, 'N': 1.43}


def create_composition():
    return VolatileComposition(composition=composition,
                               ultimate_analysis=ultimate_analysis)


def test_calc_element_fraction():
    vol = create_composition()
    assert sum(vol.composition.values()) == 1
    assert sum(vol.ultimate_analysis.values()) == 1
    assert 'C6H6' in vol.composition

    gas = vol.gas
    i_C6H6 = gas.species_index('C6H6')
    assert vol.calc_element_fraction('C', 'C6H6') == \
        6 * gas.atomic_weight('C') / gas.molecular_weights[i_C6H6]
    assert vol.calc_element_fraction('C', 'Char') == 1.0
    assert vol.calc_element_fraction('H', 'Char') == 0.0


def test_calc_tot_element_fraction():
    vol = create_composition()
    vol.calc_tot_element_fraction('C')

vol = create_composition()
vol.define_composition_empirical(CO=0)

print 'CH4=', vol.heat_of_reaction_species('CH4') / 1e6
print 'Solid=', vol.heat_of_reaction_species('Solid') / 1e6
print 'N2=', vol.heat_of_reaction_species('N2') / 1e6
print 'Vol=', vol.heat_of_volatiles() / 1e6
print 'Heat pyro=', vol.heat_of_pyrolysis(heat_coal=28e6) / 1e6
