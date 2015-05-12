""" Database for all thermo physical data """

from Models import BalancedComposition

R = 8.314 # [kJ/(kmol*K)]

MolWeights = { #g/mol
    'Oxygen':15.999,
    'Carbon':12.011,
    'Hydrogen':1.008,
    'Nitrogen':14.007,
    'CO': 28.010,
    'CO2': 44.011,
    'H2O': 18.016,
}

EnthOfForm = { # in [kJ/kmol]
    'CO':  -110541.0, # Turns p. 622
    'CO2': -393546.0, # Turns p. 623
    'H2O': -241845.0, # Turns p. 632
}

EnthOfFormKG = {name: value/MolWeights[name]
        for name, value in EnthOfForm.iteritems()}

CoresProd = {
    'Carbon':'CO',
    'Hydrogen':'H2O',
}

class Coal(object):
    """ Class to hold all coal properties """

    def __init__(self, input_dict):
        self.pa = BalancedComposition(input_dict["Proximate Analysis"])
        self.pa_dry = self.pa.remove_elems_rebalance(['Moisture'])
        self.pa_daf = self.pa.remove_elems_rebalance(['Moisture', 'Ash'])
        # NOTE by default we neglect any sulphur content
        self.ua_wS = BalancedComposition(input_dict["Ultimate Analysis"])
        self.ua = self.ua_wS.remove_elems_rebalance(['Sulphur'])
        self.ua_vm = self.ua.remove_elem_mass_rebalance('Carbon',
            self.pa_daf['Fixed Carbon'])
        self.hhv = input_dict["hhv"]
        self.MW_PS = input_dict.get("MW_PS",100)


class PostulateSubstance(object):

    def __init__(self, coal, qFactor=1.0):
        self.coal = coal
        self.q = qFactor

    def VolatileCompositionMass(self):
        """ m_species/m_tot
            where:
                m_tot = m_h + m_o + m_vc_cur

            The difficulty is to know the carbon content of the
            volatile yield.
                m_c_ua = m_fc_prox + m_vc_prox
                m_c_ua = m_fc_cur + m_vc_cur

                m_vc_cur = m_c_ua - m_fc_prox/q_factor
        """
        #TODO use BalancedComposition methods
        ua = self.coal.ua
        pa = self.coal.pa_daf
        carbon = (ua['Carbon'] - pa['Fixed Carbon']
                 + (self.q-1) * pa['Volatile Matter'])
        oxygen = ua['Oxygen']
        hydrogen = ua['Hydrogen']
        nitrogen = ua['Nitrogen']
        tot = carbon + oxygen + hydrogen + nitrogen
        assert tot < 100.1
        return {'Carbon': carbon/tot,
                'Hydrogen': hydrogen/tot,
                'Oxygen': oxygen/tot,
                'Nitrogen': nitrogen/tot
                }

    def VolatileCompositionMol(self):
        """ the species composition of the fuel in mol per mol fuel """
        molar_mass_vm = self.coal.MW_PS
        comp_mass = self.VolatileCompositionMass()
        # species_mass_fraction*mw = moles of species per kg fuel
        # and we normalise that by multiplying with molar_mass_vm
        # to get moles of species per mol fuel
        return {elem: comp_mass[elem]/mw*molar_mass_vm
                for elem, mw in MolWeights.iteritems()
                if elem in comp_mass.keys()}

    def ProductCompositionMass(self, partialOx=True):
        pass

    def ProductCompositionMol(self, partialOx=True):
        """ compute the molar composition of the products per mol fuel

            assumes that each mol of carbon produces a mol CO or CO2
                         each mol of hydrogen produces half a mol H20
        """
        # first we get the molar composition of the volatile matter,
        # with that we can compute the product compostion per mol vm
        comp = self.VolatileCompositionMol()
        CO   =     comp['Carbon']
        H2O  = 0.5*comp['Hydrogen']
        N2   =     comp['Nitrogen']
        # if we assume partial oxidaten CO is produced instead of
        # CO2, but only the product name changes, stoichiometric
        # constants are the same. This is important since enthalpy
        # of formation is computed assuming full oxidation
        prodName = ('CO' if partialOx else 'CO2')
        return {prodName: CO, 'H2O': H2O, 'N2': N2}

    def EnthalpyOfFormation(self):
        """ Computes the enthalpy of formation of the volatile matter

            h_react = LHV + sum(n_i h_prod_i) with n beeing stoich factor
        """
        vol_comp = self.ProductCompositionMol(partialOx=False)
        H_products = 0.0
        # for every element in the volatile composition we get the
        # corresponding product and its enthapy of formation kJ/kmol
        LHV = self.coal.hhv
        molar_mass_vm = self.coal.MW_PS
        for name, mol in vol_comp.iteritems():
            h_prod = EnthOfForm.get(name, 0.0)
            H_products += mol*h_prod
        h_0f = (LHV*molar_mass_vm+H_products) # [kJ/kmol]
        return h_0f

    def mechanism(self, partialOx=True):
        vol = self.VolatileCompositionMol()
        prods = self.ProductCompositionMol(partialOx=partialOx)
        prodName = ('CO' if partialOx else 'CO2')
        nO2prod = (prods[prodName] if partialOx else 2*prods[prodName])
        nO2 = (nO2prod + prods['H2O'] - vol['Oxygen'] * 0.5)*0.5 #TODO correct?
        reaction = "CxHyOz + {}O2 -> {}{} + {}H2O + {}N2".format(
            nO2, prods[prodName], prodName, prods['H2O'], prods['N2'])
        return reaction
