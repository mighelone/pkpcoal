""" Database for all thermo physical data """
import matplotlib.pyplot as plt

from Models import BalancedComposition

R = 8.314 # [kJ/(kmol*K)]

CoresProd = {
    'Carbon':'CO',
    'Hydrogen':'H2O',
}


MolWeights = { #g/mol
    'Oxygen':16.0,
    'Carbon':12.011,
    'Hydrogen':1.008,
    'Nitrogen':14.007,
    'CO': 28.010,
    'CO2': 44.011,
    'H2O': 18.016,
    'OH': 17.0,
    'H2': 2.016,
    'N2': 28.013,
    'O2': 31.999,
    'CH4': 16.043,
    'C2H2': 26.0,
    'C2': 24.0,
    'Char': 12.0,
    'CH3OH': 32.04,
    'C2H5OH': 46.07,
    'C7H6O2': 122.12,
    'C6H12O6': 180.16,
    'C6H10O5': 162.0,
    'CH3OOH': 48.4, # http://pubchem.ncbi.nlm.nih.gov/compound/Methyl_hydroperoxide#section=Chemical-and-Physical-Properties
    'CH2O': 30.0,
}

EnthOfForm = { # in [kJ/kmol], gaseous
    'CO':  -110541.0, # Turns p. 622
    'CO2': -393546.0, # Turns p. 623
    'H2O': -241845.0, # Turns p. 632
    'OH': 38987.0,
    'CH4': -74831.0,
    'C2H2': 226880.0, #http://en.wikipedia.org/wiki/Acetylene
    'Char': 0,
    'CH3OH': -201000, #https://de.wikipedia.org/wiki/Methanol
    'C2H5OH': -235000, #https://en.wikipedia.org/wiki/Standard_enthalpy_change_of_formation_%28data_table%29
    'C7H6O2': -385000,
    'C6H12O6': -1271000,
    'C6H10O5': -963000,
    'CH3OOH': 0.0, # Missing
    'CH2O': -108600.0,  # http://webbook.nist.gov/cgi/cbook.cgi?ID=C50000&Mask=1
}

EnthOfVap = { #kJ/kmol
    "H2O": 44010
    }

EnthOfFormKG = {name: value/MolWeights[name]
        for name, value in EnthOfForm.iteritems()}


class Species(object):

    def __init__(self, name, molecular_weight,
            enthalpy_of_formation, components):
        self.name = name
        self.molecular_weight = molecular_weight
        self.components = components
        self.enthalpy_of_formation = enthalpy_of_formation

    def fraction(self, elem):
        """ the fraction of mass MW_el/MW_tot """
        if elem in self.components:
            return self.components[elem]/self.molecular_weight
        else:
            return 0.0

    @property
    def hf_perKg(self):
        return self.enthalpy_of_formation/self.molecular_weight

    @property
    def hf_perMol(self):
        return self.enthalpy_of_formation

    @property
    def fractions(self):
        return {elem:self.fraction(elem) for elem in self.components}

#TODO replace instantiation by parser
CO   = Species('CO', MolWeights['CO'], EnthOfForm['CO'], {'Carbon': MolWeights['Carbon'], 'Oxygen': MolWeights['Oxygen']})
OH   = Species('OH', MolWeights['OH'], EnthOfForm['OH'], {'Hydrogen': MolWeights['Hydrogen'] ,'Oxygen': MolWeights['Oxygen']})
CO2  = Species('CO2', MolWeights['CO2'], EnthOfForm['CO2'], {'Carbon': MolWeights['Carbon'], 'Oxygen': 2*MolWeights['Oxygen']})
CH4  = Species('CH4', MolWeights['CH4'], EnthOfForm['CH4'], {'Hydrogen': 4*MolWeights['Hydrogen'] ,'Carbon': MolWeights['Carbon']})
C2H2 = Species('C2H2', MolWeights['C2H2'], EnthOfForm['C2H2'], {'Hydrogen': 2*MolWeights['Hydrogen'] ,'Carbon': 2*MolWeights['Carbon']})
N2   = Species('N2', MolWeights['N2'], 0.0, {'Nitrogen': 2*MolWeights['Nitrogen']})
H2   = Species('H2', MolWeights['H2'], 0.0, {'Hydrogen': 2*MolWeights['Hydrogen']})
O2   = Species('O2', MolWeights['O2'], 0.0, {'Oxygen': 2*MolWeights['Oxygen']})
H2O  = Species('H2O', MolWeights['H2O'], EnthOfForm['H2O'], {'Hydrogen': 2*MolWeights['Hydrogen'],'Oxygen': MolWeights['Oxygen']})
C2   = Species('C2', MolWeights['C2'], 0.0, {'Carbon': 2*MolWeights['Carbon']}) #FOR DEBUG ONLY

# Methanol
CH3OH = Species('CH3OH', MolWeights['CH3OH'], EnthOfForm['CH3OH'],
    {'Hydrogen': 4*MolWeights['Hydrogen'], 'Oxygen': MolWeights['Oxygen'], 'Carbon': MolWeights['Carbon']})

# Ethanol
C2H5OH = Species('C2H5OH', MolWeights['C2H5OH'], EnthOfForm['C2H5OH'],
    {'Hydrogen': 6*MolWeights['Hydrogen'], 'Oxygen': MolWeights['Oxygen'], 'Carbon': 2*MolWeights['Carbon']})

# Benzoic Acid
C7H6O2 = Species('C7H6O2', MolWeights['C7H6O2'], EnthOfForm['C7H6O2'],
    {'Hydrogen': 6*MolWeights['Hydrogen'], 'Oxygen': 2*MolWeights['Oxygen'], 'Carbon': 7*MolWeights['Carbon']})

# Glucose
C6H12O6 = Species('C6H12O6', MolWeights['C6H12O6'], EnthOfForm['C6H12O6'],
    {'Hydrogen': 12*MolWeights['Hydrogen'], 'Oxygen': 6*MolWeights['Oxygen'], 'Carbon': 6*MolWeights['Carbon']})

# Cellulose
C6H10O5 = Species('C6H10O5', MolWeights['C6H10O5'], EnthOfForm['C6H10O5'],
    {'Hydrogen': 10*MolWeights['Hydrogen'], 'Oxygen': 5*MolWeights['Oxygen'], 'Carbon': 6*MolWeights['Carbon']})

#  Methyl hydroperoxide # Enthalpy of Formation is missing
CH3OOH = Species('CH3OOH', MolWeights['CH3OOH'], EnthOfForm['CH3OOH'],
    {'Hydrogen': 4*MolWeights['Hydrogen'], 'Oxygen': 2*MolWeights['Oxygen'], 'Carbon': 1*MolWeights['Carbon']})

# Formaldehyde
CH2O  = Species('CH2O', MolWeights['CH2O'], EnthOfForm['CH2O'],
    {'Hydrogen': 2*MolWeights['Hydrogen'], 'Oxygen': MolWeights['Oxygen'], 'Carbon': MolWeights['Carbon']})


class Coal(object):
    """ Class to hold all coal properties """

    def __init__(self, input_dict):
        self.d = input_dict
        self.pa = BalancedComposition(input_dict["Proximate Analysis"])
        self.pa_dry = self.pa.remove_elems_rebalance(['Moisture'])
        self.pa_daf = self.pa.remove_elems_rebalance(['Moisture', 'Ash'])
        # NOTE by default we neglect any sulphur content
        self.ua_wS = BalancedComposition(input_dict["Ultimate Analysis"])
        self.ua = self.ua_wS.remove_elems_rebalance(['Sulphur'])
        self.ua_vm = self.ua.remove_elem_mass_rebalance('Carbon',
            self.pa_daf['Fixed Carbon'])
        self.hhv = input_dict["hhv"] #in kJ/kg
        self.MW_PS = input_dict.get("MW_PS", 100)

    @property
    def hhv_daf(self):
        """ HHV_{daf} = HHV_{as-recieved}/(f_{FC}+f_{VM}) """
        #NOTE the factor 100.0 comes from the fact that pa is in percents
        return self.hhv*100.0/(self.pa['Fixed Carbon'] + self.pa['Volatile Matter'])


    @property
    def hhv_dulong(self):
        """ HHV_{as_recieved} approximated by the Dulong formula in [MJ/kg_{ar} """
        ua = self.ua
        return (  0.3279*ua['Carbon']
                + 1.504*(ua['Hydrogen'] - ua['Oxygen']/8.0)
                + 0.0926*ua['Sulphur']
                + 0.0497*ua['Oxygen']
                + 0.0242*ua['Nitrogen'])

    @property
    def lhv_daf(self):
        """ LHV_{daf} = HHV_{daf} - h_{latent,H2O}(MW_O+2MW_H)/(2MW_H)f_{UA,H} """
        return self.hhv_daf - self.ua['Hydrogen']/100.0*EnthOfVap['H2O']/MolWeights['H2O']/(H2O.fractions['Hydrogen'])

    def plot_composition(self):
        # TODO refactor
        f, ax = plt.subplots(1,3)
        f.set_figwidth(18)
        f.set_figheight(5.5)
        ax[0].pie(self.pa.values(),
            labels=self.pa.elems.keys(),
            autopct='%.0f%%',
            explode=[0.05 for _ in self.pa.values()],
            colors=['w' for _ in self.pa.values()])
        ax[0].set_title("Prox. Ana.")

        ax[1].pie(self.pa_daf.values(),
            labels=self.pa_daf.elems.keys(),
            autopct='%.0f%%',
            explode=[0.05 for _ in self.pa_daf.values()],
            colors=['w' for _ in self.pa_daf.values()])
        ax[1].set_title("Prox. Ana. DAF")

        ax[2].pie(self.ua.values(),
            labels=self.ua.elems.keys(),
            autopct='%.0f%%',
            explode=[0.05 for _ in self.ua.values()],
            colors=['w' for _ in self.ua.values()])
        ax[2].set_title("Ultimate. Ana. DAF")
        plt.suptitle("Elemental Compostion Coal")
        return f, ax


def composition(species, scale=1.0):
    """ returns the elemental composition for mixture of species

        Arguments: species a dictionary with {'name':massfraction}
"""
    elems = ['Carbon', 'Oxygen', 'Hydrogen', 'Nitrogen']
    comp = {elem: 0.0 for elem in elems}
    for s, fraction in species.iteritems():
        for elem, value in globals()[s].fractions.iteritems():
            comp[elem] += fraction * value * scale
    return comp

class Yield(object):
    """ Base class for fitting the composition of the pyrolysis yield
        to the preprocessor data """

    def __init__(self,
        preProcessorResults,
        targetSpecies=None,
    ):
        self.pre = preProcessorResults
        self.targetSpecies = targetSpecies
        self.target_map = {pos: name for pos, name in enumerate(targetSpecies)}
        self.composition =  [False for _ in targetSpecies]

    # ALIASES
    @property
    def VolatileCompositionMass(self):
        # TODO inheritance instead of delegation ?
        return self.pre.VolatileCompositionMass

    @property
    def VolatileCompositionMol(self):
        return self.pre.VolatileCompositionMol

    @property
    def hhv_daf(self):
        return self.pre.coal.hhv_daf

    @property
    def lhv_daf(self):
        return self.pre.coal.lhv_daf

    @property
    def q_release_char_yield(self):
        """ heat released by char conversion
            kJ/kg_yield """
        return (392.0/12.0*self.pre.pa_raw['Fixed Carbon']/100.0*1000.0)

    @property
    def lhv_yield(self):
        """ LHV_{yield} = (LHV_{daf} - f_{char} LHV_{char})/f_{yield} """
        f_yield = self.pre.ftot
        f_char =  1.0 - f_yield
        return (self.lhv_daf - self.q_release_char_yield)/(self.pre.pa_raw['Volatile Matter']/100.0)

    @property
    def hhv_yield(self):
        """ LHV_{yield} = (LHV_{daf} - f_{char} LHV_{char})/f_{yield} """
        f_yield = self.pre.ftot
        f_char =  1.0 - f_yield
        return (self.hhv_daf - self.q_release_char_yield)/(self.pre.pa_raw['Volatile Matter']/100.0)

    @property
    def lhv_yield_species(self):
        """ The heating value computed based on
            species composition and individial species
            enthalpies """
        return -(self.EnthalpyBalanceProds - self.EnthalpyOfFormation)

    def error_func(self, species_massfractions):
        """ compute the cumulated percentual error of composition

            Parameters:
                target: list of target mass fractions of elements
                coeffs: matrix of nu*MW_elm/MW_species, i = Elm, j = spec
        """
        def elem_error_perc(target_mass, mass_elem):
            return abs((target_mass - mass_elem)/target_mass)
        def mass_elm_species(elem):
            return sum([species_massfractions[i] * spec.fraction(elem)
                        for i, spec in enumerate(self.targetSpecies)])
        return sum([elem_error_perc(target_mass/100.0, mass_elm_species(elem))
                        for elem, target_mass in self.VolatileCompositionMass])

    def fit_composition(self, optimizer):
        from scipy.optimize import brute
        from collections import OrderedDict
        res = optimizer()
        od = OrderedDict([(self.target_map[i].name, x) for i, x in enumerate(res)])
        self.composition = od.values()
        return od

    def ProductCompositionMass(self, partialOx=True):
        """ Same as ProductCompositionMol but on mass basis
            Units [kg/kg_Prod]
        """
        # NOTE This could be based on ProductCompositionMol, but as long
        # it remains unclear that self.VolatileCompositionMol is correct
        # the product composition is computed separately
        comp = self.VolatileCompositionMass
        spec = ('CO' if partialOx else 'CO2')
        # each kg carbon in the yield can produce MW_CO?/MW_C kg CO2
        COx = comp['Carbon']*MolWeights[spec]/MolWeights['Carbon']
        # each kg hydrogen in the yield can produce MW_H2O/MW_H kg H2O
        H2O = comp['Hydrogen']*MolWeights['H2O']/(2.0*MolWeights['Hydrogen'])
        # The product mass of nitrogen stays unchanged
        N2 = comp['Nitrogen']
        # Creating a BalancedComposition changes the units from
        # [kg/kg_Yield] to [kg/kg_Prod]
        # but the old basis is still availible
        return BalancedComposition({spec: COx, 'H2O': H2O, 'N2': N2})

    def ProductsMass(self, partialOx=True):
        """ [kg/kg_yield] """
        return {name: comp/self.Fuel_to_Product_Mass/100.0
                for name, comp in self.ProductCompositionMass(partialOx).iteritems()
        }


    def ProductCompositionMol(self, partialOx=True):
        """ compute the molar composition of the products per mol fuel

            assumes that each mol of carbon produces a mol CO or CO2
                         each mol of hydrogen produces half a mol H20
            Units: [mol/mol_Prod]
        """
        # the molar composition is based on the product mass composition
        comp = self.ProductCompositionMass(partialOx)
        # by using a BalancedComposition the individual components are
        # divided by the total sum product moles
        return BalancedComposition(
                {spec: val/MolWeights[spec] for spec, val in comp.iteritems()})


    def ProductMols(self, partialOx=True):
        """ compute the the number of mols of
            product species per mol Yield """
        return {name: comp/MolWeights[name]*self.MW
                for name, comp in self.ProductsMass(partialOx).iteritems()
        }


    @property
    def Oxygen_from_yield(self):
        """ the amount of oxygen in products coming from the yield
        """
        comp = composition(self.ProductCompositionMass(partialOx=False))
        f = self.Fuel_to_Product_Mass
        return self.VolatileCompositionMass['Oxygen']/(comp['Oxygen']/f)

    @property
    def Fuel_to_Product_Mass(self):
        # kg_yield/kg_products
        return 100.0/self.ProductCompositionMass(partialOx=False).basis

    @property
    def FOx(self):
        """ the Fuel to Oxygen ratio on mass basis
            m_Fuel + m_O2 -> m_Prod
            FOx = m_Fuel/m_O2
            m_O2 = m_O_Products * (1-r)
            FOx = (m_Prod - m_O2)/m_O2
            TODO DOUBLE CHECK THIS
        """
        prod = composition(self.ProductCompositionMass(partialOx=False))
        f = self.Fuel_to_Product_Mass
        m_O2 = prod['Oxygen']/(f*100.0)*(1.0 - self.Oxygen_from_yield)
        return 1.0/m_O2

    @property
    def MW(self):
        """ Return the molecular weight of the yield
            ATTENTION the right composition needs to be
            set before hand """
        if False in self.composition:
            print "WARNING composition hasn't been computed yet"
            return 0.0
        return MolecularWeight(
                    self.targetSpecies,
                    self.composition
               )

    @property
    def FOx_molar(self):
        """ the Fuel to Oxygen ratio on molar basis """
        return (self.MW/MolWeights["O2"])/self.FOx

    def mechanism(self, partialOx=True):
        """ Return a string of the reaction mechanism """
        vol = self.VolatileCompositionMol
        prods = self.ProductMols(partialOx=partialOx)
        prodName = ('CO' if partialOx else 'CO2')
        nO2 = self.FOx_molar
        yield_ = "+".join([_.name for _ in self.targetSpecies])
        reaction = "({}) + {}O2 -> {}{} + {}H2O + {}N2".format(
           yield_, nO2, prods[prodName], prodName, prods['H2O'], prods['N2'])
        return reaction

    @property
    def EnthalpyBalanceProds(self):
        """ kJ/kg_yield """
        return enthalpy_balance(self.ProductsMass(partialOx=False))

    @property
    def EnthalpyBalanceProdsMol(self):
        """ kJ/kmol_yield """
        return enthalpy_balance_mol(self.ProductMols(partialOx=False))

    @property
    def EnthalpyOfFormation(self):
        """ Enthalpy of formation of the yield from
            products and given heating value [kJ/kg] """
        # if we have species and a composion
        # we can compute the the enthalpy of formation
        # from the individual species
        if (all([isinstance(s, Species)
                for s in self.targetSpecies])
            and all([isinstance(s, float)
                for s in self.composition])):
            return sum([s.hf_perKg*w for s,w in zip(self.targetSpecies, self.composition)])

        print "Warning computing the yield enthalpy of formation from energy balance instead of individual species"
        # other wise we fallback to compute the
        # enthalpy of formation from the energy balance
        LHV = self.pre.coal.hhv # kJ/kg_yield
        return LHV + self.EnthalpyBalanceProds

    @property
    def EnthalpyOfFormationMol(self):
        """ Enthalpy of formation of the yield from
            products and given heating value [kJ/kmol] """
        return self.EnthalpyOfFormation*self.MW

    @property
    def EnthalpyOfDevol(self):
        """ kJ/kg """
        return -self.lhv_yield + self.lhv_yield_species

def enthalpy_balance(composition):
    """ returns the sum of the enthalpies of formation weighted by the composition
        in kJ/kg """
    # NOTE all species that are not in EnthOfForm are assumed to have 0
    #      enthalpy of formation DANGEROUS
    return sum([val*EnthOfFormKG[spec]
                for spec, val in composition.iteritems()
                if spec in EnthOfForm])

def enthalpy_balance_mol(composition):
    """ returns the sum of the enthalpies of formation weighted by the composition
        in kJ/kmol
    """
    # NOTE all species that are not in EnthOfForm are assumed to have 0
    #      enthalpy of formation DANGEROUS
    return sum([val*EnthOfForm[spec]
                for spec, val in composition.iteritems()
                if spec in EnthOfForm])

def generatePostulateSubstanceSpecies(preProc):
    """ generate a single postulate ^substance species """
    mw = preProc.coal.MW_PS
    hf = Yield(preProc, ['CxHyOz']).EnthalpyOfFormationMol
    return Species("CxHyOz", mw, hf, preProc.VolatileCompositionMass)

def MolecularWeight(species, composition):
    """ compute the molecular weight of a given composition """
    return sum([f*s.molecular_weight for f, s in zip(composition, species)])

def PreProcFromCoalInp(coalDict, qFactor=1.0, targetSpecies=None):
    """ a factory method to create a preProc object
        directly from a coal input dict and a qFactor
        ommiting any pre-processor """
    from pkp.src.PreProc import ManualQfactor
    return ManualQfactor(Coal(coalDict), qFactor=qFactor)

def YieldFromCoalInp(coalDict, qFactor, targetSpecies):
    """ a factory method to create a yield object
        directly from a coal input dict and a qFactor
        ommiting any pre-processor """
    from pkp.src.PreProc import ManualQfactor
    preProc = ManualQfactor(Coal(coalDict), qFactor=qFactor)
    return Yield(preProc, targetSpecies)

