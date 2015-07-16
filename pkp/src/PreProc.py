from CoalThermoPhysics import MolWeights

class PreProcResult(object):
    """ Base class for the various preprocessors """

    def __init__(self, coal):
        self.coal = coal
        self.qFactor = 1.0

    @property
    def VolatileCompositionMass(self):
        """ m_species/m_yield
            where:
                m_yield = m_h + m_o + m_vc_cur

            The difficulty is to know the carbon content of the
            volatile yield.
            Units: [kg/kg_yield]
        """
        ua = self.coal.ua_vm # uncorected ua_vm
        f = (1.0 - self.qFactor)*100.0
        return ua.remove_elem_mass_rebalance('Carbon', f)

    @property
    def pa_raw(self):
        """ composition of intermediate coal """
        ftot = self.ftot*100.0
        return {'Fixed Carbon': 100.0 - ftot,
                'Volatile Matter': ftot}

    @property
    def ftot(self):
        return self.qFactor*self.coal.pa_daf['Volatile Matter']/100.0

    @property
    def VolatileCompositionMol(self):
        """ the species composition of the fuel in mol per mol yield
            Units [mol/mol_yield]
        """
        molar_mass_vm = self.coal.MW_PS
        comp_mass = self.VolatileCompositionMass
        # species_mass_fraction*mw = moles of species per kg fuel
        # and we normalise that by multiplying with molar_mass_vm
        # to get moles of species per mol fuel
        # TODO: Does this make sense? This probably doesnt sum up to one!
        #       The function has been tested but what is the effect of
        #       arbitray MW_PS. Do we need to rescale?
        return {elem: comp_mass[elem]/mw*molar_mass_vm
                for elem, mw in MolWeights.iteritems()
                if elem in comp_mass}

class ManualQfactor(PreProcResult):
    """ Intermediate class to mock a preprocessor result """

    def __init__(self, coal, qFactor):
        PreProcResult.__init__(self, coal)
        self.qFactor = qFactor
    
