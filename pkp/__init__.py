'''
The module contains the class for reading and run pyrolysis
calculations with PKP.

Contains
--------
* :class:`pkp.PKPRunner`
* :class:`pkp.ReadConfiguration`
'''
from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals
from builtins import dict
from six import string_types

from autologging import logged
try:
    import ruamel.yaml as yaml
except:
    try:
        import ruamel_yaml as yaml
    except:
        print('Loading standard yaml module ...\n'
              'Note that it might give problems with'
              ' exponetial decoding')
        import yaml
import json
import os
import numpy as np
import pandas as pd
import cantera

from pkp.cpd import CPD
from pkp.polimi import Polimi
from pkp.biopolimi import BioPolimi
from pkp.detailed_model import M_elements, hf
import pkp.evolution

import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
try:
    plt.style.use('mystyle')
except:
    plt.style.use('ggplot')

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
col_red = "#C54E6D"
col_green = "#009380"

models = ['CPD', 'Polimi', 'BioPolimi']


def clean_dict(d, tolist=False):
    '''
    Clean dictionary fron numpy and pandas object.

    Parameters
    ----------
    d: dict
        Dictionary to convert
    tolist: bool, default=False
        Convert numpy arrays and pandas dataframe to list

    Returns
    -------
    d_new: dict
        New cleaned dictionary
    '''
    d_new = {}
    for k, v in d.items():
        if isinstance(v, dict):
            d_new[k] = clean_dict(v)
        elif isinstance(v, np.ndarray):
            if tolist:
                d_new[k] = v.tolist()
        elif isinstance(v, pd.DataFrame):
            if tolist:
                d_new[k] = v.values.tolist()
        else:
            d_new[k] = v
    return d_new


@logged
class ReadConfiguration(pkp.detailed_model.DetailedModel):
    '''
    Read configuration file for PKP
    '''

    def __init__(self, yml):
        '''
        Parameters
        ----------
        yml: str, unicode, dict
            Input dict or yaml file containing the configuration for
            run PKP. See :ref:`input-file-label`.
        '''
        if isinstance(yml, string_types):
            with open(yml, 'r') as f:
                yml_input = yaml.load(f)
        elif isinstance(yml, dict):
            yml_input = yml
        else:
            raise ValueError('Define yml as file name or dictionary')

        # coal settings
        coal_settings = yml_input['Coal']
        # Solver settings

        self.__log.debug(
            'Pressure in yml %s',
            yml_input['operating_conditions']['pressure'])
        super(ReadConfiguration, self).__init__(
            proximate_analysis=coal_settings['proximate_analysis'],
            ultimate_analysis=coal_settings['ultimate_analysis'],
            pressure=yml_input['operating_conditions'][
                'pressure'] * 101325,
            name=coal_settings['name'])
        self.__log.debug('Pressure setted %s', self.pressure)

        self.operating_conditions = yml_input['operating_conditions']

        # convert HHV from MJ/kg to J/kg
        self.HHV = coal_settings['HHV'] * 1e6
        self.rho_dry = coal_settings['rho_dry']

        # Solver settings
        [setattr(self, model, yml_input[model])
         for model in models]

    @property
    def operating_conditions(self):
        return self._operating_conditions

    @operating_conditions.setter
    def operating_conditions(self, value):
        self._operating_conditions = value


@logged
class PKPRunner(ReadConfiguration):
    '''
    PKP Runner manager class. It uses configuration in the *yaml* file
    to run multiple simulations of coal pyrolysis using different
    *detailed models* and fitting their results with
    :ref:`empmodels-label`.
    '''
    models = models

    def run(self, results_dir=None, n_p=1):
        '''
        Run detailed models and fit them.

        Parameters
        ----------
        results_dir: str, default=None
            Directory where results are stored. If None is specified is
            used the directory from where PKP is launched.
        np: int, default=1
            Number of processors for evolution fitting

        Returns
        -------
        run_results: dict
            Dictionary with the results of the detailed models.
        fit_results: dict
            Dictionary with the results of the calibration of the
            empirical models.
        '''
        results_dir = self.set_results_dir(results_dir)
        run_results = {}

        # define a information structure for the coal properties
        # TODO improve this part

        run_results['coal'] = {
            'name': self.name,
            'ultimate_analysis': self.ultimate_analysis,
            'proximate_analysis': self.proximate_analysis,
            'proximate_analysis_daf': self.proximate_analysis_daf,
            'HHV': self.HHV,
            'rho_dry': self.rho_dry,
            'operating_conditions': self.operating_conditions
        }
        fit_results = {}
        for model in self.models:
            model_settings = getattr(self, model)
            if model_settings['active']:
                self.__log.info('Run model %s', model)
                results = self._run_model(model=model,
                                          results_dir=results_dir)
                if results:
                    run_results[model] = results
                    self.__log.debug('Finish run %s %s', model,
                                     results.keys())
                else:
                    self.__log.warning('No results for %s', model)
                if model_settings['fit']:
                    self.__log.info('Start fit of %s model', model)
                    fit_results[model] = self._fit_model(
                        model, model_settings['fit'], n_p, results,
                        results_dir)

        json_fit = os.path.join(
            results_dir, '{name}-fitreport.json'.format(name=self.name))
        self.__log.debug('Export fit report to %s', json_fit)
        with open(json_fit, 'w') as f:
            json.dump(clean_dict(fit_results), f, indent=4)
        return run_results, fit_results

    def _fit_model(self, model, model_settings, n_p, results,
                   results_dir):
        '''
        Run fitting of the given model

        Parameters
        ----------
        model: str, unicode
            Name of the detailed model to fit
        model_settings: dict
            Dictionary settings of the fitting
        n_p: int
            Number of processor for multiprocessing evolution
        results: dict
            Dictionary results of the detailed model
        results_dir: str, unicode
            Name of the output directory for storing results

        Returns
        -------
        fit_results: dict
            Contains results of fitting

        '''
        fit_results = {}
        # loop over the fitting runs, fit0, fit1, etc.
        for fitname, fit in model_settings.items():
            if fit.get('active', True):
                self.__log.info(
                    'Fit %s model with %s', model, fit['model'])
                target_conditions = {
                    run: {'t': np.array(res.index),
                          'y': np.array(res[fit['species']])}
                    for run, res in results.items()}
                fit_dict = {'model': model,
                            'fit': fitname,
                            'species': fit['species']}
                fit_results[fitname] = self._evolution(
                    target_conditions, fit_dict,
                    fit, results_dir, n_p)
                fit_results[fitname]['species'] = fit['species']
                fit_results[fitname]['model'] = fit['model']
        return fit_results

    @staticmethod
    def set_results_dir(results_dir):
        if results_dir is None:
            results_dir = os.getcwd()
        return results_dir

    def _run_model(self, model, results_dir):
        '''
        Run simulations for the given model

        Parameters
        ----------
        model: str
            Name of the detailed model. Note that it should be the same
            of a defined class
        results_dir: str
            Path of results

        Returns
        -------
        results: dict
        '''
        model_settings = getattr(self, model)
        self.__log.debug('Model %s active %s', model,
                         model_settings['active'])
        if model_settings['active']:
            results = {}
            self.__log.debug('Run %s',
                             self.operating_conditions['runs'])
            vol_composition = pd.DataFrame()
            for n in range(
                    self.operating_conditions['runs']):
                res = self._run_single(model, model_settings, n,
                                       results_dir)
                results['run{}'.format(n)] = res

                # add last row to vol_composition
                vol_composition = vol_composition.append(
                    res.tail(1), ignore_index=True)
                self.__log.debug('Finish run %s', results.keys())

                # plot results
                self._plot_results(model, n, res, results_dir)

            # add index to vol_composition dataframe
            vol_composition.index = [
                'run{}'.format(n)
                for n in range(self.operating_conditions['runs'])]
            final_yield = '{name}-{model}-finalyields.csv'.format(
                name=self.name, model=model)
            self.__log.debug('Export vol_composition to csv %s',
                             final_yield)
            vol_composition.to_csv(
                os.path.join(results_dir, final_yield),
                index=True)
        else:
            results = None
        return results

    def _run_single(self, model, model_settings, n, results_dir):
        '''
        Run a single simulation for the given detailed model.

        Parameters
        ----------
        model: str, unicode
            Model name
        model_settings: dict
            Settings of the given model
        n: int
            Run number
        results_dir: str, unicode
            Store results directory

        Returns
        -------
        res: pd.DataFrame
            Results datafram
        '''
        self.__log.debug(
            'Initialize run %s for %s', n, model)
        if model == 'Polimi' and 'reference' in model_settings:
            self.__log.debug('Use reference coal for Polimi %s',
                             model_settings['reference'])
            run = globals()[model].reference_coal(
                ref_coal=model_settings['reference'],
                proximate_analysis=self.proximate_analysis,
                pressure=self.pressure
            )
            self.__log.debug(
                'Polimi coal composition is set to %s', run.composition)
        else:
            self.__log.debug('Initialize detailed model %s',
                             model)
            run = globals()[model](
                ultimate_analysis=self.ultimate_analysis,
                proximate_analysis=self.proximate_analysis,
                pressure=self.pressure,
                #name='{}-{}-Run{}'.format(self.name, model, n)
                name=self.name
            )
        run.basename = '{name}-{model}-run{run}'.format(
            name=self.name, model=model, run=n)
        self.__log.debug('Set basename %s', run.basename)
        run.path = results_dir
        self.__log.debug('Set path to: %s', run.path)
        run.set_parameters(**model_settings)
        self.__log.debug('Set property run %s for %s', n,
                         model)
        run.operating_conditions = (
            self.operating_conditions['run{}'.format(n)])
        self.__log.debug('Run %s for %s', n, model)
        res = run.run()
        return res

    def _plot_results(self, model, n, res, results_dir):
        fig, ax = plt.subplots()
        for sp in ['tar', 'light_gas', 'char', 'solid',
                   'volatiles']:
            if sp in res:
                ax.plot(res.index, res[sp], label=sp)
        ax.set_xlabel('Time, s')
        ax.set_ylabel('Yield, daf')
        ax.legend(loc='best', frameon=False)
        ax.set_title('Run{} Model {}'.format(n, model))
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_position(('outward', 20))
        ax.spines['left'].set_position(('outward', 20))
        ax.set_ylim([0, 1])
        # ax.spines['left'].set_color(col_right)
        # ax.spines['left'].set_color(col_right)
        ax1 = ax.twinx()
        ax1.plot(res.index, res['T'],
                 label='T', color=col_green)
        ax1.spines['top'].set_visible(False)
        ax1.spines['left'].set_visible(False)
        ax1.spines['bottom'].set_position(('outward', 20))
        ax1.spines['right'].set_position(('outward', 20))
        ax1.spines['right'].set_color(col_green)
        ax1.tick_params(axis='y', colors=col_green)
        ax1.set_ylabel('Other scale', color=col_green)
        ax1.set_ylabel('Temperature, K')
        ax1.set_ylim(
            [res['T'].min() - 100, res['T'].max() + 100])
        ax1.grid(False)
        fig_name = '{name}-{model}-run{run}.png'.format(
            name=self.name, model=model, run=n)
        self.__log.debug('Save plot to %s', fig_name)
        fig.savefig(
            os.path.join(results_dir, fig_name), bbox_inches='tight')

    def _evolution(self, target_conditions, fit_dict, fit_settings,
                   results_dir, n_p=1):
        '''
        Perform calibration fitting of the empirical model using
        results of the detailed model.

        Parameters
        ----------
        target_conditions: list
            List of target conditions for the calibration. Each entry
            of the list contains:
            `{t: array, 'y': array, operating_conditions: array}`
            `t` and `y` time and volatile yield arrays of length
            N_points.
            operating_conditions: array (2, N_cond) containing the op.
            conditions.
        fit_dict: dict
        fit_settings: dict
            Dictionary containing settings for the evolution algorithm.
        results_dir: str
            Path where results are stored
        n_p: int
            Number of processors for the evolution
        '''
        model = fit_settings['model']
        self.__log.debug('Fit with model %s', model)
        parameters_min = fit_settings['parameters_min']
        parameters_max = fit_settings['parameters_max']
        # parameters_init = fit_settings['parameters_init']
        method = fit_settings['method']
        fit_results = {}
        if method == 'evolve':
            # Define properties of evolutionary model
            npop = fit_settings['npop']
            ngen = fit_settings['ngen']
            mu = fit_settings['mu']
            lambda_ = fit_settings['lambda_']
            cxpb = fit_settings['cxpb']
            mutpb = fit_settings['mutpb']

            # Define Evolution method
            # add a binary field in input yaml for running binary
            # fitting
            Evolution = pkp.evolution.EvolutionBinary \
                if fit_settings.get('binary', False) else \
                pkp.evolution.Evolution

            # Init Evolution
            ga = Evolution(npop=npop, ngen=ngen, cxpb=cxpb, mutpb=mutpb,
                           mu=mu, lambda_=lambda_)
            self.__log.debug('Init GA %s', ga)
            ga.empirical_model = getattr(pkp.empirical_model, model)
            self.__log.debug('Set GA model %s', ga.empirical_model)

            # Define the range of parameters
            ga.parameters_range(parameters_min=parameters_min,
                                parameters_max=parameters_max)
            self.__log.debug('Set GA par range %s, %s',
                             ga._parameters_min, ga._parameters_max)

            # set target conditions
            [ga.set_target(
                t=res['t'], y=res['y'],
                operating_conditions=self.operating_conditions[run])
             for run, res in target_conditions.items()]

            # Register the DEAP toolbox and do the evolution! (Pearl
            # Jam)
            ga.register()
            best = ga.evolve(n_p=n_p, verbose=True)

            # TODO this has to be done inside the Evolution class
            fit_results['best'] = best
            self.__log.info('Best population: %s', fit_results['best'])

            # report only last iteration
            fit_results['log'] = ga.log[-1]
            fit_results['Coal'] = {
                'name': self.name,
                'ultimate_analysis': self.ultimate_analysis,
                'proximate_analysis': self.proximate_analysis,
                'proximate_analysis_daf': self.proximate_analysis_daf,
                'pressure': self.pressure
            }

            # run model and add to fit_results
            det_model, fitname = fit_dict['model'], fit_dict['fit']
            m = ga.empirical_model(fit_results['best'])
            emp_model = m.__class__.__name__
            self.__log.debug('Emp model %s', emp_model)
            # filename = '{}_{}_{}'.format(fitname, det_model, emp_model)
            filename = '{name}-{fit}-{det}-{emp}'.format(
                name=self.name, fit=fitname, det=det_model,
                emp=emp_model)

            # calculate postulate substance
            # it is calculated only for fix y0 variables

            # plot results (evolution history)
            self._plot_evolution(det_model, filename, fitname, ga,
                                 results_dir)
            # plot yield
            self._plot_yieldfit(det_model, emp_model, filename,
                                fit_dict, fit_results, fitname, m,
                                results_dir, target_conditions)
            # calc postulate species
            if 'y0' in m.parameters_names:
                y0 = fit_results['best']['y0']
            elif emp_model == 'C2SM':
                y0 = np.mean([fit_results[run]['y'][-1]
                              for run in sorted(target_conditions)])
                self.__log.debug('Average y0 for C2SM %s', y0)
            fit_results[
                'postulate_volatiles'] = self._postulate_species(y0)
            fit_results['empirical_comp'] = self._emirical_composition(
                y0, tar=0.3, CO=0.1)
        else:
            raise NotImplementedError(
                'Fit method {} not implemented!'.format(method))
        return fit_results

    def _plot_yieldfit(self, det_model, emp_model, filename, fit_dict,
                       fit_results, fitname, m, results_dir,
                       target_conditions):
        '''
        Plot comparison between the yields of the detailed model used as
        target and the yields obtained by the fitted empirical models

        Parameters
        ----------
        det_model: str
            Detailed model name
        emp_model: str
            Empirical Model
        filename: str
            Basename for file output
        fit_dict: dict
            Fit input dictionary
        fit_results: dict
            Fit results dictionary
        fitname: str:
            Name of fit
        m:
        results_dir: str
            Name of results directory
        target_conditions:

        Returns
        -------

        '''
        self.__log.debug('Plot yields')
        fig, ax = plt.subplots()
        runs = list(sorted(target_conditions))
        for i, run in enumerate(runs):
            fit_results[run] = {}
            res = target_conditions[run]
            if i == 0:
                l = '{} {}'.format(run, det_model)
            else:
                l = run
            self.__log.debug('Plot %s ', run)
            ax.plot(res['t'], res['y'], label=l, color=colors[i],
                    linestyle='solid')
            # use list for exporting files
            # fit_results[run]['t'] = res['t'].tolist()
            # fit_results[run]['y'] = res['y'].tolist()
            fit_results[run]['t'] = res['t']
            fit_results[run]['y'] = res['y']
            m.operating_conditions = self.operating_conditions[run]
            t_fit, y_fit = m.run(res['t'])
            if y_fit.ndim == 2:
                y_fit = y_fit[:, 0]
            # fit_results[run]['y_fit'] = y_fit.tolist()
            fit_results[run]['y_fit'] = y_fit
            if i == 0:
                l = '{} {}'.format(run, m.__class__.__name__)
            else:
                l = None
            ax.plot(t_fit, y_fit, color=colors[
                i], linestyle='dashed', label=l)
        ax.set_ylabel('Yield {}'.format(fit_dict['species']))
        ax.set_xlabel('t, s')
        # add an extra legend
        # http://matplotlib.org/users/legend_guide.html#multiple-legend
        ax.add_artist(plt.legend(ax.lines[::2], runs,
                                 loc='upper right', frameon=False))
        ax.legend(ax.lines[:2],
                  [det_model, m.__class__.__name__],
                  loc='lower right', frameon=False)
        ax.set_title(
            'Fit {} from {} with {} ({})'.format(
                fit_dict['species'],
                det_model,
                emp_model,
                fitname))
        fig.savefig(os.path.join(results_dir,
                                 'yield_{}.png'.format(filename)),
                    bbox_inches='tight')
        plt.close(fig)

    def _plot_evolution(self, det_model, filename, fitname, ga,
                        results_dir):
        '''
        Plot the evolution history

        Parameters
        ----------
        det_model
        filename
        fitname
        ga
        results_dir

        Returns
        -------

        '''
        color = 'black'
        color_min = 'red'
        fig, ax = plt.subplots()
        fit_min, fit_max, fit_avg, fit_std = ga.log.select(
            'min', 'max', 'avg', 'std')
        ax.plot(fit_min, label='Min', color=color_min)
        ax.plot(fit_max, label='Max', color=color)
        ax.plot(fit_avg, label='Avg', color=color,
                linestyle='dashed')
        ax.set_yscale('log')
        ax.legend(loc='best')
        ax.set_xlabel('N. generations')
        ax.set_ylabel('Fitness')
        ax.set_title(
            'Fit {} with {}: fitness evolution ({})'.format(
                det_model, ga.empirical_model.__class__.__name__,
                fitname))
        fig.savefig(os.path.join(
            results_dir,
            '{}-evolution.png'.format(filename)))
        plt.close(fig)

    def _postulate_species(self, y0, mw=200.0):
        '''
        Calculate the volatile composition for the fitted empirical
        model assuming a unique *postulate* species. The composition is
        calculated assuming that char is composed only by carbon and
        the remaining carbon and other elements goes to the postulate
        species.

        Parameters
        ----------
        y0: float
            Final volatile yield
        mw: float
            Postulate species molecular weight

        Returns
        -------
        dict:
            Dictionary containing the composition of the postulate
            species with molecular weight, enthalpy of formation and
            other information.
        '''
        assert 0 < y0 < 1, 'Define y0 between 0 and 1'
        ua_char = {el: (1 if el == 'C' else 0)
                   for el in self.ultimate_analysis}
        molecule = {el: ((val - (1 - y0) * ua_char[el]) * mw /
                         M_elements[el] / y0)
                    for el, val in self.ultimate_analysis.items()}
        molecule_name = ''.join('{}_{:4.3f} '.format(el, molecule[el])
                                for el in ['C', 'H', 'O', 'N', 'S'])

        lhv_vol = (self.lhv_daf - (1 - y0) * self.lhv_char) / y0

        nu = {
            'CO2': molecule['C'],
            'H2O': (molecule['H'] * 0.5),
            # reactant is negative
            'O2': -(-0.5 * molecule['O'] + molecule['C'] +
                    0.25 * molecule['H']),
            'SO2': molecule['S']
        }

        hf_vol = np.sum(n * hf[el]
                        for el, n in nu.items()) + lhv_vol * mw

        postulate_dict = {
            'name': molecule_name,
            'formula': molecule,
            'molecular_weight': mw,
            'hf': hf_vol,
            'y0': y0
        }

        return postulate_dict

    def _emirical_composition(self, y0, tar, CO):
        '''
        Set the empirical composition of volatiles using the method by
        NAME.
        http://www.sciencedirect.com/science/article/pii/S0255270104001916

        Parameters
        ----------
        y0: float
            Final volatile yield
        tar: mass fraction of tar in volatiles
        CO: fraction of O converted to CO
        '''
        def el_fraction(sp, el):
            '''Mass fraction of element el in species sp'''
            if sp == 'char':
                return 1.0 if el == 'C' else 0.0
            else:
                return (gas.species(sp).composition.get(el, 0) *
                        gas.atomic_weight(el) /
                        gas.molecular_weights[gas.species_index(sp)])

        def calc_remaining(comp):
            '''Remaining fraction of each elements'''
            return {el: (ua - tot_el_fraction(comp, el))
                    for el, ua in ultimate_analysis.items()}

        def tot_el_fraction(comp, element):
            '''Calc the total element fraction of the given element'''
            return np.sum([val * el_fraction(sp, element)
                           for sp, val in comp.items()])

        sum_ua = (sum(self.ultimate_analysis.values()) -
                  self.ultimate_analysis['S'])
        ultimate_analysis = {
            el: v / sum_ua
            for el, v in self.ultimate_analysis.items()
            if el != 'S'}
        self.__log.debug('Update ultimate_analysis %s',
                         ultimate_analysis)

        gas = cantera.Solution('52.xml')
        composition = {}
        composition['char'] = 1 - y0
        # assume tar as C6H6

        # assume N -> N2
        composition['N2'] = ultimate_analysis['N']
        composition['CO'] = (CO * ultimate_analysis['O'] /
                             el_fraction('CO', 'O'))
        composition['CO2'] = ((1 - CO) * ultimate_analysis['O'] /
                              el_fraction('CO2', 'O'))
        self.__log.debug('Vol composition %s', composition)

        remaining = calc_remaining(composition)
        self.__log.debug('Remaining element after CO/CO2/N2: %s',
                         remaining)

        C_in_tar = el_fraction('C6H6', 'C')
        self.__log.debug('C in TAR %s', C_in_tar)
        if C_in_tar * tar > remaining['C']:
            composition['C6H6'] = remaining['C'] / C_in_tar
            self.__log.debug('C in Tar > remaining C -> set tar: %s',
                             composition['C6H6'])
        else:
            composition['C6H6'] = tar
            self.__log.debug('Set TAR as %s', tar)

        # recalculate remaining
        remaining = calc_remaining(composition)
        self.__log.debug('Remaining element after tar: %s', remaining)

        c_to_h_mass = remaining['C'] / remaining['H']
        c_to_h_molar = c_to_h_mass * M_elements['H'] / M_elements['C']
        self.__log.debug('C/H molar: %s', c_to_h_molar)

        if 0 <= c_to_h_molar <= 0.5:
            # use C2H4 and H2
            composition['C2H4'] = remaining['C'] / el_fraction(
                'C2H4', 'C')
            self.__log.debug('C2H4 %s', composition['C2H4'])
        elif 0.5 < c_to_h_molar < 1:
            # use C6H6
            composition['C6H6'] = (composition['C6H6'] +
                                   remaining['C'] /
                                   el_fraction('C6H6', 'C'))
            self.__log.debug('Update C6H6 %s', composition['C6H6'])

        self.__log.debug('Remaining element after C: %s', remaining)
        remaining = calc_remaining(composition)

        composition['H2'] = remaining['H']
        remaining = calc_remaining(composition)
        self.__log.debug('Remaining %s', remaining)
        self.__log.debug('Final Vol composition %s', composition)

        emp_dict = {'composition': composition,
                    'heat_pyro': self.heat_of_pyrolysis(
                        composition, gas)}
        return emp_dict

    def calc_element_fraction(self, element, species):
        '''
        Calculate element fraction of the given species
        '''
        if species in ('Char', 'Solid'):
            if element == 'C':
                return 1.0
            else:
                return 0.0
        elif species == 'Other':
            return 0.0
        else:
            i = self.gas.element_index(element)
            return (self.gas.atomic_weights[i] *
                    self.gas.species(species).composition.get(element, 0) /
                    self.gas.molecular_weights[
                    self.gas.species_index(species)])

    def heat_of_volatiles(self, composition, gas):
        '''Heat of volatiles including char'''
        return sum(val * self.heat_of_reaction_species(sp, gas)
                   for sp, val in composition.items())

    def heat_of_pyrolysis(self, composition, gas):
        '''
        Heat of pyrolysis. It is defined as the total heat released
        during pyrolysis per unit of volatiles.

        Positive heat of pyrolysis means that extra heat is removed to
        fullfill the energy balance of coal.

        Same convention used for latent heat of evaporation!!!
        '''
        heat_vol = self.heat_of_volatiles(composition, gas)
        dh_pyro = self.lhv_daf - heat_vol
        return -dh_pyro / (1 - composition['char'])

    def heat_of_reaction_species(self, sp, gas):
        '''
        Calculate the heat of reaction for the species sp
        '''
        T_ref = 273
        if sp in ('N2'):
            return 0.0
        elif sp == 'char':
            return self.lhv_char
        spc = gas.species(sp)
        hf = spc.thermo.h(T_ref)
        n_o2 = 0.5 * (
            2 * spc.composition.get('C', 0) +
            0.5 * spc.composition.get('H', 0) -
            spc.composition.get('O', 0))
        n_co2 = spc.composition.get('C', 0)
        # print 'n_co2', n_co2
        n_h2o = spc.composition.get('H', 0) * 0.5
        # print 'n_h2o', n_h2o
        mw = gas.molecular_weights[gas.species_index(sp)]

        h_o2 = gas.species('O2').thermo.h(T_ref)
        h_co2 = gas.species('CO2').thermo.h(T_ref)
        h_h2o = gas.species('H2O').thermo.h(T_ref)

        heat_molar = (hf + n_o2 * h_o2 - n_co2 *
                      h_co2 - n_h2o * h_h2o)
        return heat_molar / mw
