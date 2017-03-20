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
import logging
try:
    import ruamel.yaml as yaml
except:
    try:
        import ruamel_yaml as yaml
    except:
        print('Loading standard yaml module ...\n'
              'Note that it might give problems with'
              ' exponential decoding')
        import yaml

import os
import numpy as np
import pandas as pd

# detailed models
# they must be loaded here if you want to use them!
from pkp.cpd import CPD
# from cpd import CPD
from pkp.cpd_fortran import CPD as CPDfortran

try:
    from pkp.polimi import Polimi
    from pkp.biopolimi import BioPolimi
    models = ['CPD', 'CPDfortran', 'Polimi', 'BioPolimi']
except ModuleNotFoundError:
    logger = logging.getLogger('pkp.runner')
    logger.warning(
        'Cantera not available. Polimi and BioPolimi models cannot be used!')
    models = ['CPD', 'CPDfortran']

# optimization
import pkp.evolution
import pkp.minimize

import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('bmh')

# try:
#     plt.style.use(['mystyle', 'mystyle-vega'])
# except:
#     plt.style.use('ggplot')

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
col_right = "#C54E6D"
col_left = "#009380"


# raise error from numpy
# np.seterr(all='raise')


# yaml serialization of numpy objects
# representer data for yaml


def ndarray_representer(dumper, data):
    # return dumper.represent_data(data.tolist())
    return dumper.represent_data([data.min(), data.max()])


def npfloat64_representer(dumper, data):
    return dumper.represent_data(data.tolist())


def npint_representer(dumper, data):
    return ndarray_representer(dumper, data)


def tuple_representer(dumper, data):
    return dumper.represent_data(list(data))


yaml.add_representer(np.ndarray, ndarray_representer)
yaml.add_representer(np.float64, npfloat64_representer)
yaml.add_representer(tuple, tuple_representer)


def runs_iterator(op_cond):
    '''Iterate over operating_conditions'''
    def get_run(i):
        return 'run{}'.format(i)
    i = 0
    while get_run(i) in op_cond:
        i += 1
        yield i - 1


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
                yml_input = yaml.safe_load(f)
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
        self.hhv = coal_settings['HHV'] * 1e6
        self.rho_dry = coal_settings['rho_dry']

        # Solver settings
        [setattr(self, model, yml_input[model])
         for model in models if model in yml_input]

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

    def run(self, results_dir=None, n_p=1, run_only=False):
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

        coal = {
            'name': self.name,
            'ultimate_analysis': self.ultimate_analysis,
            'proximate_analysis': self.proximate_analysis,
            'proximate_analysis_daf': self.proximate_analysis_daf,
            'HHV ar': (self.hhv / 1e6, 'MJ/kg'),
            'HHV daf': (self.hhv_daf / 1e6, 'MJ/kg'),
            'LHV daf': (self.lhv_daf / 1e6, 'MJ/kg'),
            'rho_dry': self.rho_dry,
            'operating_conditions': self.operating_conditions
        }
        run_results = {'coal': coal, 'version': pkp.__version__}
        fit_results = {'coal': coal, 'version': pkp.__version__}
        for model in self.models:
            if hasattr(self, model):
                model_settings = getattr(self, model)
                if model_settings['active']:
                    self.__log.info('Run model %s', model)
                    results = self.run_model(model=model,
                                             results_dir=results_dir)
                    self.__log.debug('Finish run %s %s',
                                     model, results.keys())
                    if results:
                        run_results[model] = results
                        if model_settings['fit'] and not run_only:
                            self.__log.info('Start fit of %s model', model)
                            fit_results[model] = self.fit_detmodel(
                                model, model_settings['fit'], n_p, results,
                                results_dir)
                    else:
                        self.__log.warning('No results for %s', model)

        yml_fit = os.path.join(
            results_dir, '{name}-fitreport.yml'.format(name=self.name))
        self.__log.debug('Export fit report to %s', yml_fit)
        with open(yml_fit, 'w') as f:
            # yaml.dump(clean_dict(fit_results), f, indent=4)
            yaml.dump(fit_results, f, indent=4)

        return run_results, fit_results

    def fit_detmodel(self, model, model_settings, n_p, results,
                     results_dir):
        '''
        Run all fitting of the given detailed model

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
        runs: int
            Number of results to include in the calibration

        Returns
        -------
        fit_results: dict
            Contains results of fitting

        '''
        # loop over the fitting runs, fit0, fit1, etc.
        fit_results = {}
        for fitname, fit in model_settings.items():
            if fit.get('active', True):
                self.__log.debug(
                    'Fit %s model with %s', model, fit['model'])
                target_conditions = {
                    run: {'t': np.array(res['t']),
                          'y': np.array(res[fit['species']])}
                    for run, res in results.items()}
                self.__log.debug('runs calibration %s',
                                 list(target_conditions.keys()))
                fit_dict = {'model': model,
                            'fit': fitname,
                            'species': fit['species']}
                fit_results[fitname] = self.fit_single(
                    results,
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

    def run_model(self, model, results_dir):
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
            # new implementation run all the cases
            # for n in range(
            #        self.operating_conditions['runs']):
            for n in runs_iterator(self.operating_conditions):
                self.__log.info('Run %s with %s model', n, model)
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
                for n in runs_iterator(self.operating_conditions)]
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
                # name='{}-{}-Run{}'.format(self.name, model, n)
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
                ax.plot(res['t'], res[sp], label=sp)
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
        ax1.plot(res['t'], res['T'],
                 label='T', color=col_right)
        ax1.spines['top'].set_visible(False)
        ax1.spines['left'].set_visible(False)
        ax1.spines['bottom'].set_position(('outward', 20))
        ax1.spines['right'].set_position(('outward', 20))
        ax1.spines['right'].set_color(col_right)
        ax1.tick_params(axis='y', colors=col_right)
        ax1.set_ylabel('Other scale', color=col_right)
        ax1.set_ylabel('Temperature, K')
        ax1.set_ylim(
            [res['T'].min() - 100, res['T'].max() + 100])
        ax1.grid(False)
        fig_name = '{name}-{model}-run{run}.png'.format(
            name=self.name, model=model, run=n)
        self.__log.debug('Save plot to %s', fig_name)
        fig.savefig(
            os.path.join(results_dir, fig_name), bbox_inches='tight')

    def fit_single(self, results, target_conditions, fit_dict,
                   fit_settings, results_dir, n_p=1):
        '''
        Perform calibration fitting of the empirical model using
        results of the detailed model.

        Parameters
        ----------
        results: dict
            Results dictionary
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
        # parameters_init = fit_settings['parameters_init']
        method = fit_settings['method']
        fit_results = {}
        det_model, fitname = fit_dict['model'], fit_dict['fit']
        emp_model = fit_settings['model']
        filename = '{name}-{fit}-{det}-{emp}'.format(
            name=self.name, fit=fitname, det=det_model,
            emp=emp_model)

        runs = self._operating_conditions['runs']
        target_conditions_used = {key: value for key, value
                                  in target_conditions.items()
                                  if int(key[3:]) < runs}
        self.__log.debug('Runs used for calibration %s',
                         list(target_conditions_used.keys()))

        if 'evolve' in method:
            self.__log.info('%s Evolution to fit %s with %s',
                            fitname, det_model, emp_model)
            # Define properties of evolutionary model
            best, ga = self.evolve(n_p, fit_results, fit_settings,
                                   target_conditions_used)
            # plot results (evolution history)
            self._plot_evolution(det_model, filename, fitname, ga,
                                 results_dir)

            # this is the initial parameters for fmin
            parameters_init = best
            emp_model_class = ga.empirical_model
        else:
            parameters_init = fit_settings['parameters_init']

        if 'min' in method:
            self.__log.info('%s Minimization to fit %s with %s',
                            fitname, det_model, emp_model)
            best, fmin = self.minimization(fit_results, fit_settings,
                                           target_conditions_used,
                                           parameters_init)
            emp_model_class = fmin.empirical_model

            # run optimized empirical model
        m = emp_model_class(best)
        self.__log.debug('Emp model %s', emp_model)

        # plot yield
        self._plot_yieldfit(det_model, emp_model, filename,
                            fit_dict, fit_results, fitname, m,
                            results_dir, target_conditions)
        # calc postulate species
        if 'y0' in m.parameters_names():
            y0 = best['y0']
        else:
            y0 = np.mean([fit_results[run]['y'][-1]
                          for run in sorted(target_conditions)])
            self.__log.debug('Average y0 %s', y0)

        tar_mean = np.mean([r.iloc[-1]['tar']
                            for r in results.values()])
        self.__log.debug('columns %s', results['run0'].columns)
        try:
            co_mean = np.mean([r.iloc[-1]['CO'] for r in results.values()])
        except KeyError:
            co_mean = 0.1
        self.__log.debug('tar mean %s', tar_mean)
        self.__log.debug('CO mean %s', co_mean)
        self.__log.debug('y0 mean %s', y0)

        fit_results[
            'postulate_volatiles'] = self.postulate_species(y0)
        fit_results['empirical_comp'] = self.empirical_composition(
            y0, tar=tar_mean, CO=co_mean)
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
            ax.plot(t_fit, y_fit, color=colors[i],
                    linestyle='dashed', label=l)
        ax.set_ylabel('Yield {}'.format(fit_dict['species']))
        ax.set_xlabel('t, s')
        # add an extra legend
        # http://matplotlib.org/users/legend_guide.html#multiple-legend
        nruns = self.operating_conditions['runs']
        runs_label = [r + '(fitted)' if i < nruns else r
                      for i, r in enumerate(runs)]
        ax.add_artist(plt.legend(ax.lines[::2], runs_label,
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
        fig.savefig(os.path.join(
            results_dir, '{}-yields.png'.format(filename)),
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

    def evolve(self, n_p, fit_results, fit_settings, target_conditions):
        model = fit_settings['model']
        self.__log.debug('Evolution fit with model %s', model)
        npop = fit_settings['npop']
        ngen = fit_settings['ngen']
        mu = fit_settings['mu']
        lambda_ = fit_settings['lambda_']
        cxpb = fit_settings['cxpb']
        mutpb = fit_settings['mutpb']
        skip = fit_settings.get('skip', 1)

        parameters_min = fit_settings['parameters_min']
        parameters_max = fit_settings['parameters_max']

        # Define Evolution method
        # add a binary field in input yaml for running binary
        # fitting
        Evolution = pkp.evolution.EvolutionBinary \
            if fit_settings.get('binary', False) else \
            pkp.evolution.Evolution

        # Init Evolution
        self.__log.debug('Set skip=%s', skip)
        ga = Evolution(npop=npop, ngen=ngen, cxpb=cxpb, mutpb=mutpb,
                       mu=mu, lambda_=lambda_, skip=skip)
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
        self.__log.debug('Best: %s', best)

        fit_results['evolve'] = {
            'best': {p: (best[p], unit)
                     for p, unit in zip(
                ga.empirical_model.parameters_names(),
                ga.empirical_model.parameters_units())
            },
            'log': ga.log[-1]}

        # report only last iteration
        self.__log.info('Best population: %s',
                        fit_results['evolve']['best'])

        return best, ga

    def minimization(self, fit_results, fit_settings, target_conditions,
                     init):
        model = fit_settings['model']
        self.__log.debug('Minimization fit with model %s', model)

        parameters_min = fit_settings['parameters_min']
        parameters_max = fit_settings['parameters_max']

        fmin = pkp.minimize.Minimization()

        self.__log.debug('Init fmin %s', fmin)
        fmin.empirical_model = getattr(pkp.empirical_model, model)
        self.__log.debug('Set fmin model %s', fmin.empirical_model)

        # Define the range of parameters
        fmin.parameters_range(parameters_min=parameters_min,
                              parameters_max=parameters_max)
        self.__log.debug('Set fmin par range %s, %s',
                         fmin._parameters_min, fmin._parameters_max)

        # set target conditions
        [fmin.set_target(
            t=res['t'], y=res['y'],
            operating_conditions=self.operating_conditions[run])
         for run, res in target_conditions.items()]

        # Register the DEAP toolbox and do the evolution! (Pearl
        # Jam)
        best = fmin.run(initial=init)
        self.__log.debug('Best: %s', best)

        fit_results['fmin'] = {
            'best': {p: (best[p], unit)
                     for p, unit in zip(
                fmin.empirical_model.parameters_names(),
                fmin.empirical_model.parameters_units())},
            'report': dict(fmin.results)
        }

        self.__log.info('Minimized value: %s',
                        fit_results['fmin']['best'])

        return best, fmin
