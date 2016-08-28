'''
PKP Pyrolysis Kinetic Preprocessor
==================================
'''
from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals

from autologging import logged
import ruamel_yaml as yaml
import os

from pkp.cpd import CPD
from pkp.polimi import Polimi
import pkp.evolution
import numpy as np

models = ['CPD', 'Polimi']


@logged
class ReadConfiguration(object):
    '''
    Read configuration file for PKP
    '''

    def __init__(self, yml):
        super(ReadConfiguration, self).__init__()
        if isinstance(yml, (str, unicode)):
            with open(yml, 'r') as f:
                yml_input = yaml.load(f)
        elif isinstance(yml, dict):
            yml_input = yml
        else:
            raise ValueError('Define yml as file name or dictionary')

        # coal settings
        coal_settings = yml_input['Coal']
        self.proximate_analysis = coal_settings['proximate_analysis']
        self.ultimate_analysis = coal_settings['ultimate_analysis']
        # convert HHV from MJ/kg to J/kg
        self.HHV = coal_settings['HHV'] * 1e6
        self.rho_dry = coal_settings['rho_dry']

        # Solver settings
        [setattr(self, model, yml_input[model])
         for model in models]

        # Solver settings
        self.operating_conditions = yml_input['operating_conditions']

        # fit settings
        self.fit_settings = yml_input['FIT']


@logged
class PKPRunner(ReadConfiguration):
    '''
    Run PKP case
    '''
    models = models

    def run(self, results_dir=None):
        results_dir = self.set_results_dir(results_dir)
        self.__log.info('Run models %s', self.models)
        run_results = {}
        fit_results = {}
        for model in self.models:
            self.__log.debug('Model %s', model)
            model_settings = getattr(self, model)
            results = self._run_model(model=model,
                                      results_dir=results_dir)
            if results:
                run_results[model] = results
                self.__log.debug('Finish run %s %s', model,
                                 results.keys())
            if model_settings['fit']:
                fit_results[model] = {}
                for fitname, fit in model_settings['fit'].iteritems():
                    self.__log.debug('Fit %s model %s', fit, model)
                    target_conditions = {
                        run: {'t': np.array(res.index),
                              'y': np.array(res[fit['species']])}
                        for run, res in results.iteritems()}
                    fit_results[model][fitname] = self._fit(
                        target_conditions, fitname, fit, results_dir)
        return run_results, fit_results

    @staticmethod
    def set_results_dir(results_dir):
        if results_dir is None:
            results_dir = os.getcwd()
        return results_dir

    def _run_model(self, model, results_dir):
        '''
        Run simulations for the given model
        '''
        self.__log.info('Run %s model', model)
        model_settings = getattr(self, model)
        self.__log.debug('Model %s active %s', model,
                         model_settings['active'])
        if model_settings['active']:
            results = {}
            self.__log.debug('Run %s',
                             self.operating_conditions['runs'])
            for n in range(
                    self.operating_conditions['runs']):
                self.__log.debug(
                    'Initialize run %s for %s', n, model)
                run = globals()[model](
                    ultimate_analysis=self.ultimate_analysis,
                    proximate_analysis=self.proximate_analysis,
                    pressure=self.operating_conditions['pressure'],
                    name='{}-Run{}'.format(model, n)
                )
                run.path = results_dir
                self.__log.debug('Set path to: %s', run.path)
                # run.set_parameters(**getattr(self.reader, model))
                run.set_parameters(**model_settings)
                self.__log.debug('Set property run %s for %s', n,
                                 model)
                run.operating_conditions = (
                    self.operating_conditions['run{}'.format(n)])
                self.__log.debug('Run %s for %s', n, model)
                res = run.run()
                results['run{}'.format(n)] = res
                self.__log.debug('Finish run %s', results.keys())
        else:
            results = None
        return results

    def _fit(self, target_conditions, fitname, fit_settings,
             results_dir):
        model = fit_settings['model']
        self.__log.debug('Fit with model %s', model)
        parameters_min = fit_settings['parameters_min']
        parameters_max = fit_settings['parameters_max']
        parameters_init = fit_settings['parameters_init']
        method = fit_settings['method']
        fit_results = {}
        if method == 'evolve':
            npop = fit_settings['npop']
            ngen = fit_settings['ngen']
            mu = fit_settings['mu']
            lambda_ = fit_settings['lambda_']
            cxpb = fit_settings['cxpb']
            mutpb = fit_settings['mutpb']

            ga = pkp.evolution.Evolution(npop=npop, ngen=ngen,
                                         cxpb=cxpb, mutpb=mutpb)
            self.__log.debug('Init GA %s', ga)
            ga.empirical_model = getattr(pkp.empirical_model, model)
            self.__log.debug('Set GA model %s', ga.empirical_model)
            ga.parameters_range(parameters_min=parameters_min,
                                parameters_max=parameters_max)

            self.__log.debug('Set GA par range %s, %s',
                             ga._parameters_min, ga._parameters_max)

            [ga.set_target(
                t=res['t'], y=res['y'],
                operating_conditions=self.operating_conditions[run])
             for run, res in target_conditions.iteritems()]

            # self.__log.debug('Op. conditions %s',
            #                 ga.operating_conditions)

            ga.register()
            fit_results['best'] = ga.evolve(mu=mu, lambda_=lambda_)
            # run model and add to fit_results
        else:
            raise NotImplementedError(
                'Fit method {} not implemented!'.format(method))
        return None
