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


@logged
class PKPRunner(ReadConfiguration):
    '''
    Run PKP case
    '''
    models = models

    def run(self, results_dir=None):
        if results_dir is None:
            results_dir = os.getcwd()
        self.__log.info('Run models %s', self.models)
        results = {}
        for model in self.models:
            res = self._run_model(model=model,
                                  results_dir=results_dir)
            if res:
                results[model] = res
        return results

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
        else:
            results = None
        return results
