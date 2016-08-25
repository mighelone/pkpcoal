'''
PKP Pyrolysis Kinetic Preprocessor
==================================
'''
from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals

import logging
import ruamel_yaml as yaml
import os

from pkp.cpd import CPD
from pkp.polimi import Polimi

models = ['CPD', 'Polimi']


class ReadConfiguration(object):
    '''
    Read configuration file for PKP
    '''

    def __init__(self, yml):
        self.logger = logging.getLogger(
            'main' + self.__class__.__name__)
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


class PKPRunner(object):
    '''
    Run PKP case
    '''
    models = models

    def __init__(self):
        self.logger = logging.getLogger(
            'main' + self.__class__.__name__)

    def read(self, yml):
        self.logger.info('Read configuration file')
        self.reader = ReadConfiguration(yml)

    def run(self, results_dir=None):
        if results_dir is None:
            results_dir = os.getcwd()
        results = {}
        reader = self.reader
        for model in self.models:
            self.logger.info('Run %s model', model)
            model_settings = getattr(reader, model)
            if model_settings['active']:
                results[model] = {}
                for n in range(
                        self.reader.operating_conditions['runs']):
                    self.logger.debug(
                        'Initialize run %s for %s', n, model)
                    run = globals()[model](
                        ultimate_analysis=reader.ultimate_analysis,
                        proximate_analysis=reader.proximate_analysis,
                        pressure=reader.operating_conditions[
                            'pressure'],
                        name='{}-Run{}'.format(model, n)
                    )
                    run.path = results_dir
                    self.logger.debug('Set property run %s for %s', n,
                                      model)
                    # run.set_parameters(**getattr(self.reader, model))
                    run.set_parameters(**model_settings)
                    run.operating_conditions = (
                        reader.operating_conditions['run{}'.format(n)])
                    self.logger.debug('Run %s for %s', n, model)
                    res = run.run()
                    results[model]['run{}'.format(n)] = res
            return results
