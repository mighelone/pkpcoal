from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals
from builtins import dict

import logging
import argparse
import shutil

import pytest
import pkp.runner
import os.path


@pytest.fixture(params=['input_CPD.yml', 'input_bio.yml'])
def input_yml(request):
    return os.path.join('test', request.param)


def test_CPD(input_yml):
    model = os.path.splitext(input_yml)[0].split('_')[1]

    # define a logger
    logger = logging.getLogger('pkp')
    logger.setLevel(logging.DEBUG)

    # define an handler
    handler_info = logging.StreamHandler()
    handler_info = logging.FileHandler(
        'test/{}.log'.format(model), mode='w')
    handler_info.setFormatter(logging.Formatter(
        fmt="%(name)s:%(funcName)s:%(message)s"))
    handler_info.setLevel(logging.DEBUG)

    # define a debug handler
    handler_debug = logging.StreamHandler()
    handler_debug.setFormatter(logging.Formatter(
        fmt="%(levelname)s:%(name)s:%(funcName)s:%(message)s"))
    handler_debug.setLevel(logging.DEBUG)

    logger.addHandler(handler_info)

    # set level of evolution
    logging.getLogger('pkp.evolution').setLevel(logging.INFO)
    logging.getLogger('algorithms').addHandler(handler_info)

    logger.warning('Start PKP run with %s', input_yml)

    res_dir = 'test/Results_{}'.format(model)
    runner = pkp.runner.PKPRunner(input_yml)
    results = runner.run(results_dir=res_dir, n_p=1, run_only=False)

    # best_results = results[1]['CPD']['fit0']['evolve']['best']
    # logger.info('Best results: %s', best_results)

    # no assert because of the random nature of the problem
    shutil.rmtree(res_dir)
    logger.warning('End run %s', input_yml)
