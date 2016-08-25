from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals

import pkp
import logging

yml_file = 'input.yml'

if __name__ == '__main__':
    logger = logging.getLogger('main')
    logger.setLevel(logging.INFO)
    ch = logging.StreamHandler()
    formatter = logging.Formatter(
        '%(name)s:%(levelname)s:%(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    logger.info('Create runner and read settings')
    runner = pkp.PKPRunner(yml_file)
    logger.info('Start run')
    results = runner.run(results_dir='test')
