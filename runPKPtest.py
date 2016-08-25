from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals

import pkp
import logging
from dictinput import settings

if __name__ == '__main__':
    logger = logging.getLogger('main')
    logger.setLevel(logging.INFO)
    ch = logging.StreamHandler()
    formatter = logging.Formatter(
        '%(name)s:%(levelname)s:%(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    logger.info('Create runner')
    runner = pkp.PKPRunner()
    logger.info('Read settings')
    runner.read(settings)
    logger.info('Start run')
    results = runner.run(results_dir='test')
