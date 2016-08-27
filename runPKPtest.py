#!/usr/bin/env python


from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals

import pkp
import logging
import sys
import argparse

# import coloredlogs
# coloredlogs.install(level='DEBUG')
logging.basicConfig(
    level=logging.INFO, stream=sys.stdout,
    format="%(levelname)s:%(name)s:%(funcName)s:%(message)s")

yml_file_default = 'input.yml'

if __name__ == '__main__':
    #logger = logging.getLogger('main')
    # logger.setLevel(logging.INFO)
    #ch = logging.StreamHandler()
    # formatter = logging.Formatter(
    #    '%(name)s:%(levelname)s:%(message)s')
    # ch.setFormatter(formatter)
    # logger.addHandler(ch)

    parser = argparse.ArgumentParser(
        description="PKP Runner"
    )
    parser.add_argument('yml_file', action="store",
                        help="YAML input file")
    parser.add_argument('-n', action="store", dest="n_p", type=int,
                        default=1, help="Number of processor (NOT WORKING)")
    parser.add_argument('-o', action="store", dest="results_dir", type=str,
                        default="Results", help="Results directory")
    parser.add_argument('-d', action="store_true", dest="debug",
                        help="Print debug messages")

    argument = parser.parse_args()
    if argument.debug:
        logging.basicConfig(level=logging.DEBUG)

    yml_file = argument.yml_file
    logging.info('Create runner and read settings')
    runner = pkp.PKPRunner(yml_file)
    logging.info('Start run')
    results = runner.run(results_dir=argument.results_dir)

    # plot results
