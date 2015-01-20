#!/usr/bin/python
"""
PKP

Usage:
      pkp-cli -h | --help
      pkp-cli generate (--json-input=<string> | --file-input=<loc>) --results-folder=<loc>
      pkp-cli fit      (--json-input=<string> | --file-input=<loc>) --fit-target=<string> --results-folder=<loc>

Options:
    -h  --help              Shows this screen
    --json-input=<string>   Foo
    --file-input=<loc>      Bar
    --results-folder=<loc>  Baz
    --fit-target=<string>   Bla

"""
import sys
import os
import platform
from docopt import docopt

import src
from pkpcli import generate
from pkpcli import fit

def main():
    arguments = docopt(__doc__)
    print arguments
    if arguments['generate']:
        generate(json_string=arguments['--json-input'], folder=arguments['--file-input'])
