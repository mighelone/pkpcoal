#!/usr/bin/python
"""
PKP

Usage:
      PKP.py -h | --help
      PKP.py generate (--json-input=<string> | --file-input=<loc>) --results-folder=<loc>
      PKP.py fit      (--json-input=<string> | --file-input=<loc>) --fit-target=<string> --results-folder=<loc>

Options:
    -h  --help              Shows this screen
    --json-input=<string>   Foo
    --file-input=<loc>      Bar
    --results-folder=<loc>  Baz
    --fit-target=<string>   Bla

"""
try:
    import PKP.src
except:
    import src
import sys
import os
import platform
from docopt import docopt

def main():
    arguments = docopt(__doc__)
    print arguments
    if arguments['generate']:
        generate(json_string=arguments['--json-input'], folder=arguments['--file-input'])
    main()
    
    
