#!/usr/bin/python
"""
pkp-cli

Usage:
      pkp-cli -h | --help
      pkp-cli generate (--json-input=<string> | --file-input=<loc>) --results-folder=<loc>
      pkp-cli generate-only (--json-input=<string> | --file-input=<loc>) --results-folder=<loc>
      pkp-cli fit-only      (--json-input=<string> | --file-input=<loc>) --fit-target=<string> --results-folder=<loc>

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
    args = docopt(__doc__)
    print args
    inp = args['--file-input']
    if args['generate'] or args['generate-only']:
        json = args['--json-input']
        pre = generate(json_string=json, folder=inp)
        if args['generate-only']:
            sys.exit(0)
        else:
            res = fit(
                folder=inp,
                results=pre,
                selectPyrolModel="constantRate"
            )
            fit_result = res.startFittingProcedure(pre)
            print fit_result._tsv
    elif args['fit-only']:
        pass
    
