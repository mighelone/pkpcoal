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
        output = args.get('--results-folder', False)
        pre = generate(json_string=json, input_file=inp, output_folder=output)
        if args['generate-only']:
            sys.exit(0)
        else:
            res = fit(
                input_file=inp,
                results=pre,
                selectPyrolModel="constantRate"
            )
            fit_result = res.startFittingProcedure(pre)
            write_tsv(output, fit_result._tsv)
            plot_and_save(output , pre[0], fit_result)
    elif args['fit-only']:
        pass

def write_tsv(path, content):
    print "writting "  + path
    target = path + '/fit'
    #if not os.path.exists(path + '/fit' )
    with open(target, 'w') as f:
        l = content.replace(' ','\t')
        f.write(l)

def plot_and_save(path, pre, model):
    import matplotlib.pyplot as plt
    time = pre['time']
    ftot_cpd = pre['ftot']
    ftot_model = model.res['ftot']
    print ftot_model
    fig, axs = plt.subplots()
    colors = ['c', 'm', 'y', 'k','b','g','r']
    print time
    axs.scatter(
        x=time,
        y=ftot_cpd,
        marker='.',
        label = 'ftot_cpd',
        color = colors[0],
    )
    axs.scatter(
         x=time,
         y=ftot_model,
         color = 'k',
         label = "ftot_constRate",
         marker='x',
    )
    axs.set_xlim([0, 0.00004])
    axs.set_ylim([0, 1])
    plt.legend()
    fig.savefig(path + '/plot.png')
