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
            res = fit(input_file=inp)
            fit_result = res.startFittingProcedure(pre)
            #write_tsv(output, fit_result._tsv)
            plot_and_save(output, pre, fit_result)
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
    time = pre[pre.keys()[0]][1]['time']
    print model.res['ftot'].fittedYield()
    fig, axs = plt.subplots()
    fig.set_size_inches(18.5, 10.5)
    colors = ['c', 'm', 'y', 'k', 'b', 'g', 'r', 'k', 'm']
    for runName, run in pre.iteritems():
        for i, (name, data) in enumerate(run[1].iterspecies()):
            axs.plot(time, data, label=name, color=colors[i], ls='--')
            axs.plot(time, model.res[name].fittedYield(),
                     color=colors[i], label = name + "Model")
    # axs.set_xlim([0, 0.00004])
    axs.set_ylim([0, 1])
    plt.legend(loc="center left",
            bbox_to_anchor=(0.95, 0.75))
    fig.savefig(path + '/plot.png')
    with open(path + '/parameters.dat', 'w') as p:
        p.write(str(model.res[model.res.keys()[0]].name))
        p.write(str(model.res[model.res.keys()[0]].paramNames))
        for name, spec in model.res.iteritems():
            p.write("\n{}: {}".format(name, spec.parameter))
