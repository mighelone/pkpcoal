from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals

# remove at the end
import sys
import logging
import pkp
import pkp.evolution
import numpy as np
import matplotlib.pyplot as plt
plt.style.use(['mystyle'])


yml_file_default = 'input.yml'

# USER parameters
run_CPD = True
npop = 100  # 40
ngen = 40
cxpb = 0.4  # 0.2
mutpb = 0.6
mu = npop
lambda_ = 50
# best=0.00298792

# best=0.00100143 npop=100, ngen=40, cxpb=0.4, mutpb=0.6
parameters_min = [1e4, 50e6, 0.4]  # A, E, y0
parameters_max = [1e10, 300e6, 0.8]


def calc_CPD():
    '''Calc CPD'''
    runner = pkp.PKPRunner('input.yml')
    runner.operating_conditions['runs'] = 1
    runner.Polimi['active'] = False
    res = runner.run(results_dir='./optdir')
    return res['CPD']['run0']


def read_CPD():
    '''Read CPD results'''
    import pandas as pd
    return pd.read_csv('optdir/CPD_CPD-Run0.csv', index_col=0)

if __name__ == '__main__':
    logging.basicConfig(
        level=logging.DEBUG, stream=sys.stdout,
        format="%(levelname)s:%(name)s:%(funcName)s:%(message)s")
    res = calc_CPD() if run_CPD else read_CPD()
    y_cpd = np.array(res['ftot'])
    t_cpd = np.array(res.index) * 1e-3

    reader = pkp.ReadConfiguration(yml_file_default)
    operating_conditions = reader.operating_conditions['run0']

    ga = pkp.evolution.Evolution(npop=npop, ngen=ngen, cxpb=cxpb,
                                 mutpb=mutpb)
    ga.empirical_model = pkp.empirical_model.SFOR
    ga.operating_conditions = operating_conditions
    ga.set_target(t=t_cpd, y=y_cpd)
    ga.parameters_range(parameters_min=parameters_min,
                        parameters_max=parameters_max)

    ga.register()
    best_parameters = ga.evolve(mu=mu, lambda_=lambda_)

    # run best model
    m = pkp.empirical_model.SFOR(best_parameters)
    m.operating_conditions = operating_conditions
    t, y = m.run(t=t_cpd)

    # plot fit results
    fig, ax = plt.subplots()
    ax.plot(t_cpd, y_cpd, label='CPD')
    ax.plot(t, y, label='SFOR')
    ax.legend()

    # plot convergence
    fig, ax = plt.subplots()
    fit_min, fit_max, fit_avg, fit_sdt = ga.log.select(
        'min', 'max', 'avg', 'std')
    ax.plot(fit_min, color='black')
    ax.plot(fit_max, color='black')
    ax.plot(fit_avg, color='black', linestyle='dashed')

    plt.show()
