from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals

from pyevolve import G1DList, GSimpleGA, Selectors  # , Statistics
from pyevolve import Initializators, Mutators, Consts, DBAdapters

import random
import numpy as np


import pkp
import pkp.empirical_model

import matplotlib.pyplot as plt
plt.style.use(['mystyle'])


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


def scale_parameters(parameters):
    par_array = np.array(parameters)
    par_scaled = parameters_min + par_array * delta_parameters
    par_scaled[0] = np.power(10.0, par_scaled[0])
    par_scaled[1] = par_scaled[1] * 1e6
    return par_scaled


def error(parameters):
    # scale parameters and convert log
    par_scaled = scale_parameters(parameters)
    m = pkp.empirical_model.SFOR(
        parameters=par_scaled)
    m.operating_conditions = operating_conditions
    _, y = m.run(t=t_cpd)
    return np.mean(np.power((y - y_cpd), 2)) * 100.


calc = False

res0 = calc_CPD() if calc else read_CPD()


y_cpd = np.array(res0['fsolid'])
t_cpd = np.array(res0.index) * 1e-3

runner = pkp.PKPRunner('input.yml')
operating_conditions = runner.operating_conditions['run0']

parameters_min = np.array([3., 100, 0.4])  # logA, E MJ/kg, y0
#parameters_max = np.array([9., 200, 0.8])
parameters_max = np.array([10., 300, 0.8])
delta_parameters = parameters_max - parameters_min

parameters_0 = [4, 150, 0.5]

NGEN, NPOP = 50, 40

if __name__ == '__main__':
    genome = G1DList.G1DList(3)
    genome.setParams(rangemin=0, rangemax=1)
    genome.initializator.set(Initializators.G1DListInitializatorReal)
    genome.mutator.set(Mutators.G1DListMutatorRealRange)
    # The evaluator function (objective function)
    genome.evaluator.set(error)
    # Genetic Algorithm Instance
    ga = GSimpleGA.GSimpleGA(genome)
    ga.setMinimax(Consts.minimaxType["minimize"])
    # set the population size
    ga.setPopulationSize(NPOP)
    # set the number of generation
    ga.setGenerations(NGEN)
    # Set the Roulette Wheel selector method, the number of
    # generations and the termination criteria
    ga.selector.set(Selectors.GRouletteWheel)
    ga.terminationCriteria.set(GSimpleGA.ConvergenceCriteria)
    ga.setMutationRate(0.40)
    ga.setCrossoverRate(1.0)

    # parallel processing
    # ga.setMultiProcessing(True)
    # Sets the DB Adapter, the resetDB flag will make the Adapter
    # recreate
    # the database and erase all data every run, you should use
    # this flag
    # just in the first time, after the pyevolve.db was created,
    # you can omit it.
    sqlite_adapter = DBAdapters.DBSQLite(
        identify="cpd", resetDB=True)
    ga.setDBAdapter(sqlite_adapter)
    # Do the evolution, with stats dump, frequency of 20 generations
    ga.evolve(freq_stats=1)
    # Gets the best individual
    best = ga.bestIndividual()
    print(best)

    # plot the results

    best_parameters = scale_parameters(best)

    print('Best population', best_parameters)

    m = pkp.empirical_model.SFOR(parameters=best_parameters)
    m.operating_conditions = operating_conditions
    _, y_sfor = m.run(t=t_cpd)

    fig, ax = plt.subplots()

    ax.plot(t_cpd, y_cpd, label='CPD')
    ax.plot(t_cpd, y_sfor, label='SFOR')

    ax.legend()

    plt.show()
