from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals

from deap import base
from deap import creator
from deap import tools

from deap.benchmarks import rosenbrock, bohachevsky, schaffer

import random

IND_SIZE = 2
NPOP = 100

CXPB, MUTPB, NGEN = 0.2, 0.2, 1000

creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
creator.create("Individual", list, fitness=creator.FitnessMin)

toolbox = base.Toolbox()
# Attribute generator

toolbox.register("attr_float", random.randrange, -10, 10)
# Structure initializers
toolbox.register("individual", tools.initRepeat,
                 creator.Individual, toolbox.attr_float, n=IND_SIZE)
toolbox.register("population", tools.initRepeat, list,
                 toolbox.individual)


# ind1 = toolbox.individual()
# ind1.fitness.values = rosenbrock(ind1)
# print(ind1)

# print(rosenbrock(ind1))

# mutant
# mutant = toolbox.clone(ind1)
# ind2, = tools.mutGaussian(mutant, mu=0.0, sigma=0.2, indpb=0.2)
# del mutant.fitness.values

# crossover
# child1, child2 = [toolbox.clone(ind) for ind in (ind1, ind2)]
# tools.cxBlend(child1, child2, 0.5)

# del child1.fitness.values
# del child2.fitness.values

# selection
# selected = tools.selBest([child1, child2], 1)
# print(child1 in selected)
# print(child2 in selected)


toolbox.register('mate', tools.cxTwoPoint)
toolbox.register('mutate', tools.mutGaussian, mu=0, sigma=1, indpb=0.2)
toolbox.register('select', tools.selTournament, tournsize=3)
toolbox.register('evaluate', rosenbrock)
#toolbox.register('evaluate', bohachevsky)
#toolbox.register('evaluate', schaffer)

if __name__ == '__main__':
    pop = toolbox.population(n=NPOP)
    fitnesses = list(map(toolbox.evaluate, pop))
    for ind, fit in zip(pop, fitnesses):
        ind.fitness.values = fit

    for g in range(NGEN):
        # Select the next generation individuals
        offspring = toolbox.select(pop, len(pop))
        # Clone the selected individuals
        offspring = list(map(toolbox.clone, offspring))

        # Apply crossover and mutation on the offspring
        for child1, child2 in zip(offspring[::2], offspring[1::2]):
            if random.random() < CXPB:
                toolbox.mate(child1, child2)
                del child1.fitness.values
                del child2.fitness.values

        for mutant in offspring:
            if random.random() < MUTPB:
                toolbox.mutate(mutant)
                del mutant.fitness.values

        # Evaluate the individuals with an invalid fitness
        invalid_ind = [ind for ind in offspring
                       if not ind.fitness.valid]
        fitnesses = map(toolbox.evaluate, invalid_ind)
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit

        pop[:] = offspring

        # Gather all the fitnesses in one list and print the stats
        fits = [ind.fitness.values[0] for ind in pop]

        length = len(pop)
        mean = sum(fits) / length
        sum2 = sum(x * x for x in fits)
        std = abs(sum2 / length - mean**2)**0.5

        # print("  Min %s" % min(fits))
        # print("  Max %s" % max(fits))
        # print("  Avg %s" % mean)
        # print("  Std %s" % std)
        # print('-- Generation %i --' % g)
        print(
            'Gen: {} - Min:{:5.4f} Max:{:5.4f}'
            ' Avg:{:5.4f} Std:{:5.4f}'.format(g, min(fits), max(fits),
                                              mean, std))
