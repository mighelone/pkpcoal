'''
CPD module
'''
from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals

import numpy as np
import os
import subprocess
import pandas as pd

import pkp.coalnew


cpd_correlation = np.array([[0.0, 0.0, 0.0, 0.0],
                            [421.957, 1301.41, 0.489809, -52.1054],
                            [-8.64692, 16.3879, -0.00981566, 1.63872],
                            [0.0463894, -0.187493, 0.000133046,
                             -0.0107548],
                            [-8.47272, -454.773, 0.155483, -1.23688],
                            [1.18173, 51.7109, -0.0243873, 0.0931937],
                            [1.15366, -10.0720, 0.00705248, -0.165673],
                            [-0.0434024, 0.0760827, 0.000219163,
                             0.00409556],
                            [0.556772, 1.36022, -0.0110498, 0.00926097],
                            [-0.00654575, -0.0313561, 0.000100939,
                             -0.0000826717]])


class CPD(pkp.coalnew.Coal):
    '''
    Class to run and store results from CPD model
    '''

    def __init__(self, ultimate_analysis, proximate_analysis):
        '''
        Parameters
        ----------
        ultimate_analysis: dict
            Ultimate analysis
        proximate_analysis: dict
            Proximate analysis
        '''
        self.ab = 2.602e15
        self.eb = 55400
        self.ebsig = 1800
        self.ac = 0.9
        self.ec = 0
        self.ag = 3.0e15
        self.eg = 69000
        self.egsig = 8100
        self.Acr = 3.0e15
        self.Ecr = 65000
        self.arad = 18.4
        self.erad = 6000
        self.fstable = 0.03
        self.an = 5.5e7
        self.en = 90000
        self.ensig = 0
        self.nmax = 20

        # numerical parameters default
        self.dt = 1e-5
        self.increment = 1
        self.dt_max = 1e-5

        super(CPD, self).__init__(proximate_analysis=proximate_analysis,
                                  ultimate_analysis=ultimate_analysis)

        # check if they are in %
        self.fcar = self.ultimate_analysis['C']
        self.fhyd = self.ultimate_analysis['H']
        self.fnit = self.ultimate_analysis['N']
        self.foxy = self.ultimate_analysis['O']
        self.VMdaf = self.proximate_analysis_daf['VM']

    def set_numerical_parameters(self, dt, increment, dt_max):
        self.dt = dt
        self.increment = increment
        self.dt_max = dt_max

    def calc_model_parameters(self, parameters=None):
        '''
        Calc parameters using Genetti correlation

        Parameters
        ----------
        parameters: dict, default=None
            Manually set parameters using a dictionary
            {'mdel': 0, 'mw': 0, 'p0': 0, 'sig': 0}
        '''
        parameters_keys = ['mdel', 'mw', 'p0', 'sig', 'c0']
        if parameters:
            [setattr(self, key, parameters[key])
             for key in self.parameters_keys]
        else:
            c = cpd_correlation.copy()
            self.c0 = (min(0.36,
                           max(0.118 * self.fcar * 100 - 10.1, 0.0)) +
                       min(0.15,
                           max(0.014 * self.foxy * 100 - 0.175, 0.0)))

            Y = (c[1] + c[2] * (self.fcar * 100.0) +
                 c[3] * (self.fcar * 100)**2 +
                 c[4] * (self.fhyd * 100) +
                 c[5] * (self.fhyd * 100)**2 +
                 c[6] * (self.foxy * 100) +
                 c[7] * (self.foxy * 100)**2 +
                 c[8] * (self.VMdaf * 100) +
                 c[9] * (self.VMdaf * 100)**2)
            [setattr(self, key, Y[i])
             for i, key in enumerate(parameters_keys[:4])]

    def write_input_files(self, fname='CPD_input'):
        def writeline(key):
            f.write('{}           !{}\n'.format(
                    getattr(self, key), key))

        def empty_lines(n=1):
            [f.write('\n') for _ in range(n)]

        cpd_inp_file = self.input_file(fname)
        with open(cpd_inp_file, 'w') as f:
            [writeline(key)
             for key in ['p0', 'c0', 'sig', 'mw', 'mdel']]
            empty_lines(1)
            [writeline(key) for key in ['fcar', 'fhyd', 'fnit', 'foxy']]
            empty_lines(3)
            [writeline(key) for key in ['ab', 'eb', 'ebsig', 'ac', 'ec',
                                        'ag', 'eg', 'egsig', 'Acr',
                                        'Ecr']]
            empty_lines(1)
            [writeline(key) for key in ['arad', 'erad', 'fstable',
                                        'an', 'en', 'ensig']]
            empty_lines(2)
            self.pressure_atm = self.pressure / 101325.0
            writeline('pressure_atm')
            empty_lines(1)
            f.write(
                '{}             '
                '! number of time points'
                ' time(ms), temp(K)\n'.format(
                    len(self.operating_conditions)))
            [f.write('{},{}\n'.format(point[0] * 1e3, point[1]))
             for point in self.operating_conditions]
            empty_lines(11)
            f.write(
                '{},{},{}'
                '    !dt(s),'
                'print increment,'
                'max dt(s)\n'.format(
                    self.dt,
                    self.increment,
                    self.dt_max))
            f.write(
                '{}    '
                '!timax (max residence time, s)\n'.format(
                    self.operating_conditions[-1, 0]))
            writeline('nmax')

        cpd_io = self.io_file(fname)
        with open(cpd_io, 'w') as f:
            f.write('{}\n'.format(cpd_inp_file))
            [f.write(
                os.path.join(self.path, fname +
                             '_{}.out\n'.format(n)))
             for n in range(1, 5)]
        return cpd_inp_file, cpd_io

    def io_file(self, fname):
        '''Return IO file '''
        return os.path.join(self.path, 'input_' + fname)

    def input_file(self, fname):
        '''return cinfiguration file'''
        return os.path.join(self.path, fname + '.inp')

    def run(self, fname='CPD_input', solver='./cpdnlg'):
        solver = os.path.abspath(solver)
        io_file = self.io_file(fname)
        with open(io_file, 'r') as f_in:
            with open('test.out', 'w') as f_out:
                with open('test.err', 'w') as f_err:
                    code_run = subprocess.call(
                        [solver, ],
                        stdin=f_in,
                        stdout=f_out,
                        stderr=f_err)
        if code_run:
            raise RuntimeError(
                'Error running CPD with {}'.format(io_file))

    def read_results(self, fname='CPD_input'):
        def read_file(n):
            return pd.read_csv(
                os.path.join(self.path, fname + '_{}.out'.format(n)),
                delimiter=r'\s+',
                escapechar='c',
                index_col=0)
        return pd.concat([read_file(n) for n in range(1, 5)], axis=1)
        # results.plot(y='ftot', kind='line', use_index=True)
        # reset index
        # r_res = r.reset_index()
        # rename column
        # r_res = r_res.rename(columns={' time(ms)': 'time(ms)'})
        # r.index.rename('Time(ms)', inplace=True)
