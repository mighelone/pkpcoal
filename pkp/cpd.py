'''
Module for running CPD simulations.

Example
-------

Initialize the :class:`pkp.cpd.CPD` class using the proximate and
ultimate analysis::

    >>> import pkp.cpd
    >>> ua = {'C': 69,
              'H': 5,
              'O': 24.7,
              'N': 0.8,
              'S': 0.5}

    >>> pa = {'FC': 45.1,
              'VM': 50.6,
              'Ash': 4.3,
              'Moist': 19.0}
    >>> m = pkp.cpd.CPD(proximate_analysis=pa,
                        ultimate_analysis=ua,
                        pressure=101325,
                        name='Pittsburg')

Set the operating conditions::

    >>> op_cond = [[0, 500],
                   [0.001, 1400],
                   [0.01, 1400]]
    >>> m.operating_conditions = op_cond

Set the output directory where store the results::

    >>> m.path = './Results'

Set the parameters for the CPD run::

    >>> m.set_parameters(dt=1e-5, increment=2, dt_max=1e-5,
                         basename='test')

and finally run the model::

    >>> results = m.run()

The results are stored in a ``pandas`` dataframe.

'''
from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals

import os
import subprocess

import pandas as pd
from autologging import logged
import numpy as np

import pkp.detailed_model
import pkp.bins
import platform


CPD_CORRELATION = np.array([[0.0, 0.0, 0.0, 0.0],
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


@logged
class CPD(pkp.detailed_model.DetailedModel):
    '''
    Class to run and store results from CPD model. The class set the
    input files required for running CPD using the external fortran
    solver. It is based on the CPD_NLG version of the code.
    http://www.et.byu.edu/~tom/cpd/cpdcodes.html

    It is possible to use the Genetti correlation for estimating NMR
    parameters, or they can directly defined if they are known.
    '''
    nmr_parameters = ['mdel', 'mw', 'p0', 'sig', 'c0']
    num_parameters = ['dt', 'increment', 'dt_max']

    ab = 2.602e15
    eb = 55400
    ebsig = 1800
    ac = 0.9
    ec = 0
    ag = 3.0e15
    eg = 69000
    egsig = 8100
    Acr = 3.0e15
    Ecr = 65000
    arad = 18.4
    erad = 6000
    fstable = 0.03
    an = 5.5e7
    en = 90000
    ensig = 0
    nmax = 20

    def __init__(self, ultimate_analysis, proximate_analysis,
                 pressure=101325, name='CPD coal'):
        '''
        Parameters
        ----------
        proximate_analysis: dict
            Proximate analysis dict i.e:
            `{'FC': 45.1, 'VM': 50.6, 'Ash': 4.3, 'Moist': 19.0}`
        ultimate_analysis: dict
            Ultimate analysis dictionary i.e:
            `{'C': 80, 'H': 8, 'O': 12, 'N': 0, 'S': 0}`
        pressure: float
            Pressure of pyrolysis process
        name: str, unicode
            Reference name of the modelled coal
        '''
        super(CPD, self).__init__(proximate_analysis=proximate_analysis,
                                  ultimate_analysis=ultimate_analysis,
                                  pressure=pressure,
                                  name=name)

        # check if they are in %
        self.fcar = self.ultimate_analysis['C']
        self.fhyd = self.ultimate_analysis['H']
        self.fnit = self.ultimate_analysis['N']
        self.foxy = self.ultimate_analysis['O']
        self.vm_daf = self.proximate_analysis_daf['VM']

        # set parameters -> this can be changed using
        # self.set_parameters
        self.solver = None
        self._set_NMR_parameters()
        self._set_numerical_parameters()

    def set_parameters(self, **kwargs):
        '''
        Set parameters for CPD calculation.
        If a parameter is not defined it is not changed.
        If None the default value is setted.

        Parameters
        ----------
        nmr_parameters: dict
            NMR parameters dictionary. Keys are:
            ['mdel', 'mw', 'p0', 'sig', 'c0']
            If None values are calculated using Genetti correlation
        dt: float
            Time step
        increment: int
            Number of time step saved
        dt_max: float
            Max. time step
        solver: str
            CPD solver path
        basename: str
            Basename for CPD output files
        '''
        # TODO set a warning for coal outside the parameters
        if 'nmr_parameters' in kwargs:
            self._set_NMR_parameters(
                nmr_parameters=kwargs['nmr_parameters'])

        num_parameters = ('dt', 'increment', 'dt_max')
        if any(p in kwargs for p in num_parameters):
            self._set_numerical_parameters(**kwargs)

        if 'solver' in kwargs:
            self.solver = kwargs['solver']

        if 'basename' in kwargs:
            self.basename = kwargs['basename']

    def get_parameters(self):
        nmr = {p: getattr(self, p)
               for p in self.nmr_parameters}
        par = {p: getattr(self, p)
               for p in self.num_parameters +
               ['basename', 'solver']}
        par['nmr_parameters'] = nmr
        return par

    def _get_basename(self):
        return super(CPD, self)._get_basename()

    def _set_basename(self, value):
        super(CPD, self)._set_basename(value)
        self._io_file = os.path.join(self.path, 'input_' +
                                     self._basename)
        self._input_file = os.path.join(self.path,
                                        self._basename + '.inp')

    basename = property(_get_basename, _set_basename,
                        doc='Basename for saving CPD results.')

    @property
    def io_file(self):
        '''
        Input/Output file for running CPD.
        '''
        return self._io_file

    @property
    def input_file(self):
        '''
        Input file defined for CPD.
        '''
        return self._input_file

    @property
    def solver(self):
        '''
        CPD solver path.
        '''
        return self._solver

    @solver.setter
    def solver(self, value):
        if value is None:
            cpd_path = os.path.dirname(pkp.bins.__file__)
            if platform.system() == 'Darwin':
                cpd_file = 'cpdnlg.x'
            elif platform.system() == 'Linux':
                cpd_file = 'cpdnlg'
            elif platform.system() == 'Windows':
                cpd_file = 'cpdnlg.exe'
            value = os.path.join(cpd_path, cpd_file)
        self._solver = os.path.abspath(value)
        self.__log.debug('Set CPD solver %s', self._solver)

    def _set_numerical_parameters(self, dt=None, increment=None,
                                  dt_max=None, **kwargs):
        '''
        Set numerical parameters
        '''
        if dt is None:
            dt = 1e-5
        if increment is None:
            increment = 1
        if dt_max is None:
            dt_max = 1e-5
        self.dt = dt
        self.increment = increment
        self.dt_max = dt_max

    def _set_NMR_parameters(self, nmr_parameters=None):
        '''
        Calc parameters using Genetti correlation

        Parameters
        ----------
        parameters: dict, default=None
            Manually set parameters using a dictionary
            {'mdel': 0, 'mw': 0, 'p0': 0, 'sig': 0}
        '''
        if nmr_parameters:
            [setattr(self, key, nmr_parameters[key])
             for key in self.nmr_parameters]
        else:
            c = CPD_CORRELATION.copy()
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
                 c[8] * (self.vm_daf * 100) +
                 c[9] * (self.vm_daf * 100)**2)
            [setattr(self, key, Y[i])
             for i, key in enumerate(self.nmr_parameters[:4])]

    def _write_input_files(self):
        '''
        Write the input files required for the CPD calculation
        '''
        def writeline(key):
            f.write('{}           !{}\n'.format(
                getattr(self, key), key))

        def empty_lines(n=1):
            [f.write('\n') for _ in range(n)]

        with open(self.input_file, 'w') as f:
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
            self.__log.debug('pressure in Pa %s', self.pressure)
            self.pressure_atm = self.pressure / 101325.0
            writeline('pressure_atm')
            self.__log.debug('write pressure in atm %s',
                             self.pressure_atm)
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

        with open(self.io_file, 'w') as f:
            f.write('{}\n'.format(self.input_file))
            [f.write(
                os.path.join(self.path, self.basename +
                             '_{}.out\n'.format(n)))
             for n in range(1, 5)]

    def run(self):
        '''
        Run CPD code.

        Returns
        -------
        results: pandas.Dataframe
        Dataframe containg the results of CPD as a function of the
        residence time.
        '''
        self._write_input_files()
        with open(self.io_file, 'r') as f_in:
            with open('test.out', 'w') as f_out:
                with open('test.err', 'w') as f_err:
                    code_run = subprocess.call(
                        [self.solver, ],
                        stdin=f_in,
                        stdout=f_out,
                        stderr=f_err)
        if code_run:
            raise RuntimeError(
                'Error running CPD with {}'.format(self.io_file))
        return self._read_results()

    def _read_results(self):
        def read_file(n):
            fname = os.path.join(self.path, self.basename +
                                 '_{}.out'.format(n))
            with open(fname, 'r') as f:
                header = f.readline()[2:].split()
            # self.__log.debug('file=%s header %s', fname, header)
            return pd.read_csv(
                fname,
                index_col=0,
                delimiter=r'\s+',
                names=header,
                comment='c')
        df = pd.concat([read_file(n)
                        for n in range(1, 5)], axis=1).iloc[:-1]
        df.index.rename('t', inplace=True)
        df.index = df.index * 1e-3
        self.__log.debug('CPD index max %s', df.index.max())

        df.rename(columns={'fsolid': 'char',
                           'ftar': 'tar',
                           'fgas': 'light_gas',
                           'ftot': 'volatiles',
                           'met': 'metaplast',
                           'temp': 'T',
                           'fh20': 'H2O',
                           'fco2': 'CO2',
                           'fch4': 'CH4',
                           'fco': 'CO',
                           'fother': 'other'}, inplace=True)
        # self.__log.debug('Columns %s', df.columns)
        # self.__log.debug('Last row %s', df.iloc[-1])

        # out_csv = os.path.join(self.path, self.basename + '.csv')
        df.to_csv(self._out_csv)
        return df
        # results.plot(y='ftot', kind='line', use_index=True)
        # reset index
        # r_res = r.reset_index()
        # rename column
        # r_res = r_res.rename(columns={' time(ms)': 'time(ms)'})
        # r.index.rename('Time(ms)', inplace=True)
