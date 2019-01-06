"""
Module for running CPD simulations using the Fortran CPD solver.

"""
from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals

import os
import subprocess

import pandas as pd
from autologging import logged
import numpy as np
import platform

from . import cpd
from . import bins


np.seterr(all="ignore")


@logged
class CPD(cpd.CPD):
    """
    Class to run and store results from CPD model. The class set the
    input files required for running CPD using the external fortran
    solver. It is based on the CPD_NLG version of the code.
    http://www.et.byu.edu/~tom/cpd/cpdcodes.html

    It is possible to use the Genetti correlation for estimating NMR
    parameters, or they can directly defined if they are known.
    """

    def _get_basename(self):
        return super(CPD, self)._get_basename()

    def _set_basename(self, value):
        super(CPD, self)._set_basename(value)
        self._io_file = os.path.join(self.path, "input_" + self._basename)
        self._input_file = os.path.join(self.path, self._basename + ".inp")

    basename = property(
        _get_basename, _set_basename, doc="Basename for saving CPD results."
    )

    @property
    def io_file(self):
        """
        Input/Output file for running CPD.
        """
        return self._io_file

    @property
    def input_file(self):
        """
        Input file defined for CPD.
        """
        return self._input_file

    @property
    def solver(self):
        """
        CPD solver path.
        """
        return self._solver

    @solver.setter
    def solver(self, value):
        if value is None:
            cpd_path = os.path.dirname(bins.__file__)
            if platform.system() == "Darwin":
                cpd_file = "cpdnlg.x"
            elif platform.system() == "Linux":
                cpd_file = "cpdnlg"
            elif platform.system() == "Windows":
                cpd_file = "cpdnlg.exe"
            value = os.path.join(cpd_path, cpd_file)
        self._solver = os.path.abspath(value)
        self.__log.debug("Set CPD solver %s", self._solver)

    def _write_input_files(self):
        """
        Write the input files required for the CPD calculation
        """

        def writeline(key):
            f.write("{}           !{}\n".format(getattr(self, key), key))

        def empty_lines(n=1):
            [f.write("\n") for _ in range(n)]

        with open(self.input_file, "w") as f:
            [writeline(key) for key in ["p0", "c0", "sig", "mw", "mdel"]]
            empty_lines(1)
            [writeline(key) for key in ["fcar", "fhyd", "fnit", "foxy"]]
            empty_lines(3)
            [
                writeline(key)
                for key in [
                    "ab",
                    "eb",
                    "ebsig",
                    "ac",
                    "ec",
                    "ag",
                    "eg",
                    "egsig",
                    "Acr",
                    "Ecr",
                ]
            ]
            empty_lines(1)
            [writeline(key) for key in ["arad", "erad", "fstable", "an", "en", "ensig"]]
            empty_lines(2)
            self.__log.debug("pressure in Pa %s", self.pressure)
            self.pressure_atm = self.pressure / 101325.0
            writeline("pressure_atm")
            self.__log.debug("write pressure in atm %s", self.pressure_atm)
            empty_lines(1)
            f.write(
                "{}             "
                "! number of time points"
                " time(ms), temp(K)\n".format(len(self.operating_conditions))
            )
            [
                f.write("{},{}\n".format(point[0] * 1e3, point[1]))
                for point in self.operating_conditions
            ]
            empty_lines(11)
            f.write(
                "{},{},{}"
                "    !dt(s),"
                "print increment,"
                "max dt(s)\n".format(self.dt, self.increment, self.dt_max)
            )
            f.write(
                "{}    "
                "!timax (max residence time, s)\n".format(
                    self.operating_conditions[-1, 0]
                )
            )
            writeline("nmax")

        with open(self.io_file, "w") as f:
            f.write("{}\n".format(self.input_file))
            [
                f.write(os.path.join(self.path, self.basename + "_{}.out\n".format(n)))
                for n in range(1, 5)
            ]

    def run(self, **kwargs):
        """
        Run CPD code.

        Returns
        -------
        results: pandas.Dataframe
            Dataframe containg the results of CPD as a function of the
            residence time.
        """
        self._write_input_files()
        with open(self.io_file, "r") as f_in:
            with open("test.out", "w") as f_out:
                with open("test.err", "w") as f_err:
                    code_run = subprocess.call(
                        [self.solver], stdin=f_in, stdout=f_out, stderr=f_err
                    )
        if code_run:
            raise RuntimeError("Error running CPD with {}".format(self.io_file))
        return self._read_results()

    def _read_results(self):
        try:
            df = pd.concat([self._read_cpd_results(n) for n in range(1, 5)], axis=1)
        except:
            raise IOError("Problems reading CPD results")
        df.reset_index(inplace=True)
        # df.index.rename('t', inplace=True)
        # df.index = df.index * 1e-3
        # self.__log.debug('CPD index max %s', df.index.max())
        df.rename(
            columns={
                "time(ms)": "t",
                "ftar": "tar",
                "fgas": "light_gas",
                "ftot": "volatiles",
                "met": "metaplast",
                "temp": "T",
                "fh20": "H2O",
                "fco2": "CO2",
                "fch4": "CH4",
                "fco": "CO",
                "fsolid": "char",
                "fother": "others",
            },
            inplace=True,
        )
        df["t"] = df["t"] * 1e-3
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

    def _read_cpd_results(self, n):
        fname = os.path.join(self.path, self.basename + "_{}.out".format(n))
        with open(fname, "r") as f:
            header = f.readline()[2:].split()
        # self.__log.debug('file=%s header %s', fname, header)
        df = pd.read_csv(
            fname, index_col=0, delimiter=r"\s+", names=header, comment="c"
        )
        if n == 1:
            return df.iloc[:-1]
        else:
            return df

    def _set_NMR_parameters(self, nmr_parameters=None):
        """
        Calc parameters using Genetti correlation

        Parameters
        ----------
        parameters: dict, default=None
            Manually set parameters using a dictionary
            {'mdel': 0, 'mw': 0, 'p0': 0, 'sig': 0}
        """
        self._set_NMR_parameters_from_correlation(nmr_parameters=nmr_parameters)
