import sys

from setuptools import setup
from setuptools.command.test import test as TestCommand

package_name = 'pkp'
# Get the version from FoamAna/version.py without importing the package
exec(compile(open(package_name + '/version.py').read(),
             package_name + 'src/version.py', 'exec'))


class PyTest(TestCommand):
    user_options = [('pytest-args=', 'a', "Arguments to pass to py.test")]

    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.pytest_args = []

    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        #import here, cause outside the eggs aren't loaded
        import pytest
        errno = pytest.main(self.pytest_args)
        sys.exit(errno)

config = {
    'author'                 : '',
    'author_email'           : '',
    'description'            : '',
    'license'                : '',
    'version'                : __version__,
    'include_package_data'   : True,
    'package_data'           : {
        '': ['cpdnlg*'],
    },
    'packages'               : [
         "pkp",
         "pkp.src",
         "pkp.bins", 
     ],
    'install_requires'       : [
         'numpy',
         'docopt',
         'scipy',
         'pyevolve',
    ],
   'name'                    : 'pkp',
   'tests_require' : ['pytest'],
   'cmdclass' :  {'test': PyTest},
   'entry_points': {'console_scripts': ['pkp-cli = pkp:main']}
}


setup(**config)
