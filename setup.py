'''
https://github.com/pypa/sampleproject/blob/master/setup.py
'''

from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='PKP',
    version='2.0.0',
    description='Pyrolysis Kinetic Preporcessor',
    long_description=long_description,
    url='',
    author='Michele Vascellari',
    author_email='Michele.Vascellari@vtc.tu-freiberg.de',
    license='',
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 4 - Beta',

        # Indicate who your project is intended for
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',

        # Pick your license as you wish (should match "license" above)
        'License :: ',

        # Specify the Python versions you support here.
        # In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or
        # both.
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.5',
    ],
    # packages=['pkp', 'pkp.bins'],
    packages=find_packages(exclude=['contrib', 'docs', 'tests']),
    packages_dir={'pkp.bins': 'pkp/bins'},
    package_data={'pkp.bins': ['cpdnlg*', 'COAL.xml', 'Biomass.xml']},
    scripts=['runPKP'],
    install_requires=[
        'numpy',
        'matplotlib',
        'scipy',
        'argparse',
        'autologging',
        'ruamel.yaml',
        'pandas',
        'deap',
        'cantera (>=2.2.0)',
        'tabulate',
        'future'
    ]
)
