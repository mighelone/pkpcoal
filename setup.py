from distutils.core import setup

setup(
    name='PKP',
    version='2.0.0',
    packages=['pkp', 'pkp.bins'],
    packages_dir={'pkp.bins': 'pkp/bins'},
    package_data={'pkp.bins': ['cpdnlg*', 'COAL.xml', 'Biomass.xml']},
    scripts=['runPKP'],
    url='',
    license='',
    author='Michele Vascellari',
    author_email='Michele.Vascellari@vtc.tu-freiberg.de',
    description='Pyrolysis Kinetic Preporcessor',
    requires=[
        'numpy',
        'matplotlib',
        'scipy',
        'argparse',
        'autologging',
        'ruamel_yaml',
        'pandas',
        'deap',
        'cantera (>=2.2.0)',
        'tabulate',
        'future'
    ]
)
