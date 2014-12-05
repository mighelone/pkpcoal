from setuptools import setup

package_name = 'PKP'
# Get the version from FoamAna/version.py without importing the package
exec(compile(open(package_name + '/version.py').read(),
             package_name + 'src/version.py', 'exec'))

config = {
    'author': '',
    'author_email': '',
    'description': '',
    'license': '',
    'version': __version__,
    'packages': ["PKP", 
                 "PKP.src",
                 "PKP.bins",
                ],
    'install_requires': [
                         'numpy',
                         'scipy',
                         'pyside',
                         'pyevolve',
                        ],
   'name': 'PKP' 
}


setup(**config)
