import sys

if sys.version_info >= (3, 6):
    ImportError = ModuleNotFoundError
else:
    ImportError = ImportError
