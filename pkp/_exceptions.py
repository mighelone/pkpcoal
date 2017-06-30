"""Module for error exception."""
import sys

if sys.version_info >= (3, 6):
    ImportError = ModuleNotFoundError
else:
    ImportError = ImportError


class PKPError(Exception):
    """Base class for PKP exception."""

    pass


class PKPCompositionError(PKPError):
    """Composition errors."""

    pass


class PKPConvertNumber(PKPError):
    """Error converting str to numbers."""

    pass


class PKPModelError(PKPError):
    """Raise error for an unknown model."""

    pass


class PKPKeyError(KeyError):
    """Raise key error."""

    pass


class PKPMethodError(KeyError):
    """Method error."""

    pass


class PKPParametersError(PKPError):
    """Raise a parameters error."""

    pass
