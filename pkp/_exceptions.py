class PKPError(Exception):
    """
    Base class for PKP exception
    """
    pass


class PKPCompositionError(PKPError):
    """
    Composition errors
    """
    pass


class PKPConvertNumber(PKPError):
    """
    Error converting str to numbers
    """
    pass


class PKPModelError(PKPError):
    """
    Raise error for an unknown model
    """
    pass


class PKPKeyError(KeyError):
    pass


class PKPMethodError(KeyError):
    pass


class PKPParametersError(PKPError):
    pass
