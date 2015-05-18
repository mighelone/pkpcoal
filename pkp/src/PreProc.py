class PreProcResult(object):
    """ Base class for the various preprocessors """

    def __init__(self, coal):
        self.coal = coal


class ManualQfactor(PreProcResult):

    def __init__(self, coal, qFactor):
        PreProcResult.__init__(self, coal)
        self.qFactor = qFactor
    
