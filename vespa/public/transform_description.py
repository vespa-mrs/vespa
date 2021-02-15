# Python modules


# 3rd party modules

             

class TransformDescription(object):
    """
    Describes the object that's passed to user pulse create transform code.
    
    This class implements __str__() so it's easy to print a nicely formatted
    version of it.
    """
    def __init__(self):
        # The Vespa version as a string, e.g. '0.1.1'
        self.vespa_version = ""
        # User defined parameters without types
        self.parameters = {}
        # Extra global parameters
        self.extra = {}
        # The rf waveform and time axis from the previous transform
        self.previous_rf = None
        self.previous_rf_xaxis = None
        # The rf gradient and time axis from the previous transform
        self.previous_gradient = None
        self.previous_grad_xaxis = None
        
    
    def __str__(self):
        return self.__unicode__()

    def __unicode__(self):
        lines = [ ]
        lines.append("--- Transform Description ---")
        lines.append("Vespa Version: %s" % self.vespa_version)

        for key, value in self.extra.items():
            lines.append("Extra - %s: %s" % (key,str(value)))
        for key, value in self.parameters.items():
            lines.append("Parameter - %s: %s" % (key,str(value)))

        
        # __unicode__() must return a Unicode object. In practice the code
        # above always generates Unicode, but we ensure it here.
        return '\n'.join(lines)





