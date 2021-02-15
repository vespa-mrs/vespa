# Python modules

# 3rd party modules
import numpy as np

# Our modules



def create_zero(length=0):
    """A convenience function that creates an SvdOutput instance populated
    with 0-filled arrays of the specified length (defaults to 0).
    """
    return SvdOutput(np.zeros(length), np.zeros(length), np.zeros(length), 
                     np.zeros(length))


class SvdOutput(object):
    """A container class that groups our HLSVD function's output (frequencies, 
    damping factors/decays, amplitudes, phases, and "in model" flags).
    Frequencies are in Hertz, damping factors are in seconds, phases are in
    degrees.

    All are 1D numpy arrays of the same length; the length consistency is
    enforced by the class.

    The first four arrays contain floats and are read-only (both the 
    attributes and the arrays themselves).
    
    The last array (in_model) differs in a few ways. It contains True/False 
    values representing whether or not values from the other arrays at that 
    index are to be used in calculations. Unlike the float arrays, it's 
    editable. Individual items can be toggled, and the attribute can be
    overwritten.

    When instantiating, in_model is optional and will be populated with False
    values if not provided.

    Side note: the 'in model' flags are often toggled directly via the GUI so 
    they are often input, despite being here in the SvdOutput class. But when 
    specified via the 'cursor span' or 'threshold' mechanisms of the GUI, 
    these flags are the calculated output of the SVD algorithm, hence they 
    reside here.

    This class implements __len__(). Calling len() on an instance returns the 
    length of the arrays.

    Python also uses __len__() for casting this class to bool. One can use
    expressions like 'if svd_output:' and they'll return False if the arrays
    are empty and True otherwise.

    This class doesn't support inflate()/deflate().
    """
    def __init__(self, frequencies, damping_factors, amplitudes, phases,
                 in_model=None):
        length = len(frequencies)

        if in_model is None:
            # Caller didn't specify this, so init all entries to False.
            in_model = np.empty(length, bool)
            in_model.fill(False)

        if (length == len(damping_factors)) and \
           (length == len(amplitudes))      and \
           (length == len(phases))          and \
           (length == len(in_model)):
            # All is well
            pass
        else:
            raise ValueError("Arrays must be of equal length")

        self._frequencies = frequencies
        self._damping_factors = damping_factors
        self._amplitudes = amplitudes
        self._phases = phases
        self.in_model = in_model

        # We set these to read-only to emphasize that these are to be replaced
        # as a group or not at all. in_model is an exception; it's editable.
        self._frequencies.flags.writeable = False
        self._damping_factors.flags.writeable = False
        self._amplitudes.flags.writeable = False
        self._phases.flags.writeable = False


    def __len__(self):
        return len(self._frequencies)


    def __str__(self):
        return self.__unicode__()


    def __unicode__(self):
        lines = [ ]
        lines.append("--- SvdOutput (%d items) ---" % len(self))
        lines.append("frequencies: %s" % self._frequencies)
        lines.append("damping_factors: %s" % self._damping_factors)
        lines.append("amplitudes: %s" % self._amplitudes)
        lines.append("phases: %s" % self._phases)
        lines.append("in_model: %s" % self.in_model)

        return '\n'.join(lines)


    @property
    def frequencies(self):
        return self._frequencies

    @property
    def damping_factors(self):
        return self._damping_factors

    @property
    def amplitudes(self):
        return self._amplitudes

    @property
    def phases(self):
        return self._phases


