# Python modules
import copy

# 3rd party modules
import numpy as np

# Our modules
from vespa.common.util.generic_spectral import create_spectrum
from vespa.common.mrs_generic_basis import GenericBasis
from vespa.common.constants import Deflate


# The _values used when UserPriorSpectrum.reset_to_default() is called.

_DEFAULT_MODEL = { 
   "check"     : [True, True, True],
   "ppm"       : [2.01, 3.01, 3.22],
   "area"      : [0.51, 0.42, 0.40],
   "phase"     : [0.00, 0.00, 0.00],
   "lw"        : [5.0,  5,    5],
   "ppm_lim"   : [0.0,  0.0,  0.0],
   "area_lim"  : [0.0,  0.0,  0.0],
   "phase_lim" : [0.0,  0.0,  0.0],
   "lw_lim"    : [0.0,  0.0,  0.0]    }

# _DEFAULT_PRIOR is a copy and reformat of the default model. It's calculated 
# once by the function below and then saved in a dict. A copy of that dict is 
# returned when someone calls UserPriorSpectrum.default_prior()

def _calculate_default_prior():
    lines = []
    for i in range(len(_DEFAULT_MODEL["check"])):
        d = {}
        for key in list(_DEFAULT_MODEL.keys()):
            d[key] = _DEFAULT_MODEL[key][i]
        lines.append(d)

    default = {}
    default["source"] = "default"
    default["source_id"] = "default"
    default["line"] = lines

    return default

_DEFAULT_PRIOR = _calculate_default_prior()



class UserPriorSpectrum(GenericBasis):
    """
    Contains a 'simple' simulated spectrum described by the user as a group of
    area, ppm and phase terms, with some amount of lineshape applied. It is 
    used in various parts of Analysis when an ideal spectral description is 
    needed for an algorithm, such as Ph0 correction via correlation.
    
    """
    XML_VERSION = "1.0.0"

    def __init__(self, attributes=None):

        super().__init__(attributes)

        # _summed caches the value returned by the summed property. We could
        # expose this as an attribute, but wrapping this in a property 
        # emphasizes that it is read-only. It's updated by a call to
        # create_spectrum().
        self._summed = None


    @property
    def summed(self):
        """Returns a numpy array containing the sum of the lines checked 
        in the widget.
        """
        return self._summed


    ##### Standard Methods and Properties #####################################

    def deflate(self, flavor=Deflate.ETREE):
        if flavor == Deflate.ETREE:
            # Make my base class do its deflate work
            return super().deflate(flavor, tag="user_prior_spectrum", version=self.XML_VERSION)

            # e = super().deflate(flavor)
            # # Alter the tag name & XML version info
            # e.tag = "user_prior_spectrum"
            # e.set("version", self.XML_VERSION)
            # return e

        elif flavor == Deflate.DICTIONARY:
            # My base class does all the work
            return super().deflate(flavor)


    def inflate(self, source):
        # Make my base class do its inflate work
        super().inflate(source)



    ##### Object Specific Methods and Properties ###########################

    def calculate_spectrum(self, dataset, acq=False, zforce=None):
        """
        Given a dataset, calculates the list of spectral lines that are
        described in the widget and returns them. It also updates the value
        returned by this object's summed property.

        This method will return a single line of zeros if there are no lines in
        the widget OR there are lines but none have been checked to include.

        """
        all_spectra = []
        checks = self.check
        areas  = self.area       
        ppms   = self.ppm        
        phases = self.phase      
        lwhz   = self.lw  
        if areas and any(checks):
            for i in range(len(areas)):
                if checks[i]:
                    # have to call one at a time so each spectral peak
                    # can have its own different line width if desired
                    apod = 1.0 / (lwhz[i] * 0.6 * np.pi)
                    spectrum = create_spectrum(areas[i],
                                               ppms[i],
                                               phases[i],
                                               dataset,
                                               10000.0,
                                               apod,
                                               acq=acq,
                                               zforce=zforce)
                    all_spectra.append(spectrum)
        else:
            zfmult   = dataset.zero_fill_multiplier
            acqdim0  = dataset.raw_dims[0]
            all_spectra = [np.zeros((acqdim0 * zfmult), complex)]

        # Update the sum.
        self._summed = np.sum(all_spectra, axis=0)
            
        return all_spectra


    def reset_to_default(self):
        """Resets model values to the defaults expressed in _DEFAULT_MODEL."""
        for key in _DEFAULT_MODEL:
            setattr(self, key, _DEFAULT_MODEL[key][:])
        

    def default_prior(self):
        """ Return copy of _DEFAULT_MODEL dict which can be used by GenericBasis class """
        return copy.deepcopy(_DEFAULT_PRIOR)

