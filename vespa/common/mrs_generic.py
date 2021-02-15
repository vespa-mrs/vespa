"""
BJS - Sept 8, 2020. This was a start at amalgamating a number of files
  and classes into one module including:
  
  util.generic_spectral
  common.mrs_generic_basis
  analysis.mrs_user_prior
  analysis.mrs_user_prior_spectrum
  
  and maybe more.
  
  Sort of ran out of steam. And vision.  And time.  Sigh.


"""




# Python modules

# 3rd party modules
import numpy as np
import xml.etree.cElementTree as ElementTree

# Our modules
import vespa.common.util.xml_ as util_xml
import vespa.common.util.time_ as util_time
from vespa.common.util.math_ import safe_exp
from vespa.common.constants import Deflate


def calculate_peakppm(fid, sw, freq, zfmult=1, zfnew=4, resppm=4.7):

    acqdim0 = len(fid)
    dim0 = acqdim0 * zfmult * zfnew     # increase zerofill for better accuracy
    hpp = sw / dim0
    apod = safe_exp(-(2.0 * np.pi * 0.6 * (np.arange(acqdim0) / sw)) ** 2, 0)   # 2Hz Gauss
    indx = np.argmax(np.fft.fft(np.pad(fid * apod, (0,dim0-acqdim0))))          # max peak of spectrum
    pkppm = ((dim0 / 2) - indx) * (hpp / freq) + resppm

    return pkppm


def create_fid(area, ppm, phase, sw, frequency, dim0,
                  ta=[1e6,], tb=[1e6,], resppm=4.7, ideal=True):
    """
    Create an 'ideal' FID via a sum of individual FIDs from lists of spectral
    peak descriptors (ie. area, ppm, phase, (opt) ta, (opt) tb values). No
    global Phase 0 or Phase 1 are applied here. The use of Ta and Tb terms here
    is generally to allow static values to only be calculated once when the
    basis set is generated.

    Args:
        area (ndarray, floats): area under each resonance peak
        ppm (ndarray, floats): peak center for each resonance peak
        phase (ndarray, floats): phase 0 for each resonance peak
        sw (float): sweep width in Hz
        frequency (float): scanner central frequency in MHz
        dim0 (int): number of spectral points in FID
        ta (ndarray, floats): default=1e6, lorentzian decay term, can be a
            single value in the array that is broadcast for all lines.
        tb (ndarray, floats): default=1e6, gaussian decay term, can be a
            single value in the array that is broadcast for all lines.
        resppm (float): default=4.7, ppm value of the center point of spectrum.
        ideal (bool): default=True, flag to create FID with (False) or without
            (True) ta and tb factors included. Assumes line shape added later.

    Usage:

      fid, fids = calculate_fid(area, ppm, phase, sw, frequency, dim0, resppm=4.69, indiv=True)

    """
    npk  = len(area)
    hpp  = sw/dim0
    arr1 = np.zeros(dim0, np.float) + 1.0
    t = (np.arange(npk * dim0, dtype='float') % dim0) * (1.0 / sw)
    t.shape = npk, dim0

    # ppm -> hz -> rads
    freq = ((dim0/2) - (frequency*(ppm-resppm)/hpp)) * hpp * 2.0 * np.pi
    # deg -> rads -> complex
    phase = np.exp(1j * phase * np.pi / 180.0)

    # set up lineshape
    if ideal:
        lshape = 1.0
    else:
        if len(ta) != npk: ta = np.repeat(ta[0], npk)
        if len(tb) != npk: tb = np.repeat(tb[0], npk)
        expo = (t / np.outer(ta, arr1)) + (t / np.outer(tb, arr1)) ** 2
        indx = np.where(expo > 50)
        if len(indx) > 0: expo[indx] = 50
        lshape = safe_exp(-expo)

    am = np.outer(np.array(area), arr1)
    fr = np.exp( 1j * np.outer(freq, arr1) * t)
    ph = np.outer(phase, arr1)

    fid = np.sum(am * fr * ph * lshape, axis=0)

    return fid


def create_spectrum(area, ppm, phase, sw, frequency, dim0, zfmult=1,
                    echopeak=0.0, ta=[1e6,], tb=[1e6,], resppm=4.7,
                    ideal=True, acq=False, zforce=None, peakppm=False):
    """
    Create a spectrum via a sum of spectral peaks. Each peak has an area, ppm,
    and phase used to create an 'ideal' FID. The set of summed FIDs is modified
    by a Voigt lineshape described by Ta and Tb parameters. No global Phase 0
    or Phase 1 is applied in here.

    Args:
        area (ndarray, floats): area under each resonance peak
        ppm (ndarray, floats): peak center for each resonance peak
        phase (ndarray, floats): phase 0 for each resonance peak
        sw (float): sweep width in Hz
        frequency (float): scanner central frequency in MHz
        dim0 (int): number of spectral points in FID
        zfmult (int): default=1, zero fill factor
        echopeak (float): default=0.0, value 0.0-1.0, 0.0 is FID, 0.5 is a
            half-echo.
        ta (ndarray, floats): default=1e6, lorentzian decay term, can be a
            single value in the array that is broadcast for all lines.
        tb (ndarray, floats): default=1e6, gaussian decay term, can be a
            single value in the array that is broadcast for all lines.
        resppm (float): default=4.7, ppm value of the center point of spectrum.
        ideal (bool): default=True, flag to create FID with (False) or without
            (True) ta and tb factors included. Assumes line shape added later.
        acq (bool, optional): default=False, determines if returned spectrum is
            sized based on current zero filling value, or on acquisition dims.
        zforce (int, optional): default=None, can be used to force returned
            spectrum to be zero filled to an arbitrary size based on integer
            value of the zforce parameter.

    Usage:

      spectrum = create_spectrum(areas, ppms, phases, dataset, ta=0.07, tb=0.07)

    """
    fid = create_fid(area, ppm, phase, sw, frequency, dim0,
                       ta=ta, tb=tb, resppm=resppm, ideal=ideal)

    pkppm = calculate_peakppm(fid, sw, frequency, zfmult=zfmult, zfnew=4, resppm=resppm) if peakppm else None

    zfmult = 1 if acq else zfmult
    zfmult = zforce if zforce is not None else zfmult

    spectrum = np.zeros(int(round(dim0 * zfmult)), np.complex)
    spectrum[0:dim0] = fid
    spectrum[0] *= 0.5
    spectrum = np.fft.fft(spectrum) / len(spectrum)

    return spectrum, pkppm



class GenericBasis(object):
    """
    This is the fundamental object that represents the data being
    manipulated in the Analysis program.

    """
    XML_VERSION = "1.0.0"

    def __init__(self, attributes=None):

        self.ppm    = None
        self.area   = None
        self.phase  = None
        self.ta     = None         # in Hz for Gaussian
        self.tb     = None

        if attributes is not None:
            self.inflate(attributes)
        else:
            self.inflate(self.default_prior())


    ##### Standard Methods and Properties #####################################

    def __str__(self):

        lines = []
        lines.append("------- {0} Object -------".format(self.__class__.__name__))
        if self.ppm:
            for i in range(len(self.ppm)):
                line = "Line = "+str(self.ppm[i])+" "+str(self.area[i])+" "+str(self.phase[i])+" "+str(self.ta[i])+" "+str(self.tb[i])
                lines.append(line)
        return '\n'.join(lines)


    def deflate(self, flavor=Deflate.ETREE):

        if flavor == Deflate.ETREE:
            e = ElementTree.Element("generic_basis", {"version" : self.XML_VERSION})
            if self.ppm is not None:   e.append(util_xml.numpy_array_to_element(self.ppm, 'ppm'))
            if self.area is not None:  e.append(util_xml.numpy_array_to_element(self.area, 'area'))
            if self.phase is not None: e.append(util_xml.numpy_array_to_element(self.phase, 'phase'))
            if self.ta is not None:    e.append(util_xml.numpy_array_to_element(self.ta, 'ta'))
            if self.tb is not None:    e.append(util_xml.numpy_array_to_element(self.tb, 'tb'))
            return e

        elif flavor == Deflate.DICTIONARY:
            return self.__dict__.copy()


    def inflate(self, source):

        if hasattr(source, "makeelement"):
            # Quacks like an ElementTree.Element
            for attribute in ("ppm",  "area", "phase", "ta", "tb"):
                item = source.find(attribute)
                if item is not None:
                    setattr(self, attribute, util_xml.element_to_numpy_array(item))

        elif hasattr(source, "keys"):
            # Quacks like a dict
            for key in list(source.keys()):
                setattr(self, key, source[key])


    ##### Object Specific Methods and Properties #####################################

    def default_prior(self):
        """
        Set by user to pre-populate the object if no attribute keyword
        is sent in during initialization.

        """
        d = {'ppm'   : [1.0,],
             'area'  : [1.0,],
             'phase' : [0.0,],
             'ta'    : [0.07,],
             'tb'    : [0.07,]  }

        return d


    def get_rows(self):
        ta = np.repeat(self.ta[0], len(self.ppm)) if len(self.ta) != len(self.ppm) else self.ta
        tb = np.repeat(self.tb[0], len(self.ppm)) if len(self.tb) != len(self.ppm) else self.tb
        return [(r[0],r[1],r[2],r[3],r[4]) for r in zip(self.ppm, self.area, self.phase, ta, tb)]


    def get_row(self, indx):
        if indx >= len(self.values_ppm) or indx < 0:
            return None
        ta = self.ta[0] if len(self.ta) != len(self.ppm) else self.ta[indx]
        tb = self.tb[0] if len(self.tb) != len(self.ppm) else self.tb[indx]
        return (self.ppm[indx], self.area[indx], self.phase[indx], ta, tb)


    def set_values(self, vals):
        # vals are a list of dicts, each dict will have a
        # value for value_ppm, value_area ... limit_lwhz

        self.ppm   = [item[0] for item in vals]
        self.area  = [item[1] for item in vals]
        self.phase = [item[2] for item in vals]
        self.ta    = [item[3] for item in vals]
        self.tb    = [item[4] for item in vals]



class GenericSpectrum(object):

    def __init__(self):

        self.bob = 10
        self.peak_ppm = []

    @property
    def peak_ppm(self):
        """List of max peak locations in PPM. It's read only."""
        return list(self.values_ppm)




#--------------------------------------------------------------------
# test code

def _test():

    test = GenericBasis()

    class_name = test.__class__.__name__
    filename = "_test_output_"+class_name+".xml"
    element = test.deflate()
    root = ElementTree.Element("_test_"+class_name, { "version" : "1.0.0" })
    util_xml.TextSubElement(root, "timestamp", util_time.now().isoformat())
    root.append(element)
    tree = ElementTree.ElementTree(root)
    tree.write(filename, "utf-8")

    tom = 10


if __name__ == '__main__':
    _test()

