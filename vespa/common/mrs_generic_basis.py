# Python modules
import copy

# 3rd party modules
import numpy as np
import xml.etree.cElementTree as ElementTree

# Our modules
import vespa.common.util.xml_ as util_xml
import vespa.common.util.time_ as util_time
from vespa.common.constants import Deflate
from vespa.common.util.generic_spectral import create_fid



# The _values used when GenericBasis.reset_to_default() is called.
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


class GenericBasis(object):
    """
    Describes a single basis function using lists of ppm, area, phase and line
    width values. The Nth element of each list corresponds to the Nth element
    in the other lists. This function is constructed by the user, so it may
    contain one or more 'metabolites' as part of its makeup. Its typical use
    is to create an ideal simulated spectrum for use in an algorithm like
    phase 0 corrections via correlation.

    NB. The 'lw' term is a float value indicating FWHM in Hz. It should be
    either a list of N terms (same length as ppm, area, etc.) OR a single term
    (still in a list) that can be broadcast across all lines in the spectrum.

    The 'xxx_lim' attributes are used to store descriptions of valid regions
    for the corresponding 'xxx' attribute.

    """
    XML_VERSION = "1.0.0"

    def __init__(self, attributes=None):

        self.source     = ''         # like 'file' or 'experiment' or 'default'
        self.source_id  = ''         # filename or uuid
        self.check      = []         # bool, shows if a line is active
        self.ppm        = []
        self.area       = []
        self.phase      = []
        self.lw         = []         # in Hz for Gaussian
        self.ppm_lim    = []
        self.area_lim   = []
        self.phase_lim  = []
        self.lw_lim     = []         # in Hz for Gaussian

        self._fid          = None     # cache for calculated array
        self._spectrum_sum = None     # cache for calculated array
        self._spectrum_all = None

        if attributes is not None:
            self.inflate(attributes)
        else:
            self.inflate(self.default_prior())


    ##### Utilities and Properties #####################################

    @property
    def peak_ppm(self):
        """List of max peak locations in PPM. It's read only."""
        return list(self.ppm)

    def get_fid(self, dataset=None):
        if dataset is not None:
            self.update(dataset)
        return self._fid

    def get_spectrum_sum(self, dataset=None):
        if dataset is not None:
            self.update(dataset)
        return self._spectrum_sum

    def get_spectrum_all(self, dataset=None):
        if dataset is not None:
            self.update(dataset)
        return self._spectrum_all


    def __str__(self):

        lines = []
        lines.append("------- {0} Object -------".format(self.__class__.__name__))
        lines.append("Source: "+str(self.source))
        lines.append("Source ID: "+str(self.source_id))
        if self.ppm:
            for i in range(len(self.ppm)):
                line  = "Line (+/-Limits) ="
                line += " "+str(self.check[i])
                line += " "+str(self.ppm[i])
                line += " ("+str(self.ppm_lim[i])+"), "
                line += " "+str(self.area[i])
                line += " ("+str(self.area_lim[i])+"), "
                line += " "+str(self.phase[i])
                line += " ("+str(self.phase_lim[i])+"), "
                line += " "+str(self.lw[i])
                line += " ("+str(self.lw_lim[i])+") "
                lines.append(line)
        return '\n'.join(lines)


    ##### Standard Methods #####################################

    def default_prior(self):
        """ Return copy of _DEFAULT_MODEL dict which can be used by GenericBasis class """
        return copy.deepcopy(_DEFAULT_PRIOR)


    def get_rows(self):
        res = []
        for i in range(len(self.ppm)):

            t = (self.check[i],  self.ppm[i], self.area[i], self.phase[i], self.lw[i],  \
                 self.ppm_lim[i], self.area_lim[i], self.phase_lim[i], self.lw_lim[i])
            res.append(t)

        return res


    def get_row(self, indx):
        if indx >= len(self.ppm) or indx < 0:
            return None
        t = (self.check[indx], self.ppm[indx], self.phase[indx], self.lw[indx], \
             self.ppm_lim[indx], self.phase_lim[indx], self.lw_lim[indx])

        return t


    def reset_to_default(self):
        """Resets model values to the defaults expressed in _DEFAULT_MODEL."""
        for key in _DEFAULT_MODEL:
            setattr(self, key, _DEFAULT_MODEL[key][:])


    def set_values(self, vals):
        # vals are a list of dicts, each dict will have a
        # value for value_ppm, value_area ... lw_lim

        if vals != []:

            self.check     = [item['check'] for item in vals]
            self.ppm       = [item['ppm'] for item in vals]
            self.area      = [item['area'] for item in vals]
            self.phase     = [item['phase'] for item in vals]
            self.lw        = [item['lw'] for item in vals]        # in Hz for Gaussian
            self.ppm_lim   = [item['ppm_lim'] for item in vals]
            self.area_lim  = [item['area_lim'] for item in vals]
            self.phase_lim = [item['phase_lim'] for item in vals]
            self.lw_lim    = [item['lw_lim'] for item in vals]         # in Hz for Gaussian


    def update(self, dataset, zfmult=None):
        """
        Given a dataset, calculates spectrum from the list of spectral lines
        described in the widget and sets the self._fid and self._spectrum
        attributes. Use the 'fid' and 'spectrum' properties to return these
        results locally.

        This method will return a single line of zeros if there are no lines in
        the widget OR there are lines but none have been checked to include.

        """
        zfmult = dataset.zero_fill_multiplier if zfmult is None else zfmult

        dim0 = dataset.raw_dims[0]
        td   = 1.0 / dataset.sw
        off  = round(1.0 * dim0 * dataset.echopeak / 100.0) * td
        xx   = np.arange(dim0) * td
        xxx  = abs(xx - off)
        tmp  = np.zeros(int(dim0 * zfmult), dtype=np.complex64)

        npk = len(self.area)
        fids = []
        specs = []
        if npk > 0 and any(self.check):
            for i in range(npk):
                if self.check[i]:

                    fid, _ = create_fid(self.area[i], self.ppm[i], self.phase[i], dataset, calc_peakppm=False)
                    fids.append(fid)

                    expo = (self.lw[i] * 0.6 * np.pi * xxx) ** 2       # only gaussian here
                    indx = np.where(expo > 50)
                    cnt = np.size(indx)
                    if cnt > 0: expo[indx] = 50
                    lshape = np.exp(-expo)

                    tmp[0:dim0] = fid * lshape
                    tmp[0] *= 0.5

                    specs.append(np.fft.fft(tmp) / len(tmp))

            self._fid          = np.sum(fids, axis=0)
            self._spectrum_sum = np.sum(specs, axis=0)
            self._spectrum_all = np.array(specs)

        else:
            self._fid          = np.zeros((dim0 * zfmult), dtype=np.complex64)
            self._spectrum_sum = np.zeros((dim0 * zfmult), dtype=np.complex64)
            self._spectrum_all = np.zeros((dim0 * zfmult), dtype=np.complex64)

        return True






    ##### Object Specific Methods and Properties #####################################

    def deflate(self, flavor=Deflate.ETREE, tag='', version=''):

        if flavor == Deflate.ETREE:
            tag = "generic_basis" if tag == '' else tag
            version = self.XML_VERSION if version == '' else version

            e = ElementTree.Element(tag, {"version" : version})

            util_xml.TextSubElement(e, "source",    self.source)
            util_xml.TextSubElement(e, "source_id", self.source_id)

            for r in zip(self.check, self.ppm, self.area, self.phase, self.lw,
                         self.ppm_lim, self.area_lim, self.phase_lim, self.lw_lim):
                line_element = ElementTree.SubElement(e, "line")
                util_xml.TextSubElement(line_element, "check",     r[0])
                util_xml.TextSubElement(line_element, "ppm",       r[1])
                util_xml.TextSubElement(line_element, "area",      r[2])
                util_xml.TextSubElement(line_element, "phase",     r[3])
                util_xml.TextSubElement(line_element, "lw",        r[4])
                util_xml.TextSubElement(line_element, "ppm_lim",   r[5])
                util_xml.TextSubElement(line_element, "area_lim",  r[6])
                util_xml.TextSubElement(line_element, "phase_lim", r[7])
                util_xml.TextSubElement(line_element, "lw_lim",    r[8])

            return e

        elif flavor == Deflate.DICTIONARY:
            d = { }
            attributes = ("source", "source_id" )

            for attribute in attributes:
                d[attribute] = getattr(self, attribute)

            line = []
            if self.ppm:
                for i in range(len(self.ppm)):
                    l = {}
                    l["check"]     = self.check[i]
                    l["ppm"]       = self.ppm[i]
                    l["area"]      = self.area[i]
                    l["phase"]     = self.phase[i]
                    l["lw"]        = self.lw[i]
                    l["ppm_lim"]   = self.ppm_lim[i]
                    l["area_lim"]  = self.area_lim[i]
                    l["phase_lim"] = self.phase_lim[i]
                    l["lw_lim"]    = self.lw_lim[i]
                    line.append(l)
            d["line"] = line

            return d


    def inflate(self, source):

        if hasattr(source, "makeelement"):
            # Quacks like an ElementTree.Element
            self.source     = source.findtext("source")
            self.source_id  = source.findtext("source_id")
            self.check      = [util_xml.BOOLEANS[val.text] for val in source.iter("check")]
            self.ppm        = [float(val.text) for val in source.iter("ppm")]
            self.area       = [float(val.text) for val in source.iter("area")]
            self.phase      = [float(val.text) for val in source.iter("phase")]
            self.lw         = [float(val.text) for val in source.iter("lw")]
            self.ppm_lim    = [float(val.text) for val in source.iter("ppm_lim")]
            self.area_lim   = [float(val.text) for val in source.iter("area_lim")]
            self.phase_lim  = [float(val.text) for val in source.iter("phase_lim")]
            self.lw_lim     = [float(val.text) for val in source.iter("lw_lim")]


        elif hasattr(source, "keys"):
            # Quacks like a dict
            self.source     = source["source"]
            self.source_id  = source["source_id"]
            self.check      = []
            self.ppm        = []
            self.area       = []
            self.phase      = []
            self.lw         = []
            self.ppm_lim    = []
            self.area_lim   = []
            self.phase_lim  = []
            self.lw_lim     = []
            if "line" in source:
                for line in source["line"]:
                    self.check.append(line["check"])
                    self.ppm.append(line["ppm"])
                    self.area.append(line["area"])
                    self.phase.append(line["phase"])
                    self.lw.append(line["lw"])
                    self.ppm_lim.append(line["ppm_lim"])
                    self.area_lim.append(line["area_lim"])
                    self.phase_lim.append(line["phase_lim"])
                    self.lw_lim.append(line["lw_lim"])


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

