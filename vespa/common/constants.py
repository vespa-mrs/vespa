# Python modules

import math

# Our modules


# DATABASE_VERSION is an int that describes the database format used by
# this code.
DATABASE_VERSION = 12

DEGREES_TO_RADIANS = math.pi / 180
RADIANS_TO_DEGREES = 180 / math.pi
MINUTES_TO_SECONDS = 60

# We encode numeric lists (also numpy arrays) in a three step process.
# First is XDR (http://en.wikipedia.org/wiki/External_Data_Representation),
# second is zlib to save space, third is base64 to make the output of
# zlib palatable to XML.
#
# NUMERIC_LIST_ENCODING = "xdr zlib base64"
#
# As of 2019-12-2 we are using Numpy 'save()' format to encode numpy arrays
# for saving to xml because it reads/writes much faster than the decode_xdr
# methods in xdrlib. Vespa can still read in the 'xdr' encoded data, but it
# will subsequently save it using 'npy' formatting.
#
# This constant used to be in vespa.common.util.xml_ module
NUMERIC_LIST_ENCODING = "npy zlib base64"



# Used in common.util.db.fetch_experiment_previews() when no standard
# (ie. existing) B0 bin range is returned for a db query for a given B0.
# In this case we assume that this B0 range is unique and set up a left
# and right B0 bin range using the keyword value +/- this constant value
# in MHz (e.g. +/- 1 MHz)
AD_HOC_B0_BIN_DRIFT = 1

class Deflate(object):
    """Constants for the internally-used methods inflate() & deflate()
    that appear on many of our objects.

    These constants are arbitrary and may change.
    However, bool(NONE) is always guaranteed to be False while
    bool(X) == True is guaranteed for all other values.
    """
    NONE = 0
    # ETREE stands for ElementTree which implies serialization to an
    # ElementTree.Element object which is trivial to turn into XML.
    ETREE = 1
    # DICTIONARY is to represent objects as dicts. We don't often
    # deflate to dicts, but we inflate from dicts in the database code and
    # in some GUI code.
    DICTIONARY = 2


class Export(object):
    """Constants for our export/import format."""
    # VERSION is the version of our XML format.
    # When we change our format at some point, it will be helpful to
    # have a version number embedded in our exported files.
    VERSION = "1.0.0"

    ROOT_ELEMENT_NAME = "vespa_export"

# These are used as default attribute values in a number of data objects
# as well as widget start values in Simulation (and maybe Analysis)
DEFAULT_ISOTOPE = "1H"
DEFAULT_B0 = 64.0
DEFAULT_PROTON_CENTER_FREQUENCY = 64.0  # [MHz]
DEFAULT_PROTON_CENTER_PPM = 4.7         # [ppm]
DEFAULT_XNUCLEI_CENTER_PPM = 0.0        # [ppm]
DEFAULT_SWEEP_WIDTH = 2000.0            # [Hz]
DEFAULT_LINEWIDTH = 3.0                 # [Hz]
DEFAULT_SPECTRAL_POINTS = 2048          # integer

RESULTS_SPACE_DIMENSIONS = 4



# Default code for a pulse sequence's binning algorithm. Shows up in
# Simulation's pulse sequence description dialog.
DEFAULT_BINNING_CODE = """import pygamma
import numpy as np

def run(sim_desc):
    area   = pygamma.DoubleVector(0)
    ppm    = pygamma.DoubleVector(0)
    phase  = pygamma.DoubleVector(0)
    field  = sim_desc.field
    nspins = sim_desc.nspins
    tolppm = sim_desc.blend_tolerance_ppm
    tolpha = sim_desc.blend_tolerance_phase
    ppmlow = sim_desc.peak_search_ppm_low
    ppmhi  = sim_desc.peak_search_ppm_high

    specfreq = field
    numberspins = nspins
    freqtol = tolppm
    phasetol = tolpha
    lowppm = ppmlow
    highppm = ppmhi

    freqs    = []
    freqout  = []
    ampout   = []
    phaseout = []
    bincount = 0
    foundone = 0
    
    PI = 3.14159265358979323846
    RAD2DEG = 180.0/PI

    mx_index = sim_desc.mx.Sort(0,-1,0)
    nlines   = sim_desc.mx.size()
    sys      = sim_desc.spin_system
    obs_iso  = sim_desc.observe_isotope

    # Scaling is based on the quantum number of each spin, the qnStates.size
    # method sums 2*qn[i]+1 and then we divide by 2*qn[observe]+1 of the
    # observe spin.  This gets us back to the 'physics' norm where values are
    # =/- 1/2, so then we multiply by 2.0 to get values of 1.0 for each spin 
    # that is being normalized here.
    
    obs_qn  = pygamma.Isotope(obs_iso).qn()
    qnscale = sys.qnStates().size()
    qnscale = qnscale / (2.0 * (2.0*obs_qn+1))

    for ii in range(nlines):
        freqs.append(-1 * sim_desc.mx.Fr(mx_index[ii])/(2.0*PI*specfreq))

    for i in range(nlines):
        freq = freqs[i]
        if (freq > lowppm) and (freq < highppm):
            val = sim_desc.mx.I(mx_index[i])
            valr = val.real()
            vali = val.imag()
            amptemp = np.sqrt(valr**2 + vali**2) / qnscale
            phasetemp = -RAD2DEG * np.angle(valr+vali*1j)

        if bincount == 0:
            freqout.append(freq)
            ampout.append(amptemp)
            phaseout.append(phasetemp)
            bincount += 1
        else:
            for k in range(bincount):
                if (freq >= freqout[k]-freqtol) and (freq <= freqout[k]+freqtol):
                    if (phasetemp >= phaseout[k]-phasetol) and (phasetemp <= phaseout[k]+phasetol):
                        ampsum      =  ampout[k]+amptemp;
                        freqout[k]  = (ampout[k]*freqout[k]  + amptemp*freq)/ampsum;
                        phaseout[k] = (ampout[k]*phaseout[k] + amptemp*phasetemp)/ampsum;
                        ampout[k]  +=  amptemp;
                        foundone = 1; 
            if foundone == 0:
                freqout.append(freq)
                ampout.append(amptemp)
                phaseout.append(phasetemp)
                bincount += 1
            foundone = 0
    
    ppm.resize(bincount)
    area.resize(bincount)
    phase.resize(bincount)

    for i in range(bincount):
        ppm[i]   = freqout[i]
        area[i]  = ampout[i]
        phase[i] = -1.0*phaseout[i]

    return (ppm, area, phase)
"""





# Default code for a pulse sequence's binning algorithm. Shows up in
# Simulation's pulse sequence description dialog.
DEFAULT_BINNING_CODE_V2 = """import pygamma
import numpy as np

def run(sim_desc):
    area   = pygamma.DoubleVector(0)
    ppm    = pygamma.DoubleVector(0)
    phase  = pygamma.DoubleVector(0)
    field  = sim_desc.field             # freq in MHz
    nspins = sim_desc.nspins
    tolppm = sim_desc.blend_tolerance_ppm
    tolpha = sim_desc.blend_tolerance_phase
    ppmlo  = sim_desc.peak_search_ppm_low
    ppmhi  = sim_desc.peak_search_ppm_high

    sys     = sim_desc.spin_system
    obs_iso = sim_desc.observe_isotope
    nlines  = sim_desc.mx.size()

    # need the quantum number of the observe isotope
    tmp = pygamma.Isotope(obs_iso)
    obs_qn = tmp.qn()

    # The original scaling factor fails when a heteronuclear spin system. We
    # have found empirically that scaling needs to be a product of the spin
    # quantum number of the non-observe nuclei minus the observe nuclei
    # quantum number. It also depends on the number of non-observe nuclei.
    # So ... we get the factor calculated below, 'qnscale', with which we
    # modify each line returned from the transition table.
    qnscale = 1.0
    for i in range(nspins):
        qnscale *= 2*sys.qn(i)+1
    qnscale = qnscale / (2.0 * (2.0*obs_qn+1))

    freqs = []
    outf  = []
    outa  = []
    outp  = []
    nbin  = 0
    found = False

    PI = 3.14159265358979323846
    RAD2DEG = 180.0/PI

    # Based on TTable1D::calc_spectra() ---------------------------------------

    indx = sim_desc.mx.Sort(0,-1,0)

    for i in range(nlines):
        freqs.append(-1 * sim_desc.mx.Fr(indx[i])/(2.0*PI*field))

    for i in range(nlines):
        freq = freqs[i]
        if (freq > ppmlo) and (freq < ppmhi):
            val = sim_desc.mx.I(indx[i])
            tmpa = np.sqrt(val.real()**2 + val.imag()**2) / qnscale  #normal
            tmpp = -RAD2DEG * np.angle(val.real()+1j*val.imag())

        if nbin == 0:
            outf.append(freq)
            outa.append(tmpa)
            outp.append(tmpp)
            nbin += 1
        else:
            for k in range(nbin):
                if (freq >= outf[k]-tolppm) and (freq <= outf[k]+tolppm):
                    if (tmpp >= outp[k]-tolpha) and (tmpp <= outp[k]+tolpha):
                        ampsum   =  outa[k]+tmpa
                        outf[k]  = (outa[k]*outf[k] + tmpa*freq)/ampsum
                        outp[k]  = (outa[k]*outp[k] + tmpa*tmpp)/ampsum
                        outa[k] +=  tmpa;
                        found = True
            if not found:
                outf.append(freq)
                outa.append(tmpa)
                outp.append(tmpp)
                nbin += 1
            found = False

    ppm.resize(nbin)
    area.resize(nbin)
    phase.resize(nbin)

    for i in range(nbin):
        ppm[i]   = outf[i]
        area[i]  = outa[i]
        phase[i] = -1.0*outp[i]

    return (ppm, area, phase)
"""

DEFAULT_BINNING_CODE_VERSION1 = """import pygamma

def run(sim_desc):
    area   = pygamma.DoubleVector(0)
    ppm    = pygamma.DoubleVector(0)
    phase  = pygamma.DoubleVector(0)
    field  = sim_desc.field
    nspins = sim_desc.nspins
    tolppm = sim_desc.blend_tolerance_ppm
    tolpha = sim_desc.blend_tolerance_phase
    ppmlow = sim_desc.peak_search_ppm_low
    ppmhi  = sim_desc.peak_search_ppm_high

    bins = sim_desc.mx.calc_spectra(ppm, area, phase, field, nspins, tolppm, tolpha, ppmlow, ppmhi)

    return (ppm, area, phase)
"""

# This is the color of the list items for frozen objects (metabs, pulse
# sequences, pulse projects, etc.). You can pass it to wx.Colour() like so:
#   frozen_color = wx.Colour(*constants.FROZEN_COLOR)
# The only reason we don't create the color directly here is to avoid creating
# a dependency on wx in this file.
FROZEN_COLOR = (236, 236, 255)


class DataTypes(object):
    """Internal representation of data type and routines to convert between
    the three type systems Vespa deals with: Python, XDR and numpy.
    """

    # Vespa has to deal with types in three contexts -- Python, XDR, and numpy.
    # They mostly overlap, but not entirely seamlessly. Here's some bumps to
    # be aware of --
    # - Python ints are only guaranteed to be at least 32 bits. On (some? all?)
    #   64-bit systems, numpy ints default to the numpy.int64 type. That means
    #   that a numpy int can't safely be represented as a Python int. One has
    #   to resort to the rarely-used Python long type to safely represent a
    #   numpy.int64 in Python. The distinction between int and long disappears
    #   in Python 3 (where ints can represent arbitrarily large values, just
    #   like longs in Python 2).
    #
    # - XDR supports 64-bit ints under the name "hyper".
    #
    # - Going in the other direction (Python long ==> numpy.int64 or hyper) is
    #   also an overflow risk since a Python long can contain an arbitrarily
    #   large value. Fortunately we never deal with large int values.
    #
    # - Python floats are "usually implemented using double in C" (per Python
    #   2.7 doc on standard types) which means 64 bits. XDR and numpy both
    #   support higher precision (128-bit) floats, calling them "quadruple"
    #   and "float128" respectively. Python's XDR module doesn't support writing
    #   128-bit floats (which makes sense, since Python's float can't represent
    #   them), so if numpy ever hands us a 128-bit float, we're in trouble.
    #
    # - XDR doesn't understand complex numbers. It's trivial to represent them
    #   into (real, imag) float pairs, but it's nevertheless something that
    #   our code needs to (and does) handle.



    # These constants are arbitrary and may change except that bool(NONE) is
    # always guaranteed to be False and bool() of any other value will always
    # be True.
    NONE        = 0
    BOOL        = 1
    BYTE        = 2
    INT32       = 3
    INT64       = 4
    FLOAT32     = 5
    FLOAT64     = 6
    COMPLEX64   = 7
    COMPLEX128  = 8

    ALL = (BOOL, BYTE, INT32, INT64, FLOAT32, FLOAT64, COMPLEX64, COMPLEX128)

    # The mapping XDR_TYPE_SIZES maps our standard data types to the
    # XDR sizes in bytes (as guaranteed by the XDR standard).
    # bool(NONE) is always guaranteed to be False and bool() of any other
    # value will always be True.
    #
    # XDR knows about lots more types than this, but these are the only ones
    # we care about. We use them to infer the # of elements contained in a
    # glop of XDR data.
    #
    # Note that XDR doesn't define a byte, so it's ambiguous what a "byte"
    # would mean in XDR terms. I defined it below just for completeness sake
    # but it might not be correct.
    #
    # ref: http://tools.ietf.org/html/rfc1014
    # ref: http://tools.ietf.org/html/rfc1832
    # ref: http://tools.ietf.org/html/rfc4506
    XDR_TYPE_SIZES = {
                        NONE            :     0,
                        BOOL            :     4,
                        BYTE            :     4,
                        INT32           :     4,
                        INT64           :     8,
                        FLOAT32         :     4,
                        FLOAT64         :     8,
                        COMPLEX64       :     8,
                        COMPLEX128      :    16,
                     }


    _EXTERNAL_TO_INTERNAL = {
        # Maps external type strings and object to internal values. External
        # strings include the many variations one can find in VASF format.
        # They also include the numpy type strings. Older code (Python & IDL)
        # for reading & writing VASF tended to strip the spaces out of values
        # read from the INI file, so all the VASF type names have to appear
        # in spaceless (e.g. "doublefloat") as well as "spaced" form.

        # These are VASF strings
        'float'                     :    FLOAT32,
        'double'                    :    FLOAT64,
        'doublefloat'               :    FLOAT64,
        'double float'              :    FLOAT64,
        'shortinteger'              :    INT32,
        'short integer'             :    INT32,
        'integer'                   :    INT32,
        'unsignedinteger'           :    INT32,
        'unsigned integer'          :    INT32,
        'integer16bit'              :    INT32,
        'integer 16bit'             :    INT32,
        'integer 16 bit'            :    INT32,
        'integer'                   :    INT32,
        'long'                      :    INT64,
        'unsignedlong'              :    INT64,
        'unsigned long'             :    INT64,
        'complexinteger8bit'        :    COMPLEX64,
        'complex integer8bit'       :    COMPLEX64,
        'complex integer 8bit'      :    COMPLEX64,
        'complex integer 8 bit'     :    COMPLEX64,
        'complexinteger16bit'       :    COMPLEX64,
        'complex integer16bit'      :    COMPLEX64,
        'complex integer 16bit'     :    COMPLEX64,
        'complex integer 16 bit'    :    COMPLEX64,
        'complexfloat'              :    COMPLEX64,
        'complex float'             :    COMPLEX64,
        'complex'                   :    COMPLEX64,
        'complexdouble'             :    COMPLEX128,
        'complex double'            :    COMPLEX128,
        'byte'                      :    BYTE,
        # These are numpy types
        # We need to handle all of the numeric types in numpy.sctypeDict.values().
        # Here's a useful piece of code that prints them:
        # print '\n'.join([str(type_) for type_ in sorted(set(numpy.sctypeDict.values()))])

        # We ignore the numpy types float128 and complex256. Python can't
        # guarantee support for them, since Python floats map to C doubles
        # which are only guaranteed a 64 bit minimum.
        "bool"                      :    BOOL,
        "character"                 :    BYTE,
        "int8"                      :    INT32,
        "uint8"                     :    INT32,
        "int16"                     :    INT32,
        "uint16"                    :    INT32,
        "int32"                     :    INT32,
        "uint32"                    :    INT32,
        "int32"                     :    INT32,
        "uint32"                    :    INT32,
        "int64"                     :    INT64,
        "uint64"                    :    INT64,
        "float32"                   :    FLOAT32,
        "float64"                   :    FLOAT64,
        "complex64"                 :    COMPLEX64,
        "complex128"                :    COMPLEX128,

        # These are Python types
        bool                        :    BOOL,
        int                         :    INT32,
        int                        :    INT64,
        float                       :    FLOAT64,
        complex                     :    COMPLEX128,
    }

    _INTERNAL_TO_NUMPY = {
        # Maps internal types to numpy type strings
        # Valid numpy type names are in numpy.sctypeDict.keys()
        BOOL              :   "bool",
        BYTE              :   "byte",
        INT32             :   "int32",
        INT64             :   "int64",
        FLOAT32           :   "float32",
        FLOAT64           :   "float64",
        COMPLEX64         :   "complex64",
        COMPLEX128        :   "complex128",
    }


    @staticmethod
    def is_complex(the_type):
        return the_type in (DataTypes.COMPLEX64, DataTypes.COMPLEX128)


    @staticmethod
    def any_type_to_internal(the_type):
        if the_type in DataTypes.ALL:
            pass
            # This is already an internal type
        else:
            # If it's a string, lower case it.
            if hasattr(the_type, "lower"):
                the_type = the_type.lower()

            if the_type in DataTypes._EXTERNAL_TO_INTERNAL:
                the_type = DataTypes._EXTERNAL_TO_INTERNAL[the_type]
            else:
                raise ValueError('Unknown type "%s"' % the_type)

        return the_type


    @staticmethod
    def any_type_to_numpy(the_type):
        the_type = DataTypes.any_type_to_internal(the_type)
        return DataTypes._INTERNAL_TO_NUMPY[the_type]


class MrsFileTypes(object):
    # These constants are arbitrary and may change.
    # However bool(NONE) is guaranteed to be False and bool() of anything
    # else will always be True.
    NONE                    =  0
    VASF                    =  1
    VASF_DATA               =  2
    VASF_PARAMETERS         =  3
    VIFF                    =  4
    DICOM_SIEMENS           =  5
    SIEMENS_RDA             =  6
    VARIAN                  =  7
    VIFF_MRS_DATA_RAW       =  8
    DICOM_SIEMENS_FIDSUM    =  9
    PHILIPS_SPAR            = 10




######################     RFPulse constants     ########################

# The next set of classes are basically enumerations.
# They include a logical enumeration e.g. MachineType.GENERAL_ELECTRIC,
# along with their database, 'db', and display, 'display'
# strings (if and when appropriate).
#
# We follow this convention for the enums below:
# bool(EnumClass.NONE) evaluates to False
# bool(EnumClass.ANY_OTHER_VALUE) evaluates to True
#
# Many (but not all!) of the enums inherit from the helper base
# class _TypeEnum.


class _TypeEnum(object):
    """This is a base class to simplify life for the enums below that have
    to implement get_type_for_value(). This implements it for them,
    provided that they set self.ALL properly.
    """

    NONE = { }
    ALL = [ ]

    @classmethod
    def get_type_for_value(self, value, source):
        """Given a value and a source (i.e. 'db' or 'display'), returns
        the type associated with it.

        For instance, this:
            FilterType.get_type_for_value('cosine', 'db')
        would return FilterType.COSINE.
        """
        the_type = self.NONE

        for type_ in self.ALL:
            if type_.get(source) == value:
                the_type = type_
                break
        return the_type


class FilterType(_TypeEnum):
    NONE = {}
    NO_SELECTION = {'display':'None-selected', 'db':'no_selection'}
    COSINE = {'display': 'Cosine', 'db':'cosine'}
    HAMMING = {'display': 'Hamming', 'db':'hamming'}

    ALL = (NO_SELECTION, COSINE, HAMMING)


class FreqProfileType(object):
    NONE = ''
    M_XY = 'm_xy'
    M_X_MINUS_Y = 'm_x_minus_y'
    M_Z = 'm_z'
    M_MINUS_Z = 'm_minus_z'

    @classmethod
    def get_type_for_usage(self, usage):
        if usage == UsageType.EXCITE:
            return self.M_XY
        elif usage == UsageType.SATURATION:
            return self.M_Z
        elif usage == UsageType.INVERSION:
            return self.M_MINUS_Z
        elif usage == UsageType.SPIN_ECHO:
            return self.M_X_MINUS_Y
        else:
            return NONE


class ImportFileFormat(_TypeEnum):
    NONE = {}
    AMPLITUDE_PHASE_DEGREES = {'display':'Amplitude-Phase(Degrees)', 'db':'amplitude_phase_degrees'}

    ALL = (AMPLITUDE_PHASE_DEGREES, )


class MachineType(_TypeEnum):
    NONE = {}
    BRUKER = {'display':'Bruker', 'db':'bruker'}
    FONAR = {'display':'Fonar','db':'fonar'}
    GENERAL_ELECTRIC = {'display':'General Electric','db':'general_electric'}
    HITACHI = {'display':'Hitachi','db':'hitachi'}
    PHILIPS = {'display':'Philips', 'db':'philips'}
    SIEMENS = {'display':'Siemens','db':'siemens'}
    SIMS = {'display':'Sims','db':'sims'}
    TOSHIBA = {'display':'Toshiba','db':'toshiba'}
    VARIAN = {'display':'Varian', 'db':'varian'}

    ALL = (BRUKER, FONAR, GENERAL_ELECTRIC, HITACHI,
           PHILIPS, SIEMENS, SIMS, TOSHIBA, VARIAN)

    @classmethod
    def get_all_display(self):
        vals = [i.get('display') for i in self.ALL]
        return vals


class NewPulseSubtype(_TypeEnum):
    NONE = {}
    GAUSSIAN = {'display':'Gaussian',}
    HYPERBOLIC_SECANT = {'display':'Hyperbolic Secant',}
    HYPERBOLIC_SECANT_TAL = {'display':'Hyperbolic Secant Tal',}
    RANDOM = {'display':'Random',}
    SINC_GAUSSIAN = {'display':'Sinc Gaussian',}
    SLR = {'display':'SLR',}

    ALL = (GAUSSIAN, HYPERBOLIC_SECANT, HYPERBOLIC_SECANT_TAL, RANDOM, SINC_GAUSSIAN, SLR)


class NonCoalescedPhaseSubtype(_TypeEnum):
    NONE = {}
    LINEAR = {'display':'Linear','db':'linear'}
    MAX = {'display':'Max', 'db':'max'}
    MIN = {'display':'Min', 'db':'min'}

    ALL = (LINEAR, MAX, MIN)


class PhaseType(_TypeEnum):
    NONE = {}
    COALESCED = {'display':'Coalesced', 'db':'coalesced'}
    NON_COALESCED = {'display':'Non-coalesced', 'db':'non_coalesced'}

    ALL = (COALESCED, NON_COALESCED)


class PlotType(object):
    NONE = ''
    WAVEFORM = 'waveform'
    FREQUENCY_PROFILE = 'frequency_profile'
    FREQUENCY_PROFILE_ABS = 'frequency_profile_abs'
    EXTENDED_PROFILE = 'extended_profile'
    EXTENDED_PROFILE_ABS = 'extended_profile_abs'
    CONTOUR = 'contour'


class PulseConvention(_TypeEnum):
    NONE = {}
    HALF_HEIGHT = {'display': "FW at Half Height", 'db':'half_height'}
    MINIMUM = {'display': "FW at Minimum", 'db':'minimum'}
    MAXIMUM = {'display': "FW at Maximum", 'db': 'maximum'}

    ALL = (HALF_HEIGHT, MINIMUM, MAXIMUM)


class ThirdPartyExportFormat(object):
    """Constants for 3rd party export for RFPulse."""
    NONE = 0

    # Format can be IDEA, Vision or generic ASCII.
    IDEA = 1
    VISION = 2
    ASCII_MAGN_PHASE = 3
    ANNOTATED_ASCII_MAGN_PHASE = 4
    IDEA_C_HEADER = 5

    # The exported field can be one of B1, B2, or G2.
    FIELD_B1 = 3
    FIELD_B2 = 4
    FIELD_G2 = 5



class SLRFilterType(_TypeEnum):
    NONE = {}
    REMEZ = {'display':'SLR (Remez)', 'db':'remez'}
    LEAST_SQUARES = {'display':'Least Squares', 'db':'least_squares'}

    ALL = (REMEZ, LEAST_SQUARES)


class StepSizeModification(_TypeEnum):
    NONE = {}
    AVERAGE = {'display':'Average', 'db':'average'}
    FIXED = {'display':'Fixed', 'db':'fixed'}
    COMPUTE = {'display':'Compute', 'db':'compute'}

    ALL = (AVERAGE, FIXED, COMPUTE)


class TransformationType(_TypeEnum):
    NONE = {}
    BASIC_INFO = {'display':'Basic Info'} # Not for db
    CREATE_GAUSSIAN = {'display':'Create Gaussian', 'db':'create_gaussian'}
    CREATE_RANDOMIZED = {'display':'Create Randomized', 'db':'create_randomized'}
    CREATE_HYPERBOLIC_SECANT = {'display':'Create Hyperbolic Secant', 'db':'create_hyperbolic_secant'}
    CREATE_SINC_GAUSSIAN = {'display':'Create Sinc Gaussian', 'db':'create_sinc_gaussian'}
    CREATE_SLR = {'display':'Create SLR', 'db':'create_slr'}
    CREATE_IMPORT = {'display':'Create Import', 'db':'create_import'}
    GRADIENT_REFOCUSING = {'display':'Gradient Refocusing','db':'gradient_refocusing'}
    INTERPOLATE_RESCALE = {'display':'Interpolate and Rescale', 'db':'interpolate_rescale'}
    OCN = {'display':'Optimal Control, Non-Selective', 'db':'optimal_control_nonselective'}
    PULSE_CONCATENATION = {'display':'Pulse Concatenation', 'db':'pulse_concatenation'}
    REMAPPING = {'display':'Remapping', 'db':'remapping'}
    ROOT_REFLECTION = {'display':'Root Reflection', 'db':'root_reflection'}
    SAMPLE_EFFECTS = {'display':'Sample Effects', 'db':'sample_effects'}

    CREATE_TYPES = (CREATE_GAUSSIAN, CREATE_RANDOMIZED,
                    CREATE_HYPERBOLIC_SECANT, CREATE_SINC_GAUSSIAN,
                    CREATE_SLR, CREATE_IMPORT, )

    # "ALL" contains the create types plus everything else...
    ALL = tuple( list(CREATE_TYPES) + [BASIC_INFO, GRADIENT_REFOCUSING,
                                      INTERPOLATE_RESCALE, OCN,
                                      PULSE_CONCATENATION, REMAPPING,
                                      ROOT_REFLECTION, SAMPLE_EFFECTS,] )


class UsageType(_TypeEnum):
    NONE = {}
    EXCITE = {'display': 'Excite', 'db':'excite'}
    SATURATION = {'display':'Saturation', 'db':'saturation'}
    INVERSION = {'display':'Inversion', 'db':'inversion'}
    SPIN_ECHO = {'display':'Spin Echo', 'db':'spin_echo'}

    ALL = (EXCITE, INVERSION, SPIN_ECHO, SATURATION)


######################     Pulse constants     ########################

DEFAULT_CALC_RESOLUTION = 5000
DEFAULT_PULSE_BANDWIDTH_TYPE = 'half_height'
DEFAULT_GYROMAGNETIC_NUCLEI = '1H'
DEFAULT_BLOCH_RANGE_VALUE = 4.0
DEFAULT_BLOCH_RANGE_UNITS = 'cm'
DEFAULT_BLOCH_OFFSET_VALUE = 0.0

# sourced gamma vals from website http://bio.groups.et.byu.net/LarmourFreqCal.phtml
# web site refs are :
# -M A Bernstein, K F King and X J Zhou (2004). Handbook of MRI Pulse Sequences. 
#   San Diego: Elsevier Academic Press. p. 960. ISBN 0-1209-2861-2.
# -R C Weast, M J Astle, ed (1982). Handbook of Chemistry and Physics. Boca 
#   Raton: CRC Press. p. E66. ISBN 0-8493-0463-6.

GAMMA_VALUES = {'1H'   : 42.576,       # MHz/T and Hz/gauss (I think)
                '13C'  : 10.705,       # NB. spin  1/2
                '17O'  :  5.7716,      # NB. spin -5/2
                '19F'  : 40.055,
                '23Na' : 11.262,
                '31P'  : 17.235,
                '129Xe': 11.777,}      # NB. spin  1/2

