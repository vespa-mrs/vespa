"""A home for constants."""

# Python modules
import collections

# 3rd party modules

# Our modules


########           Table of Contents         ##########
#
# This file is divided into several sections. First are some some
# miscellaneous constants.
#
# Next are some science-y constants, like min/max apodization values,
# spatial processing options, etc.
#
# Last are Fit (general) and Voigt (specific) options which get their 
# own section because there's so many of them.


##################   Miscellaneous constants

try:
    import hlsvdpro
except ImportError:
    HLSVDPRO_AVAILABLE = False
else:
    HLSVDPRO_AVAILABLE = True
    # Some 64-bit versions of Linux have issues with the hlsvdpro libraries we
    # use for water filtering in Analysis Spectral tab. For those users we
    # provide a native python based solution. Just uncomment the line below to
    # force the native python method to be used.
    HLSVDPRO_AVAILABLE = False

# HLSVD_METHOD = 'hlsvdpro'       # to fix Linux issue comment out this line
# # HLSVD_METHOD = 'hlsvdpropy'    #  and uncomment this line

# PyWavelets is an optional dependency. There's a couple of places in the
# code that need to know whether or not it's installed. Here we make that
# determination so that any code which cares can just check this constant.
try:
    import pywt
except ImportError:
    PYWAVELETS_AVAILABLE = False
else:
    PYWAVELETS_AVAILABLE = True


# class UserButton(object):
#     """Choices available for the User Defined Function button on the Spectral
#     tab bottom. These are used as the string value for the button.
#     """
#
#     # These constants are arbitrary and may change.
#     AUTOPHASE  = 'Do Automatic Phasing'
#     AREAOUTPUT = 'Output Area Value'



##################   Science-y constants


# class AmplitudeMultiplier(object):
#     """ Amplitude multiplier constants """
#     MIN = 0
#     MAX = 1e12
#     MULTIPLIER = 1.1


class Apodization(object):
    """ Apodization constants """
    # MIN =   0
    # MAX = 100
    # INCREMENT = 0.5

    # These constants are arbitrary and may change.
    # However bool(NONE) is guaranteed to be False while
    # bool(X) == True is guaranteed for all other values.
    NONE       = ''
    GAUSSIAN   = 'gaussian'
    LORENTZIAN = 'lorentzian'

    # Items for the spectral processing options dropdown
    choices = collections.OrderedDict(((NONE ,       "None"),
                                       (GAUSSIAN ,   "Gaussian"),
                                       (LORENTZIAN , "Lorentzian"),
                                      ))

# class DcOffset(object):
#     """ DC offset constants """
#     MIN = -1e5
#     MAX =  1e5
#     INCREMENT = 0.5
#
#
# class FrequencyShift(object):
#     """ Frequency shift constants """
#     MIN = -1e4
#     MAX =  1e4
#     INCREMENT = 0.5
#
#
# class Phase_0(object):
#     """ First order phase constants """
#     MIN = 0
#     MAX = 360
#     INCREMENT = 0.01
#
#
# class Phase_1(object):
#     """ First order phase constants """
#     MIN = -1e4
#     MAX =  1e4
#     INCREMENT = 10.0
#
#     MIN_PIVOT = -1000
#     MAX_PIVOT =  1000
#     INCREMENT_PIVOT = 0.5


class SpatialFilter(object):
    """ Spatial filter constants """
    # These constants are arbitrary and may change.
    # However bool(NONE) is guaranteed to be False while
    # bool(X) == True is guaranteed for all other values.
    NONE = ''
    HAMMING = 'hamming'
    EXPONENTIAL = 'exponential'
    GAUSSIAN = 'gaussian'

    # Items for the spatial processing options dropdown
    choices = collections.OrderedDict(((NONE ,        "None"),
                                       (HAMMING ,     "Hamming"),
                                       (EXPONENTIAL , "Exponential"),
                                       (GAUSSIAN ,    "Gaussian"),
                                      ))

class SpatialTranspose(object):
    """ Spatial transposition constants """
    # These constants are arbitrary and may change.
    # However bool(NONE) is guaranteed to be False while
    # bool(X) == True is guaranteed for all other values.
    NONE = ''
    TRANSPOSE_XY  = 'x_to_y'
    TRANSPOSE_XZ  = 'x_to_z'
    TRANSPOSE_YZ  = 'y_to_z'
    TRANSPOSE_XYZ = 'x_to_y_to_z_to_x'

    # Items for the spatial processing options
    choices = collections.OrderedDict(((NONE , "None"),
                                       (TRANSPOSE_XY , "Transpose_XY"),
                                       (TRANSPOSE_XZ , "Transpose_XZ"),
                                       (TRANSPOSE_YZ , "Transpose_YZ"),
                                       (TRANSPOSE_XYZ, "Transpose_XYZ"),
                                      ))

class SvdThreshold(object):
    """ PPM threshold constants """
    MIN = -200
    MAX =  200

class SvdThresholdUnit(object):
    """ SVD Threshold Value Unit constants """

    HZ  = 'Hz'
    PPM = 'PPM'

    # Items for the svd processing options dropdown
    choices = collections.OrderedDict(((HZ ,  "Hz"),
                                       (PPM , "PPM"),
                                      ))


class WaterExtrapolation(object):
    """ Water extrapolation constants """
    # These constants are arbitrary and may change.
    # However bool(NONE) is guaranteed to be False while
    # bool(X) == True is guaranteed for all other values.
    NONE = 0
    LINEAR = 1
    AR_MODEL = 2

    MIN_POINTS = 1
    MAX_POINTS = 1000

    # Items for the spectral processing options dropdown
    choices = collections.OrderedDict(((NONE ,     "None"),
                                       (LINEAR ,   "Linear"),
                                       (AR_MODEL , "AR Model"),
                                      ))

class ZeroFillMultiplier(object):
    """
    Zero fill multiplier constants. The zero fill is not only limited
    to a specific range, but it must also be an integral power of 2.
    """
    _EXPONENT = 5
    MIN = 1
    MAX = 2 ** _EXPONENT

    # Items for the spatial processing options
    choices = [ (2 ** i, str(2 ** i)) for i in range(0, _EXPONENT + 1) ]
    choices = collections.OrderedDict(choices)



##################   Preprocess Summed Constants

class PrepFidsumDataExcludeDisplayMethod(object):
    """
    Select what to display in the GUI plot to reflect what has changed over
    the acquisition of all the FIDs
    
    """
    # These constants are arbitrary and may change.
    INTEGRAL_ENTIRE  = 'integral_entire'
    INTEGRAL_INITIAL = 'integral_initial'
    PEAK_SHIFT_HZ    = 'peak_shift_hz'
    PEAK_PHASE_DEG   = 'peak_phase_deg'

    # Items for the spectral processing options dropdown
    choices = collections.OrderedDict(((INTEGRAL_ENTIRE,  "Integral (entire)"),
                                       (INTEGRAL_INITIAL, "Integral (initial 20%)"),
                                       (PEAK_SHIFT_HZ,    "Peak Shift [Hz]"),
                                       (PEAK_PHASE_DEG,   "Peak Phase [deg]"),
                                      ))


class PrepWbnaa(object):
    """
    Summed FID raw tab widget constants for both data processing
    and data correction routines.
    """
    MIN_APODIZE  = -100.0
    MAX_APODIZE  =  100.0
    STEP_APODIZE =  0.5

    MIN_LEFT = 0
    MAX_LEFT = 100

    MIN_PEAK  =  -1000.0
    MAX_PEAK  =   1000.0
    STEP_PEAK = 0.5

    MIN_PHASE  =  -360.0
    MAX_PHASE  =   360.0
    STEP_PHASE = 1.0

    MIN_PHASE1  =  -10000.0
    MAX_PHASE1  =   10000.0
    STEP_PHASE1 = 50.0

    MIN_CENTER  = -1.0
    MAX_CENTER  =  8.0
    STEP_CENTER = 0.01

    MIN_WIDTH =   0.0
    MAX_WIDTH = 100.0
    STEP_WIDTH = 0.01

    MIN_LEFT_B0 = 0
    MAX_LEFT_B0 = 100

    MIN_LEFT_PHASE0 = 0
    MAX_LEFT_PHASE0 = 100

    MIN_REF_LINE_WIDTH = 2
    MAX_REF_LINE_WIDTH = 100

    MAX_CONSTANT_PH0 = 360
    MIN_CONSTANT_PH0 = -360


class WbnaaRefSpectrumSource(object):
    """ Ref spectrum source for Phase0 correction constants """

    # These constants are arbitrary and may change.
    AVERAGE_ALL_FIDS = 'average_all_fids'
    SINGLET_CENTERED_IN_RANGE = 'singlet_centered_in_range'

    # Items for the spectral processing options dropdown
    choices = collections.OrderedDict(((AVERAGE_ALL_FIDS ,   "Average of All FIDs"),
                                       (SINGLET_CENTERED_IN_RANGE , "Singlet Centered in Range"),
                                      ))

class PrepMegalaser(object):
    """
    Summed FID raw tab widget constants for both data processing
    and data correction routines.
    """
    MIN_APODIZE  = -100.0
    MAX_APODIZE  =  100.0
    STEP_APODIZE =  0.5

    MIN_LEFT = 0
    MAX_LEFT = 100

    MIN_PEAK  =  -1000.0
    MAX_PEAK  =   1000.0
    STEP_PEAK = 0.5

    MIN_PHASE  =  -360.0
    MAX_PHASE  =   360.0
    STEP_PHASE = 1.0

    MIN_PHASE1  =  -10000.0
    MAX_PHASE1  =   10000.0
    STEP_PHASE1 = 100.0

    MIN_CENTER  = -1.0
    MAX_CENTER  =  8.0
    STEP_CENTER = 0.01

    MIN_WIDTH =   0.0
    MAX_WIDTH = 100.0
    STEP_WIDTH = 0.01

    MIN_LEFT_B0 = 0
    MAX_LEFT_B0 = 100

    MIN_LEFT_PHASE0 = 0
    MAX_LEFT_PHASE0 = 100

    MIN_REF_LINE_WIDTH = 2
    MAX_REF_LINE_WIDTH = 100

    MAX_CONSTANT_PH0 = 360
    MIN_CONSTANT_PH0 = -360


class MegalaserRefSpectrumSource(object):
    """ Ref spectrum source for Phase0 correction constants """

    # These constants are arbitrary and may change.
    AVERAGE_ALL_FIDS = 'average_all_fids'
    SINGLET_CENTERED_IN_RANGE = 'singlet_centered_in_range'

    # Items for the spectral processing options dropdown
    choices = collections.OrderedDict(((AVERAGE_ALL_FIDS ,   "Average of All FIDs"),
                                       (SINGLET_CENTERED_IN_RANGE , "Singlet Centered in Range"),
                                      ))

###################  Dataset Constants

class Dataset(object):
    """ Dataset constants """
    SCALE_INCREMENT = 1.1


##################   Voigt constants

class VoigtDynMetabolite(object):
    """ Metabolites fixed T2 value range and default """
    AREA_SCALE_MIN       =   0.00001
    AREA_SCALE_MAX       = 100000.0
    AREA_SCALE_INCR      = 1.0
    AREA_SCALE_MULT      = 1.25
    AREA_SCALE_DIGITS    = 5

    SEARCH_CENTER_INCR   = 0.1
    SEARCH_CENTER_DIGITS = 3

    SEARCH_WIDTH_MIN      =   0.001
    SEARCH_WIDTH_MAX      = 100.0
    SEARCH_WIDTH_INCR     =   0.01
    SEARCH_WIDTH_DIGITS   = 3

    SEARCH_PHASE0_MIN     = -360.0
    SEARCH_PHASE0_MAX     = 360.0
    SEARCH_PHASE0_INCR    = 5
    SEARCH_PHASE0_DIGITS  = 1

class VoigtDefaultFixedT2(object):
    """ Metabolites fixed T2 value range and default """
    MIN     =    0.001
    MAX     = 1000.0
    INCR    = 0.01
    DEFAULT = 1000.0
    CENTER  =    0.250
    DIGITS  = 3


##################   Giso constants

class GisoDynMetabolite(object):
    """ Metabolites fixed T2 value range and default """
    AREA_SCALE_MIN       =   0.00001
    AREA_SCALE_MAX       = 100000.0
    AREA_SCALE_INCR      = 1.0
    AREA_SCALE_MULT      = 1.25
    AREA_SCALE_DIGITS    = 5

    SEARCH_CENTER_INCR   = 0.1
    SEARCH_CENTER_DIGITS = 3

    SEARCH_WIDTH_MIN      =   0.001
    SEARCH_WIDTH_MAX      = 100.0
    SEARCH_WIDTH_INCR     =   0.01
    SEARCH_WIDTH_DIGITS   = 3

    SEARCH_PHASE0_MIN     = -360.0
    SEARCH_PHASE0_MAX     = 360.0
    SEARCH_PHASE0_INCR    = 5
    SEARCH_PHASE0_DIGITS  = 1

class GisoDefaultFixedT2(object):
    """ Metabolites fixed T2 value range and default """
    MIN     =    0.001
    MAX     = 1000.0
    INCR    = 0.01
    DEFAULT = 1000.0
    CENTER  =    0.250
    DIGITS  = 3

##################   General Fit constants

class FitPriorMaskSource(object):
    """ Source for spatial map used as prior information """
    # These constants are arbitrary and may change.
    # However bool(NONE) is guaranteed to be False while
    # bool(X) == True is guaranteed for all other values.
    NONE         = ''
    USER_DEFINED = 'user_defined'
    ALL_ON       = 'all_on'
    ALL_OFF      = 'all_off'

    # Items for the Voigt fitting tool radio buttons
    choices = collections.OrderedDict(((NONE ,         "None"),
                                       (USER_DEFINED , "User Defined"),
                                       (ALL_ON ,       "All Voxels On"),
                                       (ALL_OFF ,      "All Voxels Off"),
                                      ))

class FitPriorCalculateCombinations(object):
    """ lists available metabolite combinations that can be calculated """
    choices = ["naa+naag", "cr+pcr", "gpc+pcho", "cr2+pcr2", "glu+gln", "tau+glc"]
    

class FitLineshapeModel(object):
    """ Lineshape model """
    NONE    = ''
    VOIGT   = 'voigt'
    LORENTZ = 'lorentzian'
    GAUSS   = 'gaussian'

    # Items for the Voigt lineshape method combo
    choices = collections.OrderedDict(((VOIGT ,   "Voigt"),
                                       (LORENTZ , "Lorentzian"),
                                       (GAUSS ,   "Gaussian"),
                                      ))


class FitInitialB0ShiftMethod(object):
    """ Initial B0 Shift Method """
    MANUAL         = 'manual'
    AUTO_CORRELATE = 'auto_correlate'

    # Items for the Voigt inital B0 shift combo
    choices = collections.OrderedDict(((MANUAL ,         "Manual"),
                                       (AUTO_CORRELATE , "Auto Correlate"),
                                      ))

class FitInitialBaselineMethod(object):
    """ Initial Baseline Method """
    # These constants are arbitrary and may change.
    # However bool(NONE) is guaranteed to be False while
    # bool(X) == True is guaranteed for all other values.
    NONE          = ''
    LOWESS_FILTER = 'lowess_filter'
    SAVGOL_FILTER = 'savgol_filter'

    # Items for the Voigt inital Baseline Method combo
    choices = collections.OrderedDict(((NONE ,          "None"),
                                       (LOWESS_FILTER , "Lowess Filter"),
                                       (SAVGOL_FILTER , "Savitzky-Golay"),
                                      ))

class FitInitialSmallPeakAreas(object):
    """ Initial Small Peak Areas method """
    PEAK_SEARCH = 'peak_search'
    NAA_RATIO   = 'naa_ratio'
    CR_RATIO    = 'cr_ratio'
    MEAN_RATIO  = 'mean_ratio'
    PEAK_2PPM   = 'peak_2ppm'
    PEAK_3PPM   = 'peak_3ppm'

    # Items for the Voigt inital Baseline Method combo
    choices = collections.OrderedDict(((PEAK_SEARCH , "Peak search"),
                                       (NAA_RATIO ,   "NAA ratio"),
                                       (CR_RATIO ,    "Cr ratio"),
                                       (MEAN_RATIO ,  "Mean NAA+Cr ratio"),
                                       (PEAK_2PPM,    "2.01 PPM Peak"),
                                       (PEAK_3PPM,    "3.02 PPM Peak"),
                                      ))

class FitInitialSmallPeakFreqs(object):
    """ Initial Small Peak Freqs method """
    PEAK_SEARCH = 'peak_search'
    REF_PEAK    = 'reference_peak'

    # Items for the Voigt inital Baseline Method combo
    choices = collections.OrderedDict(((PEAK_SEARCH , "Peak sarch"),
                                       (REF_PEAK ,    "Ref peak"),
                                      ))

class FitInitialLinewidthMethod(object):
    """ Initial Linewidth Method """
    MANUAL         = 'manual'
    DECONVOLUTION  = 'deconvolution'
    AUTO_CORRELATE = 'auto_correlate'

    # Items for the Voigt inital linewidth combo
    choices = collections.OrderedDict(((MANUAL ,         "Manual"),
                                       (DECONVOLUTION ,  "Deconvolution"),
                                       (AUTO_CORRELATE , "Auto Correlate"),
                                      ))

class FitInitialPhaseMethod(object):
    """ Initial Phase Method """
    MANUAL              = 'manual'
    CORRELATION_PHASE0  = 'correlation_phase0'
    CORRELATION_PHASE01 = 'correlation_phase01'
    INTEGRATION_PHASE0  = 'integration_phase0'
    INTEGRATION_PHASE01 = 'integration_phase01'

    # Items for the Voigt inital phase combo
    choices = collections.OrderedDict(((MANUAL ,              "Manual"),
                                       (CORRELATION_PHASE0 ,  "Correlation - Phase0 only"),
                                       (CORRELATION_PHASE01 , "Correlation - Phase0+Phase1"),
                                       (INTEGRATION_PHASE0 ,  "Integration - Phase0 only"),
                                       (INTEGRATION_PHASE01 , "Integration - Phase0+Phase1"),
                                      ))


class FitBaselineMethod(object):
    """ Baseline Algorithm """
    # These constants are arbitrary and may change.
    # However bool(NONE) is guaranteed to be False while
    # bool(X) == True is guaranteed for all other values.
    NONE                    = ''
    WAVELET_FILTER_BASIC    = 'wavelet_filter_basic'
    BSPLINE_FIXED_KNOT      = 'bspline_fixed_knot'
    BSPLINE_VARIABLE_KNOT   = 'bspline_variable_knot'

    # Feel free to change the DEFAULT to anything but WAVELET_FILTER_BASIC.
    # That one doesn't work well as a default because PyWavelets is optional
    # and we don't want to set the default to be something that people
    # might not have installed.
    DEFAULT = WAVELET_FILTER_BASIC

    # Items for the Voigt baseline method combo
    choices = collections.OrderedDict(((NONE , "None"),
                                       (WAVELET_FILTER_BASIC ,  "Wavelet Filter (basic)"),
                                       (BSPLINE_FIXED_KNOT ,    "B-spline (fixed knot)"),
                                       (BSPLINE_VARIABLE_KNOT , "B-spline (variable knot)"),
                                      ))

class FitBaselineUnderestimateMethod(object):
    """ Baseline Underestimation Flavor """
    # These constants are arbitrary and may change.
    # However bool(NONE) is guaranteed to be False while
    # bool(X) == True is guaranteed for all other values.
    FIRST_ITERATION_ONLY       = 'first_iteration_only'
    MULTI_STEP_UNDERESTIMATION = 'multi_step_underestimation'

    # Feel free to change the DEFAULT to anything but WAVELET_FILTER_BASIC.
    # That one doesn't work well as a default because PyWavelets is optional
    # and we don't want to set the default to be something that people
    # might not have installed.
    DEFAULT = FIRST_ITERATION_ONLY

    # Items for the Voigt baseline method combo
    choices = collections.OrderedDict(((FIRST_ITERATION_ONLY ,       "First Iteration Only"),
                                       (MULTI_STEP_UNDERESTIMATION , "Multi-step Underestimation"),
                                      ))


class FitMacromoleculeMethod(object):
    """ Macromolecule Function """
    # These constants are arbitrary and may change.
    # However bool(NONE) is guaranteed to be False while
    # bool(X) == True is guaranteed for all other values.
    NONE                    = ''
    SINGLE_BASIS_DATASET    = 'single_basis_dataset'

    # Feel free to change the DEFAULT to anything
    DEFAULT = NONE

    # Items for the Voigt macromolecule model combo
    choices = collections.OrderedDict(((NONE , "None"),
                                       (SINGLE_BASIS_DATASET , "Single Basis Function - from Dataset"),
                                      ))


# class FitMacromoleculeMethod(object):
#     """ Macromolecular Peak Method """
#     # These constants are arbitrary and may change.
#     # However bool(NONE) is guaranteed to be False while
#     # bool(X) == True is guaranteed for all other values.
#     NONE                = ''
#     GROUPED_GROUPED     = 'grouped_grouped'
#     GROUPED_SEPARATE    = 'grouped_separate'
#     SEPARATE_GROUPED    = 'separate_grouped'
#     SEPARATE_SEPARATE   = 'separate_separate'
#
#     # Items for the Voigt method combo
#     choices = collections.OrderedDict(((NONE , "None"),
#                                        (GROUPED_GROUPED ,   "Area/PPM (grouped), Linewidths (grouped)"),
#                                        (GROUPED_SEPARATE ,  "Area/PPM (grouped), Linewidths (separate)"),
#                                        (SEPARATE_GROUPED ,  "Area/PPM (separate), Linewidths (grouped)"),
#                                        (SEPARATE_SEPARATE , "Area/PPM (separate), Linewidths (separate)"),
#                                       ))

class FitMacromoleculeMethodInitVal(object):
    """ Macromolecule Function """
    # These constants are arbitrary and may change.
    MANUAL   = 'manual'
    ONEPOINT = 'onepoint'
    REGION   = 'region'
    
    # Feel free to change the DEFAULT to anything
    DEFAULT = MANUAL

    # Items for the Voigt macromolecule single dataset initial value method
    choices = collections.OrderedDict(((MANUAL   , "Manual"),
                                       (ONEPOINT , "One Point"),
                                       (REGION   , "Region"),
                                      ))


class FitOptimizeMethod(object):
    """ Optimization Algorithm """
    # These constants are arbitrary and may change.
    # However bool(NONE) is guaranteed to be False while
    # bool(X) == True is guaranteed for all other values.
    NONE                            = ""
    CONSTRAINED_LEVENBERG_MARQUARDT = "constrained_levenburg_marquardt"
    LMFIT_DEFAULT                   = "lmfit_default"

    # Items for the Voigt fitting tool radio buttons
    choices = collections.OrderedDict(((NONE , "None"),
                                       (CONSTRAINED_LEVENBERG_MARQUARDT , "Constrained Levenburg-Marquardt"),
                                       (LMFIT_DEFAULT ,                   "LMFit Default"),
                                      ))

class FitOptimizeWeightsMethod(object):
    """ Macromolecular Linewidths Method """
    EVEN_WEIGHTING  = 'even_weighting'
    LOCAL_WEIGHTING = 'local_weighting'

    # Items for the Voigt fitting tool radio buttons
    choices = collections.OrderedDict(((EVEN_WEIGHTING ,  "Even Weighting"),
                                       (LOCAL_WEIGHTING , "Local Weighting"),
                                      ))

FIT_OPTIMIZE_EXCLUDES_FOR_SMALL_PEAKS = ['naa', 'cho', 'cr', 'h2o']

class FitAmplitudeMultiplier(object):
    """ Metabolites amplitude multiplier constants """
    MIN =    0.001
    MAX = 1000.0

class FitBaselineBsplineOrder(object):
    """ Baseline B-Spline order constants """
    MIN = 1.0
    MAX = 5.0

class FitBaselineLowessWindowSize(object):
    """ Baseline metabolites region Lowess window size (Hz) constants """
    MIN =    0.00001
    MAX = 5000.0

class FitBaselineUnderestimation(object):
    """ Baseline first pass underestimation constants """
    MIN =  -50.0
    MAX =  100.0

class FitBaselineUnderestimationLast(object):
    """ Baseline first pass underestimation constants """
    MIN =  -50.0
    MAX =  100.0

class FitLineWidth(object):
    """ Metabolites line width constants """
    MIN =    0.001
    MAX = 1000.0

class FitMacroMoleculeLines(object):
    """ Metabolites Voigt macro molecule model lines constants """
    MIN =  1
    MAX = 50

class FitMacroMoleculeAdjustment(object):
    """ Metabolites Voigt macro molecule adjustment (spinner) constants """
    MIN =   0
    MAX = 100

class FitOptimizationAmplitude(object):
    """ Fitting Voigt optimization metabolite amplitude constants """
    MIN =     1
    MAX = 10000
    
class FitOptimizationAmplitudeBoundPercentMax(object):
    """ Fitting Voigt optimization metabolite amplitude constants """
    MIN =   100.0
    MAX = 10000.0

class FitOptimizationAmplitudeBoundPercentMin(object):
    """ Fitting Voigt optimization metabolite amplitude constants """
    MIN = 1.0
    MAX = 99.0

class FitOptimizationAreaWeight(object):
    """ Fitting Voigt optimization area weight constants """
    MIN =       0.0001
    MAX = 1000000.0

class FitOptimizationAlgorithmIterations(object):
    """ Fitting Voigt optimization algorithm max iterations constants """
    MIN =     1
    MAX = 10000

class FitOptimizationConfidenceAlpha(object):
    """ Fitting Voigt optimization algorithm confidence alpha constants """
    MIN = 0.05
    MAX = 0.9999

class FitOptimizationFrequency(object):
    """ Fitting Voigt optimization metabolite frequency constants """
    MIN =     0.1
    MAX = 10000.0

class FitOptimizationGlobalIterations(object):
    """ Fitting Voigt optimization global iterations constants """
    MIN =    1
    MAX = 1000

class FitOptimizationLocalMultiplier(object):
    """ Fitting Voigt optimization "LW Local Mult" (???) constants """
    MIN =   1.0
    MAX = 100.0

class FitOptimizationPhase1(object):
    """ Fitting Voigt optimization metabolite phase 1 constants """
    MIN =    1
    MAX = 5000

class FirtOptimizationStopTolerance(object):
    """ Fitting Voigt optimization algorithm stop tolerance constants """
    MIN =      0.000000001
    MAX = 100000.0

class FitOptimizationTaTb(object):
    """ Fitting Voigt optimization Ta=Tb constants """
    MIN =   0.001
    MAX = 200.0

class FitPeakPpm(object):
    """ Metabolites peak PPM constants """
    MIN = -5000.0
    MAX =  5000.0

class FitPeakSearchRange(object):
    """ Metabolites peak search range constants """
    MIN =  0.001
    MAX = 10.0

# Default params for - ppm, area, phase, lwhz, lim_ppm, lim_area, lim_phase, lim_lwhz
VOIGT_NEW_MACROMOLECULE_ROW = [True, 1.0, 1.0, 0.0, 5.0, 1.0, 5.0, 1.0, 10.0]


##################   Watref constants

class WatrefCorrections(object):
    """ Source for water reference quantification correction """
    # These constants are arbitrary and may change.
    # However bool(NONE) is guaranteed to be False while
    # bool(X) == True is guaranteed for all other values.

    WATER_AVERAGES_DEFAULT =   0.001
    WATER_AVERAGES_MIN     =   0.001
    WATER_AVERAGES_MAX     =   0.001
    MAX = 200.0
    
    
    NONE         = ''
    USER_DEFINED = 'user_defined'
    ALL_ON       = 'all_on'
    ALL_OFF      = 'all_off'

    # Items for the Voigt fitting tool radio buttons
    choices = collections.OrderedDict(((NONE ,         "None"),
                                       (USER_DEFINED , "User Defined"),
                                       (ALL_ON ,       "All Voxels On"),
                                       (ALL_OFF ,      "All Voxels Off"),
                                      ))




