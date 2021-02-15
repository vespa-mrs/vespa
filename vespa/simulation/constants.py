# Python modules



# these are default values and max/min values for widgets in
# the Visualize tab on each Experiment Tab
SPECTRAL_POINTS_MAX = 65536             # integer
SPECTRAL_POINTS_MIN =    64             # integer
SWEEP_WIDTH_MAX     = 100000.0          # [Hz]
SWEEP_WIDTH_MIN     =    100.0          # [Hz]
LINEWIDTH_MAX       =   1000.0          # [Hz]
LINEWIDTH_MIN       =      0.000001     # [Hz]


class ThirdPartyExportTypes(object):
    # These constants are arbitrary and may change.
    # However bool(NONE) is guaranteed to be False and bool() of anything
    # else will always be True.
    NONE             = 0
    ANALYSIS_PRIOR   = 1
    MIDAS_PRIOR      = 2
    LCMODEL          = 3
    JMRUI            = 4
    GAVA             = 5
    ANALYSIS_DIRECT  = 6
