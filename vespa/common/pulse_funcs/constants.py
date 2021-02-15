# Python modules


# Constants specific to pulse_funcs


class ApodizationFilterType(object):
    NONE = 0
    COSINE = 1
    HAMMING = 2
    
class AnalyticType(object):
    NONE = 0
    GAUSSIAN = 1
    SINC_GAUSSIAN = 2
    HYPERBOLIC_SECANT = 3

class ProfileType(object):
    NONE = 0
    M_XY = 1
    M_X_MINUS_Y = 2    
    M_Z  = 3
    M_MINUS_Z = 3


# Gyromagnetic ratio of 1H - the hydrogen nucleus. (units: kHz/mT)
GAMMA1H = 42.576

# Size limit on b1 for root reflection
# b1 to polynomial will automagically truncate to this
# so calling code should check and issue a warning
b1rootlimit = 65

# Small number used for floating point comparisons.
epsilon = 0.000001  

# Very small number used for double precision comparisons.
small_epsilon = 0.00000000001  

# An even smaller number; used to represent an acceptable fractional error
# or difference from an expected or desired value.
EPS = pow(2,-52)

