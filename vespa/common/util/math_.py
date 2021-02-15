# Python modules


# 3rd party modules
import numpy as np

# Our modules
import vespa.common.util.misc as util_misc

# _SAFE_EXP_BOUND is for the function safe_exp() (q.v.). It's only used there;
# we make it a constant just so that it's easy to spot here at the top of
# the file.
_SAFE_EXP_BOUND = 60


def _float_eq(a, b, rtol=1.e-5, atol=1.e-8):
    # Implementation swiped from numpy.allclose()
    # Any errors are my fault.
    return abs(a - b) <= atol + rtol * abs(b)


def eq(a, b, rtol=1.e-5, atol=1.e-8):
    """Returns True if a and b are "equal", where "equal" means "close enough"
    as defined by the relative and absolute tolerance parameters. The
    parameters and implementation match the numpy function allclose().
    
    a and b must be floats or complex numbers. In the latter case, a.real
    is compared to b.real and the imaginary components are treated the same.

    Here's an example of how to use this to compare two lists of floats:
    all([util_math.eq(a, b) for (a, b) in zip(one_list, another_list)])

    I would have called this function "close" as in, "PI is close to 22/7",
    but it's easy to confuse with "close" as in, "close the file".
    """
    if isinstance(a, (int, float)):
        return _float_eq(a, b, rtol, atol)
    else:
        return _float_eq(a.real, b.real, rtol, atol) and \
               _float_eq(a.imag, b.imag, rtol, atol)


def safe_exp(value, bound=_SAFE_EXP_BOUND):
    """Given a value (which is typically a float array but can also be a 
    scalar), applies np.exp() to the value and returns the result. The "safe"
    aspect of the function is that it offers two methods to avoid under/overflow
    errors.

    If bound is left to its default (currently 60), then the data is cropped
    to +/- that limit. In other words, for each elements x of the array, the
    following holds true:
        -bound < x < bound

    If bound is 0, then the data is not cropped. Instead we temporarily disable
    numpy under/overflow errors while calling exp(). If the result contains
    any NaNs due to under/overflow, those NaNs are set to 0.
    """
    if bound:
        if isinstance(value, complex):
            value = complex(max(min(value.real, bound), -bound), 
                            max(min(value.imag, bound), -bound))
        else:
            # It's a numpy array
            value = np.where(value >  bound,  bound, value)
            value = np.where(value < -bound, -bound, value)

        value = np.exp(value)
    else:
        restore_these_settings = np.geterr()

        temp_settings = restore_these_settings.copy()
        temp_settings["over"] = "ignore"
        temp_settings["under"] = "ignore"

        np.seterr(**temp_settings)
        value = np.exp(value)
        np.seterr(**restore_these_settings)
        
        zeros = np.zeros_like(value)

        value = np.where(np.isnan(value), zeros, value)

    return value



def _test():
    # Unit tests for eq()
    # This is guaranteed to be as correct as numpy's allclose() from which
    # we've stolen the implementation.
    import numpy
    
    # Test data is params of a, b, rtol, atol, expected outcome. 
    # None ==> use default.
    test_data = ( (42, 42, None, None, True), 
                  (42, 42 + 1.e-6, None, None, True), 
                  (42, 42 + 1.e-9, 0, None, True), 
                  (42, 42 + 1.e-7, 0, None, False), 
                  (42, 43, None, None, False), 
                  (.01, .01 + 1.e-7, None, None, True), 
                  (.01, .01 + 1.e-7, 0, None, False), 
                  (complex(42, 42), complex(42, 42), None, None, True), 
                  (complex(42, 0), complex(42, 0), None, None, True), 
                  (complex(0, 42), complex(42, 0), None, None, False), 
                  (complex(42 + 1.e-6, 0), complex(42, 0), None, None, True), 
                  (complex(42 + 1.e-8, 0), complex(42, 0), None, None, True), 
                  (complex(42 + 1.e-5, 0), complex(42, 0), None, None, True), 
                  (complex(1.e-5, 0), complex(0, 0), None, None, False), 
                )
                
    for a, b, rtol, atol, outcome in test_data:
        kwargs = { }
        if rtol is not None:
            kwargs["rtol"] = rtol
        if atol is not None:
            kwargs["atol"] = atol
        
        assert(eq(a, b, **kwargs) == outcome)
        a = numpy.array([a])
        b = numpy.array([b])
        assert(numpy.allclose(a, b, **kwargs) == outcome)
        
