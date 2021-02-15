# Python modules

# 3rd party modules
import numpy as np




def auto_correlate(x, lag, covariance=False):
    """
    This method calculates the autocorrelation or autocovariance of a data 
    collection x as a function of the lag.

    x:    An numpy array of type integer, float or complex.

    lag:  A numpy array, in the interval [-(n-2), (n-2)], of integers to
          specify the absolute distance(s) between indexed elements of x.

    covariance: keyword, bool, if set the sample autocovariance is returned

    Reference:
        INTRODUCTION TO STATISTICAL TIME SERIES, Wayne A. Fuller, ISBN 0-471-28715-6

    """
    nx = len(x)
    if nx < 2:
        raise ValueError("x array must contain 2 or more elements.") 

    nlag = len(lag)
    corr = np.zeros(nlag, dtype=x.dtype)

    # Here we subtract the mean over the entire data series, rather than 
    # subtracting the mean from the first and last n-lag points separately. 
    # This is discussed further in Jenkins & Watts, Spectral Analysis and 
    # its Applications, 1968. In short, the use of separate means for each 
    # portion is not recommended, as it is not a satisfactory estimate 
    # when a several autocorrelations at different lags are required.

    data = x - (np.sum(x) / nx)

    # Compute autocovariance
    M = np.abs(lag)
    for k in range(nlag): 
        corr[k] = np.sum(data[0:nx - M[k]] * data[M[k]:])

    # Divide by nx for autocovariance, or by variance for autocorrelation.
    if covariance:
        corr /= nx
    else:
        corr /= np.sum(data * data)

    return corr