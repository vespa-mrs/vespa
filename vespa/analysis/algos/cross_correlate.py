# Python modules

# 3rd party modules
import numpy as np


def cross_correlate(x, y, lag, covariance=False):
    """

    This method calculates the cross correlation Pxy(lag) or cross covariance 
    Rxy(lag) of two data sets x and y as a function of the lag.

    x:    a numpy array of type integer, float or complex.
    y:    a numpy array of type integer, float or complex.
    lag:  a numpy array, in the interval [-(n-2), (n-2)], of integers that 
          gives the absolute distance(s) between indexed members of x.

    covariance: bool, flag, if set the sample cross covariance is returned

    Reference: INTRODUCTION TO STATISTICAL TIME SERIES
                Wayne A. Fuller  ISBN 0-471-28715-6

    """

    nx = len(x)
    ny = len(y)
    if nx < 2 or ny < 2:
        raise ValueError("x and y arrays must contain two or more values.") 

    if nx != len(y):
        raise ValueError("x and y arrays must be same length.")

    xd = x - np.sum(x) / nx # deviations
    yd = y - np.sum(y) / nx

    nlag = len(lag)

    corr = np.zeros(nlag, dtype=x.dtype)

    for k in range(nlag):
        # reverse the variables for negative lags.
        if lag[k] > 0:
            corr[k] = np.sum(xd[0:nx - lag[k]] * yd[lag[k]:])
        else:
            corr[k] = np.sum(yd[0:nx + lag[k]] * xd[-lag[k]:])

    # Divide by N for covariance, or divide by variance for correlation.
    if covariance:
        corr /=  nx
    else:    
        corr /= np.sqrt(np.sum(xd*xd) * np.sum(yd*yd))

    return corr