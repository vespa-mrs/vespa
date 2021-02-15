# Python modules

# 3rd party modules
import numpy as np

# Our modules
from vespa.analysis.algos.cross_correlate import cross_correlate



def b0_correction( dat, ref, magn=False, cdeg=1, nlag=None):
    """
    B0 correction algorithm number 9.  Uses correlation of data region to an ideal
    spectrum to determine integer shift of data to match.

    INPUT:
     dat: complex array, raw data
     ref: complex array, ideal spectrum
    KEYWORDS:
     CORR: correlation value array, [magn,real,imag]
     CSUM: sum of real and imag ccvals
     NLAG: number of points to lag in correlation
     CDEG: real and imag cc val mixing degree
     MAGN: flag, performs calc only on Magnitude data when set

    """
    error = 0
    corr  = [-1.0,-1.0,-1.0]
    csum  = -2.0

    ndat = len(dat)

    if ndat != len(ref):
        return 0.0

    if not nlag:
        nlag = ndat

    shft  = 0.0
    ndat2 = int(ndat/2)
    nlag2 = int(nlag/2)
    lag   = np.arange(nlag) - nlag2
    indx  = [0.0,0.0,0.0]

    tmp = cross_correlate(np.abs(dat),np.abs(ref), lag)

    corr[0] = np.max(tmp)
    indx[0] = np.argmax(tmp)
    shft = indx[0] - nlag2 + 1 # based on empirical test with Prior data sets 
    csum = 0.0

    if not magn:
        tmp     = cross_correlate(dat.real, ref.real, lag)
        corr[1] = np.max(tmp)
        indx[1] = np.argmax(tmp)
        tmp     = cross_correlate(dat.imag, ref.imag, lag)
        corr[2] = np.max(tmp)
        indx[2] = np.argmax(tmp)

        isort = np.argsort(corr)
        shft  = indx[isort[2]] - nlag2 + 1
        csum  = corr[1]**cdeg + corr[2]**cdeg

    return int(shft), csum