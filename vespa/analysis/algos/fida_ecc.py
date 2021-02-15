"""

op_ecc.m
Jamie Near, McGill University 2014.

USAGE:
[out,outw]=op_ecc(data,inw)

DESCRIPTION:
Perform an eddy current correction by estimating any non-linearity in the
phase of the water unsuppressed data in the time domain and applying the
appropriate correction to both the water suppressed and water unsuppressed
data.

INPUTS:
data   = water suppressed input data in matlab structure format.
inw    = water unsuppressed input data in matlab structure format.

OUTPUTS:
out    = Water suppressed output following eddy current correction
outw   = Water unsuppressed output following eddy current correction

"""

# Python modules

# 3rd party modules
import numpy as np

# Our modules



def op_ecc(data, chain):

    nrep, ncoil, nfid, npts = data.shape
    ds = chain._dataset
    dwell = 1.0 / ds.sw
    t  = np.arange(npts)/dwell
    t150 = (npts * dwell) * np.arange(150, dtype=np.float)/149.0

    if ncoil != 1 or nfid != 0:
        msg = 'ERROR:  Must combine receivers and averages prior to running ecc!! ABORTING!!'
        raise ValueError(msg)

    # save the phase as a vector of hard numbers.
    inph = phase(data)

    # choose the part of the phase function that is most linear
    tmin, tmax = 0.02, 0.10 # in sec

    # now fit a straight line to the linear part of the phase function
    c = np.polyfit(t[t>tmin and t<tmax], inph[t>tmin and t<tmax],1)
    p = np.polyval(c, t)

    # now fit a spline to approximate a smooth version of the phase function
    # pp=splinefit(t,inph,150)
    ff = np.interpolate.splrep(t, inph, k=3, task=-1, t=t150[1:npts-1])
    pp = np.splev(t, ff)

    # subtract line from spline -> eddy current related phase offset
    ecphase = pp - p

    ecphase = np.exp(-1j * ecphase)

    return ecphase



def phase(data):
    """
    PHASE  Computes the phase of a complex vector

      PHI=phase(G)

      G is a complex-valued array and phi is returned as its
      phase (in radians), with an effort made to keep it continuous
      over the pi-borders.
    
      L. Ljung 10-2-86
      Copyright 1986-2004 The MathWorks, Inc.
      $Revision: 1.5.4.2 $  $Date: 2004/07/31 23:24:49 $
    
    PHI = unwrap(angle(G))

    """
    if len(data.shape) > 1:
        msg = "Data not a 1D array, exiting."
        raise ValueError(msg)

    phi = np.arctan2(data.imag, data.real)
    npts = len(data)
    df = phi[0:npts-1] - phi[1:npts]

    wraps = np.where(np.abs(df)>3.5)

    for i in wraps[0]:
        if i != 0:
            phi = phi + 2 * np.pi * np.sign(df[i]) * np.r_[np.zeros(i+1),np.ones(npts-i-1)]

    return phi


#------------------------------------------------------------------------------
# test code

def _test():
    import numpy as np
    from matplotlib import pyplot as plt

    npts = 256

    y = 4 * 360.0 * np.arange(npts, dtype=np.float)/(npts-1)
    y = y % 360.0
    y = np.exp(1j * (y * np.pi / 180.0))

    yy = phase(y)

    plt.plot(y)
    plt.plot(yy, 'b', linewidth=2)
    plt.show()

    bob = 10
    bob += 1


if __name__ == '__main__':
    """ Just testing the phase() method here """
    _test()



