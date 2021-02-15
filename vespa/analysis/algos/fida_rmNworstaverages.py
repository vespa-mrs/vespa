# Python modules

# 3rd party modules
import numpy as np

# Our modules


def op_rmNworstaverages(data, n, sw):
    """
    USAGE:
     badAverages, metric = op_rmNworstaverages(data, sw, n)

    DESCRIPTION:
     Remove motion corrupted averages from a dataset. The N most badly motion
     corrupted averages are discarded. Bad averages are identified by
     calculating a 'likeness' metric for each average in the Frequency domain.
     This is done by subtracting each spectra from the median spectrum, and
     then calculating the summed square of each average's difference spectrum.

    Inputs:
        data (ndarray): dims=(navg,npts) complex k-space data
        sw (float): sweep width in Hz
        n (int): number of indices for bad averages to return

    Outputs:
        badAverages (list, int): list of indices indicating the averages that
            were 'bad' in this data array. Can be empty.
        metric (list, float): list of unlikeness metrics corresponding to all
            input averages.

    Derived from FID-A project: module - op_rmNworstaverages.m
    Jamie Near, McGill University 2014.

    """
    nfid, npts = data.shape
    t = np.arange(npts) / sw
    x = np.arange(nfid)

    if nfid == 1:
        raise ValueError('ERROR:  Averaging has already been performed!  Aborting!')

    # subtract all spectra from the median spectrum and sum the squared difference values

    infilt = data.copy() * np.exp(-(t * np.pi * 10.0))          # op_filter(data, 10)
    infilt = np.fft.fftshift(np.fft.ifft(infilt, axis=1), axes=1)
    inmed  = np.squeeze(np.median(infilt.real, axis=0) + 1j * np.median(infilt.imag, axis=0))  # take median

    metric = np.sum( (infilt.real - inmed.real)**2, axis=1)

    # metric = np.zeros([nfid, ], dtype=np.float)
    # for k in range(nfid):
    #     metric[k]=np.sum( (infilt[k,:].real - inmed.real)**2)       # freq domain

    # z-transform the metric so it's centered about zero, with stdev of 1.0
    avg   = np.mean(metric)
    stdev = np.std(metric)
    zmetric = (metric-avg)/stdev
    pfit = np.polyfit(x, zmetric, 2)

    # sort the zmetric array to find the n highest values:
    delta = zmetric - np.polyval(pfit,x)         # FIXME bjs - should this be abs() here to find greatest + or - distance from median??
    badAverages = np.argsort(delta)[::-1]        # sort descending order
    badAverages = np.array(badAverages[0:n], dtype=np.int)

    return badAverages, metric