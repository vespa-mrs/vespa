# Python modules

# 3rd party modules
import numpy as np

# Our modules


def op_rmbadaverages( data, sw, nsd='3', domain='t'):
    """
    USAGE:
     badAverages, metric = op_rmbadaverages(data, sw, nsd=nsd, domain=domain)

    DESCRIPTION:
     Removes motion corrupted averages from a dataset containing multiple
     averages.  Bad averages are identified by calculating a 'likeness' metric
     for each average by subtracting each average from the median average,
     and then calculating the summed square of real part of the difference.
     Likeness metrics greater than 'nsd' above the mean are discarded.

    Inputs:
        data (ndarray, complex): shape=(navg,npts) complex k-space data
        sw (float): sweep width in Hz
        nsd (float): def=3.0, number of standard deviations to use a rejection threshold
        domain (string): def='t', domain ('t' or 'f') in which to perform calculations

    Outputs:
        badAverages (list, int): list of indices indicating the averages that
            were 'bad' in this data array. Can be empty.
        metric (list, float): list of unlikeness metrics corresponding to all
            input averages.

    Derived from FID-A project: module - op_rmbadaverages.m
    Jamie Near, McGill University 2014.

    """
    nfid, npts = data.shape
    raw = data.copy()
    t = np.arange(npts) / sw
    x = np.arange(nfid)

    if nfid == 1:
        raise ValueError('ERROR:  Averaging has already been performed!  Aborting!')

    # algorithm can be applied to Frequency or Time domain data

    if domain in ['t', 'T']:
        tmax = 0.4                                         # from fida_run_pressproc_auto.m
        trange = np.logical_and(t >= 0, t <= tmax)
    elif domain in ['f', 'F']:
        raw *= np.exp(-(t * np.pi * 10.0))                 # 10 Hz lorentz - from fida_run_pressproc_auto.m
        raw = np.fft.fftshift(np.fft.ifft(raw, axis=1), axes=1)

    inmed = np.squeeze(np.median(raw.real, axis=0) + 1j*np.median(raw.imag, axis=0))   # calculate median
    if domain == 't' or domain == 'T':
        metric = np.sum((raw[:,trange].real - inmed[trange].real)**2, axis=1)
    elif domain=='f' or domain=='F':
        metric = np.sum((raw.real - inmed.real) ** 2, axis=1)

    # z-transform the metric so it's centered about zero, with stdev of 1.0

    zmetric = (metric-np.mean(metric))/np.std(metric)
    pfit = np.polyfit(x, zmetric, 2)
    pval = np.polyval(pfit,x)

    # mask FIDs with metrics > nsd standard deviations from the mean metric value
    
    mask = zmetric > pval+nsd
    badAverages = np.array(x[mask], dtype=np.int16)      # these are the corrupted indices

    return badAverages, metric