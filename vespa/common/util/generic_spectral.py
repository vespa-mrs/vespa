# Python modules

import math

# 3rd party modules
import numpy as np

# Our modules
from vespa.common.util.math_ import safe_exp




def apodize(t, width, type_):
    """
    Given data to apodize, width, and apodization type, returns the line
    shape.

    t in this instance is the digitization time course for the spectral data
    to be apodized, namely: 0, 1*dwelltime, 2*dwelltime, ...

    The type can be 'G' for Gaussian apodization or 'L' for Lorentzian.

    Only the first letter of the type is considered and it is case insensitive,
    so 'G', 'g', 'gaussian', 'GAUSSIAN', 'Geronimo', 'glop', 'GGG', etc. all
    imply Gaussian.
    """
    if type_:
        # Grab the first letter and uppercase it.
        type_ = type_.upper()[0]

    if type_ == 'G':
        # Gaussian
        lineshape = -(width * math.pi * 0.6 * t) ** 2
        lineshape = safe_exp(lineshape, 0)
    elif type_ == 'L':
        # Lorentzian
        lineshape = -(width * math.pi * t)
        lineshape = safe_exp(lineshape, 0)
    else:
        raise ValueError("""Bad apodization type "{0}" """.format(type_))

    return lineshape


def create_fid(areas, ppms, phases, dataset, calc_peakppm=True):
    """
    Create an 'ideal' FID via a sum of individual FIDs from lists of spectral
    peak descriptors. Each peak descriptor has an area, ppm, and phase value.
    No global Phase 0 or Phase 1 are applied here.

    Args:
        areas (array-like, floats): area under each resonance peak
        ppms (array-like, floats): peak center for each resonance peak
        phases (array-like, floats): phase 0 for each resonance peak
        dataset (object):  Dataset object from mrs_dataset module. Provides
            data parameters on which the simultated spectrum is based, such as
            sweep width, zerofill, spectral points, etc.
        calc_peakppm (bool, optional): default=True, flag, determine if ppm
            value for the max peak location of the FFT of the summed FID is
            calculated using a higher zero fill factor for more precision.

    Usage:

      fid = calculate_fid(areas, ppms, phases, dataset, calc_peakppm=False)

    """
    acqdim0 = dataset.raw_dims[0]
    td      = 1.0/dataset.sw
    arr1    = np.zeros(acqdim0,float) + 1.0
    freqs   = dataset.ppm2hz(np.array(ppms)) * 2.0 * np.pi
    phases  = np.array(phases) * np.pi / 180.0
    npk     = np.size(freqs)

    xx = (np.arange(npk*acqdim0, dtype='float') % acqdim0) * td
    xx.shape = npk,acqdim0

    am  = np.outer(np.array(areas),arr1)
    fr  = np.outer(freqs,arr1)
    ph  = np.outer(phases,arr1)

    tmp = am * np.exp(1j * fr * xx) * np.exp(1j * ph)
    tmp = np.sum(tmp, axis=0)

    pkppm = None

    if calc_peakppm:

        zfnew = 4.0
        dim0  = int(round(acqdim0 * dataset.zero_fill_multiplier * zfnew))
        acqsw = dataset.sw

        xx  = np.arange(acqdim0)/acqsw
        dat = np.zeros(dim0, complex)

        #dat[0:acqdim0] = tmp * apodize(xx, 2.0, 'Gaussian') # 1Hz Gaussian Lineshape
        dat[0:acqdim0] = tmp * safe_exp(-(2.0*math.pi*0.6*xx)**2, 0)  # 1Hz Gaussian Lineshape

        dat  = np.fft.fft(dat) / len(dat)
        r    = max(dat.real)
        indx = np.where(r == dat.real)[0][0]

        hpp   = dataset.sw / dim0
        pkppm = ((dim0/2) - indx) * (hpp / dataset.frequency) + dataset.resppm

    return tmp, pkppm


def create_spectrum(areas, ppms, phases, dataset, ta=0.07, tb=0.07, acq=False, zforce=None):
    """
    Create a spectrum via a sum of spectral peaks. Each peak has an area, ppm,
    and phase used to create an 'ideal' FID. The set of summed FIDs is modified
    by a Voigt lineshape described by Ta and Tb parameters. No global Phase 0
    or Phase 1 is applied in here.

    Args:
        areas (array-like, floats): area under each resonance peak
        ppms (array-like, floats): peak center for each resonance peak
        phases (array-like, floats): phase 0 for each resonance peak
        dataset (object):  Dataset object from mrs_dataset module. Provides
            data parameters on which the simultated spectrum is based, such as
            sweep width, zerofill, spectral points, etc.
        ta (float, optional): default=0.07, Voigt lineshape Lorentzian parameter
        tb (float, optional): default=0.07, Voigt lineshape Gaussian parameter
        acq (bool, optional): default=False, determines if returned spectrum is
            sized based on current zero filling value, or on acquisition dims.
        zform (int, optional): default=None, can be used to force returned
            spectrum to be zero filled to an arbitrary size based on integer
            value of the zforce parameter.

    Usage:

      spectrum = create_spectrum(areas, ppms, phases, dataset, ta=0.07, tb=0.07)

    """
    basis, pkppm = create_fid(areas, ppms, phases, dataset, calc_peakppm=False)

    if zforce is None:
        zfmult = dataset.zero_fill_multiplier if not acq else 1.0
    else:
        zfmult = zforce

    acqdim0 = dataset.raw_dims[0]

    # create lineshape
    td    = 1.0/dataset.sw
    off   = round(1.0 * acqdim0 * dataset.echopeak/100.0) * td
    xx    = np.arange(acqdim0) * td
    xxx   = abs(xx-off)
    decay = xx/ta + (xxx/tb)**2
    indx  = np.where(decay > 50)
    cnt   = np.size(indx)
    if cnt > 0:
        decay[indx] = 50
    lshape  = np.exp(-decay)
    
    # create spectrum by applying lineshape and doing iFFT and lineshape 
    tmp = np.zeros(int(round(acqdim0 * zfmult)), complex)
    tmp[0:acqdim0] = basis * lshape
    tmp[0] *= 0.5

    final_spectrum = np.fft.fft(tmp) / len(tmp)

    return final_spectrum


def calculate_peakppm_standalone(fid, sw, freq, zfmult=1, zfnew=4, resppm=4.7):

    acqdim0 = len(fid)
    dim0 = acqdim0 * zfmult * zfnew     # increase zerofill for better accuracy
    hpp = sw / dim0
    apod = safe_exp(-(2.0 * np.pi * 0.6 * (np.arange(acqdim0) / sw)) ** 2, 0)   # 2Hz Gauss
    indx = np.argmax(np.fft.fft(np.pad(fid * apod, (0,dim0-acqdim0))))          # max peak of spectrum
    pkppm = ((dim0 / 2) - indx) * (hpp / freq) + resppm

    return pkppm


#def create_fid_standalone(area, ppm, phase, sw, frequency, dim0, ta=[1e6,], tb=[1e6,], resppm=4.7, ideal=True):
#    """
#
#    bjs-start  Not found in search all vespa
#
#    Create an 'ideal' FID via a sum of individual FIDs from lists of spectral
#    peak descriptors (ie. area, ppm, phase, (opt) ta, (opt) tb values). No
#    global Phase 0 or Phase 1 are applied here. The use of Ta and Tb terms here
#    is generally to allow static values to only be calculated once when the
#    basis set is generated.
#
#    Args:
#        area (ndarray, floats): area under each resonance peak
#        ppm (ndarray, floats): peak center for each resonance peak
#        phase (ndarray, floats): phase 0 for each resonance peak
#        sw (float): sweep width in Hz
#        frequency (float): scanner central frequency in MHz
#        dim0 (int): number of spectral points in FID
#        ta (ndarray, floats): default=1e6, lorentzian decay term, can be a
#            single value in the array that is broadcast for all lines.
#        tb (ndarray, floats): default=1e6, gaussian decay term, can be a
#            single value in the array that is broadcast for all lines.
#        resppm (float): default=4.7, ppm value of the center point of spectrum.
#        ideal (bool): default=True, flag to create FID with (False) or without
#            (True) ta and tb factors included. Assumes line shape added later.
#
#    Usage:
#
#      fid, fids = calculate_fid(area, ppm, phase, sw, frequency, dim0, resppm=4.69, indiv=True)
#
#    """
#    npk  = len(area)
#    hpp  = sw/dim0
#    arr1 = np.zeros(dim0, np.float) + 1.0
#    t = (np.arange(npk * dim0, dtype='float') % dim0) * (1.0 / sw)
#    t.shape = npk, dim0
#
#    # ppm -> hz -> rads
#    freq = ((dim0/2) - (frequency*(ppm-resppm)/hpp)) * hpp * 2.0 * np.pi
#    # deg -> rads -> complex
#    phase = np.exp(1j * phase * np.pi / 180.0)
#
#    # set up lineshape
#    if ideal:
#        lshape = 1.0
#    else:
#        if len(ta) != npk: ta = np.repeat(ta[0], npk)
#        if len(tb) != npk: tb = np.repeat(tb[0], npk)
#        expo = (t / np.outer(ta, arr1)) + (t / np.outer(tb, arr1)) ** 2
#        indx = np.where(expo > 50)
#        if len(indx) > 0: expo[indx] = 50
#        lshape = safe_exp(-expo)
#
#    am = np.outer(np.array(area), arr1)
#    fr = np.exp( 1j * np.outer(freq, arr1) * t)
#    ph = np.outer(phase, arr1)
#
#    fid = np.sum(am * fr * ph * lshape, axis=0)
#
#    return fid


#def create_spectrum_standalone(area, ppm, phase, sw, frequency, dim0, zfmult=1,
#                    echopeak=0.0, ta=[1e6,], tb=[1e6,], resppm=4.7,
#                    ideal=True, acq=False, zforce=None, peakppm=False):
#    """
#    Create a spectrum via a sum of spectral peaks. Each peak has an area, ppm,
#    and phase used to create an 'ideal' FID. The set of summed FIDs is modified
#    by a Voigt lineshape described by Ta and Tb parameters. No global Phase 0
#    or Phase 1 is applied in here.
#
#    Args:
#        area (ndarray, floats): area under each resonance peak
#        ppm (ndarray, floats): peak center for each resonance peak
#        phase (ndarray, floats): phase 0 for each resonance peak
#        sw (float): sweep width in Hz
#        frequency (float): scanner central frequency in MHz
#        dim0 (int): number of spectral points in FID
#        zfmult (int): default=1, zero fill factor
#        echopeak (float): default=0.0, value 0.0-1.0, 0.0 is FID, 0.5 is a
#            half-echo.
#        ta (ndarray, floats): default=1e6, lorentzian decay term, can be a
#            single value in the array that is broadcast for all lines.
#        tb (ndarray, floats): default=1e6, gaussian decay term, can be a
#            single value in the array that is broadcast for all lines.
#        resppm (float): default=4.7, ppm value of the center point of spectrum.
#        ideal (bool): default=True, flag to create FID with (False) or without
#            (True) ta and tb factors included. Assumes line shape added later.
#        acq (bool, optional): default=False, determines if returned spectrum is
#            sized based on current zero filling value, or on acquisition dims.
#        zforce (int, optional): default=None, can be used to force returned
#            spectrum to be zero filled to an arbitrary size based on integer
#            value of the zforce parameter.
#
#    Usage:
#
#      spectrum, pkppm = create_spectrum(areas, ppms, phases, dataset, ta=0.07, tb=0.07)
#
#    """
#    fid = create_fid(area, ppm, phase, sw, frequency, dim0,
#                       ta=ta, tb=tb, resppm=resppm, ideal=ideal)
#
#    pkppm = calculate_peakppm_standalone(fid, sw, frequency, zfmult=zfmult, zfnew=4, resppm=resppm) if peakppm else None
#
#    zfmult = 1 if acq else zfmult
#    zfmult = zforce if zforce is not None else zfmult
#
#    spectrum = np.zeros(int(round(dim0 * zfmult)), np.complex)
#    spectrum[0:dim0] = fid
#    spectrum[0] *= 0.5
#    spectrum = np.fft.fft(spectrum) / len(spectrum)
#
#    return spectrum, pkppm


def chop(data):
    """
    Alternates points between positive and negative
    That is, even points are multiplied by 1
    and odd points are multiplied by -1
    """
    return data * ((((np.arange(len(data)) + 1) % 2) * 2) - 1)


#def clip(value, minimum, maximum=None):
#    """ Returns value such that minimum <= value <= maximum """
#    if value < minimum:
#        value = minimum
#    if (maximum is not None) and (value > maximum):
#        value = maximum
#
#    return value


def full_width_half_max(data, minfloor=False, positive=False):
    """
    Finds maximum in data then searches left and right to find
    the point at half the amplitude of the max point.  Adds
    indices together to get a Full Width Half Max value.
    """

    if positive:
        data = np.abs(data)
    else:
        data = data.copy()

    if minfloor:
        floor = np.min(data)
        floor = floor if floor>0 else 0
    else:
        floor = 0.0


    data = data.real
    peak = np.max(data)
    peak_index = np.argmax(data)
    half_peak = peak - ((peak-floor)/2.0)

    # Find the half max to the right of the peak
    right_index = peak_index
    while (data[right_index] >= half_peak) and (right_index < len(data) - 1):
        right_index += 1

    # Find the half max to the left of the peak
    left_index = peak_index
    while (data[left_index] >= half_peak) and left_index:
        left_index -= 1

    return right_index - left_index


#def get_exponential_filter(length, center, line_broaden, sweep_width):
#    # Get exponential function
#    echo_at = round(length*center)
#    echo_at = np.where(echo_at > 0, echo_at, 0)
#
#    if sweep_width == 0:
#        damper = line_broaden / (length - echo_at)
#    else:
#        damper = np.pi * line_broaden / sweep_width
#
#    index = abs(np.arange(length) - echo_at)
#    expo = index * damper
#
#    return  safe_exp(-expo)
#
#
#def get_gaussian_filter(length, center, line_broaden, sweep_width):
#    # Get Gaussian function
#    echo_at = round(length*center)
#    echo_at = np.where(echo_at > 0, echo_at, 0)
#
#    if sweep_width == 0:
#        npts = float(length - echo_at)
#        damper = line_broaden / npts**2.
#    else:
#        damper = (0.6 * np.pi * line_broaden / sweep_width)**2
#
#    tindex = abs(np.arange(length) - echo_at)
#    expo = tindex**2 * damper
#
#    return safe_exp(-expo)


#def hamming2(N, echo_position):
#    # To obtain Hamming function with value=1.0 at point=N/2+1
#    n1 = int(N * echo_position)
#    n1 = np.where(n1 > 0, n1, 0)
#    n2 = N - n1
#
#    if n1 == 0:
#        filter_length = 2*n2-1
#        ret2 = np.hamming(filter_length + 1)[0:filter_length]    # right half
#        return ret2[N-1:2*N-1]
#    elif n1 == N:
#        filter_length = 2*n1-1
#        ret1 = np.hamming(filter_length + 1)[0:filter_length]    # left half
#        return ret1[0:n1-1]
#    else:
#        filter_length = 2*n1+1
#        ret1 = np.hamming(filter_length + 1)[0:filter_length]    # left half
#        filter_length = 2*n2+1
#        ret2 = np.hamming(filter_length + 1)[0:filter_length]    # right half
#        return [ret1[0:n1],ret2[n2+1:2*n2-1]]


def voigt_width( ta, tb, dataset, hzres=0.1 ):

    if hzres < 0.01:
        hzres = 0.1

    zfmult = dataset.zero_fill_multiplier

    # Calc zerofill needed for < hzres Hz per point
    hpp  = dataset.raw_hpp / zfmult
    mult = 1
    while hpp >= hzres:
       hpp  = hpp  / 2.0
       mult = mult * 2.0

    npts    = dataset.raw_dims[0]
    nptszf1 = npts * zfmult * mult
    nptszf2 = npts * zfmult
    f1      = np.zeros(int(nptszf1), complex)
    f2      = np.zeros(int(nptszf2), complex)

    td      = 1.0/dataset.sw
    off     = round(1.0 * npts * dataset.echopeak/100.0) * td
    t       = np.arange(npts) * td              # for T2, always decreases with time
    ta      = ta if ta > 0.00001 else 0.00001       # Ta > 0
    tb      = tb if tb > 0.00001 else 0.00001       # Tb > 0

    # setup Lineshape
    expo    = t/ta + (t/tb)**2
    lshape  = safe_exp(-expo)
    lshape  = lshape + (1j*lshape*0)

    f1[0:npts] = lshape
    f1[0] *= 0.5
    f1 = np.fft.fft(f1) / nptszf1
    f1 = np.roll(f1, int(nptszf1 / 2))

    f2[0:npts] = lshape
    f2[0] *= 0.5
    f2 = np.fft.fft(f2) / nptszf2

    widpts = full_width_half_max(f1.real)
    peak_width = widpts * dataset.sw / nptszf1

    peak_norm = np.max(f2.real)
    peak_norm = 1.0/peak_norm

    return peak_width, peak_norm
