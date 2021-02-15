"""This module contains a 'pure Python' implementation of the version 2.x 
HLSVDPRO package. 

The HLSVDPROPY package computes a 'sum of lorentzians' model for the complex 
'signals' data passed in via a 'black box' state space approach. We are using 
the scipy.linalg SVD libraries, rather than the PROPACK Fortran libraries used
in the HLSVDPRO package, but the algorithm is otherwise similarly based on the 
algorithm in:

Laudadio T, et.al. "Improved Lanczos algorithms for blackbox MRS data 
quantitation", Journal of Magnetic Resonance, Volume 157, p.292-297, 2002

Functions:
    hlsvd(data, nsv_sought, dwell_time, sparse=False) -> 6-tuple
    hlsvdpro(data, nsv_sought, m=None, sparse=True) -> 8-tuple
    convert_hlsvd_result(result, dwell)
    create_hlsvd_fids(result, npts, dwell, sum_results=False, convert=True)
    get_testdata()
    
Example:
    $ python hlvsd.py
    
    Running this module from the site-packages directory will run the internal 
    example. It requires matplotlib be already installed.

"""      

# Python modules
from __future__ import division
import math

# 3rd party modules
import numpy as np
import scipy.linalg
import scipy.sparse.linalg
import scipy.linalg.lapack as lapack

# Our modules



def hlsvd(data, nsv_sought, dwell_time, sparse=False):
    """
    This calls HLSVDPRO version 2.x code, but simulates the hlsvd.hlsvd()
    call from HLSVDPRO version 1.0.x to maintain the API. See doc string
    below for hlsvdpro() method

    Args:
        data (ndarray): an iterable of complex numbers.

        nsv_sought (int): The number of singular values sought. The function
            will return a maximum of this many singular values.

        dwell_time (float): Dwell time in milliseconds.

    Returns:
        tuple: a 6-tuple containing -
            (int), number of singular values found (nsv_found <= nsv_sought)
            (ndarray, floats) the singular values
            (ndarray, floats) the frequencies (in kilohertz)
            (ndarray, floats) the damping factors (in milliseconds?)
            (ndarray, floats) the amplitudes (in arbitrary units)
            (ndarray, floats) the phases (in degrees)

            Each list's length == nsv_found. The five lists are correlated
            (element N of one list is associated with element N in the other
            four lists) and are sorted by singular value with the largest
            (strongest signal) first. The frequencies and damping factors
            HAVE been adjusted by the dwell time, and phases converted to
            degrees for compatibility with HLSVDPRO version 1.x API.

    """
    m = len(data) // 2
    r = hlsvdpro(data, nsv_sought, m=m, sparse=sparse)
    r = convert_hlsvd_result(r, dwell_time)

    nsv_found, singular_values, frequencies, damping_factors, amplitudes, phases = r[0:6]

    return (nsv_found, singular_values, frequencies, damping_factors, amplitudes, phases)


    
def hlsvdpro(data, nsv_sought, m=None, sparse=False):
    """A pure Python implementation of the HLSVDPRO version 2.x package.
    
    Computes a 'sum of lorentzians' model for the complex 'data' passed in 
    using the scipy.linalg SVD libraries and based on the algorithm in:
      
    Laudadio T, et.al. "Improved Lanczos algorithms for blackbox MRS data 
    quantitation", Journal of Magnetic Resonance, Volume 157, p.292-297, 2002

    Args:
        data (list, ndarray): an iterable of complex numbers.

        nsv_sought (int): The number of singular values sought. The function 
            will return a maximum of this many singular values.

        m (int): (optional) default=len(data)/2, Use to set the size of
            the Hankel matrix used to compute the singular values. Hankel 
            matrix shape is (L+1,L) where L = len(data)-m-1

        sparse (bool): (optional) default True. If set to True, the  
            scipy.sparse.linalg.svds() is used to calculate singular values and 
            nsv_sought is passed in as a parameter. If False, scipy.linalg.svd()
            is used to calculate the singular values, and nsv_sought is used to 
            truncate the results returned.

    Returns:
        tuple: an 8-tuple containing -
            (int), number of singular values found (nsv_found <= nsv_sought)
            (ndarray, floats) the singular values
            (ndarray, floats) the frequencies (in kilohertz)
            (ndarray, floats) the damping factors 
            (ndarray, floats) the amplitudes (in arbitrary units)
            (ndarray, floats) the phases (in radians)
            (ndarray) The top nsv_found left singular vectors, shape=(L+1,3)
            (ndarray) The top nsv_found right singular vectors, shape=(3,L)

            For the first five ndarrays, length == nsv_found and element N of 
            one list is associated with element N in the other four lists. They
            are sorted by singular value with the largest (strongest signal) 
            first.  The frequencies and damping factors have NOT been adjusted 
            by the dwell time, nor the phases converted to degrees.

    """
    mode = 'f' # or 'b'
    xx = data
    k = nsv_sought
    n = len(xx)
    m = int(n/2) if m is None else m
    
    l = n - m - 1

    if mode == "f":
        x = scipy.linalg.hankel(xx[:l + 1], xx[l:])
    else:
        # for backward LP we need to make the hankel matrix:
        #    x_N-1 x_N-2 ... x_N-M-1
        #    x_N-2 x_N-3 ... x_N-M-2
        #      ...
        #    x_M   x_M-1 ... x_0

        x = scipy.linalg.hankel(xx[:m - 1:-1], xx[m::-1])
    
    if sparse:
        u, s, vh = scipy.sparse.linalg.svds(x, k=k)
    else:
        u, s, vh = scipy.linalg.svd(x, full_matrices=False)

    k = min(k,len(s))               # number of singular values found
        
    uk = np.mat(u[:, :k])           # trucated U matrix of rank K
    ub = uk[:-1]                    # Uk with bottom row removed
    ut = uk[1:]                     # Uk with top row removed

    zp, resid, rank, ss = scipy.linalg.lstsq(ub, ut)

    #----------------------------------------------------------------
    # Diagonalize Z' (=hx), yields 'signal'poles' aka 'roots'
    #   Eigenvalues are returned unordered. I sort them here to make
    #   it easier to compare them with output from the Fortran code.

    roots = scipy.linalg.eigvals(zp)
    roots = np.array(sorted(roots, reverse=True))

    #----------------------------------------------------------------
    #  Calculate dampings (damp) and frequencies (freq) from roots

    dampings = np.log(np.abs(roots))
    frequencies = np.arctan2(roots.imag, roots.real) / (math.pi * 2)

    #----------------------------------------------------------------
    #  Calculate complex-valued amplitudes, using the pseudoinverse
    #    of the Lrow*kfit Vandermonde matrix zeta.
    #
    # FIXME this code, like the Fortran, ignores possible errors
    #    reported in the "info" return value.

    zeta = np.vander(roots, N=len(data), increasing=True).T
    v1, x1, s1, rank, _, info = lapack.zgelss(zeta, data)

    #----------------------------------------------------------------
    # Discard uneeded values of x
    x1 = x1[:k]
    amplitudes = np.abs(x1)
    phases = np.arctan2(x1.imag, x1.real)

    return k, s[::-1], frequencies, dampings, amplitudes, phases, u, vh



def convert_hlsvd_result(result, dwell):
    """
    Use dwell time to convert output from hlsvdpro() method to standard units.
    Frequencies convert to [kHz] and damping factors to [ms]. Phases convert
    to [degrees]. Singular values and row and column matrices are maintained
    at their same values and output tuple locations.

    """
    nsv_found, sing_vals, frequencies, damping, ampls, phases = result[0:6]

    damping = np.array([1.0 / df for df in damping])
    damping = np.array([df * dwell for df in damping])
    frequencies = np.array([frequency / dwell for frequency in frequencies])
    phases = np.array([phase * (180.0/3.1415926) for phase in phases])

    return (nsv_found, sing_vals, frequencies, damping, ampls, phases, result[6], result[7])



def create_hlsvd_fids(result, npts, dwell, sum_results=False, convert=True):
    """ 
    This is a convenience function for creating time domain FID arrays from
    the results of the hlsvd() method. You can either send results in 
    directly from hlsvd() (the assumed condition) or run the conversion
    routine convert_hlsvd_results() manually yourself and send them in, but
    in the latter case set convert=False.
    
    See the self-test example at the bottom of this module for how to
    process the fids into spectra and compare in matplotlib.

    Args:
        result (tuple): output from hlsvd() method, an 8-tuple.

        npts (int): The number of points in the created fids.

        dwell (float): Dwell time in milliseconds for created fids.

        sum_results (bool): (optional) default False. If True, all fid lines
            will be summed up to one fid before being returned. If false, a
            (len(result[2]),npts) ndarray is returned

        convert (bool): (optional) default True. If True, the results parameter
            will be run through the convert_hlsvd_results() method before
            being used. If False, it will assume that the user has already
            converted the freq, damp and phase values.

    Returns:
        ndarray: either a (npts,) ndarray (if sum_results=True) or a
            (len(result[2]),npts) ndarray (if sum_results=False) with time
            domain fids.

    """
    if convert: result = convert_hlsvd_result(result, dwell)

    freqs, damps, areas, phase = result[2:6]

    fids = np.zeros((len(freqs), npts), dtype=np.complex)
    t = np.arange(npts) * dwell
    k = 1j * 2 * np.pi

    for i, damp in enumerate(damps):
        if damp:
            # hack for very small exp() values, temp disable numpy's error
            # report, it silently generate NaNs, then change those NaNs to 0.0
            old_settings = np.seterr(all='ignore')
            line = areas[i] * np.exp((t / damp) + k * (freqs[i] * t + phase[i] / 360.0))
            zeros = np.zeros_like(line)
            fids[i,:] = np.where(np.isnan(line), zeros, line)
            np.seterr(**old_settings)
        else:
            fids[i,:] = fids[i,:] * 0

    if sum_results: result = np.sum(fids, axis=0)

    return result


def get_testdata():
    """
    This is a convenience function for accessing an internal test data array. 
    It is part of a dict of information called TESTDATA. However, the 'data'
    entry is encoded. This method decodes that data into a 1024 point numpy
    array of complex128 data.
    
    See the self-test example at the bottom of this module for and example of
    how to process this data into a spectrum and fitted data and compare in 
    matplotlib.

    Args:
        None

    Returns:
        ndarray: 1024 complex128 numbers representative of a short TE single
            voxel time domain FID data set.
    """
    import io
    import zlib
    import base64

    return np.load(io.BytesIO(zlib.decompress(base64.b64decode(TESTDATA['data']))))

#------------------------------------------------------------------------------
# test and helper functions below

def _example():
    """
    This code can be run at the command line to test the internal methods for
    a real MRS SVS data set by running the following command:

        $ python hlsvd.py

    """
    data = get_testdata()
    indat = TESTDATA
    dwell = indat['step_size']
    
    sing0 = np.array(indat['sing0'])
    freq0 = np.array(indat['freq0'])
    ampl0 = np.array(indat['ampl0'])
    damp0 = np.array(indat['damp0'])
    phas0 = np.array(indat['phas0'])
    k = indat['n_singular_values']

    r = hlsvdpro(data, k, m=None, sparse=True)
    c = convert_hlsvd_result(r, dwell)

    nsv_found, sigma, freq, damp, ampl, phas, u, v = c

    udot = abs(np.dot(u.conjugate().T, u))
    vdot = abs(np.dot(v, v.conjugate().T))
    udot[udot < 1e-6] = 0.0
    vdot[vdot < 1e-3] = 0.0

    # print the results -------------------------------------------------------

    np.set_printoptions(suppress=True, precision=6)

    print('Singular Values hlsvdpropy  = ', sigma)
    print('Singular Values actual      = ', sing0[:k])
    print('')
    print('SingVal Diffs (hlsvdpropy vs actual ) = ', sigma - sing0[:k])
    print('Freqs   Diffs (hlsvdpropy vs actual ) = ', freq - freq0[:k])
    print('Damps   Diffs (hlsvdpropy vs actual ) = ', damp - damp0[:k])
    print('Ampls   Diffs (hlsvdpropy vs actual ) = ', ampl - ampl0[:k])
    print('Phase   Diffs (hlsvdpropy vs actual ) = ', phas - phas0[:k])
    print('')
    print('  np.clip(abs(np.dot(u.conj().T, u))) = \n',udot,'\n')
    print('  np.clip(abs(np.dot(v, v.conj().T))) = \n',vdot,'\n')

    # plot the results --------------------------------------------------------

    import matplotlib.pyplot as plt

    lines = create_hlsvd_fids(r, len(data), indat['step_size'], sum_results=True)
    datt = data.flatten()
    fitt = lines.flatten()

    plt.subplot(1, 2, 1)
    plt.plot(datt.real, color='r')
    plt.plot(fitt.real, color='b')
    plt.plot((datt-fitt).real, color='g')
    plt.title('HLSVDProPy fit overlay Data and Residual [Time]')

    datt[0] *= 0.5
    fitt[0] *= 0.5

    plt.subplot(1, 2, 2)
    plt.plot(np.fft.fft(_chop((datt))).real, color='r')
    plt.plot(np.fft.fft(_chop((fitt))).real, color='b')
    plt.plot(np.fft.fft(_chop((datt-fitt))).real, color='g')
    plt.title('HLSVDProPy fit overlay Data and Residual [Freq]')

    plt.show()
    

def _chop(data):
    return data * ((((np.arange(len(data)) + 1) % 2) * 2) - 1)

TESTDATA = {"step_size": 0.256,
            "phase": -115,
            "data": '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',
            "n_singular_values": 20,
            "sing0": [87694.18789057, 25020.31327661, 22847.44495583, 14031.88636191, 12594.34744491, 10820.16406146, 7169.92483586,  5507.71838323, 3691.66743458,  3354.61392512,  3109.67596131, 2435.45510563, 2328.06467652, 1933.81716136, 1811.50599303, 1649.60571037, 1487.51138929, 1340.96036432, 1327.15685121, 1203.24821706],
            "freq0": [3.82828265e-04, 3.94339386e-02, 3.60897262e-03, 6.45623490e-02, 4.85209932e-02, 7.26164973e-02, 8.97837749e-02,  9.41598483e-02, 1.05887952e-01,  1.30243574e-01, -1.34508469e-04, 1.41866902e-01, 5.92026751e-02, 1.63986620e-01, 1.70866772e-01, 1.54505847e-01, 2.10843967e-01, 2.43686511e-01, 2.59173617e-01, -1.70390382e-01],
            "ampl0": [756.50689836, 12.19490399, 492.50851247, 6.55643675, 64.46771301, 35.74552336, 6.60024637, 84.10198358, 82.53885535,  34.00806486, 763.3322855 ,   6.67073619, 365.81338303, 6.96895208, 140.63777114, 230.4402509, 146.46052794, 101.6278103, 9.89415129, 132.60918895],
            "damp0": [-79.51577959, -173.4202973, -55.07720671, -226.14709214, -68.24682132, -134.9270493, -157.77043494, -88.23514852, -94.40561606,  -70.92586358,   -9.77240489, -259.41737822, -11.05556204, -285.85730821, -96.05003221, -12.46085305, -10.2556549, -22.12471209, -84.12678776, -3.9940325],
            "phas0": [-57.32433736, 7.507611969999999, 39.53854928000001, -251.60476041, 4.410037430000003, 3.772625689999998, -172.18178104, -1.6056844400000045, 5.811510560000002, 17.188066969999994, 32.55942451999999, -8.932631229999998, 15.669816249999997, -24.45394605, -1.7509222700000038, 14.432275589999989, -27.217299890000007, -9.31918906, -34.798774910000006, -20.57325659],
            }


if __name__ == "__main__":
    _example()

