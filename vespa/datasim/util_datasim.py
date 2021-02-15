# Python modules

# 3rd party modules
import numpy as np

# Our modules
import vespa.common.util.math_ as util_math


def calculate_linewidth_old(ta, tb, npts=2000):
    """
    creates an on resonance FID with the Voigt (Ta/Tb) decay envelope
    specified, then FFTs to spectral domain and finds full width at
    half max linewidth of the single peak.

    """
    xx = np.arange(npts) / 100.0
    val = -(xx / ta + (xx / tb) ** 2)
    lineshape = util_math.safe_exp(val)
    peak = np.fft.fft(lineshape) / len(lineshape)
    peak /= peak[0]

    # Calculate FWHM
    for i in range(round(npts * 0.5)):
        if peak[i] < 0.5:
            return i / 10.0
    return i / 10.0


def calc_lw_fft(ta, tb, max_lw=4096, min_res=0.1):
    """
    create on resonance FID with Voigt (Ta/Tb) decay envelope, FFT and measure
    full width at half max linewidth of the single peak.

    max_lw in [Hz] is used to set sweep width
    min_res in [Hz[ is used to find next power of 2 > max_lw / mix_res

    If actual Linewidth is wider than max_lw, then max_lw is returned

    NB. bjs, Found that it is important to have a wide max_lw, so that the
      wings don't wrap around (for wide peaks) and shift the 0.5 value that
      the FWHM search is looking for. Typically will get a wider FWHM if
      max_lw < 32768 or so ... :|

    """
    sw = max_lw
    npts_thresh = max_lw / min_res
    npts = 2
    while npts <= npts_thresh: npts *=2
    hpp = sw / npts

    xx = np.arange(npts) / sw
    val = -(xx / ta + (xx / tb) ** 2)
    lshape = util_math.safe_exp(val)
    peak = np.fft.fft(lshape) / npts
    peak /= peak[0]

    # Calculate FWHM - this measures half the peak width
    indx = np.where(peak <= 0.5)
    fwhm = max_lw if indx[0].size==0 else indx[0][0] * hpp * 2

    return fwhm


def calc_lw(ta, tb, max_lw=4096, min_res=0.1):
    """
    create on resonance FID with Voigt (Ta/Tb) decay envelope, FFT and measure
    full width at half max linewidth of the single peak.

    uses 'closed form' for calc ref 1996 Grivet

    """


    # Voigt from Grivet 1996
    #
    # Gauss FWHM
    # s(t) = exp(-at**2)
    #   so in Vespa terms of Tb -(t/Tb)**2 -> a = 1/Tb**2
    #
    # Hg = (2/pi)*sqrt(a * ln(2)) = (2/pi)*sqrt( (1/Tb**2) * ln(2))
    #
    # Lorentz FWHM
    # s(t) = exp(-t/T2*)
    #
    # Hz = 1/(pi * T2*) -> Vespa terms = 1.0/(np.pi*Ta)
    #
    # Voigt FWHM
    #
    # Hv = 0.5 * [Hz + sqrt(Hz*Hz + 4 * Hg*Hg)]

    Hz = 1.0/(np.pi*ta)
    Hg = (2.0 / np.pi) * np.sqrt((1.0/tb**2) * np.log(2))      # from grivet

    fwhm_grivet = 0.5 * (Hz + np.sqrt(Hz*Hz + 4 * Hg*Hg))


    return fwhm_grivet


def make_prior_basis(amp0, ppm0, pha0, datasim=None):
    """
    Sums all the ASCII lines for a given metabolite into a complex FID at the
    resolution of the current values for the data in the FINFO structure.
    Outputs a complex array containing an FID represetation of the summed lines 
    of the metabolite, created according to the SW, AcqPts, etc values currently in
    the FINFO structure.

    ppm0 in [ppm]
    pha0 in [deg] (not radians)

    """
    sw = datasim.sw
    dim0 = datasim.dims[0]

    td   = 1.0 / sw
    arr1 = np.zeros(dim0) + 1.0
    fre  = datasim.ppm2hz(ppm0) * 2.0 * np.pi
    pha  = pha0 * np.pi / 180.0
    npk  = np.size(ppm0)
    xx   = (np.arange(dim0*npk).reshape(npk,dim0) % dim0) * td
    
    am   = np.outer(amp0,arr1)
    fr   = np.outer(fre,arr1)
    ph   = np.outer(pha,arr1)
    
    val1 = 1j * fr * xx
    val2 = 1j * ph
    tmp  = am * np.exp(val1) * np.exp(val2)
    
    if np.ndim(tmp) > 1:
        tmp = np.sum(tmp, axis=0)

    return tmp
    


def generic_spectrum(amp0, ppm0, pha0,  ta=0.07,
                                        tb=0.07,
                                        datasim=None):
    """
    Creates a spectrum from input parameters

    ppm0 in [ppm]
    pha0 in [deg] (not radians)
    ta in [sec]
    tb in [sec]

    """
    sw = datasim.sw
    dim0 = datasim.dims[0]

    basis = make_prior_basis(amp0, ppm0, pha0, datasim=datasim)
    
    # Create Lineshape ---
    td = 1.0 / sw
    off = dim0 * 0 / td             # FIXME bjs, check this divide by 'td'... replace 0 with echopk, if necessary
    xx = np.arange(dim0) * td
    xxx = abs(xx-off)
    
    decay = xx/ta + (xxx/tb)**2
    lshape = util_math.safe_exp(-decay)
    
    # Create Spectrum with FFT and LineShape ---
    tmp = np.zeros(dim0, complex)
    tmp[0:dim0] = basis * lshape
    
    return tmp

