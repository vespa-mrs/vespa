# Python modules

# 3rd party modules
import numpy as np

# Our modules
import vespa.analysis.algos.lowess as lowess
import vespa.common.constants             as common_constants
import vespa.common.minf_parabolic_info   as minf
import vespa.common.util.ppm              as util_ppm
import vespa.common.util.generic_spectral as util_spectral

from vespa.analysis.algos.b0_correction  import b0_correction
from vespa.analysis.algos.auto_correlate import auto_correlate
from vespa.analysis.algos.savitzky_golay import savitzky_golay
from vespa.analysis.algos.optimize_phase import optimize_phase0_correlation
from vespa.analysis.algos.optimize_phase import optimize_phase1_correlation
from vespa.analysis.algos.optimize_phase import optimize_phase0_integration
from vespa.analysis.algos.optimize_phase import optimize_phase1_integration

from vespa.analysis.constants import FitLineshapeModel
from vespa.analysis.constants import FitInitialSmallPeakFreqs
from vespa.analysis.constants import FitInitialSmallPeakAreas
from vespa.analysis.constants import FitInitialBaselineMethod
from vespa.analysis.constants import FitInitialPhaseMethod
from vespa.analysis.constants import FitOptimizeWeightsMethod
from vespa.analysis.constants import FitInitialB0ShiftMethod
from vespa.analysis.constants import FitInitialLinewidthMethod
from vespa.analysis.constants import FitMacromoleculeMethod
from vespa.analysis.constants import FitMacromoleculeMethodInitVal

DTOR = common_constants.DEGREES_TO_RADIANS
RTOD = common_constants.RADIANS_TO_DEGREES




#########################################################################
#
#                           Public Methods
#
#########################################################################

def find_initial_values(chain):
    """
    General call to determine initial values for optimization parameters.  Based
    on info in the optimization control structure chain (typically TE and field 
    strength), a procedure call is made to calculate these values.

    chain: ptr to optimization control structure

     Allows differents assumptions/procedures to be used to set
     the starting values for the ampl/lw-t2s/freq of the DB line
     models based on the field, nucleus, acqseq, and whatever
     else we come up with

    """
    block   = chain._block
    set     = chain._block.set
    dataset = chain._dataset
    
    prior   = set.prior
    metinfo = dataset.user_prior.metinfo

    if dataset.frequency < 250.0:
        field = "low"
        te = "short" if (dataset.seqte < 90.0) else "long"

    else:
        field = "high" 
        te = ""

    _initvals_1h(chain, field, te)


def set_weight_array(chain, lwidth=None, wtmult=None, wtmax=None):
    """
    Creates a weight array to be used in the optimization based on setting in
    the chain structure.  

    chain: ptr to optimization control structure

    """

    block   = chain._block
    set     = chain._block.set
    dataset = chain._dataset
    
    prior   = set.prior
    metinfo = dataset.user_prior.metinfo

    abbr = [metinfo.get_abbreviation(item) for item in set.prior_list]

    dim0    = dataset.spectral_dims[0]

    if not lwidth:
        lwidth = set.initial_linewidth_value
    else: 
        lwidth = float(lwidth) if lwidth > 0.1 else 0.1

    if not wtmult:
        wtmult = set.optimize_weights_width_factor  
    else: 
        wtmult = float(wtmult) if wtmult>0.0001 else 0.001

    if not wtmax:
        wtmax  = dim0-1  
    else: 
        wtmax  = float(wtmax)

    wtarr = np.zeros(dim0, float)

    if set.optimize_weights_method == FitOptimizeWeightsMethod.EVEN_WEIGHTING:
        wtarr = wtarr + 1.0  

    elif set.optimize_weights_method == FitOptimizeWeightsMethod.LOCAL_WEIGHTING:

        lw = lwidth / dataset.spectral_hpp   # in points
        centers = chain.peakpts

        wid = lw * wtmult
        wid = wid if lw<wtmax else wtmax

        for ctr in chain.peakpts:
            cs  = int(np.where(round(ctr-wid)>0, round(ctr-wid), 0))
            cs  = int(np.where(cs<dim0, cs, dim0))
            ce  = int(np.where(round(ctr+wid)>0, round(ctr+wid), 0))
            ce  = int(np.where(ce<dim0, ce, dim0))
            wtarr[cs:ce] = 1.0  

        # set small pk weight scale higher if needed len(chain.peakpts)
        if set.optimize_weights_small_peak_factor != 1.0:

            ws = np.clip(int(np.round(dataset.ppm2pts(14.0))),0,dim0)
            we = np.clip(int(np.round(dataset.ppm2pts(1.25, dataset))),0,dim0)
            wtarr[ws:we] = wtarr[ws:we] * set.optimize_weights_small_peak_factor

            if 'lac' in abbr:
                ws = np.clip(int(np.round(dataset.ppm2pts(1.45))),0,dim0)
                we = np.clip(int(np.round(dataset.ppm2pts(1.25))),0,dim0)
                wtarr[ws:we] = 1.0  

            if 'naa' in abbr:
                ws = np.clip(int(np.round(dataset.ppm2pts(2.12))),0,dim0)
                we = np.clip(int(np.round(dataset.ppm2pts(1.85))),0,dim0)
                wtarr[ws:we] = 1.0  

            if 'cr' in abbr or 'cho' in abbr:
                ws = np.clip(int(np.round(dataset.ppm2pts(3.30))),0,dim0)
                we = np.clip(int(np.round(dataset.ppm2pts(2.85))),0,dim0)
                wtarr[ws:we] = 1.0  

        # Set and filter the weights
        indx0 = np.where(wtarr == 0.0)[0]
        if np.size(indx0) != 0: 
            wtarr[indx0] = 1.0 / set.optimize_weights_scale_factor        

        # set pks in water suppression low
        if set.optimize_weights_water_flag:
            ws = np.clip(int(np.round(dataset.ppm2pts(set.optimize_weights_water_end))),0,dim0)
            we = np.clip(int(np.round(dataset.ppm2pts(set.optimize_weights_water_start))),0,dim0)
            wtarr[ws:we] = 1.0 / set.optimize_weights_scale_factor

        # set pks in lipid area low
        if set.optimize_weights_lipid_flag == 1:
            ws = np.clip(int(np.round(dataset.ppm2pts(set.optimize_weights_lipid_end))),0,dim0)
            we = np.clip(int(np.round(dataset.ppm2pts(set.optimize_weights_lipid_start))),0,dim0)
            wtarr[ws:we] = 1.0 / set.optimize_weights_scale_factor

        wtarr = wtarr / max(wtarr)

    return wtarr


#########################################################################
#
#                   Private/Internal Use Only Methods
#
#########################################################################

def _apply_phase01(data, ph0, ph1, piv, chain):
    """
    Applies zero and first order phase to a given complex array

    INPUT:

     data: complex array
     ph0: zero order phase to apply to data, in degrees
     ph1: first order phase to apply to data, in degrees
     piv: first order phase pivot, in ppm units
     chain: ptr to optimization control structure

    KEYWORDS:

     phase: returns the complex array with zero/first order phase terms

    ph0 and ph1 in Degrees
    piv in ppm
    data must be complex

    """

    block   = chain._block
    set     = chain._block.set
    dataset = chain._dataset
    
    prior   = set.prior
    metinfo = dataset.user_prior.metinfo

    dsiz = data.shape
    pts  = data.shape[0]

    if 'complex' not in data.dtype.name:
        return data

    if ph0 == 0.0 and ph1 == 0.0:
        return data, 1

    pivot  = dataset.ppm2pts(piv)
    phase0 = ph0 * DTOR
    phase1 = ph1 * DTOR * (np.arange(pts)-pivot)/pts
    phase  = np.exp(1j*(phase0+phase1))

    #FIXME bjs - this code could benefit from numpy broadcasting?

    if data.ndim == 1:
        data = phase*data

    elif data.ndim == 2:
        for i in np.arange(data.shape[1]):
            data[:,i] = phase*data[:,i]

    elif data.ndim == 3:
        for j in np.arange(data.shape[2]):
            for i in np.arange(data.shape[1]):
                data[:,i,j] = phase*data[:,i,j]

    elif data.ndim == 4:
        for k in np.arange(data.shape[3]):
            for j in np.arange(data.shape[2]):
                for i in np.arange(data.shape[1]):
                    data[:,i,j,k] = phase*data[:,i,j,k]
    return data, phase



def _baseline_estimation_lowess(data, chain, widmul=1.0):
    """
    Initial baseline smoothing/estimate via Lowess filter wrap to C/Fortran DLL call

    INPUT:
     dat:   array, data
     chain:  ptr to optimization control structure

    KEYWORDS:
     widmul: linewidth multiplier

    """
    block   = chain._block
    set     = chain._block.set
    dataset = chain._dataset
    
    prior   = set.prior
    metinfo = dataset.user_prior.metinfo

    dat = data.copy()

    if dataset.nucleus == '1H': #and np.round(chain.frequency) == 64:
        # Special rule if we have Cho/Cre
        # - keeps calcs from being made from Cho/Cr peak values
        cho = int(dataset.ppm2pts(3.35))
        cre = int(dataset.ppm2pts(2.85))

    width  = set.initial_baseline_lowess_width
    delta  = set.initial_baseline_lowess_delta

    if set.prior_list:
        # De-emphasize regions containing known peaks
        # pkloc - the location of the known peaks (in pts) in data array
        # pkwid - the width about the known peaks to de-emphasize

        prior_peak_ppm = np.array(set.prior_peak_ppm)    # these are in ppm here
        pkloc = dataset.ppm2pts(prior_peak_ppm) # these are in pts here

        pkwid = int(round(set.initial_baseline_lowess_ignore_width * widmul / dataset.spectral_hpp))

#        # new pkloc source - auto_prior ppm array ... if any
#        if set.initial_baseline_lowess_peak_source == 1:
#    
#            pkloc = nd.array(list(dataset.user_prior.values_ppm))
#            pkloc = np.round(dataset.ppm2pts(pkloc))

        npts = len(dat)
        npks = len(pkloc)

        for i in range(npks):

            ym  = int(pkloc[i])
            ys  = int(round((ym - int(pkwid/2))))
            ye  = int(round((ym + int(pkwid/2))))
            ys  = ys if ys>0 else 0
            ye  = ye if ye<npts else npts
            ny  = ye - ys 
            tmp = np.arange(ny)

            lhs = ys - 1
            rhs = ye + 1
            lhs = lhs if lhs>0 else 0
            rhs = rhs if rhs<npts else npts

            if dataset.nucleus == '1H': #AND ROUND(dataset.frequency)) EQ 64:

                # Special rule if we have Cho/Cre
                # - keeps calcs from being made from Cho/Cr peak values

                if (lhs < cre) and (lhs > cho):
                    lhs = cho
                if (rhs < cre) and (rhs > cho):
                    rhs = cre

            lhs1 = lhs - 2
            lhs1 = lhs1 if lhs1>0 else 0
            rhs1 = rhs+2
            rhs1 = rhs1 if rhs1<npts else npts

            ms  = np.median(dat[lhs1:lhs])
            me  = np.median(dat[rhs:rhs1])

            a   = (me-ms)/(ny+4)

            tmp = tmp*a + ms

            dat[ys:ye] = tmp


        if dataset.nucleus == '1H' and np.round(dataset.frequency) > 250:
            # Special rule if we have high field
            lhs = int(np.round(dataset.ppm2pts(4.12)))
            rhs = int(np.round(dataset.ppm2pts(3.4)))

            lhs = lhs - 1
            rhs = rhs + 1
            lhs = lhs if lhs>0 else 0
            rhs = rhs if rhs<npts else npts

            ny  = rhs - lhs + 1
            tmp = np.arange(ny)

            lhs1 = lhs - 2
            lhs1 = lhs1 if lhs1>0 else 0
            rhs1 = rhs+2
            rhs1 = rhs1 if rhs1<npts else npts

            ms  = np.median(dat[lhs1:lhs])
            me  = np.median(dat[rhs:rhs1])
            a   = (me-ms)/(ny+4)
            tmp = tmp*a + ms

            dat[lhs:rhs] = tmp

    width = width / dataset.sw            # convert from Hz to % points
    base  = lowess.lowess(dat, frac=width, delta=delta)

    return base


def _baseline_estimation_savgol(data, chain, widmul=1.0):
    """
    Initial baseline smoothing/estimate via Savitzgy-Golay filter 

    INPUT:
     dat:   array, data
     chain:  ptr to optimization control structure

    KEYWORDS:
     widmul: linewidth multiplier

    """
    block   = chain._block
    set     = chain._block.set
    dataset = chain._dataset
    
    prior   = set.prior
    metinfo = dataset.user_prior.metinfo

    dat = data.copy()

    if dataset.nucleus == '1H': #and np.round(dataset.frequency) == 64:
        # Special rule if we have Cho/Cre
        # - keeps calcs from being made from Cho/Cr peak values
        cho = int(dataset.ppm2pts(3.35))
        cre = int(dataset.ppm2pts(2.85))

    width  = set.initial_baseline_lowess_width
    delta  = set.initial_baseline_lowess_delta

    if set.prior_list:
        # De-emphasize regions containing known peaks
        # pkloc - the location of the known peaks (in pts) in data array
        # pkwid - the width about the known peaks to de-emphasize

        pkloc = np.array(set.prior_peak_ppm)    # these are in ppm here
        pkloc = dataset.ppm2pts(pkloc) # these are in pts here

        pkwid = int(round(set.initial_baseline_lowess_ignore_width * widmul / dataset.spectral_hpp))

#        # new pkloc source - auto_prior ppm array ... if any
#        if set.initial_baseline_lowess_peak_source == 1:
#    
#            pkloc = nd.array(list(dataset.user_prior.values_ppm))
#            pkloc = np.round(dataset.ppm2pts(pkloc))

        npts = len(dat)
        npks = len(pkloc)

        for i in range(npks):

            ym  = int(pkloc[i])
            ys  = int(round((ym - int(pkwid/2))))
            ye  = int(round((ym + int(pkwid/2))))
            ys  = ys if ys>0 else 0
            ye  = ye if ye<npts else npts
            ny  = ye - ys 
            tmp = np.arange(ny)

            lhs = ys - 1
            rhs = ye + 1
            lhs = lhs if lhs>0 else 0
            rhs = rhs if rhs<npts else npts

            if dataset.nucleus == '1H': #AND ROUND(dataset.frequency)) EQ 64:

                # Special rule if we have Cho/Cre
                # - keeps calcs from being made from Cho/Cr peak values

                if (lhs < cre) and (lhs > cho):
                    lhs = cho
                if (rhs < cre) and (rhs > cho):
                    rhs = cre

            lhs1 = lhs - 2
            lhs1 = lhs1 if lhs1>0 else 0
            rhs1 = rhs+2
            rhs1 = rhs1 if rhs1<npts else npts
            ms   = np.median(dat[lhs1:lhs])
            me   = np.median(dat[rhs:rhs1])
            a    = (me-ms)/(ny+4)
            tmp  = tmp*a + ms
            dat[ys:ye] = tmp

        if dataset.nucleus == '1H' and np.round(dataset.frequency) > 250:
            # Special rule if we have high field
            #
            # This draws a straight line from one side of this region to the
            # other because there is such a mix of narrow and broad peaks
            lhs = int(np.round(dataset.ppm2pts(4.21)))
            rhs = int(np.round(dataset.ppm2pts(3.325)))

            lhs = lhs - 1
            rhs = rhs + 1
            lhs = lhs if lhs>0 else 0
            rhs = rhs if rhs<npts else npts

            ny  = rhs - lhs
            tmp = np.arange(ny)

            lhs1 = lhs - 2
            lhs1 = lhs1 if lhs1>0 else 0
            rhs1 = rhs+2
            rhs1 = rhs1 if rhs1<npts else npts

            ms  = np.median(dat[lhs1:lhs])
            me  = np.median(dat[rhs:rhs1])
            a   = (me-ms)/(ny+4)
            tmp = tmp*a + ms

            dat[lhs:rhs] = tmp

    width = int(width / dataset.spectral_hpp)            # convert from Hz to % points
    if width % 2 != 1:
        width += 1

    base  = savitzky_golay( dat, width, 2, deriv=0)

    return base    



def _find_lw(dat, ds, de, magnitude=False):
    """
    Calculate the full width half max linewidth for a given region of data.
    Calculation takes the region, auto-correlates it with itself to get a
    correlation plot for the lineshape(s) in the region, and calcs the FWHM
    of the tallest peak.

    dat: complex array of data to search
    ds: starting point
    de: ending point

    magnitude   bool, flag, performs calc only on magnitude data

    """
    # dependencies: util_spectral.full_width_half_max, a_correlate

    if ds > de: ds, de = de, ds

    pts = de-ds+1
    lag = np.arange(pts) - int(pts/2)

    a   = auto_correlate(np.abs(dat[ds:de]), lag)
    lwa = util_spectral.full_width_half_max(a)
    if magnitude:
        return lwa  

    b   = auto_correlate(dat[ds:de].real, lag)
    c   = auto_correlate(dat[ds:de].imag, lag)
    lwb = util_spectral.full_width_half_max(b)
    lwc = util_spectral.full_width_half_max(c)
    lw = min([lwa, lwb, lwc])

    return lw


def _find_lw2(dat0, ds0, de0, chain, minfloor=0, positive=False):
    # FIXME PS - This needs some documentation, and the docstring below
    # ain't it. I can say that the chain param is at least sometimes 
    # a ChainVoigt object.
    """

    dat:
    ds:
    de:
    chain:

    minfloor:
    positive:

    """
    block   = chain._block
    set     = chain._block.set
    dataset = chain._dataset
    
    prior   = set.prior
    metinfo = dataset.user_prior.metinfo

    ds  = ds0 if ds0<de0 else de0
    de  = ds0 if ds0>de0 else de0
    pts = de-ds+1
    if pts < 5:
        ds = (ds-5) if (ds-5)>0 else 0
        de = (de+5) if (de+5)<dataset.spectral_dims[0] else dataset.spectral_dims[0]

    pkloc = np.array(list(dataset.user_prior.basis.ppm))
    pkamp = np.array(list(dataset.user_prior.basis.area))

    pkloc = np.round(dataset.ppm2pts(pkloc))
    pkloc = pkloc.astype(np.int16)

    delta = np.zeros(dataset.spectral_dims[0])
    delta[pkloc] = pkamp
    delta = delta[ds:de]
    delta = delta + 1j*delta*0.0

    iend = (pts-2) if (pts-2)>1 else 1

    bob  = np.convolve(dat0[ds:de],delta[0:iend],'same')

#    # Example of how we might calc a phase for the data from this
#
#    ctr  = np.argmax(np.abs(bob))
#    ph0  = np.angle(bob[ctr])
#    tmp  = dat0[ds:de] * np.exp(1j*(-ph0))

    lw = util_spectral.full_width_half_max(bob.real, minfloor=minfloor, positive=positive)

    return lw



def _find_peakfreq(dat, start, end, base=None, fixed_peak=False, orig=False):

    """
    =========
    Arguments
    =========
    **dat:**   [array][complex] numpy array with data from which to extract peaks
    **start:** [float] start point of search region in PPM
    **end:**   [float] end point of search region in PPM
    **base:**  [keyword][float] if included, this value is subtracted from the peak amplitude estimate

    ===========
    Description
    ===========
    Searches data in a region for a peak, returns subvoxel linear estimate
    of peak location (in pts), and amplitude. Output is two element array,
    peak amplitude and location (in points)

    ======
    Syntax
    ======
    ::

      peak_pts, peak_amp = find_peakfreq(dat, start, end, base=0)

    """
    if not base:
        fixed_peak=False

    npts   = len(dat)
    dat    = abs(dat)
    start  = start if start>0 else 0
    end    = end   if end<npts else npts

    pkmax = np.max(dat[start:end])
    pkmid = start + np.argmax(dat[start:end])              # middle - max peak

    if not fixed_peak:

        pklft = (pkmid - 1) if (pkmid - 1) > 0 else 0        # one pt left
        pkrgt = (pkmid + 1) if (pkmid + 1) < npts else npts  # one pt right

        difl = dat[pkmid] - dat[pklft]
        difr = dat[pkmid] - dat[pkrgt]

        denom = difr if difr > difl else difl
        part = 0.5*(difl - difr)/denom
        freq = pkmid + part

    else:
        freq = pkmid

    peak = np.abs(pkmax - base)

    # if orig set to True, then use data in orig keyword to return a
    # value for pkneg, else use the values in dat to return pkneg value

    pkneg = 0
    if orig is not False:
        if orig[pkmid] < 0: 
            pkneg=1
    else:
        if dat[pkmid] < 0:
            pkneg=1

    return [peak, freq, pkneg]


def _get_metinfo(abbr0, chain, defconc=0.3, defprot=1.0, deffudge=1.0, defwidth=0.17):
    """
    Values saved in tab delineated METABOLITE INFO section of the processing
    parameter file are used to guess initial values.  This procedure finds
    these values from the arrays stored in chain, given the text string
    abbreviation of the metabolite name

    If no values exist, conc = 0.3, prot = 1.0, fudge = 1.0
    Unless DEFCON, DEFPRO, or DEFFUD are set, then these default
    values get returned

    abbr0: string, abbreviation of the metabolite requested
    chain: Ptr to optimization control structure

    DEFCON: default concentration value
    DEFPRO: default proton number
    DEFFUD: default amplitude fudge multiplier value

    OUTPUT:

    conc: concentration, in mM
    prot: number of protons
    fudge: amplitude estimate fudge multiplier


    # Values saved in tab delineated METABOLITE INFO section of the processing
    # parameter file are used to guess initial values.  This procedure finds
    # these values from the arrays stored in chain, given the text string
    # abbreviation of the metabolite name
    #
    # If no values exist, conc = 0.3, prot = 1.0, fudge = 1.0
    # Unless DEFCON, DEFPRO, or DEFFUD are set, then these default
    # values get returned

    """

    block   = chain._block
    set     = chain._block.set
    dataset = chain._dataset
    
    prior   = set.prior
    metinfo = dataset.user_prior.metinfo

    all_abbr = [item.lower() for item in metinfo.abbreviations]
    abbr = abbr0.lower()
    if abbr in all_abbr:
        indx = all_abbr.index(abbr)
        conc  = metinfo.concentrations[indx]
        prot  = metinfo.spins[indx]
    else:
        conc  = defconc
        prot  = defprot

    return conc, prot


def _inital_estimate_baseline(data, chain):
    """
    Generates an estimate of baseline signal using a lowess filtered version of
    the data in which known peak contributions have been de-emphasized. 
    Returns a copy of the data with the base estimate subtracted out

    data: complex array, spectral data
    chain: ptr to optimization control structure

    """
    block   = chain._block
    set     = chain._block.set
    dataset = chain._dataset
    
    prior   = set.prior
    metinfo = dataset.user_prior.metinfo
    
    dat   = data.copy()

    if set.initial_baseline_method == FitInitialBaselineMethod.LOWESS_FILTER:
        base1 = _baseline_estimation_lowess(dat.real, chain, widmul=1.0)
        base2 = _baseline_estimation_lowess(dat.imag, chain, widmul=1.4) # 1.4 for dispersion lineshape

    elif set.initial_baseline_method == FitInitialBaselineMethod.SAVGOL_FILTER:
        base1 = _baseline_estimation_savgol(dat.real, chain, widmul=1.0)
        base2 = _baseline_estimation_savgol(dat.imag, chain, widmul=1.4) # 1.4 for dispersion lineshape

    base = base1+1j*base2
    dat  = data.copy() - base

    chain.init_baseline = base

    return dat


def _inital_estimate_phase(data, chain):
    """
    Estimates zero and or first order phase, using either a correlation or 
    integration auto optimization method.

    data: complex array of data to be phased
    chain: ptr to optimization control structure
    pha: phased copy of data

    """
    block   = chain._block
    set     = chain._block.set
    dataset = chain._dataset
    
    prior   = set.prior
    metinfo = dataset.user_prior.metinfo

    pha  = data.copy()
    ddat = data.copy()

    if set.initial_phase_method == FitInitialPhaseMethod.MANUAL:
        # manually set phase
        phres  = [chain.phase0,chain.phase1]
        phased, phase = _apply_phase01(ddat, 
                                      chain.phase0, 
                                      chain.phase1, 
                                      dataset.user_prior.auto_phase1_pivot,
                                      chain)

    elif set.initial_phase_method == FitInitialPhaseMethod.CORRELATION_PHASE0:
        # maximize correlation of phase 0 of real spectrum to ideal spectrum
        cdeg  = 5
        pstr  = dataset.user_prior.auto_b0_range_start  
        pend  = dataset.user_prior.auto_b0_range_end
        phres, phase = _optimize_phase01_correlation(ddat, 
                                          dataset.user_prior.basis.get_spectrum_sum(dataset),  
                                          chain, 
                                          zero=True,  
                                          one=False,   
                                          b0str=pstr, 
                                          b0end=pend, 
                                          cdeg=cdeg)

    elif set.initial_phase_method == FitInitialPhaseMethod.CORRELATION_PHASE01:
        # maximize correlation of phase 0/1 of real spectrum to ideal spectrum
        cdeg  = 5
        pstr  = dataset.user_prior.auto_b0_range_start
        pend  = dataset.user_prior.auto_b0_range_end
        phres, phase = _optimize_phase01_correlation(ddat, 
                                          dataset.user_prior.basis.get_spectrum_sum(dataset),  
                                          chain, 
                                          zero=True,  
                                          one=True,   
                                          b0str=pstr, 
                                          b0end=pend, 
                                          cdeg=cdeg)

    elif set.initial_phase_method == FitInitialPhaseMethod.INTEGRATION_PHASE0:
        # maximize integral of phase 0 of real spectrum to ideal spectrum 
        phres, phase = _optimize_phase01_integration(ddat, chain, zero=True, one=False)

    elif set.initial_phase_method == FitInitialPhaseMethod.INTEGRATION_PHASE01:
        # maximize integral of phase 0/1 of real spectrum to ideal spectrum 
        phres, phase = _optimize_phase01_integration(ddat, chain, zero=True, one=True)

    return np.array(phres), phase


def _check_array_range(istart, iend, npts):
    """
    Takes two indices, checks that they are in the correct order, ie. start
    is smaller than end, checks that start>=0 and end<=npts, and checks 
    that start!=end.

    """
    istart = int(istart if istart<iend else iend)
    iend   = int(istart if istart>iend else iend)

    istart = istart if istart>0   else 0
    iend   = iend   if iend<=npts else npts
    if istart == iend:    
        # ensure that istart and iend are not the same
        if istart > 0:
            istart = istart-1
        else:
            iend = iend+1

    return istart, iend



def _initvals_1h(chain, field, te):
    """
    Initial value procedure.

    Calculates from the spectral data starting values for Amplitude, Frequency,
    Phase0, Phase1, Ta and Tb (and occasionally, starting baseline estimates).
    Sets these values in the chain control structure.

    chain: ptr to optimization control structure
    field: 'low' (64/128 MHz) or 'high' (above 250 Mhz)
    te: 'long' or 'short' or ''. When field is 'high', te is ignored. 
    """
    block   = chain._block
    set     = chain._block.set
    dataset = chain._dataset
    
    prior   = set.prior
    metinfo = dataset.user_prior.metinfo

    nmet    = len(set.prior_list)

    if (field == "low") and (te == "short"):
        MED_PPM = 0.17
    else:
        MED_PPM = 0.20

    if field == "low":
        BIG_PPM = 0.25
    else:
        BIG_PPM = 0.30


    dat   = chain.data.copy()
    dim0  = dataset.spectral_dims[0]
    pivot = dataset.user_prior.auto_phase1_pivot

    #--------------------------------------------------------
    # Constant phase 1 removal due to FID delay

    constph1 = set.initial_phase1_fid_constant
    if constph1 != 0.0:
        dat, _  = _apply_phase01(dat, 0.0, constph1, pivot, chain)


    #--------------------------------------------------------
    # B0 Shift section 

    if set.initial_b0_shift_method == FitInitialB0ShiftMethod.MANUAL:
        shft = chain.init_b0

    elif set.initial_b0_shift_method == FitInitialB0ShiftMethod.AUTO_CORRELATE:
        ps0  = dataset.ppm2pts(dataset.user_prior.auto_b0_range_start)
        pe0  = dataset.ppm2pts(dataset.user_prior.auto_b0_range_end)
        ps   = int(ps0 if ps0<pe0 else pe0)
        pe   = int(ps0 if ps0>pe0 else pe0)
        ref  = dataset.user_prior_summed_spectrum.flatten()
        ref  = ref[ps:pe]
        bdat = dat[ps:pe]
        shft, _       = b0_correction(bdat,ref)
        chain.init_b0 = shft * dataset.spectral_hpp     # in Hz

    shft = int(round(chain.init_b0 / dataset.spectral_hpp))     # in pts
    dat  = np.roll(dat, shft)


    #---------------------------------------------------------
    # Phase and Baseline estimate section
    #
    # - the phase 0/1 values and the phased data are returned here
    # - this phased data may then have a baseline estimated for it,
    #   and if so, then the stored estimate has to be back phased
    #   for proper display in the plot window

    phres, phase = _inital_estimate_phase(dat.copy(), chain)   # initial phase guess

    ddat = dat * phase
    
    mmol_ddat = ddat.copy()

    chain.init_ph0 = -phres[0] * DTOR
    chain.init_ph1 = -phres[1] - constph1 

    if set.initial_baseline_method:
        # subtract out the baseline estimate
        ddat = _inital_estimate_baseline(ddat, chain)  
        chain.init_baseline, _ = _apply_phase01(chain.init_baseline,
                                             -phres[0],
                                             -phres[1]-constph1,
                                              pivot, 
                                              chain)


    #---------------------------------------------------------
    # Linewidth section

    ptmax = int(np.round(dataset.ppm2pts(set.initial_linewidth_range_start))) # NB the min/max switch here
    ptmin = int(np.round(dataset.ppm2pts(set.initial_linewidth_range_end)))

    if set.initial_linewidth_method == FitInitialLinewidthMethod.MANUAL:
        lw = set.initial_linewidth_value
    elif set.initial_linewidth_method == FitInitialLinewidthMethod.DECONVOLUTION:
        lw = _find_lw( ddat, ptmin, ptmax)
        lw = lw * dataset.spectral_hpp * np.abs(set.initial_linewidth_fudge)
    elif set.initial_linewidth_method == FitInitialLinewidthMethod.AUTO_CORRELATE:
        lw = _find_lw2(ddat, ptmin, ptmax, chain, minfloor=True)   
        lw = lw * dataset.spectral_hpp * np.abs(set.initial_linewidth_fudge)

    set.initial_linewidth_value = lw if lw>2.0 else 2.0

    lw_local  = [set.initial_linewidth_value] * nmet
    lw_scales = [1.0,] * nmet
    
    lw_scaled = [val*scale for val,scale in zip(lw_local,lw_scales)]

    if set.lineshape_model == FitLineshapeModel.VOIGT:
        lw_time = []
        for val in lw_scaled:
            res = _calc_height2area_ratio( val, chain )            # Free LorGauss constraints
            lw_time.append(res[0])
        chain.init_ta = lw_time
        chain.init_tb = lw_time
    elif set.lineshape_model == FitLineshapeModel.LORENTZ:
        lw_time = []
        for val in lw_scaled:
            res = _calc_height2area_ratio( val, chain, tb=1000.0 ) # Pure Lorentz - user sets Gauss, def=1000 sec
            lw_time.append(res[0])
        chain.init_ta = lw_time
        chain.init_tb = [1000.0,] * nmet
    elif set.lineshape_model == FitLineshapeModel.GAUSS:
        # if all fixed T2 values are large then we are doing a 'pure' 
        # gaussian model rather than a physiological Voigt model with T2
        # values that approximate metabolite T2 values. 
        if all(np.array(set.prior_fix_t2) >= 1.0):
            ctr = 1000.0
        else:
            ctr = chain.fix_t2_center    # typically 0.250 [sec]
        
        lw_time = []
        for val in lw_scaled:
            res = _calc_height2area_ratio( val, chain, ta=ctr ) # Pure Gauss - user sets Lorentz, def=0.250 sec
            lw_time.append(res[0])

        # Note. for a physiologic Gaussian lineshape, we allow the user to 
        # set a fixed T2 value for each metabolite. In the model call 
        # (chain_fit_voigt.internal_lorgauss()) we create a lineshape array 
        # that takes each fixed value into account.
        # BUT! we are still passing in a Tb parameter here, and though that 
        # parameter will be tightly constrained, it should still be in the 
        # range of the fixed physiologic params choosen by the user. At the 
        # moment, we have choosen to just set this parameter to 0.250 sec, 
        # which should be a reasonable average of 1H metabolite T2 values 
        # in the brain. It is then allowed to bop plus or minus 0.001 as set 
        # in constraints in the algorithm call below.  
        # In the fitting function, we adjust each fixed value by the delta 
        # from 0.250 so that the pder will fluctuate as the parameter 
        # changes and not confuse the poor optimization routine. In reality,
        # the 0.250 is never used, just the delta amount of the actual a[nmet*2]
        # parameter from that value.
        chain.init_ta = [chain.fix_t2_center,] * nmet
        chain.init_tb = lw_time
    else:
        lw_time = []
        for val in lw_scaled:
            res = _calc_height2area_ratio( val, chain )            # Free LorGauss constraints
            lw_time.append(res[0])
        chain.init_ta = lw_time
        chain.init_tb = lw_time



    #--------------------------------------------------------
    # Peak Areas and PPMs section
    #
    # Set if initial values are taken from real or abs(real) data

    ndat    = len(ddat)
    foff    = 0.0
    chcrflg = 0
    naoff   = 0.0
    croff   = 0.0
    naamp   = None
    cramp   = None
    xsm     = int(np.round(dataset.ppm2pts(0.03, rel=True)))
    small   = int(np.round(dataset.ppm2pts(0.10, rel=True)))
    med     = int(np.round(dataset.ppm2pts(MED_PPM, rel=True)))
    big     = int(np.round(dataset.ppm2pts(BIG_PPM, rel=True)))
    pkppms  = dataset.ppm2pts(np.array(set.prior_peak_ppm))     # these are in pts here

    noramp = res[1] # NB. this comes from Lineshape sections, be careful it does not change
    fampl  = []
    ffreq  = []

    if set.initial_lac_method == 1:
        large_metabs = ['naa','naa+naag','cr','cr2','cr+pcr','pcr+cr','cho','gpc+pcho','pcho+gpc','pcho','gpc','h2o','pk0','oth']
    else:
        large_metabs = ['naa','naa+naag','cr','cr2','cr+pcr','pcr+cr','cho','gpc+pcho','pcho+gpc','pcho','gpc','lac','h2o','pk0','oth']

    if (field == "low") and (te == "short"):
        large_metabs.append('s-ino')

    orig = ddat.copy()

    # 'Anonymous' refernce peaks. In some cases we want to use the full NAA
    #  or CR peak height to estimate starting areas for an assumed concentration
    #  and then set all basis metabs via ratio to that reference area. BUT we
    #  can't use actual NAA or CR measurements made below because these need 
    #  to be set to some smaller value such as when CR = Cr+PCr.  In this case
    #  We use the peak_2(3)ppm_amp variables to separate the ref area from the
    #  measured areas for the basis metabolites.
    #
    # The peak_2(3)ppm_area variables can be calculated in a couple of 
    # locations. First, in the respective NAA or CR calculations if they are a
    # part of the basis set. Second, recalculated in the KO-filter section if
    # that option is selected. Finally, if peak_2(3)ppm_amp is still None going
    # into the ratio calculations, then it is explicitly calculated (same as in
    # the NAA or CR section) at that point. 

    peak_2ppm_amp = None
    peak_3ppm_amp = None


    #---------------------------------------------------
    # Loop 1 set large metabolite values based on peak picking

    for i in range(nmet):

        pkppm  = pkppms[i]
        fudge  = set.prior_area_scale[i]
        width  = int(np.round(chain._dataset.ppm2pts(set.prior_search_ppm[i], rel=True)))
        dbflag = set.prior_db_ppm[i]
        pkph0  = set.prior_search_ph0[i]

            # apply indiv metabolite ph0 prior to peak search,
            # used mostly for UTE data with phase 1 issues
        if pkph0 != 0.0:
            ddat, _  = _apply_phase01(dat, pkph0, 0.0, pivot, chain)

            # convert data to abs(real) if desired.
        ddat = orig.real
        if set.initial_peak_search_abs:
            orig = orig.real
            ddat = np.abs(ddat)

            # theoretical center of pk to be measured
        ctr  = int(round(pkppm))

            # using abbreviations let us clump various metabolite spellings
            # into one algorithm for finding peak and frequency values
        tabbr = metinfo.get_abbreviation(set.prior_list[i].lower())

        if 'oth' in tabbr:
            tabbr = 'oth'
        if 'pk0' in tabbr:
            tabbr = 'pk0'
        if 'glu' in tabbr:
            tabbr = 'glu'  # for glu+gln

        if tabbr in ['cho','gpc+pcho','pcho+gpc','pcho','gpc']:
            conc, prot = _get_metinfo( 'cho', chain, defconc=3.0, defprot=9.0)
            bs  = np.round(ctr - width  )
            be  = np.round(ctr + small)
            bs, be = _check_array_range(bs, be, ndat)
            bas = np.min(ddat[bs:be])
            bas = bas if bas>0 else 0
            chcrflg = chcrflg+1
            pkval = _find_peakfreq(ddat, bs, be, base=bas)
            pkval[0] = pkval[0] * fudge / prot

        elif tabbr in ['cr','cr+pcr','pcr+cr','pcr']: 
            ccr, prot = _get_metinfo( 'cr', chain, defconc=8.0, defprot=3.0)
            bs  = np.round(ctr - small)
            be  = np.round(ctr + width  )
            bs, be = _check_array_range(bs, be, ndat)
            bas = np.min(ddat[bs:be])
            bas = bas if bas>0 else 0
            chcrflg = chcrflg+1
            pkval = _find_peakfreq(ddat, bs, be, base=bas)
            pkval[0] = pkval[0] * fudge / prot
            croff = pkval[1] - ctr
            cramp = pkval[0] * noramp
            
            peak_3ppm_amp = cramp / fudge
            
#            print ' * cr initval area = '+str(cramp)+'    peak_3ppm_amp ='+str(peak_3ppm_amp)

        elif tabbr in ['naa','naa+naag']:
            cnaa, pnaa = _get_metinfo( 'naa', chain, defconc=10.0, defprot=3.0)
            bs  = np.round(ctr - width)
            be  = np.round(ctr + width)
            bs, be = _check_array_range(bs, be, ndat)
            bas = np.min(ddat[bs:be])
            bas = bas if bas>0 else 0
            pkval = _find_peakfreq(ddat, bs, be, base=bas)
            pkval[0] = pkval[0] * fudge / pnaa
            naoff = pkval[1] - ctr
            naamp = pkval[0] * noramp
            
            peak_2ppm_amp = naamp / fudge
            
#            print ' * naa initval area = '+str(naamp)+'    peak_2ppm_amp ='+str(peak_2ppm_amp)

        elif tabbr == 'lac':
            conc, prot = _get_metinfo( 'lac', chain, defconc=2.0, defprot=2.0)
            bs  = np.round(ctr - width)
            be  = np.round(ctr + width)
            bs, be = _check_array_range(bs, be, ndat)
            bas = 0  
            bas = bas if bas>0 else 0
            pkval = _find_peakfreq(np.abs(ddat), bs, be, base=bas, fixed_peak=True)
            pkval[0] = pkval[0] * fudge / prot

        elif tabbr in ['s-ino','sino','sins']:
            conc, prot = _get_metinfo( 'sino', chain, defconc=2.0, defprot=3.0)
            bs  = np.round(ctr - width)
            be  = np.round(ctr + width)
            bs, be = _check_array_range(bs, be, ndat)
            if field == "low":
                bas = np.min(ddat[bs:be]) 
            else:
                bas = np.min(ddat[bs:be])
            bas = bas if bas>0 else 0
            pkval = _find_peakfreq(ddat, bs, be, base=bas)
            pkval[0] = pkval[0] * fudge / prot

        elif tabbr == 'h2o':
            conc, prot = _get_metinfo( 'h2o', chain, defconc=10.0, defprot=2.0)
            bs  = np.round(ctr - width)
            be  = np.round(ctr + width)
            bs, be = _check_array_range(bs, be, ndat)
            bas = np.min(ddat[bs:be])
            bas = bas if bas>0 else 0
            pkval = _find_peakfreq(ddat, bs, be, base=bas, orig=orig)
            pkval[0] = pkval[0] * fudge / prot
            if set.initial_peak_negative_flag != 0:
                if pkval[2] != 0:
                    pkval[0] = -pkval[0]

        elif tabbr == 'oth':
            conc, prot = _get_metinfo(set.prior_list[i].lower(), chain, defconc=1.0, defprot=1.0)
            bs  = np.round(ctr - width)
            be  = np.round(ctr + width)
            bas = 0
            bs, be = _check_array_range(bs, be, ndat)
            bas = bas if bas>0 else 0
            if field == "low":
                orig_param = orig
            else:
                orig_param = None
            pkval = _find_peakfreq(ddat, bs, be, base=bas, orig=orig_param)
            pkval[0] = pkval[0] * fudge / prot
            if set.initial_peak_negative_flag != 0:
                if pkval[2] != 0:
                    pkval[0] = -pkval[0]

        elif tabbr == 'pk0':
            conc, prot = _get_metinfo(set.prior_list[i].lower(), chain, defconc=1.0, defprot=1.0)
            bs  = np.round(ctr - width)
            be  = np.round(ctr + width)
            bas = 0
            bs, be = _check_array_range(bs, be, ndat)
            bas = bas if bas>0 else 0
            if field == "low":
                orig_param = orig
            else:
                orig_param = None
            pkval = _find_peakfreq(ddat, bs, be, base=bas, orig=orig)
            pkval[0] = pkval[0] * fudge / prot
            if set.initial_peak_negative_flag != 0:
                if pkval[2] != 0:
                    pkval[0] = -pkval[0]
        else:
            # there is a chance that the user will de-select n-acetylaspartate
            # in which case the naa-ratio method will not work, so here we
            # find start values using a peak_search method first, and later we
            # can replace them with naa-ratio values if valid

            conc, prot = _get_metinfo( tabbr, chain, defconc=1.0, defprot=1.0)
            bs = np.round(ctr - width)
            be = np.round(ctr + width)
            bs, be = _check_array_range(bs, be, ndat)

            bas = np.min(ddat[bs:be])
            bas = bas if bas>0 else 0
            pkval = _find_peakfreq(ddat, bs, be, base=bas)
            pkval[0] = pkval[0] * fudge / prot
            if tabbr == 'm-ino':
                pkval[1] = ctr
            if set.initial_small_peak_areas == 'peak_search':
                pkval[1] = ctr
#            else:
#                # this is for 'naa_ratio' option
#                pkval = [1.0/noramp,ctr]

        if dbflag:
            # use the database value for PPM if flag checked
            pkval[1] = ctr  

        fampl.append(pkval[0])
        ffreq.append(pkval[1])

        # this is avg of all the pk offsets
        foff = foff + (pkval[1]-ctr) / nmet


    chain.init_area = np.array(fampl) * noramp
    chain.init_freq = np.array(ffreq)

    # here we use a truncation filter to better estimate the NAA peak height
    # so we can set all small peaks via ratio to NAA. As per the paper, we have
    # to include all the metaboltites in the truncation filter iteration to 
    # account for all overlapping signals as we bootstrap to a new NAA value

    if set.initial_apply_ko_filter:

        kodata   = chain.kodata.copy()
        kopoints = set.initial_ko_points 
        kominlw  = set.initial_ko_linewidth_minimum

        dim0    = dataset.spectral_dims[0]
        acqdim0 = dataset.raw_dims[0]
        acqsw   = dataset.sw

        data = np.zeros([dim0],'complex')
        data[0:len(kodata[kopoints:])] = kodata[kopoints:]

        # numpy broadcasting deals with (N,dim0) sized data arrays
        chop_array = ((((np.arange(acqdim0) + 1) % 2) * 2) - 1)
        data[0:acqdim0] = data[0:acqdim0] * chop_array                    # deal with zfill

        data = np.fft.fft(data, n=dim0) / float(dim0)
        data = np.abs(data)

        # create a line shape to apply to basis functions 
        xx = np.arange(acqdim0)/float(acqsw)
        lshape = np.exp(-(set.initial_linewidth_value * np.pi * xx * 0.6 )**2)     # Gaussian

        # get area and frequency estimates
        height = np.zeros([nmet], 'float')
        area0  = np.zeros([nmet], 'float')
        freq0  = np.zeros([nmet], 'float')
        area   = np.zeros([nmet], 'float')
        basis0 = np.zeros([nmet,acqdim0], 'complex')

        # get start peak height and broaden basis set
        for i in range(nmet):

            ctr = int(round(pkppms[i]))
            bs  = np.round(ctr - small)
            be  = np.round(ctr + small)
            bs, be = _check_array_range(bs, be, len(data))
            val = _find_peakfreq(data, bs, be, base=0)

            area0[i] = val[0]
            freq0[i] = val[1]

            basis0[i,:] = chain.basis_mets[i,:] * lshape

        # bootstrap from kofilter[0] to kofilter[N] for each metabolite height
        for j in range(4):

            basis = basis0.copy()
            if j == 0:
                area = area0.copy()

            for k in range(nmet):
                basis[k,:] = area[k] * basis[k,:]

            basis = np.sum(basis,axis=0)

            tmp = np.zeros([dim0],'complex')
            tmp[0:len(basis[kopoints:])] = basis[kopoints:]
            tmp = np.fft.fft(tmp) / float(dim0)
            tmp = np.abs(tmp)

            for k in range(nmet):
                bs  = np.round(freq0[k] - small)
                be  = np.round(freq0[k] + small)
                bs, be = _check_array_range(bs, be, len(tmp))
                val = _find_peakfreq(tmp, bs, be, base=0)

                height[k] = val[0]
                area[k]  *=  area0[k]/val[0]

        # done with estimate, now replace as needed in chain init values  
        # - NB. we only replace NAA ampl value and NAA freq offset
        # - we also adjust the Cr and Cho freq offsets too if in model
        for i in range(nmet):
            tabbr = metinfo.get_abbreviation(set.prior_list[i].lower())
            fudge = set.prior_area_scale[i]
            ctr   = int(round(pkppms[i]))
            if tabbr in ['naa','naa+naag']:
                large_metabs = ['naa','naa+naag']
                orig = chain.init_area[i]
                chain.init_area[i] = area[i] * fudge
                chain.init_freq[i] = freq0[i]
                naoff = freq0[i] - ctr
                naamp = area[i] * fudge
                peak_2ppm_amp = area[i]
                
#                print ' *ko peak_2ppm_amp initval area = '+str(peak_2ppm_amp)
                
#                print "Name = "+tabbr+"  area = "+str(area[i])+"   height = "+str(height[i])+"  peak_pick = "+str(orig)
#                if True:
#                    name = set.prior_list[i]
#                    print "Name = "+name+"  area = "+str(area[i])+"   height = "+str(height[i])+"  peak_pick = "+str(orig)
            elif tabbr in ['cr', 'cr+pcr', 'pcr+cr']:
                if set.initial_small_peak_areas == FitInitialSmallPeakAreas.NAA_RATIO:
                    chain.init_freq[i] = ctr
                elif set.initial_small_peak_areas == FitInitialSmallPeakAreas.CR_RATIO:
                    chain.init_freq[i] = freq0[i]
                elif set.initial_small_peak_areas == FitInitialSmallPeakAreas.MEAN_RATIO:
                    chain.init_freq[i] = chain.init_freq[i]
                cramp = area[i] * fudge
                croff = freq0[i] - ctr
                peak_3ppm_amp = area[i]
                
#                print ' *ko peak_3ppm_amp initval area = '+str(peak_3ppm_amp)
                
            elif tabbr in ['cho', 'gpc+pcho','pcho+gpc', 'pcho']:
                chain.init_freq[i] = ctr


    #---------------------------------------------------
    # Explicit calc of 'anonymous' ref peak - if needed

    #if peak_2ppm_amp is None and set.initial_small_peak_areas == 'peak_2ppm' :
    if set.initial_small_peak_areas == 'peak_2ppm' :
        # may be calculated above in the NAA peak section if in basis
        _, peak_2ppm_prot = _get_metinfo( 'peak_2ppm', chain, defconc=10.0, defprot=3.0)
        ctr = dataset.ppm2pts(2.01)     # in pts here
        bs, be = _check_array_range(np.round(ctr-width), np.round(ctr+width), ndat)
        bas = np.min(ddat[bs:be])
        bas = bas if bas>0 else 0
        pkval = _find_peakfreq(ddat, bs, be, base=bas)
        pkval[0] = pkval[0] / peak_2ppm_prot
        peak_2ppm_amp = pkval[0] * noramp
        
#        print ' *** peak_2ppm_amp initval area = '+str(peak_2ppm_amp)

    if peak_3ppm_amp is None and set.initial_small_peak_areas == 'peak_3ppm' : 
        # may be calculated above in the CR peak section if in basis
        _, peak_3ppm_prot = _get_metinfo( 'peak_3ppm', chain, defconc=8.0, defprot=3.0)
        ctr = dataset.ppm2pts(3.02)     # in pts here
        bs, be = _check_array_range(np.round(ctr-small), np.round(ctr+width), ndat)
        bas = np.min(ddat[bs:be])
        bas = bas if bas>0 else 0
        pkval = _find_peakfreq(ddat, bs, be, base=bas)
        pkval[0] = pkval[0] / peak_3ppm_prot
        peak_3ppm_amp = pkval[0] * noramp
        
#        print ' *** peak_3ppm_amp initval area = '+str(peak_3ppm_amp)

    #---------------------------------------------------
    # Loop 2 sets Small Metab values via NAA conc ratio

    reference_peak_offset = 0.0    
    if set.initial_small_peak_freqs == FitInitialSmallPeakFreqs.REF_PEAK:
        if set.initial_apply_ko_filter:
            reference_peak_offset = naoff       # same value for NAA or NAA+CR
            if set.initial_small_peak_areas == FitInitialSmallPeakAreas.CR_RATIO:
                reference_peak_offset = croff
        else:
            reference_peak_offset = (naoff + croff)/2.0  # general offset - no ko-filter on

    if naamp and (set.initial_small_peak_areas == FitInitialSmallPeakAreas.NAA_RATIO):
        large_metabs = ['naa','naa+naag']
        for i in range(nmet): 

            fudge  = set.prior_area_scale[i]
            metstr = metinfo.get_abbreviation(set.prior_list[i].lower())
            if 'oth'  in metstr: metstr = 'oth'
            if 'pk0'  in metstr: metstr = 'pk0'
            if 'glu'  in metstr: metstr = 'glu'
            if 'sins' in metstr: metstr = 'sino'

            conc, prot = _get_metinfo( metstr, chain, defconc=5.0, defprot=3.0)

            if metstr not in large_metabs:
                mult   = naamp * fudge * (conc/cnaa) 
                chain.init_area[i]  = mult
                chain.init_freq[i] += reference_peak_offset

    elif cramp and (set.initial_small_peak_areas == FitInitialSmallPeakAreas.CR_RATIO):
        large_metabs = ['cr', 'cr+pcr', 'pcr+cr']
        for i in range(nmet): 

            fudge  = set.prior_area_scale[i]
            metstr = metinfo.get_abbreviation(set.prior_list[i].lower())
            if 'oth'  in metstr: metstr = 'oth'
            if 'pk0'  in metstr: metstr = 'pk0'
            if 'glu'  in metstr: metstr = 'glu'
            if 'sins' in metstr: metstr = 'sino'

            conc, prot = _get_metinfo( metstr, chain, defconc=5.0, defprot=3.0)

            if metstr not in large_metabs:
                mult   = cramp * fudge * (conc/ccr)  
                chain.init_area[i]  = mult
                chain.init_freq[i] += reference_peak_offset
                
    elif naamp and cramp and (set.initial_small_peak_areas == FitInitialSmallPeakAreas.MEAN_RATIO):
        large_metabs = ['naa','naa+naag', 'cr', 'cr+pcr', 'pcr+cr']
        
        ratio_data  = naamp/cramp
        ratio_ideal = cnaa/ccr
        data2ideal  = ratio_data / ratio_ideal 
        
        # if data2ideal > 1.0 then NAA area is likely OK, but if it's < 1.0 
        # then maybe NAA is abnormally low vs Cr and I could adjust NAA up a bit
        # to ensure that other metabs are not underestimated due to low NAA area
        # - should make sure that NAA is not bumped up too, though
        #
        # empirically 0.9 seems to be a good level
        amp_ref = naamp
        conc_ref = cnaa
        if data2ideal < 0.9:
            val = 0.9 / data2ideal
            norm_fudge = 1.7 if val > 1.7 else val
        else:
            norm_fudge = 1.0
            
#        print "util_initial_values - norm_fudge = "+str(norm_fudge)
        
        for i in range(nmet): 

            fudge  = set.prior_area_scale[i]
            metstr = metinfo.get_abbreviation(set.prior_list[i].lower())
            if 'oth'  in metstr: metstr = 'oth'
            if 'pk0'  in metstr: metstr = 'pk0'
            if 'glu'  in metstr: metstr = 'glu'
            if 'sins' in metstr: metstr = 'sino'

            conc, prot = _get_metinfo( metstr, chain, defconc=5.0, defprot=3.0)

            if metstr not in large_metabs:
                mult   = amp_ref * fudge * norm_fudge * (conc/conc_ref)  
                chain.init_area[i]  = mult
                chain.init_freq[i] += reference_peak_offset                

    elif set.initial_small_peak_areas == FitInitialSmallPeakAreas.PEAK_2PPM:
        # no list of large_metabs to skip here because all peaks reference 
        # their areas off of the entire ref peak area at the 2 or 3 ppm peak
        
        peak_2ppm_conc, _ = _get_metinfo( 'peak_2ppm', chain, defconc=10.0, defprot=3.0)
        for i in range(nmet): 

            fudge  = set.prior_area_scale[i]
            metstr = metinfo.get_abbreviation(set.prior_list[i].lower())
            if 'oth'  in metstr: metstr = 'oth'
            if 'pk0'  in metstr: metstr = 'pk0'

            conc, prot = _get_metinfo( metstr, chain, defconc=5.0, defprot=1.0)

            mult = peak_2ppm_amp * fudge * (conc/peak_2ppm_conc)  
            chain.init_area[i]  = mult
            chain.init_freq[i] += reference_peak_offset
            
    elif set.initial_small_peak_areas == FitInitialSmallPeakAreas.PEAK_3PPM: 
        # no list of large_metabs to skip here because all peaks reference 
        # their areas off of the entire ref peak area at the 2 or 3 ppm peak
        
        peak_3ppm_conc, _ = _get_metinfo( 'peak_3ppm', chain, defconc=8.0, defprot=3.0)
        for i in range(nmet): 

            fudge  = set.prior_area_scale[i]
            metstr = metinfo.get_abbreviation(set.prior_list[i].lower())
            if 'oth'  in metstr: metstr = 'oth'
            if 'pk0'  in metstr: metstr = 'pk0'

            conc, prot = _get_metinfo( metstr, chain, defconc=5.0, defprot=1.0)

            mult = peak_3ppm_amp * fudge * (conc/peak_3ppm_conc)   
            chain.init_area[i]  = mult
            chain.init_freq[i] += reference_peak_offset            

    # End Areas and PPMs
    #---------------------------------------------------------


    #---------------------------------------------------------
    # cho/cr separation section

    if set.initial_cr_cho_separation and (chcrflg >= 2):

        cho    = int(np.round(dataset.ppm2pts(3.2)))
        cre    = int(np.round(dataset.ppm2pts(3.0)))
        ppm01  = int(np.round(dataset.ppm2pts(0.1, rel=True)))
        ppm02  = int(np.round(dataset.ppm2pts(0.2, rel=True)))
        delcho = 0.0
        delcho = 0.0

        prior_list = set.prior_list
        abbr_list = [metinfo.get_abbreviation(item).lower() for item in set.prior_list]

        if 'cho' in abbr_list and 'cr' in abbr_list:
            # FIXME bjs - right now, fails silently - fix in GUI?
            chox   = abbr_list.index('cho')
            fcho   = chain.init_freq[chox]
            delcho = np.abs(fcho - cho)
            crex   = abbr_list.index('cr')
            fcre   = chain.init_freq[crex]
            delcre = np.abs(fcre - cre)
    
            if np.abs(fcho-fcre) < ppm01:               # too close
                if delcho <= delcre:                    # reset Cre
                    chain.init_freq[crex] = fcho + ppm02
                else:                                   # reset Cho
                    chain.init_freq[chox] = fcre - ppm02



    #---------------------------------------------------------
    # calculations for Mmol basis - single dataset - if needed
    
    chain.mmol_area = 1.0
    chain.mmol_fre  = chain.init_b0       # in Hz

    if set.macromol_model == FitMacromoleculeMethod.SINGLE_BASIS_DATASET:
        # in this method we assume that the basis function dataset provided
        # is 'correct' in that no freq shift or phase is needed and that the
        # real part of the data is an absorption spectra
        
        # default values
        mmol_fre  = 0.001       # in Hz - this is all we have for and algorithm now
        mmol_area = 1.0
        
        if set.macromol_single_basis_dataset_initval == FitMacromoleculeMethodInitVal.MANUAL:
            mmol_area = set.macromol_single_basis_dataset_start_area

        else:
            if chain.basis_mmol is not None:
                # get widget values for algorithm
                ppm_start = set.macromol_single_basis_dataset_ppm_start
                ppm_end   = set.macromol_single_basis_dataset_ppm_end 
                ppm_start = int(dataset.ppm2pts(ppm_start, acq=False))
                ppm_end   = int(dataset.ppm2pts(ppm_end  , acq=False))
                
                # create spectrum for mmol basis function
                npts    = dataset.raw_dims[0]
                zfmult  = dataset.zero_fill_multiplier
                nptszf  = int(round(npts * zfmult))
                mf      = np.zeros((int(nptszf),),complex)
                mdat    = chain.basis_mmol.copy() * ((((np.arange(npts)+1)%2)*2)-1)
                mdat 
                mf[0:npts] = mdat
                mf[0]      = mf[0] / 2.0            
                mf = np.fft.fft(mf)/nptszf
        
                # algorithm - calc ratio at point or range avg for data/basis
                if set.macromol_single_basis_dataset_initval == FitMacromoleculeMethodInitVal.ONEPOINT:
                    mmol_area = np.real(mmol_ddat[ppm_start]) / np.real(mf[ppm_start])
                    mmol_area *= set.macromol_single_basis_dataset_start_fudge
                    
                elif set.macromol_single_basis_dataset_initval == FitMacromoleculeMethodInitVal.REGION:
                    mmol_area = np.real(np.mean(mmol_ddat[ppm_end:ppm_start])) / np.real(np.mean(mf[ppm_end:ppm_start]))
                    mmol_area *= set.macromol_single_basis_dataset_start_fudge
    
        mmol_area_max = mmol_area * set.macromol_single_basis_dataset_limit_max
        mmol_area_min = mmol_area * set.macromol_single_basis_dataset_limit_min
        mmol_fre_max  = mmol_fre + set.optimize_bounds_range_ppm     # in Hz    default to metabolite limits for now
        mmol_fre_min  = mmol_fre - set.optimize_bounds_range_ppm     
    
        if mmol_area    == 0.0: mmol_area    = 0.01
        if mmol_fre     == 0.0: mmol_fre     = 0.01
        if mmol_fre_max == 0.0: mmol_fre_max = 0.01
        if mmol_fre_min == 0.0: mmol_fre_min = 0.01

        chain.mmol_area     = mmol_area
        chain.mmol_area_min = mmol_area_max
        chain.mmol_area_max = mmol_area_min
        chain.mmol_fre      = mmol_fre * 2.0 * np.pi       # convert to radians
        chain.mmol_fre_max  = mmol_fre_max * 2.0 * np.pi   # convert to radians
        chain.mmol_fre_min  = mmol_fre_min * 2.0 * np.pi   # convert to radians            

    #---------------------------------------------------------
    # store freq as radians
    chain.init_freq = dataset.pts2hz(chain.init_freq, rel=True) * 2.0 * np.pi    # in radians




def _optimize_phase01_correlation(data, modfn, chain, 
                                    zero=False, 
                                    one=False,   
                                    nlag=3, 
                                    cdeg=1, 
                                    b0str=3.4, 
                                    b0end=2.8  ):
    """
    Main call for correlation auto phase algorithm.  Returns the phase
    in deg at which the real part is maxed as a two element float array 
    with optimal zero and first order phases.

    INPUT:
      data:   complex data to be phased
      modfn:  complex model function against which corrections will be made
      chain:   optimization control structure

    KEYWORDS:
      zero:   flag, set to optimize zero order phase
      one:    flag, set to optimize first order phase
      nlag:   number of points to lag in the correlation
      cdeg:   degree of mixing for calculation CSUMM VALUE in correlation
      b0str:  in ppm, start point for B0 correction region
      b0end:  in ppm, end point for B0 correction region

    OUTPUT
      phased: returns the input DATA with 0/1 phase applied


    """
    block   = chain._block
    set     = chain._block.set
    dataset = chain._dataset
    
    prior   = set.prior
    metinfo = dataset.user_prior.metinfo

    if len(modfn.shape)>1: 
        modfn = modfn.flatten()
    cdeg = np.abs(cdeg)
    if nlag <= 0: 
        nlag = 1

    dim0 = len(data)
    ph0  = 0.001
    ph1  = 0.001
    res  = [0.0,0.0]
    pts2 = [0.0,0.0]
    pts4 = [0.0,0.0,0.0,0.0]

    # B0 Correction Section!  Must have good alignment or else.

    pstr = int(np.round(dataset.ppm2pts(b0str)))
    pend = int(np.round(dataset.ppm2pts(b0end)))
    ps = np.clip(pstr,0,dim0)
    pe = np.clip(pend,0,dim0)
    if ps > pe:
        ps,pe = pe,ps

    if ps == 0:
        pe = pe+3 
        pe = np.clip(pe,0,dim0)
    elif ps == dim0:
        pe = dim0
        ps = dim0-3
        ps = ps if ps>0 else 0
    else:
        ps = ps-1
        ps = np.clip(ps,0,dim0)
        pe = pe+1
        pe = np.clip(pe,0,dim0)

    dat  = data[ps:pe].copy()
    ref  = modfn[ps:pe].copy()

    shft, csum = b0_correction(dat,ref)

    dat  = np.roll(data.copy(), shft)

    # Set Correlation Limits

    str0 = dataset.user_prior.auto_phase0_range_start
    end0 = dataset.user_prior.auto_phase0_range_end
    str1 = dataset.user_prior.auto_phase1_range_start
    end1 = dataset.user_prior.auto_phase1_range_end
    max0 = str0 if str0>end0 else end0
    min0 = str0 if str0<end0 else end0
    max1 = str1 if str1>end1 else end1
    min1 = str1 if str1<end1 else end1

    pts4[0] = int(dataset.ppm2pts(max0))
    pts4[1] = int(dataset.ppm2pts(min0))
    pts4[2] = int(dataset.ppm2pts(max1))
    pts4[3] = int(dataset.ppm2pts(min1))

    # Phasing section

    if zero and not one:        # zero order phase only

        ph0, _ = optimize_phase0_correlation(dat, modfn, pts4, nlag, cdeg)
        phase  = np.exp(1j * ph0 * DTOR)
        res    = [ph0,ph1]

    if one and not zero:        # first order phase only

        pivot  = dataset.ppm2pts(dataset.user_prior.auto_phase1_pivot)
        ph1, _ = optimize_phase1_correlation(dat, modfn, pts4, pivot, nlag, cdeg)
        phase  = np.exp(1j * ph1 * DTOR * (np.arange(dim0)-pivot)/dim0)
        res    = [ph0,ph1]

    if zero and one:            # full phase

        pivot = dataset.ppm2pts(dataset.user_prior.auto_phase1_pivot)

        ph01, csum01 = optimize_phase0_correlation(dat.copy(), modfn, pts4, nlag, cdeg)
        phase0 = np.exp(1j * ph01 * DTOR)
        dat1  = phase0 * dat.copy()

        ph11, csum11  = optimize_phase1_correlation(dat1.copy(), modfn, pts4, pivot, nlag, cdeg)
        phase1 = np.exp(ph11 * DTOR * (np.arange(dim0)-pivot)/dim0)

        phase_dat1 = phase0 * phase1
        # ----------

        pts2[0] = pts4[2]
        pts2[1] = pts4[3]

        ph02, csum02 = optimize_phase0_correlation(dat.copy(), modfn, pts2, nlag, cdeg)
        phase0 = np.exp(1j * ph02 * DTOR)
        dat2   = phase0 * dat.copy()

        ph12, csum12 = optimize_phase1_correlation(dat2.copy(), modfn, pts4, pivot, nlag, cdeg)
        phase1 = np.exp(1j * ph12 * DTOR * (np.arange(dim0)-pivot)/dim0)

        phase_dat2 = phase0 * phase1

        #----------
        if csum11 > csum12:
            phase  = phase_dat1
            res = [ph01,ph11]
        else:
            phase  = phase_dat2
            res = [ph02,ph12]

    return res, phase    



def _optimize_phase01_integration(data, chain, zero=False, one=False):
    """
    Main call for real integration auto phase algorithm.  Returns the phase
    in deg at which the real part is maxed. Outputs a two element float 
    array with optimal zero and first order phases in degrees.

    INPUT:
      data:
      chain: pointer to optimization control structure

    KEYWORDS:
      ZERO: flag, set to optimize zero order phase
      ONE:  flag, set to optimize first order phase

    """
    block   = chain._block
    set     = chain._block.set
    dataset = chain._dataset
    
    prior   = set.prior
    metinfo = dataset.user_prior.metinfo

    dim0 = dataset.spectral_dims[0]

    ph0 = 0.001
    ph1 = 0.001
    phase0 = 1
    phase1 = 1

    str0 = dataset.user_prior.auto_phase0_range_start
    end0 = dataset.user_prior.auto_phase0_range_end
    str1 = dataset.user_prior.auto_phase1_range_start
    end1 = dataset.user_prior.auto_phase1_range_end
    max0 = str0 if str0>end0 else end0
    min0 = str0 if str0<end0 else end0
    max1 = str1 if str1>end1 else end1
    min1 = str1 if str1<end1 else end1

    dat = data.copy()
    weight = np.zeros(dataset.spectral_dims[0])

    if zero:        # zero order phase part

        ps = int(np.round(dataset.ppm2pts(max0)))
        pe = int(np.round(dataset.ppm2pts(min0)))
        weight[ps:pe] = 1.0
        ph0 = optimize_phase0_integration(data.copy(), weight)

        phase0 = np.exp(1j * ph0 * DTOR)
        phased = phase0 * data.copy()


    if one:         # first order phase part

        ps = int(np.round(dataset.ppm2pts(max1)))
        pe = int(np.round(dataset.ppm2pts(min1)))
        pivot = int(np.round(dataset.ppm2pts(dataset.user_prior.auto_phase1_pivot)))
        weight = weight * 0.0
        weight[ps:pe] = 1.0

        ph1 = optimize_phase1_integration(phased.copy(), weight, pivot)

        phase1 = np.exp(1j * ph1 * DTOR * (np.arange(dim0)-pivot) / dim0)
        phased = phase1 * phased

    phase = phase0 * phase1
    res = [ph0,ph1]

    return res, phase    


def _height2area_function(val, info):
    """
    This is the minimization function used by minf_parabolic_info in the
    _calc_height2area_ratio() call. The val parameter is the decay value
    for which we need to calculate a FWHM line width. Because we are minimizing
    in this optimization, we subtract the calculated value from the original
    line width values (in Hz) and take the absolute value.

    """

    ta = val if info["ta"] == -1 else info["ta"]
    tb = val if info["tb"] == -1 else info["tb"]

    width_hz, peak = util_spectral.voigt_width(ta, tb, info["chain"])

    info["peak"] = peak

    return np.abs(info["orig_lw"] - width_hz)



def _calc_height2area_ratio(lw, chain, ta=-1.0, tb=-1.0 ):
    """
    We know the value of the full width half max line width in Hz that we have
    in our data, and want to find the Ta and Tb values that yield this.

    This function uses the minf_parabolic_info routine to optimze Ta and Tb
    to values between 0.005 and 0.5, however either of the two parameters can
    also be set to constant values by setting the TA and TB keywords to this
    function to the constant value desired. This way we can calculate Pure 
    Gauss or Pure Lorentz lineshape starting values as well as Voigt/LorGauss
    line shape values.

    The optimization calls the fitt_height2area_function() to determine the 
    minimization function. As part of that call, we calculate the height of the
    peak for Ta and Tb, which is stored in the info.peak parameter. This is 
    used to provide a normalization value for peak height to peak area 
    conversion on return of this function.

     lw - float, linewidth in Hz
     chain - pointer to control structure
     ta - keyword, float, constant value for Ta in the optimization
     tb - keyword, float, constant value for Tb in the optimization

    """
    info = { 'ta':ta, 'tb':tb, 'orig_lw':lw, 'chain':chain._dataset, 'peak':-1.0 }

    # Call parabolic interpolation, Brent's method
    # 1-d minimization routine

    val_a = 0.005  # lower bound
    val_b = 0.06   # "some" point in the middle
    val_c = 0.5    # upper bound

    finalval, maxit = minf.minf_parabolic_info( val_a, val_b, val_c,
                                                _height2area_function, 
                                                info ) 
    return [finalval, info["peak"]]

