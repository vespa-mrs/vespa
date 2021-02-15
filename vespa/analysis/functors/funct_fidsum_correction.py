# Python imports

# 3rd party imports
import numpy as np
from scipy.optimize import least_squares

# Vespa imports
from vespa.analysis.algos.b0_correction import b0_correction
from vespa.analysis.algos.optimize_phase import optimize_phase0_correlation

from vespa.analysis.algos.suspect_frequency_correction import spectral_registration, rats


PI = np.pi
DTOR = np.pi/180.0

CORRECTION_MENU_ITEMS = [ 'Manual',
                          'Optimized Search (vespa)',
                          'Correlation (vespa)',
                          'Spectral Registration (suspect)',
                          'RATS (suspect)']
 

def ppm2pts(ppm, frequency, resppm, sw, npts):
    return int((npts / 2) - (frequency * (ppm - resppm) / (sw/npts)))


def correction_optimized_search_vespa(raw, chain):
    """ 
    FID frequency and phase correction method - Vespa
    
    Corrects individual FID data from an MRS data set for frequency and phase
    drifts through time

    Inputs:
        raw (ndarray, complex): shape=(1, 1, nfids, npts) dimensions, raw kspace data
        chain (object): block settings for algorithm
    Outputs:
        corr_freq (ndarray, float): Frequency shift values in [Hz] in numpy
            array of length 'nfids'.
        corr_phas (ndarray, float): Phase0 correction values in [deg] in numpy
            array of length 'nfids'.
    """
    set = chain._block.set
    ds  = chain._dataset

    _, ncoil, nfids, npts = raw.shape

    #------------------------------------------------------
    # Global inits and one time calculations

    corr_freq = np.zeros([nfids,], dtype=np.float)
    corr_phas = np.zeros([nfids,], dtype=np.float)

    zfmult = 4      # larger zfmult here improves peak shift accuracy
    raw_dim0 = ds.raw_shape[-1]
    raw_hpp  = ds.sw / raw_dim0
    fid_dim0 = raw_dim0 * zfmult
    fid_hpp  = ds.sw / fid_dim0

    # temporary arrays
    xx = np.arange(raw_dim0) / ds.sw
    search = np.zeros((fid_dim0,),complex)
    tmp_summed = np.zeros((raw_dim0,), dtype=np.complex)

    # one time calculations
    if set.chop_data:
        chop = ((((np.arange(raw_dim0) + 1) % 2) * 2) - 1)
    else:
        chop = 1.0

    #----------------------------------------------------------------
    # Shift Calculation

    # convert algorithm values from PPM to points
    ctr = set.vespa_reference_peak_center
    wid = set.vespa_peak_search_width
    search_str = ctr + wid
    search_end = ctr - wid
    refpt      = ppm2pts(ctr, ds.frequency, ds.resppm, ds.sw, fid_dim0)
    search_str = ppm2pts(search_str, ds.frequency, ds.resppm, ds.sw, fid_dim0)
    search_end = ppm2pts(search_end, ds.frequency, ds.resppm, ds.sw, fid_dim0)

    for i in range(nfids):
        # Calculate shift with zf x 4 over user-set range on magnitude data
        time = raw[0, 0, i, :].copy()
        search *= 0.0
        search[0:raw_dim0] = time * chop
        tmp = np.abs(np.fft.fft(search))
        imax = tmp[search_str:search_end].argmax()
        corr_freq[i] = (refpt - (search_str + imax)) * fid_hpp      # pts2hz

    # Now shift all averages (FIDs) and sum for phase0 correction step
    if set.vespa_target_method == 'Average all':
        fstart, fend = 0, nfids
    elif set.vespa_target_method == 'Avg first 4':
        fstart, fend = 0, 4
    elif set.vespa_target_method == 'Avg first 10%':
        fstart, fend = 0, int(nfids*0.1)
    elif set.vespa_target_method == 'Avg first 25%':
        fstart, fend = 0, int(nfids*0.25)
    elif set.vespa_target_method == 'Avg first 50%':
        fstart, fend = 0, int(nfids*0.5)
    elif set.vespa_target_method == 'Avg mid 10%':
        fstart, fend = int(nfids*0.45), int(nfids*0.55)
    elif set.vespa_target_method == 'Avg mid 30%':
        fstart, fend = int(nfids*0.35), int(nfids*0.65)
    elif set.vespa_target_method == 'User define':
        fstart, fend = 0, nfids

    if fend == fstart: fend = fend+1

    for i in range(fstart, fend):
        tmp_summed += raw[0, 0, i, :].copy() * np.exp(1j * 2.0 * np.pi * corr_freq[i] * xx)

    #----------------------------------------------------------------
    # Phase0 Calculation
    # - ref spectrum created from shift corrected, summed FID
    # - this loop only does phase0 correction

    ph0_str = ppm2pts(set.vespa_phase0_range_start, ds.frequency, ds.resppm, ds.sw, raw_dim0)
    ph0_end = ppm2pts(set.vespa_phase0_range_end,   ds.frequency, ds.resppm, ds.sw, raw_dim0)

    ref = tmp_summed * chop
    ref[0] *= 0.5
    ref = (np.fft.fft(ref) / len(ref))
    ref /= nfids  # scale for comparison to single FID (not nfids_excluded here)
    ref = ref[ph0_str:ph0_end]

    for i in range(nfids):
        # get freq shifted fid
        dat = raw[0, 0, i, :].copy() * np.exp(1j * 2.0 * np.pi * corr_freq[i] * xx)
        dat[0] *= 0.5
        dat = (np.fft.fft(dat * chop) / len(dat))
        dat = dat[ph0_str:ph0_end]

        def fph0(x):
            """ Inline function for ph0 optimize in FIDs - rad in [radians] """
            return np.array([np.sum((dat*np.exp(1j*x[0]) - ref).real ** 2),])

        r = least_squares(fph0, np.array([0.1*PI,]), bounds=([-PI,], [PI*2,]), method='trf')
        phdeg = r.x * 180.0 / np.pi

        corr_phas[i] = phdeg-360 if phdeg > 180.0 else phdeg      # remove (some) wrap around

    return corr_freq, corr_phas


def correction_correlation_vespa(raw, chain):

    set = chain._block.set
    ds  = chain._dataset
    prior = set.vespa_preprocess_prior

    _, ncoil, nfids, npts = raw.shape

    #------------------------------------------------------
    # Global inits and one time calculations

    corr_freq = np.zeros([nfids,], dtype=np.float)
    corr_phas = np.zeros([nfids,], dtype=np.float)

    raw_dim0 = ds.raw_shape[-1]
    raw_hpp  = ds.sw / raw_dim0

    # temporary arrays
    chop = ((((np.arange(raw_dim0) + 1) % 2) * 2) - 1) if set.chop_data else 1.0

    #----------------------------------------------------------------
    # B0 Shift and Phase0 Calculations

    nlag = 3  # default from fit_voigt
    cdeg = 5  # default from fit_voigt
    newzf = 4

    pe = int(np.round(ds.ppm2pts(prior.auto_b0_range_start))) * newzf      # ps and pe reversed here due to PPM
    ps = int(np.round(ds.ppm2pts(prior.auto_b0_range_end))) * newzf
    set.vespa_preprocess_prior.basis.update(ds, zfmult=newzf)
    ref = set.vespa_preprocess_prior.basis.get_spectrum_sum().flatten().copy()  # not sending in dataset to keep zfmult set in line above!
    bref = ref[ps:pe]

    min0, max0 = int(ds.ppm2pts(prior.auto_phase0_range_start)), int(ds.ppm2pts(prior.auto_phase0_range_end))
    min1, max1 = int(ds.ppm2pts(2.75)), int(ds.ppm2pts(3.6))  # prior.auto_phase1_range_start, prior.auto_phase1_range_end
    pts4 = [max0 * newzf, min0 * newzf, max1 * newzf, min1 * newzf]

    tmp = np.zeros([ds.raw_shape[-1]*newzf,])

    for i in range(nfids):

        tmp[0:ds.raw_shape[-1]] = raw[0, 0, i, :].copy() * chop
        dat = np.fft.fft(tmp)

        # Calculate B0 shift
        shft, _ = b0_correction(dat[ps:pe], bref)
        corr_freq[i] = shft * (raw_hpp/newzf)  # in Hz
        dat = np.roll(dat, int(round(shft)))

        # Phase estimate
        # - maximize correlation of phase 0 of real spectrum to ideal spectrum
        ph0, _ = optimize_phase0_correlation(dat, ref, pts4, nlag, cdeg)

        corr_phas[i] = ph0

    return corr_freq, corr_phas



def correction_spectral_registration_suspect(raw, chain):
    """ 
    FID frequency and phase correction method - Time Domain from Suspect, but like FID-A
    
    Corrects individual FID data from an MRS data set for frequency and phase
    drifts through time

    Inputs:
        raw (ndarray, complex): shape=(1, 1, nfids, npts) dimensions, raw kspace data
        chain (object): block settings for algorithm
    Outputs:
        corr_freq (ndarray, float): Frequency shift values in [Hz] in numpy
            array of length 'nfids'.
        corr_phas (ndarray, float): Phase0 correction values in [deg] in numpy
            array of length 'nfids'.
   
    """
    set = chain.set
    ds  = chain._dataset

    _, ncoil, nfids, npts = raw.shape

    opt_str = set.suspect_optimization_range_start
    opt_end = set.suspect_optimization_range_end
    fr_init = set.suspect_initial_guess_freq
    ph_init = set.suspect_initial_guess_phase * DTOR

    corr_freq = np.zeros([nfids,], dtype=np.float)
    corr_phas = np.zeros([nfids,], dtype=np.float)

    if set.suspect_target_method == 'Average all':
        fstart, fend = 0, nfids
    elif set.suspect_target_method == 'Avg first 4':
        fstart, fend = 0, 4
    elif set.suspect_target_method == 'Avg first 10%':
        fstart, fend = 0, int(nfids * 0.1)
    elif set.suspect_target_method == 'Avg first 25%':
        fstart, fend = 0, int(nfids * 0.25)
    elif set.suspect_target_method == 'Avg first 50%':
        fstart, fend = 0, int(nfids * 0.5)
    elif set.suspect_target_method == 'Avg mid 10%':
        fstart, fend = int(nfids * 0.45), int(nfids * 0.55)
    elif set.suspect_target_method == 'Avg mid 30%':
        fstart, fend = int(nfids * 0.35), int(nfids * 0.65)
    elif set.suspect_target_method == 'User define':
        fstart, fend = 0, nfids

    if fend == fstart: fend = fend + 1

    fids = np.squeeze(raw.copy())
    targ = np.squeeze(raw[0,0,fstart:fend,:].copy())
    targ = np.sum(targ, axis=0) / (fend-fstart)

    for i in range(fids.shape[0]):
        fid = np.squeeze(fids[i])

        fs, ps = spectral_registration( fid, targ, chain,
                                        initial_guess=(fr_init,ph_init),
                                        frequency_range=(opt_str, opt_end) )

        corr_freq[i] = -fs
        psd = -ps * (1.0/DTOR)
        corr_phas[i] = psd - 360 if psd > 180.0 else psd

    return corr_freq, corr_phas



def correction_rats_suspect(raw, chain):
    """ 
    FID frequency and phase correction method - RATS from Suspect
    
    Corrects individual FID data from an MRS data set for frequency and phase
    drifts through time

    Inputs:
        raw (ndarray, complex): shape=(1, 1, nfids, npts) dimensions, raw kspace data
        chain (object): block settings for algorithm
    Outputs:
        corr_freq (ndarray, float): Frequency shift values in [Hz] in numpy
            array of length 'nfids'.
        corr_phas (ndarray, float): Phase0 correction values in [deg] in numpy
            array of length 'nfids'.
    
    """
    set = chain.set
    ds  = chain._dataset

    _, ncoil, nfids, npts = raw.shape

    corr_freq = np.zeros([nfids,], dtype=np.float)
    corr_phas = np.zeros([nfids,], dtype=np.float)

    opt_str = set.rats_optimization_range_start
    opt_end = set.rats_optimization_range_end
    fr_init = set.rats_initial_guess_freq
    ph_init = set.rats_initial_guess_phase * DTOR
    bord    = set.rats_baseline_order

    if set.rats_target_method == 'Average all':
        fstart, fend = 0, nfids
    elif set.rats_target_method == 'Avg first 4':
        fstart, fend = 0, 4
    elif set.rats_target_method == 'Avg first 10%':
        fstart, fend = 0, int(nfids * 0.1)
    elif set.rats_target_method == 'Avg first 25%':
        fstart, fend = 0, int(nfids * 0.25)
    elif set.rats_target_method == 'Avg first 50%':
        fstart, fend = 0, int(nfids * 0.5)
    elif set.rats_target_method == 'Avg mid 10%':
        fstart, fend = int(nfids * 0.45), int(nfids * 0.55)
    elif set.rats_target_method == 'Avg mid 30%':
        fstart, fend = int(nfids * 0.35), int(nfids * 0.65)
    elif set.rats_target_method == 'User define':
        fstart, fend = 0, nfids

    if fend == fstart: fend = fend + 1

    fids = np.squeeze(raw.copy())
    targ = np.squeeze(raw[0,0,fstart:fend,:].copy())
    targ = np.sum(targ, axis=0) / (fend-fstart)

    for i in range(fids.shape[0]):
        fid = np.squeeze(fids[i])

        fs, ps = rats( fid, targ, chain, baseline_order=bord,
                       initial_guess=(fr_init,ph_init),
                       frequency_range=(opt_str, opt_end) )

        corr_freq[i] = -fs
        psd = -ps * (1.0/DTOR)
        corr_phas[i] = psd - 360 if psd > 180.0 else psd

    return corr_freq, corr_phas



#------------------------------------------------------------------------------
# test code

def _test():

    import glob
    import numpy as np
    from matplotlib import pyplot as plt
    from vespa.analysis.util_file_import import get_datasets_cli
    from vespa.analysis.mrs_user_prior import UserPrior

    fpath = r"D:\Users\bsoher\code\repository_svn\sample_data\siemens_dicom_export_fids\*.ima"
    fnames = glob.glob(fpath)

    indices = [item.split('.') for item in fnames]
    indices = [int(item[4]) for item in indices]

    bob = [[fname, indx] for fname,indx in zip(fnames, indices)]
    bob.sort(key=lambda bob: bob[1])

    fnames = [item[0] for item in bob]


    ds = get_datasets_cli(fnames, 'import_siemens_dicom_fidsum', None)[0]
    met = ds.blocks['raw'].data

    class Set(object):
        def __init__(self):

            self.vespa_reference_peak_center        = 2.01
            self.vespa_peak_search_width            = 0.2
            self.vespa_phase0_range_start           = 3.5
            self.vespa_phase0_range_end             = 0.5
            self.vespa_target_method = 'Avg first 4'

            self.vespa_preprocess_prior = UserPrior()

            self.suspect_initial_guess_freq         = 0.1   # in Hz
            self.suspect_initial_guess_phase        = 0.1   # in deg here, but suspect wants rads
            self.suspect_optimization_range_start   = 3.3   # empirical
            self.suspect_optimization_range_end     = 2.7   # empirical
            self.suspect_target_method = 'Avg first 4'

            self.rats_initial_guess_freq            = 0.1   # in Hz
            self.rats_initial_guess_phase           = 0.1   # in deg here, but suspect wants rads
            self.rats_optimization_range_start      = 3.3   # empirical
            self.rats_optimization_range_end        = 2.7   # empirical
            self.rats_baseline_order                = 2     # from suspect
            self.rats_target_method = 'Avg first 4'

            self.global_left_shift                  = 0         # was fid_left_shift
            self.global_phase0                      = 0.0       # was phase0_global_const
            self.global_phase1                      = 0.0
            self.global_gaussian_apodization        = 0.0       # was gaussian_apodization
            self.chop_data                          = True
            self.zero_phase1                        = True
            self.apply_peak_shift                   = True
            self.apply_phase0                       = True


    class Block(object):
        def __init__(self):
            self.set = Set()

    class Chain(object):
        def __init__(self, raw, dataset):
            self.data = raw
            self._dataset = dataset
            self.voxel = [0,0,0]
            self._block = Block()

    chain = Chain(met, ds)

    prior = chain._block.set.vespa_preprocess_prior
    prior.basis.area  = [0.51, 0.42, 0.40]
    prior.basis.ppm   = [2.01, 3.01, 3.22]
    prior.basis.phase = [ 0.0,  0.0,  0.0]
    prior.basis.lw    = [ 5.0,  5.0,  5.0]
    prior.basis.update(ds, zfmult=4)

    prior.auto_b0_range_end       = 3.4
    prior.auto_b0_range_start     = 1.7
    prior.auto_phase0_range_end   = 2.25
    prior.auto_phase0_range_start = 1.85

    raw = chain.data.copy()

#    freq, phas = correction_optimized_search_vespa(raw, chain)
    freq, phas = correction_correlation_vespa(raw, chain)

    fig, axs = plt.subplots(3)
    axs[0].plot(freq)
    axs[1].plot(phas)
    axs[2].plot(prior.basis.get_spectrum_sum(ds).flatten())
    plt.show()


    bob = 10
    bob += 1



if __name__ == '__main__':
    """ Just testing the phase() method here """
    _test()