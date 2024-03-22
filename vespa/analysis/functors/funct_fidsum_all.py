# Python modules

# 3rd party modules
import numpy as np

# Our modules
import vespa.analysis.functors.funct_fidsum_coil_combine as funct_combine
import vespa.analysis.functors.funct_fidsum_exclude      as funct_exclude
import vespa.analysis.functors.funct_fidsum_correction   as funct_correct

from vespa.common.util.math_ import safe_exp
from vespa.common.constants import DEGREES_TO_RADIANS as DTOR



def do_processing_all(chain):
    run_calculations(chain)
    apply_corrections(chain)
    apply_manual_adjustments(chain)
    return


def do_processing_correct_adjust(chain):
    apply_corrections(chain)
    apply_manual_adjustments(chain)
    return


def do_processing_adjust(chain):
    apply_manual_adjustments(chain)
    return


def run_calculations(chain):
    """
    Because we are bumping the zero fill factor here by a factor of 4 to 
    better find the peak max, we can not use the standard util_ppm 
    functions to calculate some conversions. So long hand here for us.
    
    """
    set   = chain.set
    ds    = chain._dataset

    _, ncoil, nfids, npts = chain.raw.shape

    #--------------------------------------------------------------------------
    # Coil Combine Section
    #
    # - do not combine if only 1 coil, no need

    if ncoil > 1:      # number of coils

        if set.coil_combine_method=='Siemens':
            r = funct_combine.siemens(chain.raw)
        elif set.coil_combine_method=='CMRR':
            r = funct_combine.cmrr_standard(chain.raw)
        elif set.coil_combine_method=='CMRR-Sequential':
            r = funct_combine.cmrr_sequential(chain.raw)
        elif set.coil_combine_method=='CMRR-Hybrid':
            r = funct_combine.cmrr_hybrid(chain.raw)
        elif set.coil_combine_method=='SVD (suspect)':
            r = funct_combine.svd_suspect(chain.raw)
        elif set.coil_combine_method=='External Dataset':
            r = funct_combine.external_dataset(chain)
        elif set.coil_combine_method=='External Dataset with Offset':
            r = funct_combine.external_dataset_with_offset(chain)
    else:
        # single coil, copy first channel only
        r = normalize_shape(chain.raw[0,0,:,:].copy()), None, None

    chain.raw_combined         = r[0]
    chain.coil_combine_weights = r[1]
    chain.coil_combine_phases  = r[2]

    #--------------------------------------------------------------------------
    # FID Exclusion Section
    #
    # - This section is always run (if not 'Manual'), but results not always
    #     applied depending on whether the 'apply_data_exclusion' is set.
    # - We only keep a list of indices to exclude rather than removing any data.

    data = chain.raw_combined.copy()

    if set.exclusion_input_adjust:
        data = apply_left_shift(data, set.global_left_shift)
        data = apply_apod_gauss(data, set.global_gaussian_apodization, ds.sw, npts)
        data = apply_phase0(data, set.global_phase0)

    if not set.exclude_method == 'Manual':

        # these return numpy arrays
        if set.exclude_method == 'Remove Bad Averages (fid-a)':
            exclude_indices = funct_exclude.exclude_remove_bad_averages_fida(data, chain)
        elif set.exclude_method == 'Remove N Worst Averages (fid-a)':
            exclude_indices = funct_exclude.exclude_remove_n_worst_averages_fida(data, chain)

        if exclude_indices == []:
            chain.exclude_indices = []
        else:
            chain.exclude_indices = exclude_indices.tolist()

    chain.results_indices = np.delete(np.arange(nfids, dtype=np.int16), chain.exclude_indices)

    #--------------------------------------------------------------------------
    # FID Correction Section
    #
    # - exclude data, if activated
    # - massage input data with phase0 and apodize, if activated
    # - only reset results if we are recalculating, not if 'Manual'

    data = chain.raw_combined.copy()
    results_indices = np.arange(nfids, dtype=np.int16)       # start with all indices
    if set.apply_data_exclusion:
        results_indices = chain.results_indices            # keep certain indices if exclude True
        data = np.delete(data, chain.exclude_indices, axis=2)

    if set.correction_input_adjust:
        data = apply_left_shift(data, set.global_left_shift)
        data = apply_apod_gauss(data, set.global_gaussian_apodization, ds.sw, npts)
        data = apply_phase0(data, set.global_phase0)

    if set.correction_method != 'Manual':

        chain.frequency_shift *= 0.0        # only reset if we are recalculating
        chain.phase_0         *= 0.0

        if set.correction_method == 'Optimized Search (vespa)':
            r = funct_correct.correction_optimized_search_vespa(data, chain)
        elif set.correction_method == 'Correlation (vespa)':
            r = funct_correct.correction_correlation_vespa(data, chain)
        elif set.correction_method == 'Spectral Registration (suspect)':
            r = funct_correct.correction_spectral_registration_suspect(data, chain)
        elif set.correction_method == 'RATS (suspect)':
            r = funct_correct.correction_rats_suspect(data, chain)

        chain.frequency_shift[results_indices] = r[0]
        chain.phase_0[results_indices]         = r[1]

    return


def apply_corrections(chain):
    """ Apply phase0 and freq shift correction to individual FID data """

    ds   = chain._dataset
    set  = chain.set
    time = chain.raw_combined.copy()
    xx   = np.arange(time.shape[3]) / ds.sw

    _, ncoil, nfids, npts = chain.raw.shape


    if chain.do_freq_raw:
        chain.freq_raw = chain.raw_combined.copy()

    if set.global_left_shift != 0:
        time = apply_left_shift(time, set.global_left_shift)
        if chain.do_freq_raw:
            chain.freq_raw = apply_left_shift(chain.freq_raw, set.global_left_shift)

    for i in range(nfids):

        if set.apply_peak_shift:
            time[0,0,i,:] *= np.exp(1j * chain.frequency_shift[i] * 2.0 * np.pi * xx)
        if set.apply_phase0:
            time[0,0,i,:] *= np.exp(1j * chain.phase_0[i] * DTOR)

    chain.raw_corrected = time


def apply_manual_adjustments(chain):
    """ apply manual adjustments and FFT to spectrum """

    set   = chain.set
    ds    = chain._dataset

    # One time calculations ---------------------------------------------------

    _, ncoil, nfids, npts = chain.raw.shape

    ph0  = np.exp(1j * set.global_phase0 * DTOR)
    chop = ((((np.arange(npts) + 1) % 2) * 2) - 1) if set.chop_data else 1.0

    if set.global_phase1 != 0.0:
        piv = np.round(ds.ppm2pts(ds.phase_1_pivot, acq=True))
        ff = (np.arange(npts, dtype=float) - piv) / npts
        phase1 = np.exp(1j * (set.global_phase1 * DTOR * ff))
    else:
        phase1 = None

    # Start data adjustments --------------------------------------------------

    # Final summed FID - no apodization
    # - no divide by npts here because that is done in Spectral Tab
    data = chain.raw_corrected.copy() * ph0
    if phase1 is not None:
        data[:,:,:,0] *= 0.5
        data = np.fft.fft(data * chop, axis=3) * phase1
        data = np.fft.ifft(data) * chop
        data[:,:,:,0] *= 2.0

    # Display of summed spectrum - there is apodization if not 0.0
    # - ph0, ph1, apod, and scale are same as for Spectral Tab
    disp = chain.raw_corrected.copy()
    disp = apply_apod_gauss(disp, set.global_gaussian_apodization, ds.sw, npts)
    disp[:,:,:,0] *= 0.5
    disp = np.fft.fft(disp * chop * ph0, axis=3) / npts
    if phase1 is not None: disp *= phase1

    disp_current = disp[0,0,chain.voxel,:].copy()
    scale = nfids - len(chain.exclude_indices) if set.apply_data_exclusion else nfids
    disp_current *= scale

    if chain.do_freq_raw:               # should only have to do this once
        tmp = chain.freq_raw.copy()
        tmp = apply_apod_gauss(tmp, set.global_gaussian_apodization, ds.sw, npts)
        tmp[:,:,:,0] *= 0.5
        tmp = np.fft.fft(tmp * chop * ph0, axis=3) / npts
        if phase1 is not None: tmp *= phase1
        chain.freq_raw = tmp

    chain.time_adjusted = data.copy()
    chain.freq_adjusted = disp.copy()
    chain.freq_current  = normalize_shape(disp_current)

    if set.apply_data_exclusion:
        if len(chain.exclude_indices) > 0:
            disp = np.delete(disp, chain.exclude_indices, axis=2)
            data = np.delete(data, chain.exclude_indices, axis=2)

    chain.time_summed = normalize_shape(np.sum(data, axis=2))
    chain.freq_summed = normalize_shape(np.sum(disp, axis=2))


#--------------------------------------------------------------------
# Helper methods

def apply_left_shift(data, pts):
    """ shift data N points left and zero N points at end of array  """
    if pts > 0:
        data = np.roll(data, -pts, axis=3)
        data[:,:,:,-pts:] *= 0.0
    return data

def apply_apod_gauss(data, apod, sw, npts):
    """ modify data with gaussian envelope with apod [Hz] FWHM """
    return data * safe_exp(-((np.arange(npts)/sw)*(0.6*np.pi*apod))**2)

def apply_phase0(data, ph0):
    """ add ph0 [deg] const phase to the data """
    return data * np.exp(1j * ph0 * DTOR)

def normalize_shape(data):
    while len(data.shape) < 4:
        data.shape = [1, ] + list(data.shape)
    return data



