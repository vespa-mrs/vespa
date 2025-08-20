# Python modules

# 3rd party modules
import numpy as np

# Our modules
import vespa.analysis.constants as constants
import vespa.analysis.svd_output as svd_output_module
import vespa.analysis.functors.funct_ecc as funct_ecc
import vespa.analysis.functors.funct_water_filter as funct_water
from vespa.common.util.math_ import safe_exp

if constants.HLSVDPRO_AVAILABLE:
    import hlsvdpro
else:
    import vespa.common.hlsvdpropy as hlsvdpro




def apodization(chain):
    set = chain._block.set
    t = np.arange(chain.raw_dim0) / chain.sw
    old_settings = np.seterr(all='ignore')      # avoid underflow warnings
    if set.apodization == 'gaussian':
        chain.data = chain.data * safe_exp(-(set.apodization_width * np.pi * 0.6 * t) ** 2, 0)
    elif set.apodization == 'lorentzian':
        chain.data = chain.data * safe_exp(-(set.apodization_width * np.pi * t), 0)
    np.seterr(**old_settings)


def chop(chain):
    # numpy broadcasting deals with (N,dim0) sized data arrays
    chain.data = chain.data * ((((np.arange(chain.raw_dim0) + 1) % 2) * 2) - 1)


def fft(chain):
    set = chain._block.set
    # need to half the value of first point of FID data to get proper area
    if chain.is_fid:
        if len(chain.data.shape) > 1:
            chain.data[:,0] *= 0.5
        else:
            chain.data[0] *= 0.5

    # calculate if zerofilling is required
    dim0 = chain.raw_dim0 * set.zero_fill_multiplier

    # apply FFT if needed
    if set.fft:
        chain.data = np.fft.fft(chain.data, n=dim0) / float(dim0)
    else:
        if set.zero_fill_multiplier > 1:
            temp = np.zeros(dim0, 'complex')
            temp[0:chain.raw_dim0] = chain.data
            chain.data = temp


def flip_spectral_axis(chain):
    set = chain._block.set

    if len(chain.data.shape)>1:
        for i in range(chain.data.shape[0]):
            chain.data[i,:] = chain.data[i,::-1].copy()
    else:
        chain.data = chain.data[::-1].copy()


def frequency_shift(chain):
    set = chain._block.set

    # seterr() avoids underflow error
    old_err_state = np.seterr(all='ignore')
    t = np.arange(chain.raw_dim0) / chain.sw
    chain.data = chain.data * np.exp(1j * 2.0 * np.pi * chain.frequency_shift * t)
    np.seterr(**old_err_state)


def left_shift(chain):
    """ shift data array N points left and zero N points at end """
    set = chain._block.set
    chain.data = np.roll(chain.data, -set.left_shift_value)
    chain.data[-set.left_shift_value:] = chain.data[0]*0.0


def svd_filter(chain):
    set = chain._block.set

    if sum(chain.pre_roll.real):

        # Calculate dwell time in (ms)
        dwell_time = 1000.0 / chain.sw
        
        nsv_sought = chain.nssv

        svd_data = chain.pre_roll.copy()

        if chain.do_fit:
            # Recompute HLSVD lib simulation

            signals = chain.pre_roll[:chain.ndp]

            if constants.HLSVDPRO_AVAILABLE:
                results = hlsvdpro.hlsvd(signals, nsv_sought, dwell_time)
            else:
                results = hlsvdpro.hlsvd(signals, nsv_sought, dwell_time, sparse=False)

            # if constants.HLSVD_METHOD == 'hlsvdpro':
            #     # the FORTRAN library option
            #     results = hlsvdpro.hlsvd(signals, nsv_sought, dwell_time)
            # else:
            #     # the scipy native python option - about 5x slower
            #     results = hlsvdpropy.hlsvd(signals, nsv_sought, dwell_time)

            nsv_found = results[0]

            # The second element of the return tuple is a list of singular
            # values that we don't care about. We don't care about any trailing
            # return values either
            frequencies, damping_factors, amplitudes, phases = results[2:6]

            if nsv_found:
                # hlsvd() returns lists; we want numpy arrays.
                frequencies     = np.array(frequencies)
                damping_factors = np.array(damping_factors)
                amplitudes      = np.array(amplitudes)
                phases          = np.array(phases)
            else:
                # When the HLSVD module finds nothing (due to bugs,
                # algorithmic weakness, sun spots, etc.), we dummy up fake
                # results that are all 0.0 to make downstream code happy.
                nsv_found = nsv_sought
                frequencies     = np.zeros(nsv_sought)
                damping_factors = np.zeros(nsv_sought)
                amplitudes      = np.zeros(nsv_sought)
                phases          = np.zeros(nsv_sought)


            svd_output = svd_output_module.SvdOutput(frequencies,
                                                     damping_factors,
                                                     amplitudes,
                                                     phases)
            chain.svd_output = svd_output
            chain.do_fit = False

            # PS - The code below writes the HLSVD input and output
            # data to a file in your home folder. It's great for
            # generating test data with which to exercise HLSVD outside
            # of Vespa, but not useful otherwise.
            # DON'T EVER COMMIT THIS FILE WITH THIS CODE ACTIVE!
            # It will destroy performance and cause Vespa users to wail
            # and gnash their teeth.
            #
            # import vespa.common.util.xml_ as util_xml
            # import os.path
            #
            # input_ = {"dwell_time" : dwell_time,
            #           "nsv_sought" : nsv_sought,
            #           "signals" : signals
            #          }
            #
            # output = { "nsv_found" : results[0],
            #            "singular_values" : results[1],
            #            "frequencies" : results[2],
            #            "damping_factors" : results[3],
            #            "amplitudes" : results[4],
            #            "phases" : results[5],
            #          }
            #
            # d = {"input" : input_, "output" : output }
            #
            # filename = chain.dataset.dataset_filename
            # if not filename:
            #     # Construct a filename from the raw filename.
            #     raw = chain.dataset.blocks["raw"]
            #     filename = raw.data_source
            #
            # path, filename = os.path.split(filename)
            # filename = os.path.splitext(filename)[0] + ".xml"
            #
            # home = os.path.expanduser("~")
            # filename = os.path.join(home, filename)
            #
            # util_xml._dict_to_file(d, filename)


        # create the fids for each element in the HLSVD model
        #
        # Note. Here we do NOT deal with any frequency shifts that are applied
        # to the raw data. This is dealt with when the water filter itself is
        # applied. At that point the individual HLSVD fids can be sorted to
        # determine which is below the threshold, or not, or maybe has been
        # manually selected in the SVD Filter sub-tab.
        svd_fids = _create_hlsvd_fids(chain.svd_output.frequencies,
                                      chain.svd_output.damping_factors,
                                      chain.svd_output.amplitudes,
                                      chain.svd_output.phases,
                                      chain.raw_dim0,
                                      len(chain.svd_output),
                                      0.0,
                                      dwell_time)

        chain.svd_fids_all = svd_fids             # previously chain.time_fids



        # CALCULATE AUTO SELECT for FIDs ---------------------------
        #
        #  Apply rules here so it is part of pipline not the GUI 

        fshift = chain.svd_output.frequencies + (chain.frequency_shift / 1000.0)
        ppms   = chain.resppm - (fshift * 1000.0 / chain.frequency)

        if set.svd_apply_threshold:
            chain.svd_output.in_model.fill(False)
            
        if set.svd_apply_threshold:
            if set.svd_threshold_unit == 'Hz':
                width = set.svd_threshold / 1000.0
                for i,freq in enumerate(fshift):
                    if freq <= width:
                        chain.svd_output.in_model[i] = True
            else:
                width = set.svd_threshold 
                for i,ppm in enumerate(ppms):
                    if ppm >= width:
                        chain.svd_output.in_model[i] = True
                
        if set.svd_exclude_lipid:
            for i,ppm in enumerate(ppms):
                if ppm >= set.svd_exclude_lipid_end and ppm <= set.svd_exclude_lipid_start:
                    chain.svd_output.in_model[i] = True


        # here we calculate and save svd_fids_checked
        in_model = np.where(chain.svd_output.in_model)[0]
        cnt = len(in_model)
        if cnt > 0:
            if cnt > 1:
                sum_fids = np.sum(chain.svd_fids_all[in_model,:],axis=0)
            else:
                sum_fids = chain.svd_fids_all[in_model[0],:]
        else:
            # no terms on in the filter
            sum_fids = np.zeros(chain.raw_dim0, complex)

        chain.svd_fids_checked = sum_fids              # previously sum_time_fids


        # FREQUENCY_SHIFT --------------------

        if chain.frequency_shift != 0.0:

            # seterr() avoids underflow error
            old_err_state = np.seterr(all='ignore')
            t = np.arange(chain.raw_dim0) / chain.sw
            phroll = np.exp(1j * 2.0 * np.pi * chain.frequency_shift * t)
            svd_fids = svd_fids * phroll
            svd_data = svd_data * phroll 
            chain.svd_fids_checked = chain.svd_fids_checked * phroll
            np.seterr(**old_err_state)

        # APODIZE ----------------------------
        
        if set.apodization:

            t = np.arange(chain.raw_dim0) / chain.sw
            old_settings = np.seterr(all='ignore')  # avoid underflow warnings
            if set.apodization == 'gaussian':
                lineshape = safe_exp(-(set.apodization_width * np.pi * 0.6 * t) ** 2, 0)
            elif set.apodization == 'lorentzian':
                lineshape = safe_exp(-(set.apodization_width * np.pi * t), 0)
            svd_fids = svd_fids * lineshape
            svd_data = svd_data * lineshape
            np.seterr(**old_settings)

        # CHOP ------------------------------

        if set.chop:
            chop = ((((np.arange(chain.raw_dim0) + 1) % 2) * 2) - 1)
            svd_fids = svd_fids * chop
            svd_data = svd_data * chop

        # ZEROFILL ---------------------------

        dim0 = chain.raw_dim0 * set.zero_fill_multiplier

        if chain.is_fid:
            if len(svd_fids.shape) > 1:
                svd_fids[:,0] *= 0.5
            else:
                svd_fids[0] *= 0.5

        svd_data[0] *= 0.5

        # FFT ---------------------------

        svd_fids = np.fft.fft(svd_fids, n=dim0) / dim0   #svd_fids.shape[-1]
        svd_data = np.fft.fft(svd_data, n=dim0) / dim0   #svd_data.shape[-1]

        # FLIP_SPECTRAL_AXIS -----------------

        if set.flip:

            if len(svd_fids.shape)>1:
                for i in range(svd_fids.shape[0]):
                    svd_fids[i,:] = svd_fids[i,::-1].copy()
            else:
                svd_fids = svd_fids[::-1].copy()

            svd_data = svd_data[::-1].copy()

        # AMPLITUDE and DC OFFSET  ---------------------------

        svd_fids = svd_fids * set.amplitude + set.dc_offset
        svd_data = svd_data * set.amplitude + set.dc_offset


        # SAVE THINGS ----------------------------

        chain.svd_data = svd_data.copy()      # previously chain.freq_svd

        in_model = np.where(chain.svd_output.in_model)[0]
        cnt = len(in_model)
        if cnt > 0:
            if cnt > 1:
                svd_fids      = svd_fids[in_model,:]
                svd_peaks_sum = np.sum(svd_fids,axis=0)
            else:
                svd_fids      = svd_fids[in_model[0],:]
                svd_peaks_sum = svd_fids.copy()
        else:
            svd_fids      = np.zeros(chain.spectral_dim0, complex)
            svd_peaks_sum = np.zeros(chain.spectral_dim0, complex)      # no terms on in the filter

        if len(svd_fids.shape) <= 1:
            svd_fids.shape = 1,svd_fids.shape[0]

        chain.svd_peaks_checked     = svd_fids.copy()       # previously chain.fids
        chain.svd_peaks_checked_sum = svd_peaks_sum         # previously chain.sum_fids


def _create_hlsvd_fids(freqs, decays, areas, phases, acqdim0, nlines, toff, dwell_time):
    """
    Construct time domain signal from the estimated parameters.

    We frequently force the exp() function here into very small values that
    raise a "FloatingPointError: underflow encountered in exp"
    Try/except is not workable as it interrupts the calculation. Instead,
    we temporarily disable numpy's error reporting, allow it to silently
    generate NaNs, and then change those NaNs to 0.0.

    A cleaner way to squelch these errors would be to call util_math.safe_exp().
    But the resulting array can trigger underflow when multiplying areas[i] by
    the result, so that call alone is not sufficient to stop underflow errors,
    hence we ignore errors for all of small code block below.

    """
    result = np.zeros((nlines, acqdim0),dtype=complex)
    t = np.arange(acqdim0) * dwell_time + toff
    K = 1j * 2 * np.pi      # K is invariant in loop below

    for i in range(nlines):
        if decays[i] == 0:
            result[i,:] = result[i,:] * 0
        else:
            old_settings = np.seterr(all='ignore')
            line_result = areas[i] * np.exp((t/decays[i]) + K * (freqs[i]*t+phases[i]/360.0))
            zeros = np.zeros_like(line_result)
            result[i,:] = np.where(np.isnan(line_result), zeros, line_result)
            np.seterr(**old_settings)

    return result



def do_processing_all(chain):
    """

    """
    set = chain._block.set

    # eddy current correction
    if chain._block.set.ecc_method != 'None':
        if chain._block.set.ecc_dataset is not None:
            funct_ecc.do_ecc_processing(chain)

    # save current data as 'pre_roll'
    chain.pre_roll = chain.data.copy()

    # left shift data
    if set.left_shift_value: left_shift(chain)

    # frequency shift data
    if chain.frequency_shift != 0.0: frequency_shift(chain)

    # calculate svd filter
    svd_filter(chain)

    # apply water filter
    if chain._block.set.water_filter_method != 'None':
        funct_water.do_water_filter_processing(chain)

    # apodize data
    if set.apodization: apodization(chain)

    # save current data as 'kodata'
    chain.kodata = chain.data.copy()

    # shift data by half
    if set.chop: chop(chain)

    fft(chain)

    # physically flip data
    if set.flip: flip_spectral_axis(chain)

    # perform scaling or dc offset
    chain.data = chain.data * set.amplitude + set.dc_offset

    # save current data as 'freq'
    chain.freq = chain.data.copy()










