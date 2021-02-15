# Python modules



# 3rd party modules
import numpy as np


# Our modules
import vespa.common.constants as common_constants
import vespa.common.minf_parabolic_info as minf
import vespa.common.util.generic_spectral as util_spectral




def funct_optimize_phase0(rad, info):
    """
    Optimization function used optimize the zero order phase for FIDs
    in a FidSum dataset. Each FID is compared to the sum of all FIDs.
    A least squared difference is calculated for a phase 0 value that
    best matches the FID to the reference absorption spectrum.
    
    INPUT:
        rad:  current phase in radians
        info: control structure, see optimize_phase0 for definition
    
    """
    phase = np.exp(1j*rad)
    dat   = info['dat'].copy() * phase
    
    istr = info['pts'][0]
    iend = info['pts'][1]
    datt = dat[istr:iend].copy()
    reff = info['ref'][istr:iend].copy()  

    diff = np.sum((datt.real - reff.real)**2)        
    
    return  diff       
   
   
def optimize_phase0(data, modfn, pts):
    """
    Returns the zero order phase in deg at which the least squares
    difference between the data and modfn  absorption spectra are 
    minimized
    
    INPUT:
     data:  array of complex data to be phased
     modfn: string array of model function to be used (obsolete)
     pts:   list, start and end point region used for leastsqr calculation
    
    """
    info = {'dat'     : data,  
            'ref'     : modfn,
            'pts'     : pts     }
    
    # Parabolic interpolation, Brent's method 1-d minimization routine
    
    phase_a = -1*np.pi   # lower bound
    phase_b = np.pi*0.1  # "some" point in the middle
    phase_c =  2*np.pi   # upper bound
    
    phaseat, maxit = minf.minf_parabolic_info( phase_a, phase_b, phase_c,
                                               funct_optimize_phase0, info)
    phdeg = phaseat*180.0/np.pi
    if phdeg > 180.0:
        phdeg = phdeg - 360.0
    
    return phdeg    


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
    info = { 'ta':ta, 'tb':tb, 'orig_lw':lw, 'chain':chain, 'peak':-1.0 }
    
    # Call parabolic interpolation, Brent's method
    # 1-d minimization routine
    
    val_a = 0.005  # lower bound
    val_b = 0.06   # "some" point in the middle
    val_c = 0.5    # upper bound
    
    finalval, maxit = minf.minf_parabolic_info( val_a, val_b, val_c,
                                                _height2area_function, 
                                                info ) 
    return [finalval, info["peak"]]


def do_processing_all(chain):
    """
    Because we are bumping the zero fill factor here by a factor of 4 to 
    better find the peak max, we can not use the standard util_ppm 
    functions to calculate some conversions.  So long hand here for us.
    
    """
    block   = chain._block
    set     = chain._block.set
    dataset = chain._dataset
    
    zfmult = 4      # larger zfmult here improves peak shift accuracy
    raw_dim0 = dataset.raw_dims[0]
    raw_hpp  = dataset.sw / raw_dim0
    fid_dim0 = raw_dim0 * zfmult
    fid_hpp  = dataset.sw / fid_dim0
    
    # reset results arrays and temporary arrays
    work = np.zeros((raw_dim0),complex)
    chain.time_summed_offset = np.zeros((raw_dim0),complex)
    chain.freq_current       = np.zeros((raw_dim0),complex)
    chain.freq_summed        = np.zeros((raw_dim0),complex)
    chain.freq_summed_offset = np.zeros((raw_dim0),complex)
    
    xx = np.arange(raw_dim0) / dataset.sw
    search = np.zeros((raw_dim0 * zfmult),complex)

    # convert algorithm values from PPM to points
    b0_start  = set.reference_peak_center + set.peak_search_width
    b0_end    = set.reference_peak_center - set.peak_search_width
    b0_ctr_pt = (fid_dim0 / 2) - (dataset.frequency * (set.reference_peak_center - dataset.resppm) / fid_hpp)
    b0_start  = int((fid_dim0 / 2) - (dataset.frequency * (b0_start - dataset.resppm) / fid_hpp))
    b0_end    = int((fid_dim0 / 2) - (dataset.frequency * (b0_end - dataset.resppm) / fid_hpp))
    
    ph0_start    = set.phase0_range_start
    ph0_end      = set.phase0_range_end
    ph0_ctr      = ph0_start - 0.5*(ph0_start - ph0_end)
    ph0_start    = int((raw_dim0 / 2) - (dataset.frequency * (ph0_start - dataset.resppm) / raw_hpp))
    ph0_end      = int((raw_dim0 / 2) - (dataset.frequency * (ph0_end   - dataset.resppm) / raw_hpp))

    # one time calculations        
    apod = util_spectral.apodize(xx, set.gaussian_apodization, 'Gaussian')
    chop  = ((((np.arange(raw_dim0) + 1) % 2) * 2) - 1)
    apod *= chop

    nfids = chain.raw.shape[2]

    # Depending on AutoCalc flags, calculate B0 and Phase0 corrections

    if set.apply_peak_shift:   # B0 corrections

        for i in range(nfids):
            time = chain.raw[0,0,i,:].copy()
            if set.fid_left_shift_b0 != 0:
                # shift fid to the left and set last points to zero
                time = np.roll(time, -set.fid_left_shift_b0) 
                time[-set.fid_left_shift_b0:] = time[0]*0.0  
                
            # Calculate peaks shift if flag set, use oversized zfmult
            # Peak search is performed over user-set range on magnitude data
            search *= 0.0
            search[0:raw_dim0] = time * apod
            search = np.fft.fft(search) 
            temp   = np.abs(search)
            imax   = temp[b0_start:b0_end].argmax()
            delta  = (b0_ctr_pt-(b0_start+imax))*fid_hpp
            block.frequency_shift[i] = delta

    if set.apply_phase0:       # Phase0 corrections

        # We need to do this in a second loop, since we need the peaks 
        # shifted before we create a reference spectrum from the summed FIDs

        # Create reference spectrum 
        
        if set.ref_spectrum_source == 'average_all_fids':
            work *= 0
            for i in range(nfids):
                time = chain.raw[0,0,i,:].copy()
                time *= np.exp(1j * 2.0 * np.pi * block.frequency_shift[i] * xx)
                work += time
            
            ref_spec     = work.copy() * apod
            ref_spec[0] *= 0.5
            ref_spec     = (np.fft.fft(ref_spec) / len(ref_spec))
            ref_spec    /= nfids            # scale for comparison to single FID
        else:
            # create reference peak at center of range
            res      = _calc_height2area_ratio( set.ref_peak_line_width, dataset )
            ref_spec = util_spectral.create_spectrum([1.0,], [ph0_ctr,], [0.0,], dataset, ta=res[0], tb=res[0])
        
        # Calc Phase0 correction
        
        for i in range(nfids):
            time = chain.raw[0,0,i,:].copy()

            if set.fid_left_shift_phase0 != 0:     # shift fid left, set last points to zero
                time = np.roll(time, -set.fid_left_shift_phase0) 
                time[-set.fid_left_shift_phase0:] = time[0]*0.0  
            
            time *= np.exp(1j * 2.0 * np.pi * block.frequency_shift[i] * xx)                
            
            time[0] *= 0.5
            tmp_freq = (np.fft.fft(time * apod) / len(time))
            phdeg = optimize_phase0(tmp_freq, ref_spec, [ph0_start, ph0_end])
            block.phase_0[i] = phdeg

    # Apply B0 and Phase0 corrections to raw data to create current and 
    # summed FID arrays with AutoCorr Left Shift values applied
    
    for i in range(nfids):

        time = chain.raw[0,0,i,:].copy()

        if set.fid_left_shift_phase0 != 0:
            time = np.roll(time, -set.fid_left_shift_phase0) 
            time[-set.fid_left_shift_phase0:] = time[0]*0.0  
        elif set.fid_left_shift_b0 != 0:
            time = np.roll(time, -set.fid_left_shift_b0) 
            time[-set.fid_left_shift_b0:] = time[0]*0.0  
        
        time *= np.exp(1j * 2.0 * np.pi * block.frequency_shift[i] * xx)                
        time *= np.exp(1j * block.phase_0[i] * common_constants.DEGREES_TO_RADIANS)

        chain.freq_summed += time
        if i == chain.voxel:
            chain.freq_current += time            
       
        
    # Calculate final summed Time and Freq arrays with constant phase offset

    for i in range(nfids):
        # Apply B0 and Phase0 corrections to raw data to create summed data
        # array for plotting with global Left Shift and Constant Phase0 
        # values applied. This is also the result for this tab

        time = chain.raw[0,0,i,:].copy()

        if set.fid_left_shift != 0:
            # shift fid to the left and set last points to zero
            time = np.roll(time, -set.fid_left_shift) 
            time[-set.fid_left_shift:] = time[0]*0.0  
        
        time *= np.exp(1j * 2.0 * np.pi * block.frequency_shift[i] * xx)                
        time *= np.exp(1j * block.phase_0[i] * common_constants.DEGREES_TO_RADIANS)
        time *= np.exp(1j * set.constant_phase0_offset * common_constants.DEGREES_TO_RADIANS)
        
        chain.freq_summed_offset += time    # for display
        chain.time_summed_offset += time    # result for tab
    
    # Last steps to prepare the display arrays for plotting in Freq domain

    chain.freq_current[0] *= 0.5 
    chain.freq_current = (np.fft.fft(chain.freq_current * apod) / raw_dim0) * nfids  # nfids for comparison plot
                
    chain.freq_summed[0] *= 0.5
    chain.freq_summed  = (np.fft.fft(chain.freq_summed * apod) / raw_dim0)

    chain.freq_summed_offset[0] *= 0.5
    chain.freq_summed_offset  = (np.fft.fft(chain.freq_summed_offset * apod) / raw_dim0)
   

