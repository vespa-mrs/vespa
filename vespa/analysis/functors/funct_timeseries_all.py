# Python modules


# 3rd party modules
import numpy as np


# Our modules
from . import funct_fidsum_coil_combine as funct_combine
import vespa.common.constants as common_constants
import vespa.common.minf_parabolic_info as minf
import vespa.common.util.ppm  as util_ppm
import vespa.common.util.generic_spectral as util_generic_spectral

from vespa.common.constants import DEGREES_TO_RADIANS as DTOR



def funct_optimize_phase0(rad, info):
    """
    Optimization function used optimize the zero order phase for FIDs
    in a  dataset. Each FID is compared to the sum of all FIDs.
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
#    print diff    
    
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


def do_processing_all(chain):
    """
    We process all voxels within this one call because it is not too time 
    consuming and the Spectral/Fit/Quant steps will need them downstream.
    
    We bump the zero fill factor here by a factor of 4 to better find the peak 
    max, so we can not use the standard util_ppm functions to calculate some 
    conversions.  So long hand here for us.
    
    """
    block   = chain._block
    set     = chain._block.set
    dataset = chain._dataset
    
    #------------------------------------------------------
    # Global inits and one time calculations
    
    zfmult = 4      # larger zfmult here improves peak shift accuracy
    raw_dim0 = dataset.raw_dims[0]
    raw_hpp  = dataset.sw / raw_dim0
    fid_dim0 = raw_dim0 * zfmult
    fid_hpp  = dataset.sw / fid_dim0

    xx       = np.arange(raw_dim0) / dataset.sw
    search   = np.zeros((raw_dim0 * zfmult),complex)
    
    # reset results arrays and temporary arrays
    chain.time_all = np.zeros(dataset.raw_dims[::-1],complex)
    chain.freq_all = np.zeros(dataset.raw_dims[::-1],complex)

    # convert algorithm values from PPM to points
    search_start = set.reference_peak_center + set.peak_search_width
    search_end   = set.reference_peak_center - set.peak_search_width
    refpt        = (fid_dim0 / 2) - (dataset.frequency * (set.reference_peak_center - dataset.resppm) / fid_hpp)
    search_start = int((fid_dim0 / 2) - (dataset.frequency * (search_start - dataset.resppm) / fid_hpp))
    search_end   = int((fid_dim0 / 2) - (dataset.frequency * (search_end - dataset.resppm) / fid_hpp))
    
    ph0_start    = set.phase0_range_start
    ph0_end      = set.phase0_range_end
    ph0_start    = int((raw_dim0 / 2) - (dataset.frequency * (ph0_start - dataset.resppm) / raw_hpp))
    ph0_end      = int((raw_dim0 / 2) - (dataset.frequency * (ph0_end   - dataset.resppm) / raw_hpp))

    # one time calculations        
    apod = util_generic_spectral.apodize(xx, set.gaussian_apodization, 'Gaussian')
    chop  = ((((np.arange(raw_dim0) + 1) % 2) * 2) - 1)
    apod *= chop
    
    global_phase0 = np.exp(1j * block.set.global_phase0 * common_constants.DEGREES_TO_RADIANS)

    nfids = chain.raw.shape[2]
    

    #------------------------------------------------------
    # Coil combination section
    #
    # - do not combine if only 1 coil, no need

    if chain.raw.shape[1] > 1:      # number of coils

        if set.coil_combine_method=='None':
            raw_combined = funct_combine.coil_combine_none(chain)
        if set.coil_combine_method=='Siemens':
            raw_combined = funct_combine.coil_combine_siemens(chain)
        if set.coil_combine_method=='CMRR':
            raw_combined = funct_combine.coil_combine_cmrr(chain)
        if set.coil_combine_method=='CMRR-Sequential':
            raw_combined = funct_combine.coil_combine_cmrr_sequential(chain)
    else:
        # single coil, copy first channel only
        raw_combined = funct_combine.coil_combine_none(chain)    



    #------------------------------------------------------
    # FID averaging section
    # 
    #  - if nfids % navgs is not 0, then extra FIDs at the end are ignored
    #  - this step changes the dimensionality of the final data object, all
    #      downstream processing will be affected  
    #

    if set.fids_to_average > 1:      # number of raw FIDs to average into one new FID
        navgs = set.fids_to_average
        nfids = int(nfids/navgs)

        new_shape = list(raw_combined.shape)
        new_shape[-2] = nfids
        time = np.zeros(new_shape, complex)
        
        for i in range(nfids):
            for j in range(navgs):
                time[0,0,i,:] += raw_combined[0,0,j+i*navgs,:]
                
        raw_combined = time
    
    #------------------------------------------------------
    # FID correction section
    #
    #  - global_phase0 and global_phase1 do not affect this algorithm

    for i in range(nfids):
        
        time = raw_combined[0,0,i,:].copy()
        
        if set.fid_left_shift != 0:
            # shift fid to the left and set last points to zero
            time = np.roll(time, -set.fid_left_shift) 
            time[-set.fid_left_shift:] = time[0]*0.0  
        
        if set.apply_peak_shift and chain.calculate_flag:
            # Calculate peaks shift if flag set, use oversized zfmult
            # Peak search is performed over user-set range on magnitude data
            search *= 0.0
            search[0:raw_dim0] = time * apod
            search = np.fft.fft(search) 
            temp   = np.abs(search)
            imax = temp[search_start:search_end].argmax()
            delta = (refpt-(search_start+imax))*fid_hpp
            block.frequency_shift[i] = delta
        
        # Phase 0 NOT calculated here because reference peak has to be 
        # calculated from summed peak-shifted data

        # Apply freq shift and phase0 corrections to the time data
        time *= np.exp(1j * 2.0 * np.pi * block.frequency_shift[i] * xx)
        time *= np.exp(1j * block.phase_0[i] * common_constants.DEGREES_TO_RADIANS)
        
        # Sum up FIDs for display, and calculate current voxel if needed
        chain.time_all[0,0,i,:] = time

    
    #---------------------------------------------------------------------
    # Calculate Phase0 optimization if flag set ON in widget. We need 
    # to do this in a second loop, since we need the peaks shifted before 
    # we create a reference spectrum from the summed FIDs
    #
    #  - global_phase0 DOES affect this algorithm, global_phase1 does not
    #  - ref spectrum that all individual FIDs are adjusted to match is created
    #      from an average of all FIDs with global_phase0 applied
    
    if set.apply_phase0 and chain.calculate_flag:
        # create reference spectrum and optimize range
        freq_all     = chain.time_all.copy() 
        freq_all     = np.sum(freq_all, axis=0) * apod
        freq_all[0] *= 0.5
        freq_all     = (np.fft.fft(freq_all) / len(freq_all))
        freq_all    /= nfids          # scale for comparison to single FID
        freq_all    *= global_phase0  # if there is a global_phase0
        ph0range     = [ph0_start, ph0_end]
        
        # reset global variable so as to fill in below with new ph0 values
        chain.time_all *= 0  
        
        for i in range(nfids):

            time = chain.raw[0,0,i,:].copy()

            if set.fid_left_shift != 0:
                # shift fid to the left and set last points to zero
                time = np.roll(time, -set.fid_left_shift) 
                time[-set.fid_left_shift:] = time[0]*0.0  
            
            time *= np.exp(1j * 2.0 * np.pi * block.frequency_shift[i] * xx)                
            
            # this is where phase 0 is optimized ...
            tmp = time.copy()

                # if global_phase0 apply here before we calculate FID independent phase 0
            tmp *= global_phase0
            
            tmp[0] *= 0.5
            tmp_freq = (np.fft.fft(tmp * apod) / len(tmp))
            phdeg = optimize_phase0(tmp_freq, freq_all, ph0range)
            block.phase_0[i] = phdeg
            
            time *= np.exp(1j * block.phase_0[i] * common_constants.DEGREES_TO_RADIANS)
            
            chain.time_all[0,0,i,:] = time


    if set.global_phase0 != 0.0:
        chain.time_all *= global_phase0

    if set.global_phase1 != 0.0:
        # move all time result into frequency domain
        time_all = chain.time_all.copy()

        # calc phase 1 
        piv    = np.round(dataset.ppm2pts(dataset.phase_1_pivot, acq=True))
        xx     = (np.arange(raw_dim0,dtype=float)-piv)/raw_dim0
        phase1 = np.exp(1j * (set.global_phase1 * DTOR * xx))
        
        for i in range(nfids):
            tmp     = time_all[0,0,i,:].copy()
            tmp[0] *= 0.5        
            tmp     = np.fft.fft(tmp * chop)
            # apply to spectral data and invert fourier transform
            tmp *= phase1
            tmp = np.fft.ifft(tmp)
            chain.time_all[0,0,i,:] = tmp * chop
    
    for i in range(nfids):    
        tmp = apod * chain.time_all[0,0,i,:].copy() 
        tmp[0] *= 0.5
        chain.freq_all[0,0,i,:] = (np.fft.fft(tmp) / len(tmp))
    

 
