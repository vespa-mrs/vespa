# Python modules

import copy

# 3rd party modules
import numpy

# Our modules
import vespa.common.constants as constants
import vespa.common.pulse_funcs.bloch_lib_matpulse as bloch_lib


def calc_all_profiles(rf_amplitude, dwell, calc_resolution):
    """
    Returns the profile, calculated using the Bloch equations. 
    
    Input values:
      rf_amplitude - The rf waveform.
      dwell - The dwell time for this waveform.
      calc_resolution - The number of points used for the Bloch Transformation(s).
      
    Returns:
      mz, mxy, xaxis, mz_ext, mxy_ext, xaxis_ext
      
    NB. The original method was called get_bloch_profile() and allowed the
        user to send in the usage type and receive back only the plots 
        associated with that. We are going to use most of the usage plots
        so here we create all the variations and pass back the 4 y-axis
        results and 2 x-axis results.
        
    NB. This used to loop four time to get these results, I've split this 
        back to four explicit calls to make it easier to see what's going on
      
    """
    khz_to_mm        = 1.0  # Doesn't matter if not being used.                
    convert_to_mm    = 0    # (0 or 1: 0 ==> False)
    resonance_offset = 0.0

    # Call Pypulse.bloch_b1 to get Cayley-Klein parameters for pulse
    # - repeat twice for regular and extended x-range

    # Get regular gradient field to do calcs over
    do_extended_flag = False
    vn,gx = bloch_lib.grad_field(khz_to_mm,
                                 do_extended_flag,
                                 convert_to_mm,
                                 resonance_offset,
                                 dwell,
                                 calc_resolution)

    [alpha,beta] = bloch_lib.bloch_b1(rf_amplitude, 
                                          khz_to_mm, 
                                          do_extended_flag, 
                                          convert_to_mm, 
                                          resonance_offset, 
                                          dwell, 
                                          calc_resolution, vn, gx)
    
    # Get extended gradient field to do calcs over
    do_extended_flag = True
    vnx,gxx = bloch_lib.grad_field(khz_to_mm,
                                   do_extended_flag,
                                   convert_to_mm,
                                   resonance_offset,
                                   dwell,
                                   calc_resolution)    

    [alphax,betax] = bloch_lib.bloch_b1(rf_amplitude, 
                                        khz_to_mm, 
                                        do_extended_flag, 
                                        convert_to_mm, 
                                        resonance_offset, 
                                        dwell, 
                                        calc_resolution, vnx, gxx)

    # Apply Cayley-Klein parameters to get rotated (pulse profile?) magnetizations

    # Set up array of z mag. vectors at each point    
    minit1 = numpy.zeros([calc_resolution, 3], numpy.float64)   
    minit2 = numpy.zeros([calc_resolution, 3], numpy.float64)          
    minit2[:,2] = 1.0       # for exc, inv, sat
    minit1[:,1] = 1.0       # for spin-echo 
    
    xaxis     = vn
    xaxis_ext = vnx     
    mxy       = bloch_lib.mag_recon_vector( alpha, beta, minit1 )
    mz        = bloch_lib.mag_recon_vector( alpha, beta, minit2 )
    mxy_ext   = bloch_lib.mag_recon_vector( alphax, betax, minit1 )
    mz_ext    = bloch_lib.mag_recon_vector( alphax, betax, minit2 )
    
    return (mz, mxy, xaxis, mz_ext, mxy_ext, xaxis_ext) 



def calc_all_profiles_b2g2(b2, g2, dwell, calc_resolution):
    """
    Returns the profile, calculated using the Bloch equations.
    input values:
      b2    - The complex rf waveform.
      g2    - gradient array, same length rf_waveform, in mT
      dwell - The dwell time for this waveform in usec.
      calc_resolution - The number of points used for the Bloch Transformation(s).
      
    returns:
      mz, mxy, xaxis, mz_ext, mxy_ext, xaxis_ext
      
    NB. This used to loop four time to get these results, I've split this 
        back to four explicit calls to make it easier to see what's going on
      
    """
    khz_to_mm = 1.0 # Doesn't matter if not being used.                
    convert_to_mm = 0 # (0 or 1: 0 ==> False)
    resonance_offset = 0.0
    
    rf_waveform = b2


    # Set up array of z mag. vectors at each point    
    minit1 = numpy.zeros([calc_resolution, 3], numpy.float64)   
    minit2 = numpy.zeros([calc_resolution, 3], numpy.float64)          
    minit2[:,2] = 1.0       # for exc, inv, sat
    minit1[:,1] = 1.0       # for spin-echo 

    f2 = numpy.zeros(len(rf_waveform), numpy.float64)

    # Call Pypulse.bloch_b1 to get Cayley-Klein parameters for pulse
    # - repeat twice for regular and extended x-range

    do_extended_flag = False
    [alpha,beta,vn] = bloch_lib.bloch_b2g2(  b2, g2, f2, 
                                              do_extended_flag, 
                                              resonance_offset, 
                                              dwell, 
                                              calc_resolution)


    do_extended_flag = True
    [alphax,betax,vnx] = bloch_lib.bloch_b2g2( b2, g2, f2,
                                                do_extended_flag, 
                                                resonance_offset, 
                                                dwell, 
                                                calc_resolution)

    # Apply Cayley-Klein parameters to get rotated (pulse profile?) magnetizations
       
    result_x  = numpy.zeros(calc_resolution, numpy.float64)     
    result_xx = numpy.zeros(calc_resolution, numpy.float64)
    outmag1   = numpy.zeros([calc_resolution,3], numpy.float64)
    outmag2   = numpy.zeros([calc_resolution,3], numpy.float64)
    outmag1x  = numpy.zeros([calc_resolution,3], numpy.float64)
    outmag2x  = numpy.zeros([calc_resolution,3], numpy.float64)
    for i in range(len(alpha)):
        result_x[i]   = vn[i]
        result_xx[i]  = vnx[i]   
        outmag1[i,:]  = bloch_lib.mag_recon( alpha[i],   beta[i],   minit1[i] )
        outmag2[i,:]  = bloch_lib.mag_recon( alpha[i],   beta[i],   minit2[i] )
        outmag1x[i,:] = bloch_lib.mag_recon( alphax[i],  betax[i],  minit1[i] )
        outmag2x[i,:] = bloch_lib.mag_recon( alphax[i],  betax[i],  minit2[i] )
        
    mz        = copy.copy(outmag2)
    mxy       = copy.copy(outmag1)
    xaxis     = copy.copy(result_x)
    mz_ext    = copy.copy(outmag2x)
    mxy_ext   = copy.copy(outmag1x)
    xaxis_ext = copy.copy(result_xx)

    return (mz, mxy, xaxis, mz_ext, mxy_ext, xaxis_ext) 


