# Python Modules.

import copy
import math

# 3rd party stuff
import numpy as np

# Local imports
import vespa.common.rfp_rf_result as rfp_rf_result
from vespa.common.transform_run_exception import TransformRunException 

from pylab import *


def run(trans_desc):
#def assaf_goia_wurst(B1max, Q, tp, Gmax, npts, f=0.9, wurst_n=16, wurst_m=4):
    r"""
    Stub Summary
    ------------------
    This is a template that you can use as a starting point to develop 
    algorithm code for both Create and Modify types of Transforms. 
    
    It has the same format needed for cut/paste insertion directly into 
    the Transform Editor. 
    
    It is also set up so you can call it directly from the command line. 
    The code in the   __name__ == '__main__'  section will create the
    'trans_desc' dictionary based on default values you provide, and 
    the resultant rf waveform, gradient and time axis will be used to 
    run a bloch simulation to create a frequency profile for display. 
    
    Create or Modify Transform Information
    ---------------------------------------------
    Listed below are the input parameters and extra parameters being passed
    into this script.

    time_steps      - int, the number of data points in the output waveform
    duration        - float, total time of pulse, in msec
    
    resolution      - int, the number of points used in calculating the 
                        frequency domain computations.
    bandwidth_convention    - int, choice
        0 for "conventional" (default, FWHM)
        1 for spectroscopy convention (MINIMUM)
        2 for filter convention (MAXIMUM)
    
    slr_pulse returns a tuple of 4 items defined below:
        rf_waveform  - ndarray, complex, real/imag values of the rf_pulse.
        rf_xaxis     - ndarray, float, same size as rf_waveform contains 
                                the corresponding x axis points.
        gradient     - ndarray, float OR None
        grad_xaxis   - ndarray, float OR None

    May raise (or throw):
    
        TransformRunException - to indicate user calculated failure of 
                                algorithm or parameters

    Notes
    -------------------------------
    from Andronesi et. al., JMR 203:283-293 (2010)
    
    We will follow the notation in that paper (and fix some typos along the way).

    """
    param = trans_desc.parameters
    extra = trans_desc.extra

    #--------------------------------------------------------------------------
    # unpack and convert parameters 
    #
    # - Add or Modify the values below as needed

    time_steps      = int(param["time_steps"])          # int
    duration        = float(param["duration"])          # float, msec
    b1_max          = float(param['b1_max'])
    adiabat_factor  = float(param['adiabat_factor'])
    gradient_max    = float(param['gradient_max'])
    grad_mod_factor = float(param['grad_mod_factor'])
    wurst_n         = float(param['wurst_n'])        # float, in meters
    wurst_m         = float(param['wurst_m'])       # float, in meters
    
    dwell_time      = (1000 * duration) / (time_steps)     # in usec
    
    # these extra items may be used in making waveform or profile
    bandwidth_convention = int(extra['pulse_bandwidth_type'])  # Choice, 0-HalfHeight, 1-Min, 2-Max    
    resolution           = int(extra['calc_resolution']) # int
    gamma                = float(extra['gamma'])         # float MHz/T -> 2pi * 100 -> Hz/gauss


    npts   = time_steps         #1024       # number of time steps
    tp     = duration           #3.5        # pulse duration
    B1max  = b1_max             #0.817      # kHz
    Q      = adiabat_factor     #0.5        # adiabaticity factor, should be >> 1
    Gmax   = gradient_max       #0.1        # max gradient, in kHz/mm (multiply by 1000/42.57=23.5 to convert to mT/m)
    f      = grad_mod_factor    #0.9        # Gradient modulation factor. The higher f, the more G "drops" in the
                                #   middle of the pulse and the less B1max is required. f=0.9 is typical.
    #wurst_n = 16        # We *constrain* the B1(t) profile and gradient profile to be WURST-like
    #wurst_m = 4         #   in their shapes


    # Use TransformRunException to check parameter ranges and alert users
    # - an example is shown below
    # - when your code returns with this exception, it fails 'gracefully'
    # - no result is changed by this exception 
    # - an error dialog shows the message passed back
    # - can be used *anywhere* in code to fail and pass info back to the user

    if time_steps < 1:
        error_msg =  "The number of time steps in the pulse must be > 0"
        raise TransformRunException(error_msg, 1)


    #--------------------------------------------------------------------------
    # transform algorithm code starts here


    # Define a "time-vector"
    dt = tp/npts
    tt = np.arange(npts)*dt
    ttNorm = 2*(tt/tp)-1          # Runs from -1 to +1

    b1 = B1max*(1 - np.abs(np.sin(np.pi/2*(ttNorm)))**wurst_n )
    g  = Gmax*((1-f) + f * np.abs(np.sin(np.pi/2*(ttNorm)))**wurst_m )

    # Calculate the instantaneous RF frequency
    x = b1**2/g*dt
    x = x.cumsum(axis=0)
    x = x - x[-1]/2 
    vrf = 1/Q*g*x

    temp = vrf*dt
    phi = 2*np.pi*temp.cumsum(axis=0)
    phi = phi - np.min(phi)

    b1c = (b1/4.2580) * 0.1 * np.exp(1j*phi)    # convert kHz to G  (4.258 kHz/Gauss) then G to mT
    gr  = g * 23.5                              # convert kHz/mm to mT/m 
    b1x = tt.copy() / 1000.0
    gx  = tt.copy() / 1000.0
    outs = {}

    outs['bw_calc'] = gr[0] * 4358.0 * (Gmax*100)    # G/cm  Hz/G  cm 
    #outs['bw_calc'] = (g[0] * 0.1) * 4358.0    # (mT/m->G/cm)  Hz/G  cm 
    print("BW = "+str(outs['bw_calc']))

    
    return b1c, b1x, gr, gx, outs
    
    
#------------------------------------------------------------

if __name__ == '__main__':
    
    import vespa.public.transform_description as transform_description
    import vespa.common.constants as constants
    import pylab as pyl
    import vespa.common.pulse_funcs.bloch_lib_hargreaves as bloch
    
    #----------------------------------------------------------------
    # list your default values for parameters and extra variables
    
    extra = {}
    extra['calc_resolution']      = 2000
    extra['pulse_bandwidth_type'] = 0
    extra['gamma']                = 42.576  # 1H gamma
    
    params = {}

    params["time_steps"]        = 1024      # int, 
    params["duration"]          = 3.5       # float, msec
    params['b1_max']            = 0.817     # kHz
    params['adiabat_factor']    = 0.5       # adiabaticity factor, should be >> 1
    params['gradient_max']      = 0.1       # max gradient, in kHz/mm (multiply by 1000/42.57=23.5 to convert to mT/m)
    params['grad_mod_factor']   = 0.9       # Gradient modulation factor. The higher f, the more G "drops" in the
                                            #   middle of the pulse and the less B1max is required. f=0.9 is typical.
    params['wurst_n']           = 16        # We *constrain* the B1(t) profile and gradient profile to be WURST-like
    params['wurst_m']           = 6         #   in their shapes


    trans_desc = transform_description.TransformDescription()
    trans_desc.extra                = extra
    trans_desc.parameters           = params

    # Run the algorithm

    rfy, rfx, gry, grx, outs = run(trans_desc)
    
    # Store results into object
    
    result = rfp_rf_result.RfResults()
    result.rf_waveform = np.array(rfy)  # this is in mT
    result.rf_xaxis    = np.array(rfx)  # this is in sec  
    if gry is not None:
        gry = np.array(gry)            # this is in mT/m
    if grx is not None:
        grx = np.array(grx)
    result.gradient    = gry
    result.grad_xaxis  = grx            # this is in sec

    # Call Bloch simulation on rf waveform
    # - need this dict to use Vespa bloch lib calls

    bloch_inputs = {}
    bloch_inputs['calc_resolution'] = extra['calc_resolution'] 
    bloch_inputs['gamma']             = extra['gamma'] 
    
    bloch_inputs['bloch_range_value'] = 25.0    # unit set below
    bloch_inputs['bloch_range_units'] = 'kHz'   # 'cm' or 'kHz'
    bloch_inputs['bloch_offset_value'] = 0.0    # in hz

    result.update_profiles(bloch_inputs)
        
    #----------------------------------------------------------------
    # Everything below here is for displaying results
    
    # use get_profile() to get complex magnetization to display
    # - other UsageTypes include SATURATION, INVERSION, SPINE_ECHO

    dat1y, dat1x = result.get_profile(constants.UsageType.EXCITE, False, False)
    dat2y, dat2x = result.get_profile(constants.UsageType.SATURATION, False, False)

    fig = pyl.figure()

    axes1 = fig.add_subplot(411)
    axes1.set_xlabel('Time [ms]')
    axes1.set_ylabel('[mTesla]')
    axes1.plot(rfx*1000, np.real(rfy), rfx*1000, np.imag(rfy))

    axes2 = fig.add_subplot(412)
    axes2.set_xlabel('Time [ms]')
    axes2.set_ylabel('[G/cm]')
    axes2.plot(rfx*1000, gry)
    
    axes3 = fig.add_subplot(413)
    axes3.set_xlabel('Position [cm]')
    axes3.set_ylabel('Magn |Mxy|')
    axes3.plot(dat1x, np.abs(dat1y), dat1x, np.real(dat1y), dat1x, np.imag(dat1y) )

    axes4 = fig.add_subplot(414)
    axes4.set_xlabel('Position [cm]')
    axes4.set_ylabel('Magn Mz')
    axes4.plot(dat2x, np.real(dat2y))

    pyl.show()
   
    bob = 10
    
    
