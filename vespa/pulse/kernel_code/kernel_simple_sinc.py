# Python Modules.

import copy
import math

# 3rd party stuff
import numpy as np

# Local imports
import vespa.common.rfp_rf_result as rfp_rf_result
from vespa.common.transform_run_exception import TransformRunException


# Gyromagnetic ratio of 1H - the hydrogen nucleus. (units: kHz/mT)
GAMMA1H = 42.576    



def run(trans_desc):
    """
    Simple Sinc Summary
    ---------------------
    I started with the kernel_stub.py template as a starting point. I added
    two input parameters ('ncycles' and 'sfilter') and forced them to be
    specific types. I added code to calculate the rf_waveform (rf_y) output 
    parameter and did not change the rf_xaxis (rf_x), gradient (gr_y) or 
    gradient_xaxis (gr_x) parameters already caculated in the stub.
    
    Using this template, I get the format I need to cut/paste the code
    directly into the Transform Editor. But, I keep command line 
    convenience, too, in case my debugging is more complicated than the
    Vespa GUI allows. 
    
    You can call this kernel_code module directly from the command line. 
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
    
    
    slr_pulse returns a tuple of 4 items defined below:
        rf_waveform  - ndarray, complex, real/imag values of the rf_pulse in mT.
        rf_xaxis     - ndarray, float, same size as rf_waveform contains 
                                the corresponding x axis time points in seconds.
        gradient     - ndarray, float OR None, units are mT/m
        grad_xaxis   - ndarray, float OR None (usually same as rf_axis ...)

    May raise (or throw):
    
        TransformRunException - to indicate user calculated failure of 
                                algorithm or parameters

    Notes
    -------------------------------

    This is the same as the 'hsinc' function I use to test my bloch_lib_hargreaves
    module.  It returns a sinc function of length npts, with ncycles sinc-cycles. 
    This yields a time-bandwidth value of 4 * ncycles

    """
    param = trans_desc.parameters
    extra = trans_desc.extra

    #--------------------------------------------------------------------------
    # unpack and convert parameters 
    #
    # - Add or Modify the values below as needed
    # - Don't forget to add them to __main__ as well
    
    time_steps      = int(param["time_steps"])      # int
    duration        = float(param["duration"])      # float, msec
    ncycles         = int(param["ncycles"])         # int, positive
    sfilter         = param["filter"]               # string, 'hamming' or nothing
    
    dwell_time      = (duration * 0.001) / (time_steps)     # in sec

    # these extra items may be used in making waveform or profile
    bandwidth_convention = int(extra['pulse_bandwidth_type'])  # Choice, 0-HalfHeight, 1-Min, 2-Max    
    resolution           = int(extra['calc_resolution']) # int
    gamma                = float(extra['gamma'])         # float MHz/T -> 2pi * 100 -> Hz/gauss
    

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

    # it is best if the waveform array is complex, even if the waveform is real only
    rf_waveform = np.zeros([time_steps,], dtype=np.complex)
    
    npts = time_steps - 8       # 8 is for the padding
    
    t = np.arange(npts) - (npts/2.0)
    t = t / (npts/2.0)
    val = 2*np.pi*ncycles*t + 0.00001
    rf_y = np.sin(val) / val
    if sfilter == 'hamming':
        rf_y = rf_y * 4 * ncycles * (0.54 + 0.46*np.cos(np.pi*t)) / npts
    
    rf_y = 0.05 * rf_y/np.max(rf_y)                                         # normalize and set to T here
    rf_y = np.concatenate((np.zeros(4), rf_y, np.zeros(4)))                 # pad edges
    rf_x = np.arange(time_steps) * dwell_time                               # in sec
    
    gr_y = np.concatenate((np.zeros(4), 2.5*np.ones(npts), np.zeros(4)))    # set to mT/m here
    gr_x = rf_x.copy()
    
    # end transform algorithm
    #--------------------------------------------------------------------------

    #--------------------------------------------------------------------------
    # fill an output dictionary, or leave as None

    outputs = None

    #--------------------------------------------------------------------------
    # ensure you are returning ndarrays or None

    rf_y = np.array(rf_y)
    rf_x = np.array(rf_x)
    gr_y = np.array(gr_y)
    gr_x = np.array(gr_x)


    return rf_y, rf_x, gr_y, gr_x, outputs
    

    
    
#------------------------------------------------------------

if __name__ == '__main__':
    
    import vespa.public.transform_description as transform_description
    import vespa.common.constants as constants
    import pylab as pyl 
    
    #----------------------------------------------------------------
    # list your default values for parameters and extra variables
    
    npts  = 128
    dwell = 0.02            # msec
    durat = npts * dwell    # msec
    
    extra = {}
    extra['calc_resolution']      = 2000    # not used in this example
    extra['pulse_bandwidth_type'] = 0       # not used in this example
    extra['gamma']                = 42.576  # 1H gamma
    
    params = {}
    params['time_steps']    = npts       
    params['duration']      = durat
    params['ncycles']       = 6
    params['filter']        = 'hamming'
    
    trans_desc = transform_description.TransformDescription()
    trans_desc.extra                = extra
    trans_desc.parameters           = params

    # Run the algorithm

    rfy, rfx, gry, grx, outs = run(trans_desc)
    
    # Store results into object
    
    result = rfp_rf_result.RfResults()
    result.rf_waveform = np.array(rfy)
    result.rf_xaxis    = np.array(rfx)
    if gry is not None:
        gry = np.array(gry)
    if grx is not None:
        grx = np.array(grx)
    result.gradient    = gry
    result.grad_xaxis  = grx

    # Call Bloch simulation on rf waveform
    # - need this dict to use Vespa bloch lib calls

    bloch_inputs = {}
    bloch_inputs['calc_resolution'] = extra['calc_resolution'] 
    bloch_inputs['gamma']             = extra['gamma'] 
    bloch_inputs['bloch_range_value'] = 25.0    # unit set below
    bloch_inputs['bloch_range_units'] = 'cm'    # 'cm' or 'khz'
    bloch_inputs['bloch_offset_value'] = 0.0    # in hz

    result.update_profiles(bloch_inputs)
    
    #----------------------------------------------------------------
    # Everything below here is for displaying results
    
    # use get_profile() to get complex magnetization to display
    # - other UsageTypes include SATURATION, INVERSION, SPINE_ECHO
    
    dat1y, dat1x = result.get_profile(constants.UsageType.EXCITE, False, False)
    dat2y, dat2x = result.get_profile(constants.UsageType.SATURATION, False, False)

    fig = pyl.figure()

    axes1 = fig.add_subplot(311)
    axes1.set_xlabel('Time [ms]')
    axes1.set_ylabel('[mTesla]')
    axes1.plot(rfx*1000, rfy)
    
    axes2 = fig.add_subplot(312)
    axes2.set_xlabel('Position [cm]')
    axes2.set_ylabel('Magn |Mxy|')
    axes2.plot(dat1x, np.abs(dat1y), dat1x, np.real(dat1y), dat1x, np.imag(dat1y) )

    axes3 = fig.add_subplot(313)
    axes3.set_xlabel('Position [cm]')
    axes3.set_ylabel('Magn Mz')
    axes3.plot(dat2x, np.real(dat2y))
    
    pyl.show()
    
    bob = 10
    
    
