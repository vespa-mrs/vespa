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

    <ADD COMMENTS HERE ABOUT YOUR CODE, PLEASE>

    """
    param = trans_desc.parameters
    extra = trans_desc.extra

    #--------------------------------------------------------------------------
    # unpack and convert parameters 
    #
    # - Add or Modify the values below as needed

    time_steps      = int(param["time_steps"])      # int
    duration        = float(param["duration"])      # float, msec
    
    dwell_time      = (1000 * duration) / (time_steps)     # in usec
    
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

    rf_waveform = np.zeros([time_steps,], dtype=np.complex)
    rf_xaxis    = np.arange([time_steps,]) * dwell_time
    
    gradient    = None
    grad_xaxis  = None
    
    # end transform algorithm
    #--------------------------------------------------------------------------

    #--------------------------------------------------------------------------
    # fill an output dictionary, or leave as None

    outputs = None

    #--------------------------------------------------------------------------
    # ensure you are returning ndarrays or None

    rf_y = np.array(rf_waveform)
    rf_x = np.array(rf_xaxis)
    gr_y = gradient
    gr_x = grad_xaxis


    return rf_y, rf_x, gr_y, gr_x, outputs
    

    
    
#------------------------------------------------------------

if __name__ == '__main__':
    
    import vespa.public.transform_description as transform_description
    import vespa.common.constants as constants
    import pylab as pyl
    
    #----------------------------------------------------------------
    # list your default values for parameters and extra variables
    
    extra = {}
    extra['calc_resolution']      = 2000
    extra['pulse_bandwidth_type'] = 0
    extra['gamma']                = 42.576  # 1H gamma
    
    params = {}
    params['time_steps']    = 250       
    params['duration']      = 8.0
    
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
    bloch_inputs['bloch_range_value'] = 12.0    # unit set below
    bloch_inputs['bloch_range_units'] = 'kHz'    # 'cm' or 'khz'
    bloch_inputs['bloch_offset_value'] = 0.0    # in hz

    result.update_profiles(bloch_inputs)    
    #----------------------------------------------------------------
    # Everything below here is for displaying results
    
    # use get_profile() to get complex magnetization to display
    # - other UsageTypes include SATURATION, INVERSION, SPINE_ECHO
    
    pry, prx = result.get_profile(constants.UsageType.EXCITE, False, False)

    fig = pyl.figure()

    axes1 = fig.add_subplot(211)
    axes1.set_xlabel('Time [ms]')
    axes1.set_ylabel('[mTesla]')
    axes1.plot(rfx*1000, rfy)
    
    axes2 = fig.add_subplot(212)
    dat = abs(pry)
    pad = dat.max() * 0.05    
    axes2.set_autoscaley_on(False)
    axes2.set_ylim([-pad,dat.max()+pad])
    axes2.set_xlabel('Position [cm]')
    axes2.set_ylabel('Magnetization |Mz|')
    axes2.plot(prx, abs(pry))

#     axes3 = fig.add_subplot(313)
#     dat = abs(gry)
#     pad = dat.max() * 0.05    
#     axes3.set_autoscaley_on(False)
#     axes3.set_ylim([-pad,dat.max()+pad])
#     axes3.set_xlabel('Position [cm]')
#     axes3.set_ylabel('Gradient')
#     axes3.plot(grx, abs(gry))
    
    pyl.show()
    
    bob = 10
    
    
