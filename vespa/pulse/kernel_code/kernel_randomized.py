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
    Script for pulse generation. Adjusts some input parameters to be in the 
    correct units.

    time_steps      - int, the number of data points in the output waveform
    duration        - float, total time of pulse, in msec
    
    resolution      - the calculation resolution for frequency domain 
                      computations.
    bandwidth_convention:
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
    
    Given the number of points expected in the array, returns a complex
    array of microTesla values with each real/imag value randomized in the 
    interval (-1, 1).
        
    See the following articles:

    """
    param = trans_desc.parameters
    extra = trans_desc.extra

    #--------------------------------------------------------------------------
    # unpack and convert parameters if needed

    time_steps      = int(param["time_steps"])      # int
    duration        = float(param["duration"])      # float, msec
    
    dwell_time      = (1000 * duration) / (time_steps)     # in usec
    
    # these extra items are used in making profile
    bandwidth_convention = int(extra['pulse_bandwidth_type'])   # Choice, 0-HalfHeight, 1-Min, 2-Max    
    resolution           = int(extra['calc_resolution'])        # int
    gamma                = float(extra['gamma'])         # float MHz/T -> 2pi * 100 -> Hz/gauss


    #--------------------------------------------------------------------------
    # pulse creation code starts here

    rand = np.random.rand
    y = (2 * rand(time_steps) - 1) + (2j * rand(time_steps) - 1j)
    
    # We want values in microTesla
    rf_waveform = y / complex(1000, 1000)


    # The call to analytic_pulse.randomized() only returned the "y"
    # component of the profile. We generate the x_axis here.
    dwell = dwell_time / 1000000.0
    rf_xaxis = np.arange(time_steps) * dwell

    #--------------------------------------------------------------------------
    # fill and output dictionary, if needed

    outputs = None

    #--------------------------------------------------------------------------
    # ensure you are returning ndarrays or None

    rf_y = np.array(rf_waveform)
    rf_x = np.array(rf_xaxis)
    gr_y = None
    gr_x = None


    return rf_y, rf_x, gr_y, gr_x, outputs
    

    
    
#------------------------------------------------------------

if __name__ == '__main__':
    
    import vespa.public.transform_description as transform_description
    import vespa.common.constants as constants
    import pylab as pyl
    
    extra = {}
    extra['calc_resolution'] = 2000
    extra['pulse_bandwidth_type'] = 0
    extra['gamma']                = 42.576  # 1H gamma    
    
    params = {}
    params['time_steps']    = 250       
    params['duration']      = 8.0
    
    trans_desc = transform_description.TransformDescription()
    trans_desc.extra                = extra
    trans_desc.parameters           = params

    rfy, rfx, gry, grx, outs = run(trans_desc)
    
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
    
    dat1y, dat1x = result.get_profile(constants.UsageType.EXCITE, False, False)
    dat2y, dat2x = result.get_profile(constants.UsageType.SATURATION, False, False)

    fig = pyl.figure()

    axes1 = fig.add_subplot(311)
    axes1.set_xlabel('Time [ms]')
    axes1.set_ylabel('[mTesla]')
    axes1.plot(rfx*1000, rfy)
    
    axes2 = fig.add_subplot(312)
    axes2.set_xlabel('Position [cm]')
    axes2.set_ylabel('Magnetization |Mxz|')
    axes2.plot(dat1x, np.abs(dat1y))  #, dat1x, np.real(dat1y), dat1x, np.imag(dat1y) )

    axes3 = fig.add_subplot(313)
    axes3.set_xlabel('Position [cm]')
    axes3.set_ylabel('Magnetization Mz')
    axes3.plot(dat2x, np.real(dat2y))

    pyl.show()
    
    bob = 10
    
    
