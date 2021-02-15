# Python Modules.

import copy
import math

# 3rd party stuff
import numpy as np
from scipy.interpolate import interp1d

# Local imports
import vespa.common.rfp_rf_result as rfp_rf_result
from vespa.common.transform_run_exception import TransformRunException


# Gyromagnetic ratio of 1H - the hydrogen nucleus. (units: kHz/mT)
# GAMMA1H = 42.576

# Size limit on b1 for root reflection
# b1 to polynomial will automagically truncate to this
# so calling code should check and issue a warning
b1rootlimit = 65

# Small number used for floating point comparisons.
epsilon = 0.000001  

# Very small number used for double precision comparisons.
small_epsilon = 0.00000000001  

# An even smaller number; used to represent an acceptable fractional error
# or difference from an expected or desired value.
EPS = pow(2,-52)  



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

    Interpolation routine does a two sided 
    interpolation, first using the values from one
    side of the dwell_time range, and then the other,
    and averaging the two. 
    
    This type of interpolation is appropriate for a 
    histogram type graph.
    
    Attempts to cover the current duration.
    Assumes xvals are equally spaced.


    """
    param = trans_desc.parameters
    extra = trans_desc.extra

    previous_rf        = copy.deepcopy(trans_desc.previous_rf)
    previous_rf_xaxis  = copy.deepcopy(trans_desc.previous_rf_xaxis)
    previous_dwell     = previous_rf_xaxis[1] - previous_rf_xaxis[0]
    
    #--------------------------------------------------------------------------
    # unpack and convert parameters if needed

    do_interpolate  = int(param["do_interpolate"])      # int
    new_dwell_time  = float(param["new_dwell_time"])    # float, usec

    do_rescaling    = int(param["do_rescaling"])        # int
    angle           = float(param["new_angle"])         # float, deg
    
    # these extra items may be used in making waveform or profile
    bandwidth_convention = int(extra['pulse_bandwidth_type'])  # Choice, 0-HalfHeight, 1-Min, 2-Max    
    resolution           = int(extra['calc_resolution']) # int
    gamma                = float(extra['gamma'])         # float MHz/T -> 2pi * 100 -> Hz/gauss


    #--------------------------------------------------------------------------
    # pulse creation code starts here

    if do_interpolate or do_rescaling:
    
        if do_interpolate:                
            y, x = interpolate_dwell(previous_rf, previous_rf_xaxis, new_dwell_time/1000000.0)                
            previous_rf = y
            previous_rf_xaxis = x
            
        if do_rescaling:
            dwell_time = previous_rf_xaxis[1] - previous_rf_xaxis[0]
            dwell_time *= 1e6   # need this in usec
            y = rescale(previous_rf, angle, dwell_time, gamma)
            previous_rf = y
                  
    #--------------------------------------------------------------------------
    # fill and output dictionary, if needed

    outputs = {}
    outputs['previous_dwell_time'] = previous_dwell

    #--------------------------------------------------------------------------
    # ensure you are returning ndarrays or None

    rf_y = np.array(previous_rf)
    rf_x = np.array(previous_rf_xaxis)
    gr_y = None
    gr_x = None

    return rf_y, rf_x, gr_y, gr_x, outputs
 

def interpolate_dwell(yvals, xvals, new_dwell):
    """
    This interpolation routine does a two sided 
    interpolation, first using the values from one
    side of the dwell_time range, and then the other,
    and averaging the two. 
    
    This type of interpolation is appropriate for a 
    histogram type graph.
    
    Attempts to cover the current duration.
    Assumes xvals are equally spaced.
    
    May raise PulseFuncException.
    """
    
    length = len(xvals)
    if length < 2:
        error_msg =  "Cannot do a meaningful linear interpolation with less than two points"
        raise TransformRunException(error_msg, 1)
    
    if not (new_dwell > 0):
        error_msg =  "new dwell time must be greater than zero"
        raise TransformRunException(error_msg, 1)       
    
    if not (xvals[1] > xvals[0]):
        error_msg =  "xvals[1] must be greater than xvals[0]"
        raise TransformRunException(error_msg, 1)    
    
    old_dwell = (xvals[1] - xvals[0])
    duration = xvals[length-1] + old_dwell
    
    time_step = xvals[0]+new_dwell
    # make sure at least 2 points.
    x = [xvals[0], time_step]
    while time_step + new_dwell < duration - small_epsilon:
        time_step += new_dwell
        x.append(time_step)  
        
    new_length = len(x)
    
    x = np.array(x)
    
    # Filling will be done at the very end of our waveform.
    # If we add more points we will have to fill - and it makes
    # the most sense to fill at the end with the value at the end.
    izzy1 = interp1d(xvals, yvals, kind='linear', bounds_error=False, fill_value=yvals[length-1])
    
    # Tried to do this using np.fliplr
    # but received this error:
    # raise ValueError, "Input must be >= 2-d."

    yvals_flipped = _fliplr1d(yvals)
    
    izzy2 = interp1d(xvals, yvals_flipped, kind='linear', bounds_error=False, fill_value=yvals_flipped[length-1])
    
    y1 = izzy1(x)
    y2 = izzy2(x)
    
    y = (y1 + _fliplr1d(y2))/2
    
    return y, x

   
def _fliplr1d(y_in):
    yvals_flip = []
    i = len(y_in)
    while i > 0:
        yvals_flip.append(y_in[i-1])
        i -= 1
    return np.array(yvals_flip)


def rescale(b1b2, angle, dwell_time, gamma):
       
    eps = EPS
    if np.max(np.abs(np.imag(b1b2))) > 10000*eps:
        # If there is any significant imaginary component    
        scale_factor = angle/((gamma) *  (dwell_time/1000)*np.sum(np.abs(b1b2))*360)
    else:
        # Only real, i.e. not significant imaginary component
        scale_factor = angle/((gamma) * (dwell_time/1000)*np.sum(b1b2)*360)
    
    return b1b2 * scale_factor


    
    
#------------------------------------------------------------

if __name__ == '__main__':
    
    import vespa.public.transform_description as transform_description
    import vespa.common.constants as constants
    import pylab as pyl
    
    # create previous RF waveform - simple hamming-sinc
    
    npts      = 128     # int
    duration  = 1.0     # float, msec
    ncycles   = 6
    
    dwell_time = 0.001 * float(duration)/float(npts)
    
    t = (np.arange(npts) - (npts/2.0)) / (npts/2.0)   # time, symmetric about zero
    y = 2 * np.pi * ncycles * t + 0.00001   # 0.0001 to avoid div by zero
    b1 = np.sin(y) / y
    
    # apply a hamming filter
    b1 = b1 * 4 * ncycles * (0.54 + 0.46*np.cos(np.pi*t)) / npts
    
    b1x = np.arange(npts) * dwell_time      # time from zero in sec
    g1x = np.arange(npts) * dwell_time      # time from zero in sec
    
    # scale empirically to 90 degree tip
    b1 = 0.1266667 * b1/np.max(b1)
    gr = np.ones(npts)

    b1 = np.array(b1)
    if not np.iscomplexobj(b1):
        b1r = b1.real.copy()
        b1i = b1r.copy() * 0
        b1 = b1r + 1j*b1i

    # set up transform description params    
    
    extra = {}
    extra['calc_resolution'] = 2000
    extra['pulse_bandwidth_type'] = 0
    extra['gamma']                = 42.576  # 1H gamma
    
    params = {}
    params['do_interpolate']    = 0       
    params['new_dwell_time']    = 15.625
    params['do_rescaling']      = 1       
    params['new_angle']         = 45.0
    
    trans_desc = transform_description.TransformDescription()
    trans_desc.extra            = extra
    trans_desc.parameters       = params

    trans_desc.previous_rf      = b1
    trans_desc.previous_rf_xaxis = b1x

    rfy, rfx, gry, grx, outs = run(trans_desc)
    
    # result for current state RF waveform
    
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
#    result.update_profiles(extra['calc_resolution'])
    
    pry, prx = result.get_profile(constants.UsageType.EXCITE, False, False)

    # result for previous RF waveform

    result0 = rfp_rf_result.RfResults()
    result0.rf_waveform = np.array(b1)
    result0.rf_xaxis    = np.array(b1x)
    if gr is not None:
        gr = np.array(gr)
    if g1x is not None:
        g1x = np.array(g1x)
    result0.gradient    = gr
    result0.grad_xaxis  = g1x

    result0.update_profiles(extra['calc_resolution'])

    pr0y, pr0x = result0.get_profile(constants.UsageType.EXCITE, False, False)
    
    # plot some results

    fig = pyl.figure()

    axes1 = fig.add_subplot(411)
    axes1.set_xlabel('Time [ms]')
    axes1.set_ylabel('[mTesla]')
    axes1.plot(b1x*1000, b1)
    
    axes2 = fig.add_subplot(412)
    dat = abs(pr0y)
    pad = dat.max() * 0.05    
    axes2.set_autoscaley_on(False)
    axes2.set_ylim([-pad,dat.max()+pad])
    axes2.set_xlabel('Position [cm]')
    axes2.set_ylabel('Magnetization |Mz|')
    axes2.plot(pr0x, abs(pr0y))

    axes3 = fig.add_subplot(413)
    axes3.set_xlabel('Time [ms]')
    axes3.set_ylabel('[mTesla]')
    axes3.plot(rfx*1000, rfy)
    
    axes4 = fig.add_subplot(414)
    dat = abs(pry)
    pad = dat.max() * 0.05    
    axes4.set_autoscaley_on(False)
    axes4.set_ylim([-pad,dat.max()+pad])
    axes4.set_xlabel('Position [cm]')
    axes4.set_ylabel('Magnetization |Mz|')
    axes4.plot(prx, abs(pry))

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
    
    
