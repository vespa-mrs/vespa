# Python Modules.

import copy
import math

# 3rd party stuff
import numpy as np

# Local imports
import vespa.common.util.math_ as util_math
import vespa.common.rfp_rf_result as rfp_rf_result
from vespa.common.transform_run_exception import TransformRunException


# Gyromagnetic ratio of 1H - the hydrogen nucleus. (units: kHz/mT)
GAMMA1H = 42.576    



def run(trans_desc):
    """
    Script for pulse generation. Adjusts some input parameters to be in the 
    correct units.

    time_steps      - the number of data points in the output waveform
    dwell_time      - the spacing between points in microseconds.
    
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
    ----------
    Bandwidth is set up in Hz, but we need to calc in time domain so we convert
    to time. Note that Gaussian is typically described by a standard deviation
    variable, sigma, rather than full width at half height. So we have to 
    ensure that our 'sigma^2' or 'gb' value, is solved for time at which we are
    at half height. 
    
    y = exp(-(t^2)/(2 * sigma^2))

    solving for t in this equation for t = T(1/2 height) we get a sqrt(ln(2))
    constant, which is the 0.83255461 constant below

    See the following articles:

    """
    param = trans_desc.parameters
    extra = trans_desc.extra

    #--------------------------------------------------------------------------
    # unpack and convert parameters if needed

    time_steps          = int(param["time_steps"])          # int
    duration            = float(param["duration"])          # float, msec
    tip_angle           = float(param["tip_angle"])         # float, deg
    bandwidth           = float(param["bandwidth"])         # float, Hz
    # Choice, 0-None, 1-Cosine, 2-Hamming
    filter_type         = int(param["filter_type"])         # int
    filter_application  = float(param["filter_application"])  # float, percent
    
    dwell_time      = (1000 * duration) / (time_steps)      # in usec
    
    # these extra items are used in making profile
    bandwidth_convention = int(extra['pulse_bandwidth_type'])  # Choice, 0-HalfHeight, 1-Min, 2-Max    
    resolution           = int(extra['calc_resolution']) # int
    gamma                = float(extra['gamma'])         # float MHz/T -> 2pi * 100 -> Hz/gauss


    #--------------------------------------------------------------------------
    # pulse creation code starts here

#     rf_y = analytic_pulse.gaussian3(self.parameters.time_steps,
#                                      dwell_time,
#                                      self.parameters.bandwidth,
#                                      filter_type, filter_application,
#                                      self.parameters.tip_angle)   
#    
#    def gaussian3(npts, dwell_time, bandwidth, filter_type, filter_percent, angle):
    
    npts = time_steps
    filter_percent = filter_application
    angle = tip_angle
    
    dwell = dwell_time*0.000001
    sw    = 1.0 // dwell
    gb    = 1.0 / (1.92 * bandwidth * 1000.0)   
    # this used to be a factor of 2.0, but empirically 1.92 works better
    # to get a full width half max bandwidth of 1kHz for a 90 degree 
    # pulse of 8ms duration and 250 points. This is due to the interaction
    # of a Gaussian pulse with the Bloch equation for something other than
    # a small (1-20 degrees) tip angle. For small tip angles the Fourier 
    # Transform approximation is pretty good for FWHM estimation. As we go
    # to larger tip angles, we have to be more empirical about getting the 
    # actual bandwidth we want. Which is why we plot the Freq Profile for
    # the user to look at and measure.
    
    # Determine whether npts even or odd (even = 1 or 0)
    if npts%2 == 0:
        is_even = True
        xvals = (np.arange(npts//2)+0.5) * dwell 
    else:
        is_even = False
        xvals = np.arange(npts//2 + 1) * dwell
    
    gb = (0.83255461 / gb)
    # the call to exp() tends to underflow with larger durations
    y = util_math.safe_exp(-((gb * xvals)**2 ), 0)

    # Make a "whole" pulse from half...
    if is_even:
        yw = np.hstack([y[::-1], y]) 
    else:
        z = y.copy() 
        z = np.delete(z,0)
        yw = np.hstack([z[::-1], y])
    
    yf = filter_pulse(yw, xvals, filter_type, filter_percent)    
        
    pulse_type = 1  # pf_constants.AnalyticType.GAUSSIAN
    ys = scale_pulse(yf, dwell_time, angle, pulse_type, gamma)  
    
    # Convert real array to complex array.
    rf_waveform = ys + 1j * np.zeros(len(ys))
        
    rf_xaxis = np.arange(time_steps) * dwell    # in sec

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


def scale_pulse(b1, dwell_time, angle, pulse_type, gamma):  

    # Scale b1 to requested angle.
    multiplier = 1000*(angle/360.0)
    if pulse_type == 3:  #pf_constants.AnalyticType.HYPERBOLIC_SECANT:
        sfact = np.sum(np.abs(b1))
    else:
        sfact = np.sum(b1)
    nfact = multiplier/(gamma * dwell_time * sfact) 
    return nfact*b1 
    

def filter_pulse(b1, xvals, filter_type, filter_percent):
    
    length = len(b1)
    if length%2 == 0:
        even_flag = 1  
    else:
        even_flag = 0     
    
    # seq to ~ pi/2 **    
    maxx = np.amax(xvals)
    f = (xvals*math.pi/maxx/2.0)*(filter_percent/100.0)

    # Cosine apodization
    if filter_type == 1:  #pf_constants.ApodizationFilterType.COSINE:   
        g = np.cos(f) 
    # Hamming apodization
    elif filter_type == 2:  #pf_constants.ApodizationFilterType.HAMMING:              
        g = 0.5*(1 - np.cos(math.pi+2*f)) 
    else:
        g = np.ones(len(f))

    if even_flag == 1:
        g = np.hstack([g[::-1], g]) 
    else:     
        h = g.copy() 
        h = np.delete(h,0)
        g = np.hstack([h[::-1], g]) 

    z = g*b1
    
    return z 
    
    
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
    params['time_steps']         = 250       
    params['tip_angle']          = 90.0
    params["duration"]           = 8.0
    params["bandwidth"]          = 4.0
    params["filter_type"]        = 0
    params["filter_application"] = 0.0
    
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

    bloch_inputs['bloch_range_value'] = 3.0    # unit set below
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
    axes2.set_ylabel('Magnetization |Mxy|')
    axes2.plot(dat1x, np.abs(dat1y), dat1x, np.real(dat1y), dat1x, np.imag(dat1y) )

    axes3 = fig.add_subplot(313)
    axes3.set_xlabel('Position [cm]')
    axes3.set_ylabel('Magnetization Mz')
    axes3.plot(dat2x, np.real(dat2y))
    
    pyl.show()
    
    bob = 10
    
    
