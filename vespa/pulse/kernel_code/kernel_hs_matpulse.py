# Python Modules.

import copy
import math

# 3rd party stuff
import numpy as np
import scipy.integrate as integrate
import scipy.stats.distributions as distributions

# Local imports
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
    -----------------------------

    pulse_type - HYPERBOLIC_SECANT
    time_steps - Number of points
    quality_cycles - Number of cycles before truncating pulse.
    filter_type - ApodizationFilterType: 0-None, 1-Cosine, 2-Hamming
    filter_percent - How much to apply the filter (0.0 to 100.0)
    total_rotation - Total rotation angle in degrees (used for hyperbolic-secant).
    power_n - Power of hyperbolic secant function.
    sharpness_mu - Parameter that defines sharpness of function.
    dwell_time - The dwell time in microseconds.

    """
    param = trans_desc.parameters
    extra = trans_desc.extra

    #--------------------------------------------------------------------------
    # unpack and convert parameters if needed

    time_steps          = int(param["time_steps"])          # int
    duration            = float(param["duration"])          # float, msec
    total_rotation      = float(param["total_rotation"])    # float, degrees
    quality_cycles      = float(param["quality_cycles"])    # float
    power_n             = int(param["power_n"])             # int
    sharpness_mu        = float(param["sharpness_mu"])      # float
    filter_application  = float(param["filter_application"]) # float, %
    # Choice, 0-None, 1-Cosine, 2-Hamming
    filter_type         = int(param["filter_type"])         # int
    
    dwell_time          = (1000 * duration) / (time_steps)  # in usec
    dwell               = dwell_time / 1000000.0            # in sec
    
    # these extra items are used in making profile
    bandwidth_convention = int(extra['pulse_bandwidth_convention'])  # Choice, 0-HalfHeight, 1-Min, 2-Max    
    resolution           = int(extra['calc_resolution']) # int
    gamma                = float(extra['gamma'])         # float MHz/T -> 2pi * 100 -> Hz/gauss


    #--------------------------------------------------------------------------
    # pulse creation code starts here

#     rf_y = analytic_pulse.hyperbolic_secant(time_points, cycles,
#                                 filter_type, filter_application, 
#                                 total_rotation, n, mu, dwell_time)
#
#     def hyperbolic_secant(nop, cycles, filter_type, filter_percent, angle, n, mu, dwell_time):
#         return analytic_pulse(pulse_type,nop,cycles,filter_type,filter_percent,angle,n,mu,dwell_time)
#    
#     def analytic_pulse(pulse_type, nop, cycles, filter_type, filter_percent, \
#                                        angle, n=1, mu=1, dwell_time=100):
    
    pulse_type = 3  #pf_constants.AnalyticType.HYPERBOLIC_SECANT:
    nop = time_steps
    cycles = quality_cycles
    angle = total_rotation
    n = power_n
    mu = sharpness_mu
    filter_percent = filter_application
    
    # Set local values
    
    g = np.ones(nop)
    y = g.copy()
    z = g.copy()
    
    # Determine whether nop even or odd (even = 1 or 0)
    if nop%2 == 0:
        even_flag = 1  
    else:
        even_flag = 0 
    
    # Generate symmetric sequence
    # starts at 0 for odd, at +1/2 for even
    if even_flag == 1:          
        # even_flag, starts at 1/2
        # matlab code; x = (1:nop/2) - 1/2
        x = np.arange(nop/2) + 1/2
    else:
        # odd, starts at 0
        # matlab code; x = (1:(nop+1)/2) -1
        x = np.arange((nop+1)/2)

    # HS pulse
    # sech(beta*t), where tmax = 1; so beta = cycles/2
    if n == 1:
        y = pow( math.pi*distributions.hypsecant.pdf(x*cycles/nop), complex(1, mu) ) 
    else:
        sqrtn = math.sqrt(n)
        s = (cycles/(sqrtn*nop))*x
        y = math.pi*distributions.hypsecant.pdf(s**n) 
    
        wfunct = lambda q: math.pi*distributions.hypsecant.pdf(q**n)**2
        
        w = np.zeros(len(s))
           
        for k in range(len(s)): 
            integral = integrate.quad(wfunct, 0, s[k])
            w[k] = integral[0]
            w[k] *= -mu*2.0*math.pi
             
        dphi = s[1]*w
        phi = dphi.cumsum(axis=0)
        
        y = y * np.exp(complex(0,1)*phi)
    
    # Note that beta, when n = 1, is cycles/2;
    # Gives width = mu*beta/pi = cycles*mu/pi
    
    if even_flag == 1:
        y = np.hstack([y[::-1], y]) 
    else:
        z = y.copy() 
        z = np.delete(z,0)
        y = np.hstack([z[::-1], y])
        
    # Check for apodization (Filter funct = g)
    # Then use apodization
    if not filter_type:
        z = y        
    else:
        z = filter_pulse(y, x, filter_type, filter_percent)

    z = scale_pulse(z, dwell_time, angle, pulse_type, gamma)    

    rf_waveform = z
    rf_xaxis = np.arange(time_steps) * dwell    # in sec
    
    gradient = None
    grad_xaxis = None
    
    #--------------------------------------------------------------------------
    # fill and output dictionary, if needed

    bandwidth = dwell_bandwidth_interconvert(dwell_time, sharpness_mu, quality_cycles, time_steps)
    outputs = {}
    outputs['bandwidth'] = bandwidth

    #--------------------------------------------------------------------------
    # ensure you are returning ndarrays or None

    rf_y = np.array(rf_waveform)
    rf_x = np.array(rf_xaxis)
    gr_y = gradient
    gr_x = grad_xaxis

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


def dwell_bandwidth_interconvert(dtb, sharpness_mu, cycles, time_points):        
    # Note: dwell_time in microseconds and bandwidth in kHz.
    multiplier = sharpness_mu/math.pi
    return multiplier * cycles * 1000 / (time_points * dtb)    
    
    
    
    
    
#------------------------------------------------------------

if __name__ == '__main__':
    
    import vespa.public.transform_description as transform_description
    import vespa.common.constants as constants
    import pylab as pyl
    
    extra = {}
    extra['calc_resolution'] = 2000
    extra['pulse_bandwidth_convention'] = 0
    extra['gamma']                = 42.576  # 1H gamma
    
    params = {}
    params['time_steps']            = 500       
    params['duration']              = 4.5
    params['total_rotation']        = 1240.0
    params['quality_cycles']        = 6.0
    params['power_n']               = 6
    params['sharpness_mu']          = 16.0
    params['filter_application']    = 0.0
    params['filter_type']           = 0
    
    trans_desc = transform_description.TransformDescription()
    trans_desc.extra                = extra
    trans_desc.parameters           = params

    rfy, rfx, gry, grx, outs = run(trans_desc)
    
    print("Bandwidth = "+str(outs['bandwidth']))
    
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
    bloch_inputs['bloch_range_value'] = 5.0     # unit set below
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

    axes1 = fig.add_subplot(311)
    axes1.set_xlabel('Time [ms]')
    axes1.set_ylabel('[mTesla]')
    axes1.plot(rfx*1000, np.real(rfy), rfx*1000, np.imag(rfy))
    
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
    
    
