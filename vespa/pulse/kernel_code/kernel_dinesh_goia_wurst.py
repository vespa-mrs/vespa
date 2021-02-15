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

from pylab import *

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

    pulse_type      - 'FOCI', 'GOIA-HS' or 'GOIA-WURST'
    time_steps      - Number of points
    duration        - float, [sec]
    bandwidth       - float, [kHz]
    hsn_modulation  - float, 
    grad_modulation - float, 
    grad_factor     - float, [0.0-1.0] gradient modulation factor. The higher 
                             the value the more G "drops" in the middle of the 
                             pulse and the less B1max is required. 0.9 is typical
    total_rotation  - float, [kHz]
    offset_hz       - float, [Hz] offset from central frequency
    slice_thick     - float, [cm] 

    Based on Matlab code from Dinesh Deelchand, CMRR, 14 September 2015
    % Updated 29 Sept 2015

    simulB1Goia2(pulseType,Tpms,HSn,RkHz,Gradm,fMod,B1vals,myOffset)
    %*************************************************
    % FOCI, GOIA-HS or GOIA-WURST gradient-modulated RF pulses
    % - RF and gradient waveforms are generated internally
    % - based on Andronesi et al. JMR 2010 paper 

    """
    param = trans_desc.parameters
    extra = trans_desc.extra

    #--------------------------------------------------------------------------
    # unpack and convert parameters if needed

    pulse_type          = params['pulse_type']              # string
    time_steps          = int(param["time_steps"])          # int
    duration            = float(param["duration"])          # float, [sec]
    bw                  = float(param["bandwidth"])         # float
    hsn                 = float(param["hsn_modulation"])    # float
    gmod                = float(param["grad_modulation"])   # float
    gfact               = float(param["grad_factor"])       # float
    b1rot               = float(param["total_rotation"])    # float
#    offset              = float(param["offset_hz"])         # float, %
    slice_thick         = float(param["slice_thick"])       # float, [cm]
    
    rval                = (duration * 1e3) * bw
    dwell_time          = (duration * 1e6) / time_steps     # in usec
    # gamma               = 4257.0                            # Hz/G
    
    # these extra items are used in making profile
    bandwidth_convention = int(extra['pulse_bandwidth_convention'])  # Choice, 0-HalfHeight, 1-Min, 2-Max    
    resolution           = int(extra['calc_resolution']) # int
    gamma                = float(extra['gamma'])         # float MHz/T -> 2pi * 100 -> rad*Hz/gauss
    gamma = 100.0 * gamma # correct units -> Hz/G

    #--------------------------------------------------------------------------
    # simulate GOIA/FOCI RF and grad waveforms

    if pulse_type =='FOCI' and (hsn>1):
        hsn = 1

    beta = np.arccosh(1.0/0.01)   # arcsech(x) = arccosh(1/x)   # 10% cut off
    t    = np.linspace(0, duration, time_steps) 
    tau  = 2*t/duration - 1 

    #***************************************************
    # Modulation parameters and Equations
    #***************************************************
    # Amplitude
    
    b1max = 1

    aa = np.zeros(len(tau))

    if pulse_type == 'GOIA-HS' or pulse_type =='FOCI':

        if pulse_type == 'FOCI':
            # shaping function, i.e. C-shaped, Kinchesh et al. JMR 2005
            thresh = 10     # 100-f*100
            temp2 = 1.0 / np.cosh(beta * tau)   # sech = 1/cosh
            for i in range(len(temp2)): 
                if temp2[i] > (1/thresh):
                    aa[i] = np.cosh(beta*tau[i])
                else:
                    aa[i] = thresh

            # for gradient
            aag = aa / thresh

            # reset params
            gfact = 0   # need to set OFF gradient modulation factor
            hsn   = 1   # order of B1 modulation, e.g. hsn

        else:
            aa  = 1
            aag = 1

        # Amplitude waveform
        amp = aa * b1max / np.cosh(beta * tau**hsn)          # sech = 1/cosh

        # Gradient waveform
        grad = aag * (1 - gfact / np.cosh(beta * tau**gmod))   # sech = 1/cosh

        # Frequency modulation
        feq = lambda x: (b1max / np.cosh(beta * x**hsn))**2 / ((1 - gfact / np.cosh(beta * x**gmod)))


    elif pulse_type == 'GOIA-WURST':
        
        # Amplitude waveform
        amp = b1max * (1 - np.abs(np.sin(np.pi/2*tau))**hsn)

        # Gradient waveform
        grad = ((1 - gfact) + gfact * np.abs(np.sin(np.pi/2*tau))**gmod)

        # Frequency modulation
        feq = lambda x: (b1max * (1 - np.abs(np.sin(np.pi/2*x))**hsn))**2 / (((1 - gfact) + gfact * np.abs(np.sin(np.pi/2*x))**gmod))

    else:
        print('Type of RF pulse to generate NOT defined')
        return


    #***************************************************
    # Frequency modulation
    # - cumulative integration built-in Matlab script
    # - wc is center freq where w(Tp/2)=0, i.e. freq at 
    #   middle of pulse is 0
    #***************************************************
    
    ctr = -1+len(tau)/2
    tmp = integrate.cumtrapz(feq(tau), x=tau, initial=0);
    tmp = tmp - tmp[ctr]

    w1 = grad * tmp
    A  = rval / (2*duration)
    w  = w1 / np.max(w1) * A


    #***************************************************
    # Phase modulation
    # - using cumulative trapz integration
    # - shift phase to have zero value at bottom
    #***************************************************
    
    phas = 2 * np.pi * integrate.cumtrapz(w, t, initial=0);
    phas = phas - np.min(phas)      # in rad here
    

    #***************************************************
    # Amplitude waveform - normalized
    amp_norm = amp / np.max(amp)

    
    #***************************************************
    # Complex RF waveform 
    re = amp_norm * np.cos(phas)
    im = amp_norm * np.sin(phas)
    b1 = re + 1j*im

    # Complex RF waveform 
    # - normalized so sum(rf)=1
    # - rescaled based on user-defined B1/Hz
    
    mypulse       = b1 / np.sum(np.real(b1))      # normalize sum(b1)=1
    mypulseG      = mypulse / (2 * np.pi * gamma * duration/len(mypulse))
    mypulseHz     = mypulseG * gamma                # convert to Hz
    mypulseHzNorm = mypulseHz / np.max(np.abs(mypulseHz))

    # gradient waveform
    g   = (rval/duration)/(gamma*slice_thick)
    gx  = grad * g * 10

    # RF pulse at total rotation value
    mypulseHzreq = mypulseHzNorm * b1rot

    # convert B1 to Gauss
    mypulseGnew = 0.1 * mypulseHzreq / gamma

    rf_waveform = mypulseGnew
    rf_xaxis    = np.arange(time_steps) * 0.00001    # in sec
    
    gradient    = gx
    grad_xaxis  = np.arange(time_steps) * 0.00001    # in sec
    
    #--------------------------------------------------------------------------
    # fill and output dictionary, if needed

    outputs = {}

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
    import vespa.common.pulse_funcs.bloch_lib_hargreaves as bloch

    pulse_type          = 'GOIA-HS'
    time_steps          = 400
    duration            = 0.004         # [sec]    prev Tp = duration x 1e3 in msec
    bandwidth           = 10.0
    hsn_modulation      = 8.0
    grad_modulation     = 4.0
    grad_factor         = 0.9
    total_rotation      = 575.0         # ~ 15uT I think
    slice_thick         = 2.0           # slice thickness in cm
    
    dp                  = np.linspace(-3,3,time_steps)
    t                   = np.linspace(0, duration, time_steps) 
    offset              = 0.0
    
    extra = {}
    extra['calc_resolution'] = 2000
    extra['pulse_bandwidth_convention'] = 0
    extra['gamma']                = 42.576  # 1H gamma
    
    params = {}
    params['pulse_type']            = 'GOIA-HS'   # GOIA-WURST, FOCI       
    params['time_steps']            = time_steps       
    params['duration']              = duration
    params['bandwidth']             = bandwidth
    params['hsn_modulation']        = hsn_modulation
    params['grad_modulation']       = grad_modulation
    params['grad_factor']           = grad_factor
    params['total_rotation']        = total_rotation
    params['slice_thick']           = slice_thick
    
    
    trans_desc = transform_description.TransformDescription()
    trans_desc.extra                = extra
    trans_desc.parameters           = params

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
    
    
