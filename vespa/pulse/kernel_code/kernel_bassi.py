# Python Modules.

import copy
import math

# 3rd party stuff
import numpy as np

# Local imports
import vespa.common.rfp_rf_result as rfp_rf_result
from vespa.common.transform_run_exception import TransformRunException

from pylab import *


PI = np.pi



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

    Derived from Warnking paper MRM 52:1190-1199 (2004)

    Note. this code was rewritten to be all in one function (ie. no separate
    aBt or bBt functs so as to make port to C/C++ easier

    vref [uT] = reference voltage for scaling, default for Siemens 1ms hard 180 pulse = 11.7uT

    Pulse Shape Parameters and Variables and Their Units and Definitions
    ----------------------------------------------------------------------------------------------
    Variable        Units       Description                             Definition
    ----------------------------------------------------------------------------------------------
    a, a0, aB       rad         Amplitude parameter                     Eqs. [1] and [2]
    alpha           deg         Effective flip angle                    alpha = acos(1-2Pe)
    AM              rad/s       RF amplitude modulation function        Eq. [1]
    b, b0, bB       rad         Bandwidth parameter                     Eqs. [1] and [2], b = pi*mu
    beta            rad         Truncation factor of driving function   g(+/-Tp/2)=+/-2*beta
    c               1           Normalized off-resonance                Eqs. [3] and [5]
    f0              1           Maximal relative bandwidth increase     Eq. [15]
    FM              rad/s       RF frequency modulation function        Eq. [2]
    g               rad         HS driving function                     Eqs. [1] and [2]
    gdot            rad/s       Derivative of the driving function
    gamma           rad/s*T     Gyromagnetic ratio                      2.675222exp8 rad/s*T
    G               T/m         Slice-select gradient                   Eq. [18]
    kappa           1           BASSI amplitude scaling parameter       Eqs. [14] and [15]
    L               1           relative RF energy, as function of Pc   Eq. 9
    mu              1           Bandwidth parameter, notation from (14) mu = beta/pi
    Pe              1           Population inversion                    Pe=(1-cos(alpha))/2
    Tp              s           Pulse duration
    x0              m           Center of inversion / saturation band
    deltax          m           Width of inversion / saturation slab

    The population Pe inversion corresponds to an effective flip angle (i.e., the 
    flip angle of a plane rotation that would produce the same longitudinal 
    magnetization), of alpha = arccos(1-2Pe).


    Pulse parameters, as well as peak B1, peak gradient amplitude (Gmax) and width of the frequency 
    sweep (BW), for the pulses compared in the simulations and phantom experiments
    -----------------------------------------------------------------------------------------------
    Pulse               Tp      alpha       beta    b0      f0      Peak B1     Gmax       BW
                       (ms)     (deg)       (rad)                    (uT)       mT/m)     (kHz)
    -----------------------------------------------------------------------------------------------
    HS                  8.1     90          5.3     157     1       23.0        4.9         20.8
    FOCI                8.1     90          5.3     157     3.2     23.0        15.7        66.7
    BASSI(kappa=2)      8.1     90          5.3     168     3       23.1        15.7        66.7
    VERSE-HS            8.1     90          5.3     503             22.8        19.0        81.0

    HS                  10.24   170         5.3     37.7    1       22.9        0.9         4.0
    FOCI                10.24   175         6.3     22.8    22      23.0        14.7        62.7
    BASSI(kappa=2)      10.24   175         6.3     22.8    22      23.0        14.7        62.7
    BASSI(kappa=1.6)    10.24   175         6.3     22.8    22      23.0        14.7        62.7
    BASSI(GOIA)         10.24   175         6.3     22.8    22      23.0        14.7        62.7
    VERSE-HS            10.24   175         5.3     251             23.0        23.6        100.5

    """
    param = trans_desc.parameters
    extra = trans_desc.extra

    #--------------------------------------------------------------------------
    # unpack and convert parameters 
    #
    # - Add or Modify the values below as needed

    time_steps      = int(param["time_steps"])          # int
    duration        = float(param["duration"])          # float, msec
    flip_angle      = float(param['tip_angle'])
    trunc_factor    = float(param['trunc_factor'])
    bw_factor       = float(param['bw_factor'])
    max_relative_bw = float(param['max_relative_bw'])
    slab_width      = float(param['slab_width'])        # float, in meters
    slab_center     = float(param['slab_center'])       # float, in meters
    ampl_scale      = float(param['ampl_scale'])        # float, 2 or 1.6
    
    dwell_time      = (1000 * duration) / (time_steps)     # in usec
    
    # these extra items may be used in making waveform or profile
    bandwidth_convention = int(extra['pulse_bandwidth_type'])  # Choice, 0-HalfHeight, 1-Min, 2-Max    
    resolution           = int(extra['calc_resolution']) # int
    gamma                = float(extra['gamma'])         # float MHz/T
    
    # convert gamma to required units -> 1e8*rad/s*T
    gamma = 2.0 * PI * gamma / 100.0

    # note on gamma units MHz/T -> 2pi -> rad*MHz/T = rad*1e6*Hz/T = rad*1e6/(sec*T)
    #                                               = rad*1e2*Hz/gauss 
    #  - so, MHz/T -> 2pi * 100 -> rad*Hz/gauss 
    #   -> 42.576 MHz/T -> 26751.3 Hz/gauss
    #   -> 42.576 T/(rad MHz) -> 2.67513 1e8 rad/(sec T)

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

    npts   = time_steps         #1024       # 8100 with dwell 0.000001
    dwell  = dwell_time * 1e-6  #0.000002
    Tp     = npts * dwell
    alpha  = flip_angle         #90.0
    beta   = trunc_factor       #5.3
    b0     = bw_factor          #15.0  #168.0
    f0     = max_relative_bw    #3.0
    deltax = slab_width         #0.04
    x0     = slab_center        #0.0
    kappa  = ampl_scale         #2.0

    # Output Units: Tesla, Hz, T/m, deg
    amt,  fmt,  g,  phit  = pulse_bassi_warnking(npts, dwell, 
                                                 alpha, beta, 
                                                 b0, f0, 
                                                 kappa=2.0, 
                                                 deltax=deltax, 
                                                 x0=x0,
                                                 gamma=gamma*1e8)

    amt = amt * 1000.0      # convert T to mT
    g   = g   * 1000.0      # convert T/m to mT/m
    
    rf_waveform = amt * np.exp(1j*phit)    
    rf_xaxis    = np.arange(time_steps) * dwell
    
    gradient    = g 
    grad_xaxis  = np.arange(time_steps) * dwell
    
    # end transform algorithm
    #--------------------------------------------------------------------------

    #--------------------------------------------------------------------------
    # fill an output dictionary, or leave as None

    outputs = {}

    bw = (g[0] * 0.1) * 4358.0 * (deltax*100)    # (mT/m->G/cm)  Hz/G  cm 
    outputs['bw'] = bw
    
#    wid = bw * 10 / (g[0] * 4358)       # this is in cm
#    outputs['wid'] = wid

    #--------------------------------------------------------------------------
    # ensure you are returning ndarrays or None

    rf_y = np.array(rf_waveform)
    rf_x = np.array(rf_xaxis)
    gr_y = gradient
    gr_x = grad_xaxis

    return rf_y, rf_x, gr_y, gr_x, outputs
    
    

def pulse_bassi_warnking(npts, dwell, alpha, beta, b0, f0, 
                         vref=11.7, b1ref=11.7, kappa=2.0, 
                         x0=0.0, deltax=0.02, gamma=267513000.0):
    r"""

    See comments above for description of this pulse algorithm
    
    gamma - defaults to 1H value
    
    """
    # gamma = 267522200.       # 2.675222 x 10^8 rad/s*T

    Tp = npts * dwell
    t  = np.arange(npts) * dwell
    t  = t - (Tp/2.0)

    Pe = 0.5*(1-np.cos(alpha*PI/180.0))      # from comment above
    
    # Equation 9
    L = (1.0/(2*PI))*np.log(1.0/(np.sqrt(1-Pe)))

    # Equation 15,  bB(t) = bBt
    bBt = np.power(np.power((np.power(np.cosh(2*beta*t/Tp),kappa)*(b0-(PI*L)) + (PI*L)),-2) + np.power(f0*b0,-2), -0.5)
    bB0 = np.power(np.power((np.power(np.cosh(2*beta*0/Tp),kappa)*(b0-(PI*L)) + (PI*L)),-2) + np.power(f0*b0,-2), -0.5)
    
    # Equation 14
    aBt = np.power((bBt-PI*L)/(b0-PI*L),1.0/kappa) * np.sqrt(b0*b0 - 4*np.power(np.arccosh(np.cosh(b0/2)*np.sqrt(1-Pe)),2))
    aB0 = np.power((bB0-PI*L)/(b0-PI*L),1.0/kappa) * np.sqrt(b0*b0 - 4*np.power(np.arccosh(np.cosh(b0/2)*np.sqrt(1-Pe)),2))
    
    # Equation 16
    #
    # Acc'd to paper, if we plug in vref and b1ref output is in Volts
    # For now, we set vref=b1ref thus removing that term from the equation
    # and it returns in T (I think)
    
    amt = (vref/b1ref)*(2*aB0*beta/(PI*gamma*Tp)) * (aBt.copy()/aB0)*(1.0/np.cosh(2*beta*t/Tp)) # for sech()
        
    # Equation 18
    #
    # Acc'd to paper this returns in Tesla/meter
    
    g = bBt.copy()*4*beta /(PI*gamma*deltax*Tp)
    
    # Equation 17
    # 
    # first we create the function, then integrate under it for each point
    fmt  = ((2*bBt.copy()*beta)/(PI*Tp))*((-2*x0/deltax) + np.tanh(2*beta*t/Tp))
    phit = np.cumsum(fmt*dwell)
    phit = phit % (PI*2)
    
    return amt, fmt, g, phit



def b_acosh(x):
    return np.log(x + np.sqrt(np.power(x,2) - 1))


_PI     = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348;
_TWO_PI = 6.2831853071795864769252867665590057683943387987502116419498891846156328125724179972560696;

# Floating-point modulo
# The result (the remainder) has same sign as the divisor.
# Similar to matlab's mod(); Not similar to fmod() -   Mod(-3,4)= 1   fmod(-3,4)= -3
def Mod(x,y):

    if 0.0 == y:
        return x

    m = x - y * np.floor(x/y)

    # handle boundary cases resulted from floating-point cut off:
    if y > 0:                   # modulo range: [0..y)
        if m >= y:              # Mod(-1e-16             , 360.    ): m= 360.
            return 0
        if m <0:
            if (y+m) == y:
                return 0        # just in case...
            else:
                return y+m      # Mod(106.81415022205296 , _TWO_PI ): m= -1.421e-14 
    else:                       # modulo range: (y..0]
        if m <= y:              # Mod(1e-16              , -360.   ): m= -360.
            return 0
        if m > 0:
            if (y+m) == y:
                return 0    # just in case...
            else:
                return y+m  # Mod(-106.81415022205296, -_TWO_PI): m= 1.421e-14 
    return m


# wrap [rad] angle to [-PI..PI)
def WrapPosNegPI(fAng):
    return Mod(fAng + _PI, _TWO_PI) - _PI

# wrap [rad] angle to [0..TWO_PI)
def WrapTwoPI(fAng):
    return Mod(fAng, _TWO_PI)

# wrap [deg] angle to [-180..180)
def WrapPosNeg180(fAng):
    return Mod(fAng + 180.0, 360.0) - 180.0

# wrap [deg] angle to [0..360)
def Wrap360(fAng):
    return Mod(fAng ,360.0)


    
    
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

    params["time_steps"]         = 1024      # int, 
    params["duration"]           = 4.096     # float, msec
    params['tip_angle']          = 90.0      # float, alpha
    params['trunc_factor']       = 5.3       # float, beta
    params['bw_factor']          = 15.0      # float, b0
    params['max_relative_bw']    = 3.0       # float, f0
    params['slab_width']         = 0.04      # float, deltax, in meters
    params['slab_center']        = 0.0       # float, x0, in meters
    params['ampl_scale']         = 2.0       # float, kappa, 2 or 1.6

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
        gry = np.array(gry)             # this is in mT/m
    if grx is not None:
        grx = np.array(grx)
    result.gradient    = gry
    result.grad_xaxis  = grx            # this is in sec

    # Call Bloch simulation on rf waveform
    # - need this dict to use Vespa bloch lib calls

    bloch_inputs = {}
    bloch_inputs['calc_resolution'] = extra['calc_resolution'] 
    bloch_inputs['gamma']             = extra['gamma'] 
    bloch_inputs['bloch_range_value'] = 12.0    # unit set below
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
    
    
