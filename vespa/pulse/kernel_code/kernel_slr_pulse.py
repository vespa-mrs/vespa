# Python Modules.

import copy
import math

# 3rd party stuff
import numpy as np
import scipy as sp
import scipy.signal as sps

# Local imports
import vespa.common.rfp_rf_result as rfp_rf_result
from vespa.common.transform_run_exception import TransformRunException


# Gyromagnetic ratio of 1H - the hydrogen nucleus. (units: kHz/mT)
GAMMA1H = 42.576    

# Small number used for floating point comparisons.
epsilon = 0.000001  

# Very small number used for double precision comparisons.
small_epsilon = 0.00000000001  

# An even smaller number; used to represent an acceptable fractional error
# or difference from an expected or desired value.
EPS = pow(2,-52)


def run(trans_desc):
    """
    Adapted from Matpulse - mplcalcc.m: 
    
    Script for pulse generation. Adjusts some input parameters to be in the 
    correct units.

    time_steps      - int, the number of data points in the output waveform
    duration        - float, total time of pulse, in msec
    tip_angle       - float, the desired tip angle, in radians
    pass_ripple     - float, the desired pass band ripple, as a percent
    reject_ripple   - float, the desired reject band ripple, as a percent
    phase_type      - int, the non-coalesced phase type:
                        1 for linear, 2 for Max, 3 for Min
    bandwidth       - float, the bandwidth for the pulse in kilohertz
    is_single_band  - bool, True for single band; False for 
                      dual band.
    separation      - float, the separation between the bands in kilohertz
    use_remez       - bool, True remez filter, False least-squares filter
    n_zero_pad      - int, number of zero pad pts at the begin/end of the pulse
    
    resolution      - int, the calculation resolution for frequency domain 
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
    See the following articles:

    J. Pauly, P. Le Roux, D. Nishimura, A. Macovski, 
    'Parmater Relations for the Shinnar-Le Roux Selective Excitation Pulse Design Algorithm',
    IEEE Transactions on Medical Imaging, Vol. 10, No. 1, pp 53-65, March 1991.

    G.B. Matson, 'An Integrated Program for Amplitude-Modulated RF Pulse Generation  
    andRe-mapping with Shaped Gradients', Magnetic Resonace Imaging, 
    Vol. 12, No. 8, pp 1205-1225, 1994'''

    """
    
    param = trans_desc.parameters
    extra = trans_desc.extra

    #--------------------------------------------------------------------------
    # parse parameters and convert if needed

    time_steps      = int(param["time_steps"])      # int
    duration        = float(param["duration"])      # float, msec
    tip_angle       = float(param["tip_angle"])     # float, deg 
    pass_ripple     = float(param["pass_ripple"])   # float, deg
    reject_ripple   = float(param["reject_ripple"]) # float, deg
    bandwidth       = float(param["bandwidth"])     # float, deg
    is_dual_band    = int(param["is_dual_band"])  # bool
    use_lsq         = int(param["use_lsq"])       # bool
    separation      = float(param["separation"])    # float, Hz

    phase_type = int(param["phase_type"])+1 # Choice, 1-Linear (refocused), 2-MinPhase, 3-MaxPhase 
    
    bandwidth_convention = int(extra['pulse_bandwidth_type']) # Choice, 0-HalfHeight, 1-Min, 2-Max     
    resolution           = int(extra['calc_resolution']) # int
    gamma                = float(extra['gamma'])         # float MHz/T -> 2pi * 100 -> Hz/gauss
    n_zero_pad           = int(extra["zero_padding"])    # int

    is_single_band = 1 - is_dual_band     # need this so default can be 0 choice
    use_remez      = 1 - use_lsq        # need this so default can be 0 choice

    # Pulse angle category (90, 180, shallow tip)
    # - Tip angles are currently limited as shown below. 
    # - Kim's paper has shown that we can use tip angles outside of this range, 
    #   but we haven't written the code yet.
    if not (tip_angle == 90) or (tip_angle == 180) or ((tip_angle > 0) and (tip_angle <=30)):
        # FIXME bjs Implement handling of non-standard tip angles, a la Kim's paper
        error_msg =  "The tip angle must be 90, 180 or <= 30."
        raise TransformRunException(error_msg, 1)

    tip_angle = tip_angle * np.pi/180.0     # convert deg to radians 

    dwell_time = (1000 * duration) / (time_steps)     # in usec
    
    if is_single_band: separation = 0.0

    # Check for excessive bandwidth
    if is_single_band == True:
        mplttbw = bandwidth
    else:        
        mplttbw = separation/2 + bandwidth

    if (0.8 * mplttbw) > (1000.0 / dwell_time) :
        error_msg = 'Error: Excessive pulse bandwidth'
        raise TransformRunException(error_msg, 1)



# DID CODE INSPECTION - THIS CALL DOES NOT SEEM TO AFFECT
# THE OUTCOME OF slr_pulse.py
#
#    mpltran = 0
#    mpltran = TB.transition_band(time_steps, resolution, dwell_time, duration, tip_angle, 
#                                 pass_ripple, reject_ripple, phase_type, bandwidth, is_single_band, 
#                                 separation, n_zero_pad, bandwidth_convention)


    # Pulse Setup, starting with initial length (for z-padding)
    if n_zero_pad > 0:
        # Reduced for z-padding
        nofp = time_steps - 2 * n_zero_pad
    else:
        nofp = time_steps

    if phase_type == 1:
        # Refocused or Linear pulse

        mplfr, mplmag, mplwts, mplcrrctf = linear_setup(nofp, dwell_time, tip_angle, pass_ripple, \
                                                        reject_ripple, bandwidth, is_single_band, \
                                                        separation, bandwidth_convention)

    else: 
        # Max or min phase

        mplfr, mplmag, mplwts, mplcrrctf = max_min_setup(nofp, dwell_time, tip_angle, pass_ripple, \
                                                         reject_ripple, bandwidth, is_single_band, \
                                                         separation, bandwidth_convention) 

    # ... test for monotonicity and for values between zero and one.
    #
    # Specific Bandwidth Tests - per Jerry Matson
    error_msg = ''
    if bandwidth_convention == 0: 
        # default bandwidth convention (width at mid height)
        if is_single_band == True:
            if mplfr[1] < 0.0:
                error_msg = 'Bandwidth Error: Increase Bandwidth or Increase Time for Pulse'
            if mplfr[2] > 1.0:
                error_msg = 'Bandwidth Error: Decrease Bandwidth or Increase Time for Pulse'
        else:    # if dual band
            if mplfr[1] < 0.0:
                error_msg = 'Bandwidth Error: Decrease Bandwidth or Increase Time for Pulse'
            if mplfr[4] > 1.0:
                error_msg = 'Bandwidth Error: Decrease Bandwidth or Increase Time for Pulse'
          
    elif bandwidth_convention == 1:  
        # specroscopic bandwidth convention (width at base of pulse)
        # :FIXME: implement this ... and the one just below it.
        raise NotImplementedError('Oops, forgot this case')        
    else:
        # filter bandwidth convention (width at top of pulse)
        raise NotImplementedError('Oops, forgot this case')   

    if error_msg:
        raise TransformRunException(error_msg, 1)

    # General Check
    xold = 0.0
    for x in mplfr:
        if x < xold or x > 1.0:
            error_msg = 'Values in fr need to be monotonic increasing, and between 0 and 1'
            raise TransformRunException(error_msg, 1)
        else:
            xold = x        

    # Calculate the expansion factor
    expansion_factor = 1/(2*mplfr[2])
    if not is_single_band:
        expansion_factor = 1/(2*mplfr[4])       # Dual band pulse

    # Guarantee that expansion_factor >= 1 
    if expansion_factor < 1.0:
        expansion_factor = 1.0


    # Calculate polynomials
    #
    # Note: The reason min_max_filter gets passed the "is_single_band" flag and 
    # the linear_filter does not is that min_max_filter does some pre-calculation
    # adjustments to try to reduce or eliminate the dc offset in the reject band or
    # region. This step is not done in the linear filter case.

    if phase_type == 1:    # Linear pulse
        # mplfr  - pass and stop
        # mplmag - magnetization
        # mplwts - weights        
        
        mpgbb, mpgaa, mplw = linear_filter(mplfr, mplmag, mplwts, nofp,
                                           use_remez, resolution)
    else:
        mpgbb, mpgaa, mplw = max_min_filter(mplfr, mplmag, mplwts, nofp,
                                             use_remez, resolution, is_single_band)

    # Calculatel b1:
    #  
    #  mmgaa and mmgbb - returned from linear_filter or max_min_filter above.
    mpgb1 = inverse_slr(nofp, dwell_time, mpgaa, mpgbb, gamma)

    # mplcrrctf returned from setup routines above.
    if mplcrrctf != 0:
        mpgb1 = mplcrrctf * mpgb1

    # Add zero padding if asked for:
    if n_zero_pad > 0:
        mpgb1 = np.hstack([np.zeros(n_zero_pad), mpgb1, np.zeros(n_zero_pad)])

    # Reverse for max phase pulse
    if phase_type == 2:
        mpgb1 = mpgb1[::-1]  # Max phase pulse 

    waveform = mpgb1

    temparr = []
    pulse_time = duration/1000.0
    for ix in range(time_steps):
        temparr.append(ix*pulse_time/(time_steps))    

    waveform_xaxis = temparr

    rf_y = np.array(waveform)
    rf_x = np.array(waveform_xaxis)
    
    
    rf_waveform = np.array(rf_y)
    rf_xaxis    = np.array(rf_x)
    gradient    = None
    grad_xaxis  = None


    return rf_waveform, rf_xaxis, gradient, grad_xaxis, None
    


def linear_setup(nlength, dwell, tip_angle, pass_ripple, reject_ripple, bandwidth, 
                 is_single_band, separation, bandwidth_convention):

    '''
    Setup script for linear (or refocused) rf pulses. 

    Input values: 
       nlength - number of points.
       dwell - dwell time in microseconds
       tip_angle - in radians        
       pass_ripple - bandpass ripple
       reject_ripple - bandreject ripple
       bandwidth - in kHz        
       is_single_band - band type; single or dual.        
       separation - in kHz
       bandwidth_convention
    
    May raise (or throw) PulseFuncException
    
    For detailed information on algorithms, see:
    K.J. Lee, 'General parameter relations for the Shinnar-Le Roux pulse design algorithm',
    J. Magnetic Resonance, 186, pp 272-258, 2007, 
    and:
    J. Pauly, P. Le Roux, D. Nishimura, A. Macovski, 
    'Parmater Relations for the Shinnar-Le Roux Selective Excitation Pulse Design Algorithm',
    IEEE Transactions on Medical Imaging, Vol. 10, No. 1, pp 53-65, March 1991.
    '''
                
    # Adapted from Matlab: MPLRB1C.M

    if not (pass_ripple > 0 and reject_ripple > 0):
        error_msg = 'Error: passband and reject band ripple must both be greater than zero'
        raise TransformRunException(error_msg, 1)

    # Set desired pass_ripple (bandpass) and reject_ripple (bandreject); 
    sigm1=pass_ripple/100 
    sigm2=reject_ripple/100

    # Initialize corrctf for later correction if needed
    corrctf = 0

    # Set local variable names
    na = nlength
    ang = tip_angle

    if is_single_band == True:
        bw = bandwidth        
    else:           
        # Dual band pulse
        # Get narrowest band.
        bw = separation
        if bandwidth < separation: 
            bw = bandwidth
        
    # pulse length (ms)
    ms = dwell*(na)/1000      
       
    # Time-bandwidth product
    tb = ms*bw
        
    # D(l), wts, and transb calculated using Lee
    sigi = sigm1
    sigo = sigm2
    phi = ang

    is180 = False
    # See if angle is either 180 or 90.
    if abs(tip_angle - np.pi) < epsilon:                 
        # SE pulse
        sigm1 = (1/(2*sigi))*(4*(1-sigi/2) - 4*math.sqrt(1-sigi))
        sigm2 = (1+sigm1)*math.sqrt(sigo) 
        is180 = True
        
    elif abs(tip_angle - np.pi/2.0) < epsilon:             
        # 90 pulse
        n = 0
        f = 0.1
        
        # sigm1 calc
        while abs(f) > 0.00001 and n < 100:
            sigin = abs( (2*math.sqrt(1 - math.pow((1+sigm1),2) * math.pow(math.sin(phi/2),2)) \
                                         *(1+sigm1)*math.sin(phi/2)-math.sin(phi))/math.sin(phi) )
            f = (sigi-sigin)/(2*sigi)
            sigm1 = abs(sigm1*(1 + f))
            n = n+1
        
        # sigm2 calc
        n = 0
        f = 0.1
        while abs(f) > 0.00001 and n < 100:
            sigin = abs(2*math.sqrt(1 - math.pow(sigm2,2) * math.pow(math.sin(phi/2),2))\
                                                           *sigm2*math.sin(phi/2)/math.sin(phi))
            f = (sigo-sigin)/(2*sigo)
            sigm2 = abs(sigm2*(1 + f))
            n = n+1

    # Re-adjust sigmas for dual band pulses
    if not is_single_band:             
        sigm1 = 0.5*sigm1
        sigm2 = 2*sigm2

    # Calculate transition band (only approximate for dual band pulses)
    l1 = math.log10(sigm1)
    l2 = math.log10(sigm2)
    a1 = 5.309e-3
    a2 = 7.114e-2
    a3 = -4.761e-1
    a4 = -2.66e-3
    a5 = -5.941e-1
    a6 = -4.278e-1

    d = (a1 * pow(l1,2) + a2*l1 + a3)*l2 + (a4 * pow(l1,2) + a5*l1 + a6)

    wt1 = 1
    wt2 = sigm1/sigm2
    wid = d/tb
    transb = wid*bw
    nq = na/(2*ms)
    magn = math.sin(ang/2)
        
    # SE pulse        
    if is180:        
        if (magn + sigm1) > 1:    
            # B1 pulse needs correction
            magn = 1-sigm1
            corrctf = magn/(1-sigm1)
        

    # Calculate the bands (fp and fs)
    if is_single_band:
        if bandwidth_convention == 0:             
            # Usual convention (FWHM)
            fp = (bandwidth/2 - transb/2)/nq
            fs = (bandwidth/2 + transb/2)/nq
        elif bandwidth_convention == 1:            
            # Minimum
            fp = (bandwidth/2 - transb)/nq
            fs = (bandwidth/2) / nq
        elif bandwidth_convention == 2:            
            # Maximum
            fp = (bandwidth/2) / nq
            fs = (bandwidth/2 + transb)/nq

        fr = [0, fp, fs, 1]
        mag = [magn, magn, 0, 0]
        wts = [wt1, wt2]
        
    else:         
        
        # Dual band pulse
        if bandwidth_convention == 0:
            fp = (separation/2 - transb/2)/nq
            fs = (separation/2 + transb/2)/nq
        elif bandwidth_convention == 1:
            fp = (separation/2 - transb)/nq
            fs = (separation/2) / nq
        elif bandwidth_convention == 2:
            fp = (separation/2) / nq
            fs = (separation/2 + transb)/nq        
                
        if bandwidth_convention == 0 :
            fss = fs+(bandwidth-transb)/nq
            fpp = fss+(transb)/nq
        elif bandwidth_convention == 1:
            fss = fs+(separation/nq)
            fpp = fss+(transb)/nq
        elif bandwidth_convention == 2:
            fss = fp+(bandwidth-transb)/nq
            fpp = fss+(transb)/nq

        fr = [0, fp, fs, fss, fpp, 1]
        mag = [0, 0, magn, magn, 0, 0]
        wts = [wt2, wt1, wt2]
    
    
    rval = (fr, mag, wts, corrctf)
    
    return rval


def max_min_setup(nlength, dwell, tip_angle, pass_ripple, reject_ripple, bandwidth, 
                  is_single_band, separation, bandwidth_convention):
    '''
    Setup script for RF pulse for max & min phase pulses
    
    '''
    if not (pass_ripple > 0 and reject_ripple > 0):
        error_msg = 'Error: passband and reject band ripple must both be greater than zero'
        raise TransformRunException(error_msg, 1)    

    # Bandpass ripple (%) 
    sigm1 = pass_ripple/100.0 
    
    # Bandreject ripple (%)
    sigm2 = reject_ripple/100.0

    # Initialize corrctf for later correction if needed
    corrctf = 0

    # Set local variable names
    na = nlength  
    ang = tip_angle
        
    if is_single_band == True:
        bw = bandwidth
    else:
        bw = separation
        # Get narrowest band for dual band pulses 
        if bandwidth < separation:
            bw = bandwidth

    # pulse length (ms)
    ms = dwell*(na)/1000.0
        
    # Time-bandwidth product
    tb = ms*bw            
        
    # D(l), wts, and transb calculated

    is180 = False
    if abs(tip_angle - np.pi) < epsilon:            
        # 180 pulse
        sigm1 = sigm1/8.0 
        sigm2 = math.sqrt(sigm2/2.0)
        is180 = True
    elif abs(tip_angle - np.pi/2.0) < epsilon:            
        # 90 pulse
        sigm1 = sigm1/2.0 
        sigm2 = math.sqrt(sigm2)

    # add'l fudges (JM). 
    # Notes, from conversaton with JM on 04/19/11: Apparently this
    # fudge (and the 0.3 "fudge" 21 lines down) is here to reduce,
    # or zero out, the dc offset in the rejection band.
    # Per JM, these "fudges" should be calculate dynamically and 
    # not statically to have the best effect. 
    sigm1 = 2.0*sigm1 
    sigm2 = (sigm2**2.0)/2.0 

    # Re-adjust sigmas for dual band pulses
    if not is_single_band:             
        # Dual band pulse
        sigm1 = 0.5*sigm1 
        sigm2 = 2.0*sigm2        

    # Calculate transition band (only approximate for dual band pulses)
    l1 = math.log10(sigm1) 
    l2 = math.log10(sigm2)
    a1 = 5.309e-3 
    a2 = 7.114e-2 
    a3 = -4.761e-1
    a4 = -2.66e-3 
    a5 = -5.941e-1 
    a6 = -4.278e-1

    # Note: 0.3 is fudge for min/max phase filters
    d = 0.3*( a1*pow(l1,2) + a2*l1 + a3 )*l2 + ( a4*pow(l1,2) + a5*l1 + a6 )

    wt1 = 1.0
    wt2 = sigm1/sigm2
    wid = d/tb 
    transb = float(wid)*bw 
    nq = na/(2.0*ms)
    magn = math.pow(math.sin(ang/2.0), 2.0) 
        
    if is180:         
        # SE pulse
        if (magn + sigm1) > 1.0:    
            # B1 pulse needs correction
            magn = 1.0-sigm1 
            corrctf = magn/(1.0-sigm1)


    if is_single_band:

        # Calculate the bands (fp and fs)
        if bandwidth_convention == 0:                     
            # Usual convention
            fp = (bandwidth/2.0 - transb/2.0)/nq  
            fs = (bandwidth/2.0 + transb/2.0)/nq
        elif bandwidth_convention == 1:                    
            # Minimum
            fp = (bandwidth/2.0 - transb)/nq 
            fs = (bandwidth/2.0)/nq
        elif bandwidth_convention == 2:                    
            # Maximum
            fp = (bandwidth/2.0)/nq  
            fs = (bandwidth/2.0 + transb)/nq

        fr = [0.0, fp, fs, 1.0]
        mag = [magn, magn, 0.0, 0.0]
        wts = [wt1, wt2]
        
    else:
        
        # Calculate the bands (fp and fs)
        if bandwidth_convention == 0:                     
            # Usual conven.
            fp = (separation/2.0 - transb/2.0)/nq  
            fs = (separation/2.0 + transb/2.0)/nq
        elif bandwidth_convention == 1:                    
            # Minimum
            fp = (separation/2.0 - transb)/nq 
            fs = (separation/2.0)/nq
        elif bandwidth_convention == 2:                    
            # Maximum
            fp = (separation/2.0)/nq  
            fs = (separation/2.0 + transb)/nq        
        
        if bandwidth_convention == 0:
            fss = fs + (bandwidth-transb)/nq  
            fpp = fss + (transb)/nq
        elif bandwidth_convention == 1:
            fss = fs + (separation/nq) 
            fpp = fss + (transb)/nq
        elif bandwidth_convention == 2: 
            fss = fp + (bandwidth-transb)/nq  
            fpp = fss + (transb)/nq

        fr = [0.0, fp, fs, fss, fpp, 1.0]
        mag = [0.0, 0.0, magn, magn, 0.0, 0.0]
        wts = [wt2, wt1, wt2]

    rval = (fr, mag, wts, corrctf)
    
    return rval


def linear_filter(fr, mag, wts, na, use_remez, mmgresol):
    """
    Function to create a (linear phase) digital filter
    
    Parameters [fp fs], mag, [wt1, wt2], and na from pulse setup
    Returns h (b-polynomial (evaluated around the unit circle), 
    a (a-polynomial evaluated around the unit circle), and w (freq)
    
    """
    # fr - pass and stop
    # mag - magnetization
    # wts - weights
    
    f = np.array(fr)
    m = np.array(mag)
    wt = np.array(wts)    

    # Check for filter function wanted
    if use_remez:
        # FYI: The matlab documentation for the remez() function
        # is now included as part of "firpm()", the parks Mclellan function.
        #
        # remez() - 
        # Calculate the filter-coefficients for the finite impulse response
        # (FIR) filter whose transfer function minimizes the maximum error
        # between the desired gain and the realized gain in the specified bands
        # using the remez exchange algorithm.
        #
        # b = sps.remez(na-1, f, m, wt)
        # default sampling frequency is 1 Hz.  
        # Put in 2 Hz in order for this to work.        
        #
        b = sps.remez(na, f, m[::2], wt, 2)
    else:        
        # Group f values into pairs.
        fx = f/2.0
        i = 0
        fin = []
        for fi in fx:
            if i%2 == 0:
                fi_0 = fi
            else:
                fin.append( (fi_0, fi) )
            i += 1
                
        # Also, take every other value out of m
        min_ = m[0::2]
               
        # NOTE: RESULTS (without weights) were NOT IN AGREEMENT 
        # WITH MATPULSE.
        #
        # ::FIXME:: Compare results with Matpulse to be sure
        # we are doing the weights correctly.
        # It looks good by visual inspection.
        #
        # b = firls(na-1, fin, m, wt)
        # Changed from na-1 to na, as we did in remez, so that it will 
        # work with an odd number of time steps / time points.
        if na%2 == 0:
            error_msg = 'An odd number of time steps is required for least squares filtering.'
            raise TransformRunException(error_msg, 1)
                
        b = firls(na, fin, min_, wt)         


    whole = 1
    
    (wr, hr) = sps.freqz(b, 1, mmgresol, whole) 
    h = hr.copy()
    w = wr.copy()
    
    # Next calculate the A(z)
    at2 = (1 - ( (np.real(h))**2 + (np.imag(h))**2) ) 

    #x = np.sqrt(at2) 
    # Need to use this version of sqrt because it will
    # work with negative numbers.
    # After fixing the remez calculation above this may no longer
    # be necessary, but probably won't hurt...
    x = np.lib.scimath.sqrt(at2)
    
    # (Filter from RCEPS.M)
    n = x.size
    xhat = np.real(sp.fftpack.ifft(np.log(x)))
    
    # Not sure if epsilon is required, but will guarantee 
    # that we get the expected result: One of n is odd, else zero.
    odd = np.fix( n%2 + epsilon )
    
    aa = (n+odd)/2 - 1
    bb = (n-odd)/2 - 1
    wm = np.vstack( [1, 2*np.ones((aa,1)), 1, np.zeros((bb,1))] )
 
    yhat = np.zeros(np.size(x), dtype=np.complex128)

    yhat[:] = np.real(np.transpose(wm) * xhat)
    
    gg = sp.fft(yhat)
    theta = np.imag(gg)
    
    y = x * np.exp(1j*theta)
    
    a = y 
    
    return (h,a,w)


def max_min_filter(frz, magz, wtz, na, use_remez, mmgresol, is_single_band):
    '''
    Function to create the (min/max phase) digital filter
    
    Parameters [fp fs], mag, [wt1, wt2], and na from pulse setup.
    Returns h (b-polynomial), a (a-polynomial), and w (freq)
    
    '''
    fr  = np.array(frz)
    mag = np.array(magz)
    wts = np.array(wtz)         
    
    # A pre-calculation will be done first to determine the dc offset and then 
    # correct for it in the real, full length calculation, so the reject band
    # is closer to zero. 
    
    # 2011.04.06 
    # Jerry believes he got this idea from Pauly and
    # his collaborators in a private conversation.
    
    nb = 2*na-1
    
    # First determine bias on shorter length pulse   
    frnn = 0 
    nn = na

    # determine pulse length
    while frnn < 0.6:        
        if is_single_band:     
            frnn = fr[2]*nb/nn      # Single band
        else:    
            frnn = fr[3]*nb/nn      # Dual band
        nn = nn/2

    # Suggested initial pulse length
    nn = 4*nn-1        

    # Don't let nn be too small
    nn = int(math.ceil(nn))
    if nn < 63:
        nn = 63

    if is_single_band:

        nbdnn = float(nb)/nn
        frn = np.hstack([ fr[0], fr[1]*nbdnn, fr[2]*nbdnn, fr[3] ])
        
        fs = frn[2] 

        # check that fs less than 0.82
        if fs > 0.82:
            error_msg = 'Error: Try a larger number of steps'
            raise TransformRunException(error_msg, 1)


        # Check for filter function wanted
        if use_remez:        
            # b-poly coefficients
            bp = sps.remez(nn, frn, mag[::2], wts, 2)
        else:            
            # Group f values into pairs.
            fx = frn/2.0
            i = 0
            fin = []
            for fi in fx:
                if i%2 == 0:
                    fi_0 = fi
                else:
                    fin.append( (fi_0, fi) )
                i += 1
                    
            # Also, Take every other value out of m
            min_ = mag[0::2]            
            
            bp = firls(nn, fin, min_, wts)

        # evaluated on unit/2 circle
        wp,hp = sps.freqz(bp, 1, mmgresol)
        
        local_resol = len(wp)
        
        selct = int(np.ceil( 1.2*fs*local_resol ) )

        # partt = [zeros(1,selct) ones(1,mmgresol-selct)]' ;
        # The conjugation method (') does not seem to be needed here since 
        # the arrays are all reals, and the subsequent array multiplication 
        # works without any problem.
        partt = np.hstack( [np.zeros((1, selct)), np.ones((1, local_resol-selct))] )

        bias = 0.93 * np.max(np.ravel((np.abs(hp))*partt))

        # Set bias for calc
        magn = mag + np.hstack([0.0, 0.0, bias, bias])

    else:     
        # Dual band pulse
        nbdnn = float(nb)/nn
        frn = np.hstack( [0.0, fr[1]*nbdnn, fr[2]*nbdnn, fr[3]*nbdnn, 1.1*fr[3]*nbdnn, 1.0] )

        fs = frn[4] 

        # check that fs is less than 0.82
        if fs > 0.82:
            error_msg = 'Error: Try a larger number of steps'
            raise TransformRunException(error_msg, 1)

        # Check for filter function wanted
        if use_remez:
            # b-poly coefficients
            bp = sps.remez(nn, frn, mag[::2], wts, 2)
        else:
            # Group f values into pairs.
            fx = frn/2.0
            i = 0
            fin = []
            for fi in fx:
                if i%2 == 0:
                    fi_0 = fi
                else:
                    fin.append( (fi_0, fi) )
                i += 1
                    
            # Also, Take every other value out of m
            min_ = mag[0::2]                   
            
            #bp = firls(nn-1,frn,mag,wts)
            bp = firls(nn, fin, min_, wts)

        # evaluated on unit/2 circle
        wp,hp = sps.freqz(bp, 1, mmgresol)
        
        local_resol = len(wp)
        
        fs = frn[1]
        
        selct = int( math.ceil(0.8 * fs * local_resol)  )
        
        partt = np.hstack([np.ones( (1,selct) ), np.zeros( (1, local_resol-selct) )])
        
        # Changed from .9 to .02 (7/15/05) (JM)
        bias = 0.02 * np.max(np.ravel(np.abs(hp) * partt))

        # Set bias for calc
        magn = mag + np.hstack([bias, bias, 0.0, 0.0, 4*bias, 4*bias])


    # Next the B(z) (h,w) (evaluated on the unit circle)

    # Check for filter function wanted
    #
    # FIXME: The calls to remez (and firls below) are using
    # nb, and not na. Why? nb = 2*na-1 so if na=16, nb=31    
    if use_remez:
        b = sps.remez(nb, fr, magn[::2], wts, 2)
    else:
        # Group f values into pairs.
        fx = fr/2.0
        i = 0
        fin = []
        for fi in fx:
            if i%2 == 0:
                fi_0 = fi
            else:
                fin.append( (fi_0, fi) )
            i += 1
                
        # Also, Take every other value out of m
        min_ = magn[0::2]    
    
        b = firls(nb, fin, min_, wts)

    whole = 1
    wr, hr = sps.freqz(b, 1, mmgresol, whole)

    w = wr.copy()
    hr = np.lib.scimath.sqrt(np.abs(hr))
    x = np.abs(hr)
    n = x.size
    xhat = np.real(sp.fftpack.ifft(np.log(x))) 
    
    odd = int(np.fix( n % 2 ))    
    
    aa = (n+odd)/2 - 1
    bb = (n-odd)/2 - 1
    wm = np.vstack([1.0, 2*np.ones((aa,1)), 1.0, np.zeros((bb,1))])
    
    yhat = np.zeros(x.size) 
    wmxf = wm.transpose() * xhat    

    yhat[:] = np.real(wmxf)
    theta = np.imag(sp.fft(yhat))
    
    y = x * (np.exp(1j*theta))
    h = y.copy()

    
    # Next the A(z)
    at2 = (1 - (np.real(h)**2 + np.imag(h)**2) )
    x = np.sqrt(at2)

    # (Filter from RCEPS.M)
    n = x.size
    xhat = np.real(sp.fftpack.ifft(np.log(x)))
    odd = int(np.fix(n % 2))
    
    aa = (n+odd)/2 - 1
    bb = (n-odd)/2 - 1
    wm = np.vstack([1, 2*np.ones((aa,1)), 1, np.zeros((bb,1))])
    
    yhat = np.zeros( x.size, dtype=np.complex128)

    yhat[:] = np.real(wm.transpose() * xhat)
    theta = np.imag(sp.fft(yhat))
    
    y = x * np.exp(1j*theta)
    a = y.copy()
    
    return (h,a,w)
    

def firls(N,f,D,wts=None):
    """
    Least-squares FIR filter routine.
    
    N -- filter length, must be odd
    f -- list of tuples of band edges
         Units of band edges are Hz with 0.5 Hz == Nyquist
         and assumed 1 Hz sampling frequency
    D -- list of desired responses, one per band
    wts -- list of weights for each of the bands
    
    May raise (or throw) PulseFuncException
    
    """
    if wts == None:
        wts = [1.0 for d in D]
        
    assert len(D) == len(f), "must have one desired response per band"

    h = np.zeros(N)    
    
    if N % 2:    
        L = (N-1)//2
        k = np.arange(L+1)
        k.shape = (1, L+1)
        j = k.T
    
        R = 0
        r = 0
        for i, (f0, f1) in enumerate(f):
            R += wts[i]*np.pi*(f1*np.sinc(2*(j-k)*f1) - f0*np.sinc(2*(j-k)*f0) + \
                               f1*np.sinc(2*(j+k)*f1) - f0*np.sinc(2*(j+k)*f0))
    
            r += wts[i]*D[i]*2*np.pi*(f1*np.sinc(2*j*f1) - f0*np.sinc(2*j*f0))
    
        a = np.dot(np.linalg.inv(R), r)
        a.shape = (-1,)

        h[:L] = a[:0:-1]/2.
        h[L] = a[0]
        h[L+1:] = a[1:]/2.
    else:
        # FIXME: Can we make this work for an even number of points?
        # It looks plausible that doing so should be possible.
        # L = N//2
  
        error_msg = 'Error: Increase (or decrease) the number of points by one'
        raise TransformRunException(error_msg, 1)

    return h    


def inverse_slr(na, mpldtmu, mpgaa, mpgbb, gamma):
    '''
    Perform the inverse SLR transformation. 
    Gives b1 for a set of a & b polynomials. 
       
    See also:
    J. Pauly, P. Le Roux, D. Nishimura, A. Macovski, 
    'Parmater Relations for the Shinnar-Le Roux Selective Excitation Pulse Design Algorithm',
    IEEE Transactions on Medical Imaging, Vol. 10, No. 1, pp 53-65, March 1991.
    '''
       
    # Adapted from Matlab, MPLPMB1F.M
       
    # Define local constants (mtesla)
    
    # gamma - user selected MHz/T
    const = (gamma/2) * (mpldtmu/1000.0) * 2.0 * np.pi
    a = mpgaa
    b = mpgbb 

    # Generate coefficients bc and ac 

    bc = sp.fftpack.ifft(b) 
    ac = sp.fftpack.ifft(a)
 
    # Preserve orig bc and ac (as 'aa' and 'bb')
    bb = bc.copy() 
    aa = ac.copy()

    # Inverse SLR Transform

    # idx = na
    # changed to na-1 to adapt to zero based indexing.
    idx = na-1
    
    etheta = np.zeros(na, dtype=complex)
    b1 = np.zeros(na, dtype=complex)
    c  = np.zeros(na, dtype=complex)
    s  = np.zeros(na, dtype=complex)

    while idx > -1:  # >= 0
    
        #etheta(idx) = bc(1)/ac(1)/(i*abs(bc(1)/ac(1)))
        # DCT 02/27/09 - convert to bc(0), etc.
        etheta[idx] = ( bc[0]/ac[0] ) / (1j * np.abs(bc[0]/ac[0]))

        b1[idx] = (1/const) * (etheta[idx]) * np.arctan( np.abs(bc[0]/ac[0]) )
                       
        c[idx] = np.cos( const * np.abs(b1[idx]) )
        s[idx] = 1j * etheta[idx] * np.sin(const * np.abs(b1[idx]))

        acu = np.delete(ac, 0)
        acl = np.delete(ac, idx)        
        bcu = np.delete(bc, 0)
        bcl = np.delete(bc, idx)
        
        qc = c[idx]        
        qs = s[idx]
        
        ac =  qc*acl + (qs.conj())*bcl    
        bc = -qs*acu + qc*bcu
      
        idx = idx-1

    # b1 = -np.imag(b1) is bad! It converts the complex array to a real array.
    # Doing as shown below preserves the values as complex numbers.
    for ii in range(len(b1)):
        b1[ii] = -np.imag(b1[ii])
    
    return b1
    
    
    
#------------------------------------------------------------

if __name__ == '__main__':
    
    import vespa.public.transform_description as transform_description
    import vespa.common.constants as constants
    import pylab as pyl
    
    extra = {}
    extra['calc_resolution'] = 2000
    extra['gamma']                = 42.576  # 1H gamma
    extra['pulse_bandwidth_type'] = 0
    extra['zero_padding']    = 0           # from machine specs
    
    params = {}
    params['time_steps']    = 250       
    params['tip_angle']     = 90.0
    params['duration']      = 8.0
    params['bandwidth']     = 1.0
    params['separation']    = 0.0
    params['is_dual_band']  = 0
    params['phase_type']    = 0     # linear
    params['use_lsq']       = 0
    params['pass_ripple']   = 2.0
    params['reject_ripple'] = 2.0
    
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
    bloch_inputs['bloch_range_value'] = 1.0    # unit set below
    bloch_inputs['bloch_range_units'] = 'khz'    # 'cm' or 'khz'
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
    
    
#     result.update_profiles(extra['calc_resolution'])
#     
#     pry, prx = result.get_profile(constants.UsageType.EXCITE, False, False)
# 
# 
#     fig = pyl.figure()
# 
#     axes1 = fig.add_subplot(211)
#     axes1.set_xlabel('Time [ms]')
#     axes1.set_ylabel('[mTesla]')
#     axes1.plot(rfx*1000, rfy)
#     
#     axes2 = fig.add_subplot(212)
#     dat = abs(pry)
#     pad = dat.max() * 0.05    
#     axes2.set_autoscaley_on(False)
#     axes2.set_ylim([-pad,dat.max()+pad])
#     axes2.set_xlabel('Position [cm]')
#     axes2.set_ylabel('Magnetization |Mz|')
#     axes2.plot(prx, abs(pry))
# 
# #     axes3 = fig.add_subplot(313)
# #     dat = abs(gry)
# #     pad = dat.max() * 0.05    
# #     axes3.set_autoscaley_on(False)
# #     axes3.set_ylim([-pad,dat.max()+pad])
# #     axes3.set_xlabel('Position [cm]')
# #     axes3.set_ylabel('Gradient')
# #     axes3.plot(grx, abs(gry))
#     
#     pyl.show()
    
    bob = 10
    
    
