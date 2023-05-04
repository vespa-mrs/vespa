# Python Modules.

import copy
import math
import cmath

# 3rd party stuff
import numpy as np

# Local imports
import vespa.common.rfp_rf_result as rfp_rf_result
from vespa.common.transform_run_exception import TransformRunException


# Gyromagnetic ratio of 1H - the hydrogen nucleus. (units: kHz/mT)
GAMMA1H = 42.576    

DEGREES_TO_RADIANS = math.pi / 180



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
    # unpack and convert parameters if needed

    file_path       =  param['file1']
    dwell_time      =   int(param["dwell_time"])    # float, usec
    max_intensity   = float(param["max_intensity"]) # float, uT
    phase_units     =   int(param["phase_units"])   # choice, 0-degs, 1-rads
    
    
    # these extra items are used in making profile
    bandwidth_convention = int(extra['pulse_bandwidth_type'])  # Choice, 0-HalfHeight, 1-Min, 2-Max    
    resolution           = int(extra['calc_resolution']) # int
    gamma                = float(extra['gamma'])         # float MHz/T -> 2pi * 100 -> Hz/gauss

    is_phase_degrees = True if phase_units == 0 else False

    #--------------------------------------------------------------------------
    # pulse creation code starts here

    rf_waveform = read_pulse(file_path, max_intensity, is_phase_degrees)
    
    dwell = dwell_time/1000000
    npts  = len(rf_waveform)

    rf_xaxis = np.arange(npts) * dwell
    
    gradient = None
    grad_axis = None

    #--------------------------------------------------------------------------
    # fill and output dictionary, if needed

    outputs = {}
    outputs['out_time_steps'] = npts
    outputs['out_duration'] = npts * dwell * 1000.0  # msec

    #--------------------------------------------------------------------------
    # ensure you are returning ndarrays or None

    rf_y = np.array(rf_waveform)
    rf_x = np.array(rf_xaxis)
    gr_y = gradient
    gr_x = grad_axis


    return rf_y, rf_x, gr_y, gr_x, outputs
    

def _crect(r, phi):
    # The useful function rect() wasn't added to the cmath module until
    # Python 2.6. We use it when it exists. Python 2.5 users get the 
    # simplified implementation (which doesn't handle inf correctly).
    if hasattr(cmath, "rect"):
        return cmath.rect(r, phi)
    else:
        return r * (math.cos(phi) + math.sin(phi)*1j)

    
def read_pulse(file_path, max_intensity=1.0, is_phase_degrees=True):
    
    """ 
    file_path - complete path of file to read, as a string.
    max_intensity - the maximum intensity to which to scale the 
                    pulse at its absolute peak, as a float
    is_phase_degrees - bool, 
    Return a numpy array containing the pulse field intensity
    as a function of time. 

    """
    try:
        f = open(file_path)        
    except IOError as ioe:
        errstr = "File I/O Error: %s" % ioe.message
        raise TransformRunException(errstr, -1)          

    rf_y = []

    for ll, line in enumerate(f.readlines()):
        line = line.strip()

        # Ignore blank lines and comments
        if (not len(line)) or line.startswith('#') or line.startswith(';'):
            continue
            
        # Split any/all whitespace; also strip that whitespace from the substrings
        substrings = line.split()

        # validate that this data line has a valid format.
        if len(substrings) != 2:
            errstr = "Only two arguments per line allowed, found %d: Error on line %d" %(len(substrings), ll)
            raise pulse_func_exception.PulseFuncException(errstr, -1)

        amplitude, phase = substrings
        
        try:
            amplitude = float(amplitude)
        except ValueError as tpe:
            errstr = "Could not convert %s into a float, on line %d" %(amplitude, ll)
            raise TransformRunException(errstr, -1)
        
        try:
            phase = float(phase)
        except ValueError as tpe:
            errstr = "Could not convert %s into a float, on line %d" %(phase, ll)
            raise TransformRunException(errstr, -1)            
        
        ph_radians = phase
        if is_phase_degrees == True:
            # Convert phase to radians.
            ph_radians = phase * DEGREES_TO_RADIANS
        c = _crect(amplitude, ph_radians)
        rf_y.append(c)
        
    rf_y = np.array(rf_y)
    
    if len(rf_y) < 2:
        errstr = "Pulse file has less than 2 points"
        raise TransformRunException(errstr, -1)             
        
    max_ = np.max(np.abs(rf_y))
    # Scale pulse to user specified maximum intensity.
    # The 1000 in the denominator is used because we internally
    # represent these pulses in milliTesla (but display as microTesla).
    rf_y = rf_y * max_intensity / (max_*1000)
    
    return rf_y


    
    
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
    params['file1']         = 'adiabatic_90_v2.txt'
    params['dwell_time']    = 64.0       
    params['max_intensity'] = 20.0
    params['phase_units']   = 0

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
    
    
