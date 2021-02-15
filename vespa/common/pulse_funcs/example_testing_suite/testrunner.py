# Python modules


# 3rd party modules
import numpy as np
import cmath
import math

# Our modules
import vespa.common.pulse_funcs.slr_pulse as slr
import vespa.common.pulse_funcs.pulse_func_exception as pulse_func_exception
import vespa.common.pulse_funcs.constants as pf_constants

# if available import pylab (from matplotlib)
try:
    import matplotlib.pylab as plt
except ImportError:
    print("Error importing matplotlib.pylab")

# Adopted from Matpulse, MPLPM.M

RADIANS_TO_DEGREES = 180/math.pi

def run_tests():

    fileout = 'slr_output.txt'
    npoints = 64
    calc_resolution = 5000
    # dwell time in microseconds 
    dwell_time = 125
    # Pulse length (duration) in milliseconds
    pulse_duration = 8.
    
    # Tip angle in radians
    tip_angle = np.pi/2.
    # Bandpass ripple percent
    pass_ripple = 1.0
    # Reject ripple percent
    reject_ripple = 1.0
    # phase type: 1=>Linear, 2=>Max, 3=>Min      
    phase_type = 1
    # Generates bandwidth in kHz
    bandwidth = 1.0    
    is_single_band = True
    # Separation is not used if is_single_band = True
    separation = 0.0
    # If use_remez is False, uses weighted least squares. 
    use_remez = True
    zero_padding = 0
    bandwidth_convention=0
    
    # Other usage types: 
    # "excite", "inversion", "saturation", "spinecho"
    usage_type = "none"
    extended_profile=False

    # NOTE: pulse_duration must equal npoints*dwell_time/1000
    assert ( abs(pulse_duration - npoints*dwell_time/1000) < pf_constants.epsilon )
    
    # NOTE: calc_resolution should be greater than 4*npoints
    assert (calc_resolution >= 4*npoints)
    
    try:
        rf_y, rf_x, profiles = slr.slr_pulse(npoints, calc_resolution, dwell_time, 
                                             pulse_duration, tip_angle, pass_ripple, 
                                             reject_ripple, phase_type, bandwidth, 
                                             is_single_band, separation, use_remez, 
                                             zero_padding, bandwidth_convention, 
                                             usage_type, extended_profile)
        
    except pulse_func_exception.PulseFuncException as pfe:
        print("\n" + "Error Generating Pulse: " + pfe.message + "\n")
        return pfe.code


    if rf_y.any():
        
        # Print data to file and get real and imaginary components of rf_y
        f = open(fileout, 'w')
        for y in rf_y.tolist():
            phase, angle = cmath.polar(y)
            str = "%s %s\n" %(phase, angle*RADIANS_TO_DEGREES)
            f.write(str)

        f.close()



if __name__ == "__main__":
    run_tests()
    
