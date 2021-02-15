# Python modules


# 3rd party modules
import numpy as np
import cmath
import math

# Our modules
from . import bloch_lib_hargreaves as bloch_h
from . import bloch_multi
import vespa.common.pulse_funcs.slr_pulse as slr
import vespa.common.pulse_funcs.pulse_func_exception as pulse_func_exception
import vespa.common.pulse_funcs.constants as pf_constants

import vespa.common.bloch_call as bloch_call

# if available import pylab (from matplotlib)
try:
    import matplotlib.pylab as plt
except ImportError:
    print("Error importing matplotlib.pylab")

from pylab import *


RADIANS_TO_DEGREES = 180/math.pi

GAMMA = 42.576
TWOPI = 6.283185

def run_tests():

    gamma = 26751.3  # 1H Hz/gauss

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
        for i in range(10):
            
            bout = bloch_call.calc_all_profiles(rf_y, dwell_time, calc_resolution)

#         mz_std    = bout[0]
#         mxy_std   = bout[1]
#         xaxis_std = bout[2]
#         mz_ext    = bout[3]
#         mxy_ext   = bout[4]
#         xaxis_ext = bout[5]     
#          
#         bob = 10
#         bob += 1
#          
#         mxy = np.zeros(len(mz_std), np.complex128) 
#         mxy = np.zeros(len(mz_std), np.complex128) 
#         for i in range(len(mz_std)):
#             mxy[i] = complex(mz_std[i][0], mz_std[i][1])
#  
#         subplot(3,1,1)
#         xlabel('Time [ms]')
#         plot(rf_x,rf_y)
#         subplot(3,1,2)
#         plot(xaxis_std, abs(mxy), xaxis_std, real(mxy))  #, xaxis_std, imag(mxy) )
#         xlabel('Position [cm]')
#         ylabel('Magnetization |Mxy|')
#         subplot(3,1,3)
#         plot(xaxis_std, mz_std[:,2])
#         xlabel('Position [cm]')
#         ylabel('Magnetization Mz')
#         show()

        #----------------------------------------
        # HARGREAVES ENGINE
        #----------------------------------------

        # Calc.Nyquist
        nq = (1./(2.*dwell_time)) * 1000.0
        
        # Array for spatial values
        vn = nq*np.arange(-calc_resolution/2,(calc_resolution/2))/float(calc_resolution/2)
        
        # in mtesla -> 10x for Gauss
        gx = 10 * 0.25 * vn / GAMMA
        
        T = dwell_time / 1000000.0
        b1 = rf_y
        g  = gx
        b1 = b1 * 10.0      #0.5*b1/np.max(b1);
        x = vn          #np.arange(-5,5,0.05)
        f = np.array([0.0,])  # np.arange(-1000.0,2000.0,200.0)
        t = np.arange(1,len(b1)+1)*T;
 
        for i in range(10):
            mx, my, mz = bloch_h.bloch(b1,g,t,1,.2,f,x,mode=0,gamma=gamma)
 
     
#        for i in range(10):
#            mx, my, mz = bloch_h.bloch(b1,g,t,1,.2,f,x,mode=0,gamma=gamma)
#
#        for i in range(10):
#            mx, my, mz = bloch_multi.bloch_multi(b1,g,t,1,.2,f,x,mode=0,gamma=gamma)

     
#         mxy = -mx + 1j*my
#         ioff = 1
#         subplot(3,1,1)
#         xlabel('Time [ms]')
#         plot(t*1000,b1)
#         subplot(3,1,2)
#         plot(x, abs(mxy[:]), x, real(mxy[:]) )      #, x, imag(mxy[:]) )
#         xlabel('Position [cm]')
#         ylabel('Magnetization |Mxy|')
#         subplot(3,1,3)
#         plot(x, mz[:])
#         xlabel('Position [cm]')
#         ylabel('Magnetization Mz')
#         show()

                
        



if __name__ == "__main__":
#    run_tests()

    import cProfile
    cProfile.run('run_tests()')
    
