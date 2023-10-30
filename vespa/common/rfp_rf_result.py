# Python modules

import copy
import xml.etree.cElementTree as ElementTree
import math

# 3rd party modules
import numpy as np

# Our modules
from vespa.common.constants import Deflate
import vespa.common.util.xml_ as util_xml
import vespa.common.util.time_ as util_time
import vespa.common.constants as constants

import vespa.common.pulse_funcs.bloch_lib_hargreaves as bloch_lib_h




class PulseFuncException(Exception):
    """
    Raised when am error occurs in one of the modules in the pulse_funcs
    directory. Includes an error string, that should be displayed to the
    user.
    """
    
    def __init__(self, error_str, error_code):
        self._message = error_str
        self.code    = error_code 


    __doc = """Getter/setter for the message associated with this exception"""
    # exceptions.Exception has a .message attribute, but it's deprecated in
    # Python 2.6 so we re-implement it ourselves.
    def __get_message(self):
        return self._message
    def __set_message(self, message):
        self._message = message

    message = property(__get_message, __set_message, doc=__doc)


class RfResults(object):
    """
    The result-set of a 'transformation' in the context of a synthesizing and
    modifying a radio frequency pulse, its gradient(s), and (raw) frequency 
    profile(s).

    """
    # The XML_VERSION enables us to change the XML output format in the future 
    XML_VERSION = "1.0.0"

    def __init__(self, attributes=None):
        """
        The attributes param to init() should contain 

        created     - A DateTime object representing the creation time.
        
        rf_waveform - list, complex numbers, units in microTesla
        rf_xaxis    - list, floats, units of time in seconds

        gradient    - list, complex numbers, units in microTesla
        grad_xaxis  - list, floats, units of time in seconds
        
        Note. display should be in milliseconds
        
        """
        self.created = util_time.now()
        
        # RF Waveform 
        self.rf_waveform = None     # [mT]  was rf.waveform
        self.rf_xaxis    = None     # [sec] was rf.waveform_x_axis

        # Gradient Waveform 
        self.gradient    = None     # [mT/m] was gradient.gradient_waveform
        self.grad_xaxis  = None     # [sec]  was gradient.time_points
 
        # Waveform and x-axis unit temporary results

        self.mz         = None  # bloch results for M starting out relaxed
        self.mxy        = None  # bloch results for M start in transverse plane
        self.xaxis      = None
        
        self.mz_ext    = None
        self.mxy_ext   = None
        self.xaxis_ext = None
        
        self.b1immunity = None
        self.last_b1immunity = None
                
        # Gradient Refocusing fractional Value
        self.grad_refocus_fraction = 0.0

        self.refocused_profile = None
        
        # An rfp_opcon_state.OpConState object. - need to create this object if used
        self.opcon_state = None
   
        # Explicit test for None necessary here.
        if attributes is not None:
            self.inflate(attributes)


    @property
    def dwell_time(self):
        '''Returns the dwell_time in microseconds'''
        if len(self.rf_xaxis) > 1:
            return 1000000 * (self.rf_xaxis[1] - self.rf_xaxis[0])
        else:
            return 0.0

    @property
    def first_gradient_value(self):
        '''returns value of first non-zero term in gradient in mT/m'''
        if self.gradient is None:
            return 10.0              # default is 10 mT/m == 1 gauss/cm
        else:
            ig = (np.nonzero(self.gradient))[0][0]  # first non-zero index
            return self.gradient[ig]                # return mT/m

      
    def gradient_refocusing(self, val, xunits, gamma):
        # This function calls a routine that may throw an exception of type
        #  PulseFuncException

        self.grad_refocus_fraction = val        # assume val if for kHz units xaxis

        profile, ax = self.get_profile(constants.UsageType.EXCITE)
        length = len(self.rf_xaxis)
        pulse_dur_ms = 1000.0*(self.rf_xaxis[length-1]) + self.dwell_time/1000.0

        gamma0 = gamma * 0.1                    # for 1H -> 4.2576 kHz/gauss
        g0 = self.first_gradient_value * 0.1    # for gauss/cm
        scale = gamma0 * g0                     # scaled for gradient

        self.refocused_profile = np.exp(1j*2.0*math.pi*val * scale * ax*pulse_dur_ms)*profile

        return val


    def gradient_refocusing_auto(self, xunits, gamma):
        # This automatically calculates a value for Gradient Refocus percentage
        # It does it for an x-axis in kHz units, regardless of what the display
        # in the GUI is set to.  The 'scale' keyword provides the conversion
        # factor for cm to kHz if needed. It's set to 1.0 otherwise, no harm no
        # foul.
        #
        # self.grad_refocus_fraction is always assumed to have been relevant to
        # an x-axis in kHz units.
        #
        # This function calls a routine that may throw an exception of type
        #  PulseFuncException

        gamma0 = gamma * 0.1                    # for 1H -> 4.2576 kHz/gauss
        g0 = self.first_gradient_value * 0.1    # for gauss/cm
        scale = gamma0 * g0                     # scaled for gradient

        profile, ax = self.get_profile(constants.UsageType.EXCITE)
        length = len(self.rf_xaxis)
        pulse_dur_ms = 1000.0 * (self.rf_xaxis[length - 1]) + self.dwell_time / 1000.0

        try:
            val = self.grad_refocus(ax, profile, pulse_dur_ms)
            val2 = self.grad_refocus_orig(ax, profile, pulse_dur_ms)
            print('orig val = '+str(val2)+'  new val = '+str(val))
            val = val / scale
        except PulseFuncException:
            val = 0.5
            self.grad_refocus_fraction = val

        self.refocused_profile = np.exp(1j*2.0*math.pi*val*ax * scale * pulse_dur_ms)*profile
        self.grad_refocus_fraction = val

        return val

    def interpolate(self, a, a0, a1, b0, b1):
        if not ((a0 <= a and a <= a1) or (a1 <= a and a <= a0)):
            error_message = 'Interpolation Failed. Target point is not between end points'
            raise PulseFuncException(error_message, 1)
        # Approximate from the closer point  
        if abs(a-a0) < abs(a1-a):
            val = b0 + (a-a0)*(b1-b0)/(a1-a0)
        else:
            val = b1 - (a1-a)*(b1-b0)/(a1-a0)

        return val

    def get_bandwidth(self, ax, profile):

        start = False
        middle = False
        end = False
        minimum = False
        restart = False
        x1 = 0
        x2 = 0
        for i in range(len(ax)):
            if abs(profile[i]) >= 0.5 and not start:
                x1 = ax[i]
                if i>0 and abs(profile[i]) != 0.5: # Then interpolate
                    x1 = self.interpolate(.5, abs(profile[i-1]), abs(profile[i]), ax[i-1], ax[i])
                start = True
            if abs(profile[i]) > .8 and not middle:
                middle = True
            if middle and not end and abs(profile[i]) <= 0.5:
                x2 = ax[i]
                if i>0 and abs(profile[i]) != 0.5: # Then interpolate
                    x2 = self.interpolate(.5, abs(profile[i-1]), abs(profile[i]), ax[i-1], ax[i])
                end = True
            if end and not minimum and abs(profile[i]) <= 0.25:
                minimum = True
            if minimum and not restart and abs(profile[i]) >= 0.5:
                error_message = 'Interpolation Failed. This appears to be a dual band pulse'
                raise PulseFuncException(error_message, 1)

        if start and middle and end:
            return x2 - x1
        else:        
            return -1


    def grad_refocus(self, xaxis0, yprofile0, pulse_dur_ms):
        """
        For single band pulses: Calculates the optimal
        refocusing gradient as a fractional real number.

        xaxis0 - The "x-axis" points of the profile - as an array (units of kHz)
        yprofile0 - The "y-axis" points of the profile (units of fractional magnetization -1 to 1)
        pulse_dur_ms - The pulse length in milliseconds.
        pulse_bw_khz - The pulse bandwidth in khz.

        BJS - just FYI, after beating my head against this, here's my wisdom ...
         The value that comes out of this is the factor that multiplies the
         waveform to minimize the area of the imag, and max the area of the
         real plot(s). And that is based on the 'standard' gradient value that
         is typically set (10mT/m, 1g/cm) if no gradient is calculated for the
         waveform algorithm.  So if I get 0.515, that means that I have to
         apply 51.5% of the area of the gradient to refocust the signal. In
         Siemens-world, the PulseAsymmetry() is set as the % from the START
         of the waveform to the 'phase-center' of the pulse. So I would have
         to set that to 1-0.515 (in this example) to properly set their pulse
         attribute.

        Oct 30, 2023 - bjs first pass refactor, mainly 1st half the code, not the actual algorithm

        """
        # Get the bandwidth from the profile
        pulse_bw_khz = self.get_bandwidth(xaxis0, yprofile0)
        if pulse_bw_khz < 0:
            error_message = 'Calculated bandwidth was less than zero'
            raise PulseFuncException(error_message, 1)           

        # Get axis range over just bandwidth
        xaxis = xaxis0.copy()
        yprofile = yprofile0.copy()
        bw_half = pulse_bw_khz/2.0

        ax1  = next(x[0] for x in enumerate(xaxis) if x[1] > -bw_half)
        ax2  = len(xaxis) - next(x[0] for x in enumerate(xaxis[::-1]) if x[1] < bw_half) - 1
        ax11 = next(x[0] for x in enumerate(xaxis) if x[1] > -2*bw_half)
        ax22 = len(xaxis) - next(x[0] for x in enumerate(xaxis[::-1]) if x[1] < 2*bw_half) - 1

        if ax1==0 or ax2==len(xaxis) or ax11==0 or ax22==len(xaxis):
            error_message = 'Unable to locate the bandwidth edges'
            raise PulseFuncException(error_message, 1)

        if (ax1 == ax2) or (ax11 == ax22):
            error_message = 'Invalid Bandwidth edges: They are the same point'
            raise PulseFuncException(error_message, 1)


        # FIRST, ENSURE THAT REFOCUS AXIS IS Y - Find angle at zero frequency
        indx = int(ax1 + math.ceil((ax2-ax1)/2.0))
        rot = np.angle(yprofile[indx])

        yprofile = np.exp(-1j*rot)*yprofile

        # INITIAL REFOCUS APPROXIMATION
        # - Use fft to get no of cycles (fre)
        # - First zero-pad the profile (rpro)

        nrpro = np.pad(yprofile[ax1:ax2+1], [ax2+1-ax1, ax2+1-ax1])
        fre = np.fft.fftshift(np.fft.fft(nrpro))

        # If two values in fre are the same in magnitude matlab actually looks 
        # at the phase and finds the one with the biggest phase, which we're 
        # not doing here. However, since we are using the first value in the
        # array that matches maxfre, fre_pos = np.nonzero(...), the phase thing 
        # is not and issue - according to Jerry.

        fre_pos = np.argmax(np.abs(fre))

        # The +1" is added at the end to purposely throw off the fit so that
        # the iterations (below) give a good curve for getting the maximum area.
        fre_n = (1.0/3.0)*(fre_pos - len(nrpro)/2 +1)    

        # Get rotation factor, Apply refocus to profile (yprofile)
        rotf = -fre_n/(pulse_bw_khz*pulse_dur_ms)
        npro = np.exp(1j*2*np.pi*rotf*xaxis*pulse_dur_ms) * yprofile

        # Get initial area under real component, over twice the bandwidth
        # Added the next line so the following one would make sense.
        aarea = np.zeros(1)
        aarea[0] = np.sum(np.real(npro[ax11-1:ax22]))

        # ITERATIVE FINE RE-PHASING ******************************
        # - select just half the bandwidth of imag data

        lnhnax = int(math.floor((ax2+1-ax1) / 4.0))
        hnpro = np.imag(npro[ax1+lnhnax-1:ax2-lnhnax]) 
        hnax = xaxis[ax1+lnhnax-1:ax2-lnhnax]
        p = np.polyfit(hnax,hnpro,1)    # Slope of line is p[0], intercept is p[1]

        nrotf = np.zeros((40))
        barea = np.zeros((40))
        mtheta = math.atan(p[0])

        nrotf[0] = -mtheta/(2.0*np.pi)
        nnpro =  np.exp(1j*2.0*np.pi*nrotf[0]*xaxis*pulse_dur_ms/8.0) * npro
        barea[0] = np.sum(np.real(nnpro[ax11-1:ax22]))

        niter = 1
        while abs(p[0]) > 0.001 and niter < 19:
            hnpro = np.imag(nnpro[ax1+lnhnax-1:ax2-lnhnax]) 
            p = np.polyfit(hnax,hnpro,1)
            mtheta = math.atan(p[0])
            nrotf[niter] = -mtheta/(2.0*np.pi)
            nnpro =  np.exp(1j*2.0*np.pi*nrotf[niter]*xaxis*pulse_dur_ms/8.0)*nnpro
            barea[niter] = np.sum(np.real(nnpro[ax11-1:ax22]))
            niter = niter + 1

        # See if convergence failed.
        if niter >= 19:
            error_message = 'Convergence Failed'
            raise PulseFuncException(error_message, 1)   

        nrotf = nrotf[0:niter]
        refocus_value = (rotf + np.sum(nrotf)/8.0)

        # Continue rotations to symetrize data
        nrotf = np.append(nrotf, nrotf[::-1] )

        for niter in range(niter,(2*niter)):
            nnpro =  np.exp(1j*2.0*np.pi*nrotf[niter]*xaxis*pulse_dur_ms/8.0)*nnpro
            barea[niter] = np.sum(np.real(nnpro[ax11-1:ax22]))

        barea = np.append(aarea, barea[0:niter])
        carea =  1000.0/barea
        nrotf = nrotf[0:niter]
        frot = rotf + np.cumsum(nrotf/8)
        frot = np.append(rotf, frot)

        # FINE TUNE FOR MAX AREA UNDER REAL PROFILE
        pa = np.polyfit(frot,carea,2)
        y = np.polyval(pa,frot)

        # This is the value considering areas
        frotn = -pa[1]/(2.0*pa[0])    

        # If p[0] > .001, or p[1]>.001, or niter > 19, iterations failed?

        if frotn < 0:
            error_message = 'Invalid Result: Fraction less than zero'
            raise PulseFuncException(error_message, 1)           

        return frotn


    def grad_refocus_orig(self, xaxis0, yprofile0, pulse_dur_ms):
        """
        For single band pulses: Calculates the optimal
        refocusing gradient as a fractional real number.

        xaxis0 - The "x-axis" points of the profile - as an array (units of kHz)
        yprofile0 - The "y-axis" points of the profile (units of fractional magnetization -1 to 1)
        pulse_dur_ms - The pulse length in milliseconds.
        pulse_bw_khz - The pulse bandwidth in khz.

        """
        # Get the bandwidth from the profile
        pulse_bw_khz = self.get_bandwidth(xaxis0, yprofile0)
        if pulse_bw_khz < 0:
            error_message = 'Calculated bandwidth was less than zero'
            raise PulseFuncException(error_message, 1)

        # Get axis range over just bandwidth
        xaxis = xaxis0.copy()
        yprofile = yprofile0.copy()
        bw_half = pulse_bw_khz/2.0

        found_ax1 = False
        found_ax2 = False
        index = 0
        for ii in xaxis:
            if ii > -bw_half and not found_ax1:
                # Want the first one, so if find one, don't reassign.
                ax1 = index
                found_ax1 = True
            if ii < bw_half:
                # Keep reassigning, since want last one.
                ax2 = index
                found_ax2 = True
            index += 1

        found_ax11 = False
        found_ax22 = False
        index = 0
        for ii in xaxis:
            if ii > -2*bw_half and not found_ax11:
                # Want the first one, so if find one, don't reassign.
                ax11 = index
                found_ax11 = True
            if ii < 2*bw_half:
                # Keep reassigning, since want last one.
                ax22 = index
                found_ax22 = True
            index += 1

        # If not (found_ax1 and found_ax2), throw an error.
        if not (found_ax1 and found_ax2 and found_ax11 and found_ax22):
            error_message = 'Unable to locate the bandwidth edges'
            raise PulseFuncException(error_message, 1)

        if (ax1 == ax2) or (ax11 == ax22):
            error_message = 'Invalid Bandwidth edges: They are the same point'
            raise PulseFuncException(error_message, 1)

        # FIRST, ENSURE THAT REFOCUS AXIS IS Y - Find angle at zero frequency
        indx = int(ax1 + math.ceil((ax2-ax1)/2.0))
        rot = np.angle(yprofile[indx])

        yprofile = np.exp(-1j*rot)*yprofile

        # INITIAL REFOCUS APPROXIMATION
        # - Use fft to get no of cycles (fre)
        # - First zero-pad the profile (rpro)
        rpro = yprofile[ax1:ax2+1]

        nrpro = np.hstack( [ np.zeros(len(rpro)), rpro, np.zeros(len(rpro))] )

        fre0 = np.fft.fft(nrpro)
        length = len(fre0)
        fre =  np.append(fre0[int(math.ceil(length/2)):], fre0[:int(math.ceil(length/2))])    # this is a roll?

        # If two values in fre are the same in magnitude matlab actually looks
        # at the phase and finds the one with the biggest phase, which we're
        # not doing here. However, since we are using the first value in the
        # array that matches maxfre, fre_pos = np.nonzero(...), the phase thing
        # is not and issue - according to Jerry.
        absfre = np.abs(fre)
        maxfre = np.max(absfre)
        fre_pos = np.nonzero(absfre==maxfre)

        # Since fre_pos is an array of arrays, (most likely with only one
        # value), this pulls out the first value.
        fre_pos = fre_pos[0][0]

        # The +1" is added at the end to purposely throw off the fit so that
        # the iterations (below) give a good curve for getting the maximum area.
        fre_n = (1.0/3.0)*(fre_pos - len(nrpro)/2 +1)

        # Get rotation factor
        rotf = -fre_n/(pulse_bw_khz*pulse_dur_ms)

        # Apply refocus to profile (yprofile)
        npro = np.exp(1j*2*math.pi*rotf*xaxis*pulse_dur_ms)*yprofile

        # Get initial area under real component, over twice the bandwidth
        # Added the next line so the following one would make sense.
        aarea = np.zeros(1)
        aarea[0] = np.sum(np.real(npro[ax11-1:ax22]))

        # ITERATIVE FINE RE-PHASING ******************************
        # - Now select just half the bandwidth
        lnhnax = int(math.floor(len(rpro)/4.0))
        hnpro = np.imag(npro[ax1+lnhnax-1:ax2-lnhnax])
        hnax = xaxis[ax1+lnhnax-1:ax2-lnhnax]

        p = np.polyfit(hnax,hnpro,1)
        # Slope of line is p[0], intercept is p[1]

        nrotf = np.zeros((40))
        barea = np.zeros((40))

        mtheta = math.atan(p[0])

        nrotf[0] = -mtheta/(2.0*math.pi)
        nnpro =  np.exp(1j*2.0*math.pi*nrotf[0]*xaxis*pulse_dur_ms/8.0) * npro
        barea[0] = np.sum(np.real(nnpro[ax11-1:ax22]))

        niter = 1
        while abs(p[0]) > .001 and niter < 19:
            hnpro = np.imag(nnpro[ax1+lnhnax-1:ax2-lnhnax])
            p = np.polyfit(hnax,hnpro,1)
            mtheta = math.atan(p[0])
            nrotf[niter] = -mtheta/(2.0*math.pi)
            nnpro =  np.exp(1j*2.0*math.pi*nrotf[niter]*xaxis*pulse_dur_ms/8.0)*nnpro
            barea[niter] = np.sum(np.real(nnpro[ax11-1:ax22]))
            niter = niter + 1

        # See if convergence failed.
        if niter >= 19:
            error_message = 'Convergence Failed'
            raise PulseFuncException(error_message, 1)

        nrotf = nrotf[0:niter]
        refocus_value = (rotf + np.sum(nrotf)/8.0)

        # Continue rotations to symetrize data
        nrotf = np.append(nrotf, nrotf[::-1] )

        for niter in range(niter,(2*niter)):
            nnpro =  np.exp(1j*2.0*math.pi*nrotf[niter]*xaxis*pulse_dur_ms/8.0)*nnpro
            barea[niter] = np.sum(np.real(nnpro[ax11-1:ax22]))

        barea = np.append(aarea, barea[0:niter])
        carea =  1000.0/barea
        nrotf = nrotf[0:niter]
        frot = rotf + np.cumsum(nrotf/8)
        frot = np.append(rotf, frot)

        # FINE TUNE FOR MAX AREA UNDER REAL PROFILE
        pa = np.polyfit(frot,carea,2)
        y = np.polyval(pa,frot)

        # This is the value considering areas
        frotn = -pa[1]/(2.0*pa[0])

        # If p[0] > .001, or p[1]>.001, or niter > 19, iterations failed?

        if frotn < 0:
            error_message = 'Invalid Result: Fraction less than zero'
            raise PulseFuncException(error_message, 1)

        return frotn


    def update_profiles(self, outputs):  
        """
        Recalculate the profile(s) using the Bloch equation
        
        Results are stored as mT and mT/m, Hargreaves lib requires G and G/cm
        - mT -> G = x10
        - mT/m -> G/cm = 10/100 = 0.1
        
        """
        igamma = outputs['gyromagnetic_nuclei']
        gamma0 = constants.GAMMA_VALUES[igamma]
        gamma = 2.0 * np.pi * 100.0 * gamma0

        resolution         = outputs['calc_resolution']
        bloch_range_value  = outputs['bloch_range_value']
        bloch_range_units  = outputs['bloch_range_units']
        bloch_offset_value = outputs['bloch_offset_value']
        
        if self.rf_waveform is not None:
            # b1 needs to be Gauss for Hargreaves library, hence the x10
            b1 = self.rf_waveform * 10.0    # convert mT to gauss
            g, _  = self.get_gradient()     # yout,x-axis
            g = g * 0.1  # convert mT/m to gauss/cm

            dwell  = self.dwell_time
            npts   = resolution

            # Array for spatial values 
            # - nyquist times +/- resolution range
            # - Convert to cm for display
            
            if bloch_range_units == 'cm':
                nq = bloch_range_value
            else:
                # bloch_range_value units are in kHz
                # convert to cm using GAMMA and gradient strength
                gam_khz_per_gauss = gamma0 * 0.1        # MHz/T -> kHz/G = 1000kHz/10000G
                g0 = self.first_gradient_value * 0.1    # first non-zero gradient value
                nq = bloch_range_value / (gam_khz_per_gauss * g0)  # kHz/gauss * gauss/cm -> kHz/cm

            x  = nq*np.arange(-npts/2,npts/2)/float(npts/2)      # always in cm
            xx = 4.0*x.copy()      # For extended range

            # this all used to be in bloch_call.py
        
            t  = np.arange(1,len(b1)+1)*dwell*0.000001  # time axis in [sec]
            f  = np.array([bloch_offset_value,])        # on resonance

            mx1, my1, mz1 = bloch_lib_h.bloch(b1,g,t,None,None,f,x, mode=0, gamma=gamma)     # no decays
            mx2, my2, mz2 = bloch_lib_h.bloch(b1,g,t,None,None,f,x, mode=0,do_mxy=True, gamma=gamma)     # no decays
            mx3, my3, mz3 = bloch_lib_h.bloch(b1,g,t,None,None,f,xx,mode=0, gamma=gamma)     # no decays
            mx4, my4, mz4 = bloch_lib_h.bloch(b1,g,t,None,None,f,xx,mode=0,do_mxy=True, gamma=gamma)     # no decays

            self.mz      = np.hstack((mx1,my1,mz1))
            self.mxy     = np.hstack((mx2,my2,mz2))
            self.mz_ext  = np.hstack((mx3,my3,mz3))
            self.mxy_ext = np.hstack((mx4,my4,mz4))

            self.xaxis     = x.copy()   # always in cm
            self.xaxis_ext = xx.copy()  # always in cm


    def get_profile(self, use_type, extended=False, absolute=False):
        """
        We can display results for extended or standard width x-range, and
        for four usage types = Excite, Inversion, Saturation, SpinEcho.
        
        This is a total of 8 profiles we can extract from our RF complex
        waveform. 
        
        This method lets user specify which one to return. Note. all of these
        plots are already pre-calculated, we are just choosing which to show.
        
        """
        if use_type == constants.UsageType.EXCITE     or \
           use_type == constants.UsageType.INVERSION  or \
           use_type == constants.UsageType.SATURATION:
            
            if extended:
                y = self.mz_ext
                x = self.xaxis_ext
            else:
                y = self.mz
                x = self.xaxis
                
        elif use_type == constants.UsageType.SPIN_ECHO:
            
            if extended:
                y = self.mxy_ext
                x = self.xaxis_ext
            else:
                y = self.mxy
                x = self.xaxis

        if use_type == constants.UsageType.EXCITE:            
            yout = y[:,0] + 1j*y[:,1]
        elif use_type == constants.UsageType.SATURATION:
            yout = y[:,2] + 0j
        elif use_type == constants.UsageType.INVERSION:
            yout = -y[:,2] + 0j
        elif use_type == constants.UsageType.SPIN_ECHO:
            yout = y[:,0] - 1j*y[:,1]
        
        if absolute:
            yout = np.abs(yout)       

        return (yout, x)


    def get_gradient(self):
        
        if self.gradient is None: 
            x = self.rf_xaxis
            y = np.ones(len(x)) * 10    # default 1 gauss/cm = 10 mT/m
            # 1 gauss/cm = 0.1 mT / 0.01 m = 10 mT/m => 10 * 10gauss/100cm
        else:
            x = self.grad_xaxis
            y = self.gradient
            
        return (y,x)


    def update_b1_immunity(self, outputs):
        """
        Recalculate the matrix of B1 immunity using the Bloch equation
        
        """
        if self.rf_waveform is not None: 

            pulse_type    =   int(outputs["pulse_type"])
            freq_center   = float(outputs["freq_center"])
            freq_range    = float(outputs["freq_range"])   # in +/- Hz
            freq_steps    =   int(outputs["freq_steps"])
            b1_strength   = float(outputs["b1_strength"])
            b1_range      = float(outputs["b1_range"])     # in +/- uT
            b1_steps      =   int(outputs["b1_steps"])
            magn_vector   =   int(outputs["magn_vector"])
            ref_b1_limits = float(outputs["ref_b1_limits"])
            
            if self.last_b1immunity is not None:
                if self.last_b1immunity["pulse_type"]    == outputs["pulse_type"]  and  \
                   self.last_b1immunity["freq_center"]   == outputs["freq_center"] and  \
                   self.last_b1immunity["freq_range"]    == outputs["freq_range"]  and  \
                   self.last_b1immunity["freq_steps"]    == outputs["freq_steps"]  and  \
                   self.last_b1immunity["b1_strength"]   == outputs["b1_strength"] and  \
                   self.last_b1immunity["b1_range"]      == outputs["b1_range"]    and  \
                   self.last_b1immunity["b1_steps"]      == outputs["b1_steps"]:
                    # all params are the same, we don't need to run bloch sims
                    return
            
            b1 = self.rf_waveform.copy()    
            if self.gradient is None:         
                g = np.ones(len(b1))        # this is 1.0 G/cm as default
            else: 
                g = self.gradient * 0.1     # convert mT/m to G/cm
            dwell  = self.dwell_time

            # units are in kHz - convert to cm using GAMMA and gradient strength
            fctr = freq_center / (4.2576 * g[0])
            frng = freq_range  / (4.2576 * g[0])

            x = np.linspace(fctr-frng, fctr+frng, freq_steps)
            f = [0.0,]  
            t = np.arange(1,len(b1)+1)*dwell*0.000001  # time axis in [sec]
            
            b1max    = np.max(np.abs(b1)) * 1000    # mT to uT
            bcenter  = b1_strength
            brange   = b1_strength * (b1_range / 100.0)
            bstr     = (b1_strength - brange) / b1max
            bend     = (b1_strength + brange) / b1max
            b1scales = np.linspace(bstr, bend, b1_steps)
            
            resx = []
            resy = []
            resz = []
            
            if pulse_type == 3:     # spin-echo simulation
                do_mxy = True
            else:                   # excite, inversion, saturation simulation
                do_mxy = False
            
            for b1scale in b1scales: 
                
                # b1 needs to be Gauss for Hargreaves library
                b1in = b1.copy() * b1scale * 10.0  # mT to G

                mx1, my1, mz1 = bloch_lib_h.bloch(b1in,g,t,None,None,f,x, mode=0, do_mxy=do_mxy)     # no decays

                resx.append(mx1)
                resy.append(my1)
                resz.append(mz1)

            resx = np.array(resx)
            resy = np.array(resy)
            resz = np.array(resz)

            self.b1immunity = np.concatenate((resx,resy,resz), axis=2)
            
            self.last_b1immunity = outputs



    def deflate(self, flavor=Deflate.ETREE):
        if flavor == Deflate.ETREE:
            
            e = ElementTree.Element("result", {"version" : self.XML_VERSION})
            
            util_xml.TextSubElement(e, "created", self.created)
            
            item = util_xml.numeric_list_to_element(self.rf_waveform, complex, "rf_waveform")
            e.append(item)

            item = util_xml.numeric_list_to_element(self.rf_xaxis, float, "rf_xaxis")
            e.append(item)

            if self.gradient is not None:
                item = util_xml.numeric_list_to_element(self.gradient, float, "gradient")
                e.append(item)

            if self.grad_xaxis is not None:
                item = util_xml.numeric_list_to_element(self.grad_xaxis, float, "grad_xaxis")
                e.append(item)

#             if self.opcon_state:
#                 e.append(self.opcon_state.deflate(flavor))
                
            return e
        
        elif flavor == Deflate.DICTIONARY:
            return copy.copy(self.__dict__)


    def inflate(self, source):
        if hasattr(source, "makeelement"):

            # Quacks like an ElementTree.Element
            self.created = source.findtext("created")
            self.created = util_time.datetime_from_iso(self.created)
            
            # Need explicit test for None below; cast to bool doesn't work as expected.
            
            e = source.find("rf_waveform")
            if e is not None:
                item = util_xml.element_to_numeric_list(e)
                self.rf_waveform = np.array(item)

            e = source.find("rf_xaxis")
            if e is not None:
                item = util_xml.element_to_numeric_list(e)
                self.rf_xaxis = np.array(item)
                
            e = source.find("gradient")
            if e is not None:
                item = util_xml.element_to_numeric_list(e)
                self.gradient = np.array(item)

            e = source.find("grad_xaxis")
            if e is not None:
                item = util_xml.element_to_numeric_list(e)
                self.grad_xaxis = np.array(item)

#            e = source.find("opcon_state")
#            if e is not None:
#                self.opcon_state = rfp_opcon_state.OpConState(e)  # FIXME - bjs need to create this object
                
        elif hasattr(source, "keys"):
            # Quacks like a dict
            for key in list(source.keys()):
                if hasattr(self, key):
                    setattr(self, key, source[key])

            if "rf" in source:
                self.rf = source["rf"]




            
           
                    
                    