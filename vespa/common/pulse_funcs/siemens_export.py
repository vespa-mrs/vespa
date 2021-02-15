# Python modules


# 3rd party modules
import numpy as np


# Our modules
import vespa.common.constants as constants




def siemens_export(pulse_uT, pulse_type, dwell_us, bw_khz, slice_thick_max, 
                   slice_thick_min, pulse_name, comment, output_format, 
                   family_name='', tip_angle= 90.0, var_base='MyPulse',
                   output_filename='MyPulse.h'):
                   
    # Arguments:
    #
    # pulse_uT        (numpy array) - B1,B2, or G2 field
    #                  for B1,B2 - complex waveform in mT
    #                  for G2    - not implemented
    # pulse_type       values specifying what pulse_uT is. Must be one of the 
    #                  constants from constants.ThirdPartyExportFormat.FIELD_XXX
    # dwell_us        (float)  - dwell time in usec
    # bw_khz          (float)  - band width in kHz
    # slice_thick_max (float)  - slice thickness maximum in mm
    # slice_thick_min (float)  - slice thickness minimum in mm
    # pulse_name      (string) - pulse name
    # comment         (string) - comment
    # output_format    what format to use. Must be either IDEA or VISION 
    #                  from constants.ThirdPartyExportFormat
    # family_name     (string) - family name for Vision format
    #
    # Returns:
    # 
    # A list of text lines in the format requested 
    
    # The sample file we have has a very consistent format of exactly 9 places
    # after the decimal and at least 1 before it, so we make sure our 
    # format matches that.
    DEFAULT_FLOAT_FORMAT = "%1.9f"

    # These guys are floats; we want them to be strings.
    slice_thick_min = DEFAULT_FLOAT_FORMAT % slice_thick_min
    slice_thick_max = DEFAULT_FLOAT_FORMAT % slice_thick_max

    # Working with RF pulse (rather than gradient)
    if pulse_type in (constants.ThirdPartyExportFormat.FIELD_B1, 
                      constants.ThirdPartyExportFormat.FIELD_B2):
        
        # Normalize to unity --------------------------------------------------
        pulse_magn      = np.abs(pulse_uT)
        max_magn        = np.max(pulse_magn)
        pulse_uT_norm   = pulse_uT *(1.0/max_magn)
        pulse_magn_norm = np.abs(pulse_uT_norm)
        
        
        # Amplitude integral --------------------------------------------------  
        #   This algorithm to calculate the amplitude integral was pulled from 
        #   Siemens documentation.  This is different than the one that Assaf
        #   sent me, and even the Siemens docs say that we may need a pulse 
        #   specific answer for and SLR or adiabatic pulse ... sigh
        
        pulse_uT_norm_sum = sum(pulse_uT_norm)
        amp_integral   = pulse_uT_norm_sum.real**2 + pulse_uT_norm_sum.imag**2
        amp_integral   = np.sqrt(amp_integral)

        abs_integral   = np.sum(pulse_magn_norm)
        power_integral = sum(np.power(pulse_magn_norm,2))

        power_integral = DEFAULT_FLOAT_FORMAT % power_integral
        amp_integral   = DEFAULT_FLOAT_FORMAT % amp_integral
        abs_integral   = DEFAULT_FLOAT_FORMAT % abs_integral
        
        # Calclulate pulse length in [ms] -------------------------------------
        duration_us = (len(pulse_uT) * dwell_us)
        duration_ms = duration_us / 1000.0
    else:
        # pulse_type == ThirdPartyExportFormat.FIELD_G2 
        raise NotImplementedError("Exporting gradient isn't supported")
        
    ref_grad = compute_reference_gradient_siemens(duration_ms, bw_khz, csa=0)    
    ref_grad = DEFAULT_FLOAT_FORMAT % ref_grad     
    
    # Array of text lines to return
    lines = [ ]
    
    if output_format == constants.ThirdPartyExportFormat.IDEA:

        datas = np.array(np.hstack((np.abs(pulse_uT_norm), np.angle(pulse_uT_norm)))).conj().T

        lines.append( ("PULSENAME", pulse_name) )
        lines.append( ("COMMENT",   comment) )
        lines.append( ("REFGRAD",   ref_grad) )
        lines.append( ("MINSLICE",  slice_thick_min) )
        lines.append( ("MAXSLICE",  slice_thick_max) )
        lines.append( ("AMPINT",    amp_integral) )
        lines.append( ("POWERINT",  power_integral) )
        lines.append( ("ABSINT", abs_integral) )

        data_arr = datas
        length_data_arr = len(data_arr)
        magn_arr  = data_arr[0:int(length_data_arr/2.0)]
        phase_arr = data_arr[int(length_data_arr/2.0+1.0)-1:length_data_arr]
        
        # Format the header lines
        lines = ["%s:\t%s" % line for line in lines]
        lines.append("\n")
        
        format = DEFAULT_FLOAT_FORMAT + "\t" + DEFAULT_FLOAT_FORMAT + "\t; (%d)"
        for i in np.arange(1., (length_data_arr/2.)+1):
            if magn_arr[int(i) - 1] < 0.:
                line = format % (-1.*magn_arr[int(i)-1], np.pi, (i-1.))
            else:
                line = format % (magn_arr[int(i)-1], phase_arr[int(i)-1], (i-1.))
            lines.append(line)

    elif output_format == constants.ThirdPartyExportFormat.VISION:
        
        datas = np.array(np.hstack((np.abs(pulse_uT_norm), np.angle(pulse_uT_norm)))).conj().T
        
        npts_text = "%d" % len(pulse_uT)
        
        lines.append( ("Begin_Entry",       "VESPA_RFPulse_Generated_Pulse") )
        lines.append( ("Entry_Type",        6) )
        lines.append( ("Pulse_Name",        pulse_name) )
        lines.append( ("Entry_Description", comment) )
        lines.append( ("Num_Points",        npts_text) )
        lines.append( ("Family_Name",       family_name) )
        lines.append( ("Slice_Thick_Min",   slice_thick_min) )
        lines.append( ("Slice_Thick_Max",   slice_thick_max) )
        lines.append( ("Reference_Grad",    ref_grad) )
        lines.append( ("Power_Integral",    power_integral) )
        lines.append( ("Ampl_Integral",     amp_integral) )
        lines.append( ("Envelope_Mode",     1) )
        
        lines = ["%s:\t%s" % line for line in lines]
        # Entry Values
        lines.append("Entry_Values:")
        for data in datas:
            lines.append(DEFAULT_FLOAT_FORMAT % data)
        lines.append("End_Entry:\t")

    elif output_format == constants.ThirdPartyExportFormat.IDEA_C_HEADER:

        # Description: exports a given pulse structure (sans gradients) to
        # a C++ include (.h) file, which can then be included in IDEA code
        # and loaded into an arbitrary pulse structure. 
        # 
        # Inputs:
        # 
        # Variable Name   Units    Description
        #  pulse           -        Input RF pulse structure having the following
        #                           fields:
        #                           pulse.tp - duration of pulse, in ms
        #                           pulse.RFamp - amplitude envelope (kHz)
        #                           pulse.RFphase - phase of pulse as func. of time
        #                           (radians)
        #                           pulse.Gx - amplitude of x-gradient (in kHz/mm)
        #                           pulse.Gy - amplitude of y-gradient (in kHz/mm)
        #                           pulse.Gz - amplitude of z-gradient (in kHz/mm)
        #                           Gx, Gy, Gz, RFamp, RFphase are all arrays with
        #                           the same number of elements.
        #  pulsename       -        Name of RF pulse 
        #  filename        -        Output filename (including directory!) [1]
        #  pulseBW         kHz      Bandwidth of pulse, in kHz. Necessary for 
        #                           calibrating reference gradient.
        #  flipAngle       rad.     Flip angle of pulse, in radians. Be sure to set
        #                           the same flip angle in the Siemens sequence, using:
        #                           myExternalPulse.setFlipAngle(flipAngle in degrees)
        #  comment         -        Comment to be added. Can be left blank (i.e.: ''),
        #                           in which case 'Just a pulse' will be used.
        #  minSlice        mm       Minimum excitable slice. [2]
        #  var_base         -        String, used to define the variables in the
        #                           include file. [3]
        # 
        #  [1] The filename will always have a .h extension (if omitted, or another
        #      is provided, a .h will be placed automatically instead).
        # 
        #  [2] A seemingly "harmless" quantity, but one that should be set by the
        #      designer for the following reason: the pulse has a certain stepsize
        #      (dwell time), which determines the overall affected bandwidth:
        #      FOVBW = 1/(dwell time)
        #      Beyond that, the spectral pattern of the pulse will replicate
        #      itself. Thus, a pulse designed to excite a 1 cm slice with a 10 cm
        #      FOV will excite a 1 mm slice with a 1 cm FOV if the gradient is 
        #      multiplied by 10. As a result, if the object is 10 cm long, 
        #      artifacts will arise by going below a 1 cm slice.
        # 
        #  [3] The .h file will have the pulse's attributes defined as static
        #      variables. For example, suppose var_base = 'myPulse'. Then
        #      static double myPulseAmp = {0.0 0.1 0.4 .... };
        #      static double myPulsePhase = {0.0 180.0 90.0 ...};  // in deg.
        #      static double myPulseAmpInt = 42.1; 
        #      static double myPulseAbsInt = ...
        #      and so forth.
        
        max_magn_khz = max_magn * 42.57    # mT * kHz/mT = kHz
        
        magn_arr  = np.abs(pulse_uT_norm)
        phase_arr = np.angle(pulse_uT_norm)
        phase_arr = np.mod(phase_arr, 2*np.pi)
        
        nsteps   = len(magn_arr)
#        duration = nsteps * dwell_us / 1000.0      # in ms, dwell is in usec
        
        abs_integral   = sum(magn_arr)
        power_integral = sum(np.power(magn_arr,2))
        
        # Compute amplitude integral. See "External RF Pulses.doc" 
        # in the "IDEA - Educational\My Notes" folder for more details.
        # amp_integral = flipAngle / (2*pi * dwellTime * maxAmp);

        ref_amp = 0.5  # Amplitude of reference pulse is 0.5 kHz
        amp_integral = (tip_angle/180.0) * (1.0/(dwell_us*0.001)) * (ref_amp/max_magn_khz);

        abs_integral   = DEFAULT_FLOAT_FORMAT % abs_integral
        power_integral = DEFAULT_FLOAT_FORMAT % power_integral
        amp_integral   = DEFAULT_FLOAT_FORMAT % amp_integral

        ref_grad = compute_reference_gradient_siemens(duration_ms, bw_khz, csa=0)
        ref_grad_khz = ref_grad * 42.57 / 1000.0

        # Create external file header
        lines.append('// =------------------=')
        lines.append('// Arbitrary pulse file')
        lines.append('// =------------------=')
        lines.append('//')
        lines.append('// General statistics (as exported): ')
        lines.append('//    Duration:     %.3f ms' % duration_ms)
        lines.append('//    Max B1:       %.3f kHz' % max_magn)
        lines.append('//    Num. Steps:   %d' % nsteps)
        lines.append('//    SAR*:         %.3f' % 0.0)     # bjs FIXME CalcSAR(magn_arr, duration_ms, 'ref'))
        lines.append('//    Ref. Grad:    %.3f mT/m ' % ref_grad)
        lines.append('//    Ref. Grad:    %.3f kHz/mm for 1H ' % ref_grad_khz)
        lines.append('//    Flip Angle:   %.1f deg. ' % tip_angle)
        lines.append('// * - Relative to a 1 ms pi-pulse.')
        lines.append('//')
        lines.append('// User comment: %s' % comment)
        lines.append('//')
        lines.append('// Hint: do not forget to include these libraries and definitions: ')
        lines.append('//    #include "MrServers\\MrMeasSrv\\SeqIF\\libRT\\libRT.h"')
        lines.append('//    #include "MrServers\\MrMeasSrv\\SeqFW\\libSSL\\SSL_local.h"')
        lines.append('//    static sSample %sPulseArray[%d];' % (var_base, nsteps))
        lines.append('')
        lines.append('float %sRefGrad = %0.4f;'  % (var_base, ref_grad))
        lines.append('float %sMinSlice = %.1f;'  % (var_base, float(slice_thick_min)))
        lines.append('float %sMaxSlice = 200.0;' % var_base)
        lines.append('float %sAmpInt = %s; // calculated for 1H'    % (var_base, amp_integral))
        lines.append('float %sPowerInt = %s;'  % (var_base, power_integral))
        lines.append('float %sAbsInt = %s;'    % (var_base, abs_integral))
        lines.append('')
        lines.append('')
        for idx in range(nsteps):
            lines.append('%sPulseArray[%d].flAbs = float(%.5f);    %sPulseArray[%d].flPha = float(%.5f);' % (var_base, idx, magn_arr[idx], var_base, idx, phase_arr[idx]))

        lines.append('')
        lines.append('')
        lines.append('/* ')
        lines.append('')
        lines.append('// The following code is an example of how to make use of this header file ')
        lines.append('// in a pulse sequence. Instructions on where to put the code have single line ')
        lines.append('// comment characters in front of them. Actual code lines are not commented out. ')
        lines.append('// Copy these sections into your file as instructed to activate them. ')
        lines.append('')
        lines.append('// 1. Put header file in the same directory as the pulse sequence. ')
        lines.append('// ')        
        lines.append('// 2. In the "global" part of the sequence (i.e. before fSEQInit, where the pulses & gradients are defined), add: ')
        lines.append('//  ')
        lines.append('// <-------------------------------------- BEGIN --------------------------------> ')
        lines.append('') 
        lines.append('static int %sNumSamples = %s;  ' % (var_base, str(nsteps)))
        lines.append('static sSample %sPulseArray[%s]; ' % (var_base, str(nsteps)))
        lines.append('')
        lines.append('static sRF_PULSE_ARB   %sPulse("%sPulse"); ' % (var_base, var_base))
        lines.append('static sFREQ_PHASE     %sPulseSet( "%sPulseSet" ); ' % (var_base, var_base))
        lines.append('static sFREQ_PHASE     %sPulseNeg( "%sPulseNeg" ); ' % (var_base, var_base))
        lines.append('')
        lines.append('// <-------------------------------------- END --------------------------------> ')
        lines.append('//  ')
        lines.append('// 3. Somewhere in the fSEQPrep function, around where you prepare the pulses, add: ')
        lines.append('// (This may not be "optimal" or the most elegant way. e.g., you may be able to  ')
        lines.append('// define SLR90NumSamples already in SLR90.h. I was not sure about that one because ')
        lines.append('// I wanted my pulse arrays to be global, i.e. static, and did not want to deal with ')
        lines.append('// dynamic memory allocation which I am never sure about) ')
        lines.append('//  ')
        lines.append('// <-------------------------------------- BEGIN --------------------------------> ')
        lines.append('')
        lines.append('#include "%s" ' % output_filename)
        lines.append('')
        lines.append('// Prepare Excitation Pulse ')
        lines.append('%sPulse.setTypeUndefined(); // Whatever. Never understood what this is good for really. ' % var_base)
        lines.append('%sPulse.setSamples(%sNumSamples);   ' % (var_base, var_base))
        lines.append('%sPulse.setDuration( %d );   // in us. Put whatever you want here ' % (var_base, duration_us))
        lines.append('%sPulse.setFlipAngle( %.1f ); // In degrees. Put whatever you want here  ' % (var_base, tip_angle))
        lines.append('%sPulse.setInitialPhase( 0.0 ); ' % var_base)
        lines.append('%sPulse.setThickness( pMrProt->spectroscopy().VoI().thickness() ); // For example, if the pulse works along the "slice" direction   ' % var_base)
        lines.append('')
        lines.append('if (!%sPulse.prepArbitrary(pMrProt,pSeqExpo, %sPulseArray, %sAmpInt)) ' % (var_base, var_base, var_base))
        lines.append('{ ')
        lines.append('    cout << "ERROR: "<< %sPulse.getNLSStatus() << endl; ' % var_base)
        lines.append('    return (false); ')
        lines.append('}; ')
        lines.append('')
        lines.append('%sPulseGrad = (10.0/%sPulse.getThickness())*%sRefGrad*(5120.0/%sPulse.getDuration()); // mT/m ' % (var_base, var_base, var_base, var_base))
        lines.append('lFrequency = specMiddleFreq; ')
        lines.append('lFrequency += (long)( 0.5 + %sPulseGrad * larmorconst * VOI.getSliceShift() ); ' % var_base)
        lines.append('')
        lines.append('%sPulseSet.setFrequency( lFrequency );   ' % var_base)                             
        lines.append('%sPulseNeg.setFrequency( 0L ); ' % var_base)
        lines.append('')
        lines.append('dPhaseExcite = - lFrequency * (360.0/1e6) * %sPulse.getDuration() * 0.5; ' % var_base)
        lines.append('%sPulseSet.setPhase( dPhaseExcite ); ' % var_base)
        lines.append('%sPulseNeg.setPhase( dPhaseExcite ); ' % var_base)
        lines.append('')   
        lines.append('// <-------------------------------------- END --------------------------------> ')
        lines.append('')
        lines.append('*/ ')
        
    # Output ends with a single blank line
    lines.append("")
    
    return lines


def compute_reference_gradient_siemens(duration_ms, bandwidth, csa=0):
    """
    Description: computes the reference gradient for exporting RF files
    to SIEMENS format, assuming the gradient level curGrad is desired.
    
    Theory: the reference gradient is defined as that gradient for which 
    a 1 cm slice is excited for a 5.12 ms pulse. Demanding the product
    Slicethickness * gamma * gradient * duration 
    to be equal in both cases (reference and current), one obtains
    gamma*refGrad*(10 mm)*(5.12 ms) = gamma*curGrad*curThickness*pulse.tp
    However, gamma*curGrad*curThickness = the pulses's bandwidth, pulseBW, so
    
    refGrad = pulseBW*pulse.tp / (gamma*(10 mm)*(5.12 ms))

    In general, the formula is,
                               (pulse_duration[ms]*pulse_bandwidth[kHz])
     Ref_grad [mT/m] =              --------------------------------
                       (Gyr[kHz/mT] * Ref_slice_thickness[m]* Ref_pulse_duration[ms])
    
    Input Variables
    
    Variables Name  Units  Description
    ------------------------------------
    duration_ms        ms     Duration of pulse
    bandwidth       kHz    Bandwidth of current pulse
    csa             kHz    Chemical shift artifact "immunity" - see below.
                           Optional, set to 0 if not present.
    
    Output Variables
    
    Variables Name  Units  Description
    ------------------------------------
    ref_grad        mT/m   Reference gradient
    
    Chemical Shift Artifact immunity: 
    Since different chemical shifts shift the excitation region, it follows
    the if we want to excite a range [-x,+x], we will not actually excite
    that range for any offset other than 0 if we calibrate our gradient
    for 0 offset. However, we CAN calibrate our gradient for 0 offset BUT
    excite a larger range [-x-dx, x+dx] such that the pulse will affect
    all chemical shifts equally. This of course comes at the price of
    exciting a larger region which might have unwanted signals. This however
    is good for:
    1. Cases in which there are not external unwanted signals.
    2. For dual-band suppression pulses, one sometimes uses the PASSBAND,
       which is also the VOI, to calibrate the pulse. If we don't want
       any spins in the VOI affected despite their varying chemical shifts
       we can grant them immunity, at the cost of pushing away the
       suppression bands - this works if, e.g., we're interested in killing
       off fat away from the VOI, so we don't care if a bit of signal comes
       from the region close to the VOI.
    To use, set CSA to the range of +-chemical shifts you want to feel the
    pulse. e.g., if you want all spins +-100 Hz from resonance to be affected
    equally within the VOI, set CSA = 0.1.
    
    """
    
    ref_slice_thickness = 0.01   # in meters
    ref_duration = 5.12          # ms
    gyromagnetic_ratio = 42.57   # kHz/milliTesla
    
    ref_grad = ((bandwidth-2*csa)*duration_ms)/(gyromagnetic_ratio*ref_slice_thickness*ref_duration)
    
    return ref_grad


def calc_sar(magn_arr, duration_ms, pulse_reference):
    """
    Calculates the (relative) SAR, in kHz^2*ms, of a pulse.
    The optional input pulseReference is used to normalize the SAR of 
    the input pulse. It can be set to one of two inputs:
    
       1. Another pulse structure.
       2. 'rect' - normalize the SAR of the input pulse by that 
          of a rectangular pulse of equal duration and amplitude 
       3. 'ref' - normalize the SAR of the input pulse by the of
          a 1 ms 0.5 kHz 180 rectangular pulse
    """
    nsteps = len(magn_arr)
    dwellTime = duration_ms/nsteps
    sar = sum(np.power(magn_arr, 2)*dwellTime)
    
#     if pulse_reference == 'rect':
#             pulse_reference = PulseCreateConst(pulse.tp, numel(pulse.RFamp), max(pulse.RFamp), 0)
#     elif pulse_reference == 'ref':
#             pulse_reference = PulseCreateConst(1, 1, 0.5, 0)
#     
#     else:
#         pulse_reference = 1.0
# 
#     nsteps_ref = len(pulse_reference.RFamp)
#     dwell_time_ref = pulse_reference.tp/nsteps_ref
#     sar_reference = sum(pulse_reference.RFamp.^2*dwell_time_ref)
#     sar = sar/sar_reference

    return sar




#  The following comments were taken from the MR-IDEA forum    
#
#  Title:  How to determine reference gradient for external RF pulse
#  Date:   4/29/2007 3:49 AM
#  Poster: Eric Li
#
# In the IDEA manual, it is said that the reference gradient ampitude 
# displayed in pulse tool is the amplitude for 10mm slice thickness.
#
# In the absence of a gradient, an rf pulse will have a bandwidth that's 
# determined by the pulse design. You can choose the boundaries of the 
# bandwidth by some criterion, ie FWHM, or some other metric according to 
# your needs (ie tip > 160 or such). 
#
# Then the gradient for any slice thickness can be calculated from the slice
# thickness and the gyromagnetic ratio: 
#  (bw) Hz/ (slice thickness) M / (gamma) Hz/mT = mT/M (gradient) 
#
# If the display of the excitation or refocusing profile in PulseTool is 
# correct, you can easily verify the reference gradient.
#
#
#
#
#  Title:  External RF-Pulses. How to implement them on the host?
#  Date:   2/25/2008 11:40 AM
#  Poster: Robert Merwa
#
# Hi Robert,
# 
# The simplest way is to rename your extrf.dat file, e.g. into 
# my_extrf_spec.dat, and then just place it into the directory %CustomerSeq%
# on the host (I believe it is C:\MedCom\MriCustomer\seq on VB15, but check 
# on this to make sure by typing %CustomerSeq% in a command shell).
# 
# Then read the RF pulses into your sequence via 
# 
# pSeqLim->setExtSrfFilename( "%CustomerSeq%/my_extrf_spec.dat" );
# 
# Yes, the host has to be rebooted once, before the new RF pulses are recognized.
# 
# Hope this helps,
# Ralf.
# 
# : I would like to work with SLR pulses generated by an external program. I 
# : have already generated the right .pta ASCII file, which was checked through 
# : import the file in the pulsetool-program. Furthermore the generated pulses 
# : had been imported in the extrf.dat database in order to work with them with 
# : the .setFamilyName instruction at ease. 
# : 
# : Now I am wondering how the extrf.dat database is updated on the 
# : scanner-host. I read something about restarting the scanner in order to get 
# : new pulseshapes into the memory. But how is that done? Should I copy the 
# : new version of the extrf.dat file somewhere on the scanner-host? Or should 
# : an external RF-Pulse only be defined and imported via 
# : pSeqLim->setExtSrfFilename, with the .pta file situated at an fixed 
# : location on the scanner?
# : 
# : Since I had generated many different pulses, the possibility to work with 
# : them through the extrf.dat libary would be much more convenient to me.
#
#
#
#  Title:  How to know B1 amplitude of pulse in microTesla ?
#  Date:   9/19/2011 11:32 AM
#  Poster: Su-Yeon Park
#
# 
# I am not sure what the question means. If it is to find out how to convert 
# the transmit voltage seen in the System transmitter page for any pulse:
# 
# 1. Reference voltage corresponds to 500Hz B1 ie: 11.75 uT.
# 2. For any pulse, peak B1 in uT = (11.75/ref voltage)* pulse peak voltage.
# 
# Only caveat is in some sequences, the actual flip angle of the refocussing 
# pulses may not be 180 deg if the maximum transmit voltage allowed for the 
# pulse is reached. 
# 
# Reference voltage corresponds to 180 deg flip for a 1ms rectangle pulse per 
# Siemens documentation. This is where the 11.75 microtesla comes from.
# 
# Pulse Peak voltage can be seen from the protocol page in the System Tab. I 
# believe that getTransmitterVoltage() returns the actual voltage used by a 
# specific pulse. 
# 
# getMaxtraV() is probably the maximum voltage allowed for the RF coil and is 
# probably the same for all pulses..
#
#
#
#  Title:  RF files (bjs edit - where to put them)
#  Date:   9/19/2011 11:32 AM
#  Poster: Michael Gach
#
# Michael,
# 
# From what I have gathered they seem to be put into %MEASCONST% in product
# sequences and %MEASDAT% in your own sequences.
# 
# MEASCONST is MedCom/MriProduct/Measurement
# MEASDAT is MedCom/MriSitedata/Measurement
# 
# You mustn't overwrite the product pulsefiles, so its best not to put files
# into MEASCONST. It might break the product sequences. So you should put it
# into MEASDAT and make sure your code expects to find it there.
# 
# i.e.
# product sequence might use 
# pSeqLim->setExtSrfFilename( "%MEASCONST%/extrf_spec.pls" );
# pSeqLim->setExtSrfFilename( "%MEASDAT%/my_extrf_spec.pls" );
# 
# On the standalone MEASCONST and MEASDAT might be in different places (don't 
# forget to copy the product sequence pulse files over to MEASCONST on your 
# standalone). These can easily be found by looking at the environment variables.
# 
# Also.
# 
# Welcome to the peculiar world of pulsetool. I assume you are building your 
# own pulses. My advice is "expect the unexpected" and after saving the .dat 
# file try opening it again.
# 
# I hope this helps, this does seem to make some kind of sense and it is what 
# I have done for a while now although much of this was home learned.
# 
# Matt
#
#
# : Hi,
# : Can anybody tell me the location where the external RF files (extrf.dat) 
# : for the pulse sequences should go in the scanner for VB15A?
# : Thanks,
#
#  Same thread as above (Michael Gach) but from Eddie Auerbach ...
#
# I put mine in %CustomerSeq%, since I couldn't think of a good reason for 
# them to be in %MEASDAT%. The product pulse database files are in there 
# because syngo maintains the capability of having a different pulse 
# database for each system type (although currently extrf.dat for every 
# system is the same as far as I know). I use the same modified pulse 
# database for all system types, so it just makes life more difficult, 
# particularly in the IDEA environment where %MEASDAT% changes with the 
# system type.
# 
# The other advantages of %CustomerSeq% are that you only have to worry 
# about one directory for all of your sequence-related files, and
# --the big advantage--  \MedCom\MriCustomer can be maintained by the 
# Backup/Restore facility in the service interface. So, if the service 
# engineer reinstalls software and does a proper backup/restore of 
# everything, your sequences will work again right away. (Unfortunately 
# this doesn't help for ICE program binaries.)