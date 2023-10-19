# Python modules

# 3rd party modules
import numpy as np

# Our modules
import vespa.common.minf_parabolic_info as minf
import vespa.common.constants as common_constants

from vespa.analysis.algos.b0_correction import b0_correction

DTOR = common_constants.DEGREES_TO_RADIANS
RTOD = common_constants.RADIANS_TO_DEGREES



def optimize_phase0_correlation(data, modfn, pts, nlag, cdeg):
    """
    Returns the zero order phase in deg at which the real part is maxed

    INPUT:
     data:  array of complex data to be phased
     modfn: string array of model function to be used (obsolete)
     pts:
     nlag:  correlation lag
     cdeg:

    KEYWORDS:
     csum: optimal correlation value for optimization

    """
    info = {'dat'     : data,  
            'ref'     : modfn,
            'pts'     : pts,
            'nlag'    : nlag, 
            'cdeg'    : cdeg,
            'csum'    : 0.0        }

    # Call parabolic interpolation, Brent's method
    # 1-d minimization routine

    phase_a = -1*np.pi   #0.0       ; lower bound
    phase_b = np.pi*0.5  #np.pi     ; "some" point in the middle
    phase_c = 2*np.pi    #np.pi*2   ; upper bound

    phaseat, maxit = minf.minf_parabolic_info( phase_a, 
                                               phase_b, 
                                               phase_c,
                                               _corr_func_phase0,
                                               info)
    ph = phaseat * RTOD
    if ph > 360.0:
        ph = ph - 360.0
    if ph < -360:
        ph = ph + 360.0

    return ph, info['csum']


def optimize_phase1_correlation(data, modfn, pts, pivot, nlag, cdeg):
    """
    Returns the first order phase in deg at which the real part is maxed

    INPUT:
      data:     array of complex data to be phased
      modfn:
      pts:
      pivot:
      nlag:
      cdeg:

    """
    info = {'dat'     : data, 
            'piv'     : pivot, 
            'ref'     : modfn,
            'pts'     : pts,
            'nlag'    : nlag, 
            'cdeg'    : cdeg,
            'csum'    : 0.0        }

    # Call parabolic interpolation, Brent's method
    # 1-d minimization routine

    phase_a = -1500.0   # lower bound
    phase_b =  0.0      # "some" point in the middle
    phase_c =  1500     # upper bound

    phaseat, maxit = minf.minf_parabolic_info( phase_a, 
                                               phase_b, 
                                               phase_c,
                                               _corr_func_phase1,
                                               info)
    return phaseat, info['csum']



def optimize_phase0_integration(data, weight):
    """
    Optimization call for real integral zero order auto phase.
    Returns the float zero order phase in degrees at which the 
    real part is maxed

    INPUT:
      Data: complex array to be phased
      Weight: weight to be applied to summed real part

    """

    info = { 'dat' : data, 'wt' : weight }

    # Call parabolic interpolation, Brent's method
    # 1-d minimization routine

    phase_a = np.pi * (-1.0)  # lower bound
    phase_b = np.pi * (0.001) # "some" point in the middle
    phase_c = np.pi * (2.0)   # upper bound

    phaseat, maxit = minf.minf_parabolic_info( phase_a, 
                                               phase_b, 
                                               phase_c, 
                                               _integ_func_phase0,
                                               info)
    ph = phaseat * RTOD
    if ph > 360.0:
        ph = ph - 360.0
    if ph < -360:
        ph = ph + 360.0

    return ph



def optimize_phase1_integration( data, weight, pivot):
    """
    Optimization call for real integral first order auto phase.
    Returns float first order phase in degrees

    INPUT:
      data:   complex array to be phased
      weight: weight to be applied to summed real part
      pivot:  first order phase pivot in ppm

    """

    info = { 'dat' : data, 'wt' : weight, 'piv' : pivot }

    # Call parabolic interpolation, Brent's method
    # 1-d minimization routine

    phase_a = -2000.0 # lower bound
    phase_b =  0.0    # "some" point in the middle
    phase_c =  2000.0 # upper bound

    phaseat, maxit = minf.minf_parabolic_info( phase_a, 
                                               phase_b, 
                                               phase_c, 
                                               _integ_func_phase1,
                                               info)
    return phaseat
    
    
    
###############################################################################    
#    
# Penalty functions for the Integral and Correlation optimization above
#
###############################################################################

def _corr_func_phase0(rad0, info):
    """
    Optimization function used in correlation auto zero order phase.
    Returns the sum of the cross correlation for the real and
    imaginary data bits with respective bits from an ideal
    reference function

    The negative sum of correlations is returned since we are
    minimizing in the function, but want to maximize the ccval

    INPUT:
        rad0: current phase in radians
        info: control structure, see _optimize_phase0_correlation for definition

    """
    phase = np.exp(1j*rad0)
    dat   = info['dat'].copy() * phase

    istr = info['pts'][0]
    iend = info['pts'][1]
    datt = dat[istr:iend].copy()
    reff = info['ref'][istr:iend].copy()  

    shft, csum = b0_correction(datt, reff, nlag=info['nlag'], cdeg=info['cdeg'])

    info['csum'] = csum

    return  -csum


def _corr_func_phase1( deg1, info):
    """
    Optimization function used in correlation auto first order phase.
    Returns the negative of the area since we are minimizing
    in the function, but want to maximize the real spectrum

    deg1: current phase in degrees
    info: control structure, see _optimize_phase1_correlation for definition

    """
    dim0 = len(info['dat'])
    dat  = info['dat'].copy()

    phase = deg1 * DTOR * (np.arange(dim0)-info['piv'])/dim0
    phase = np.exp(1j*phase)

    dat = phase * dat

    istr = info['pts'][2]
    iend = info['pts'][3]

    datt = dat[istr:iend].copy()
    reff = info['ref'][istr:iend].copy()

    shft, csum = b0_correction(datt, reff, nlag=info['nlag'], cdeg=info['cdeg'])

    info['csum'] = csum

    return  -csum

    
def _integ_func_phase0(rad0, info):
    """
    Optimization function for real integral zero order auto phase.
    Returns the negative of the area since we are minimizing
    in the function, but want to maximize the real spectrum

    INPUT:
      rad0: optimization parameter in radians
      info: control structure

    """
    phase = np.exp(1j * rad0) 
    dat   = phase * info['dat'].copy()
    res   = -1* np.sum(dat.real * info['wt'])

    return  res


def _integ_func_phase1(deg1, info):
    """
    Optimization function for real integral first order auto phase.
    Returns the negative of the area since we are minimizing
    in the function, but want to maximize the real spectrum

    INPUT:
      deg1: optimization variable in degrees
      info: control structure

    """
    dat  = info['dat'].copy()
    dim0 = len(dat)

    phase = deg1 * DTOR * (np.arange(dim0)-info['piv'])/dim0
    phase = np.exp(1j*phase)
    dat   = phase * dat
    res   = -1* np.sum(dat.real * info['wt'])

    return  res    