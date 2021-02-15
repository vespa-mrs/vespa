"""
 SHIFTOPT is called to optimize the freq shift necessary to put the
 complex data back onto resonance.

   Checked the results from Grigsby data and RefWaterSI.sid data
   Lineshapes (magnitude) stayed the same, while the FID smoothness
   increased by using a two step optimization.  Step one was to locate
   the max peak after fft and shifting.  Step two minimizes the second
   deriv over a narrow shift range (in Hz).

"""

# Python modules

# 3rd party modules
import numpy as np

# Our modules
import vespa.common.minf_parabolic_info as minf
import vespa.common.util.generic_spectral as util_generic_spectral



def optimize_b0_shift(fid_original, acqsw, ideal_width=2.0, search_range=5):
    '''
    This routine take a complex numpy array containing MRS time data (an FID
    data set) and calculates the phase roll (aka frequency shift) needed to
    center the data on the resonance frequency. 
    
    This algorithm is based on the assumption that there is a significant
    residual water signal in the overall signal, such that when this large
    water signal is shifted close to the resonance frequency, then the time
    data will be at its smoothest. Conversely, all other shifts should 
    increase the variability of the time data.
    
    Historical - this algorithm was translated from ShiftOpt code in libidl
    originally written by Brian J. Soher and translated by Jeff Steinberg.
    
    '''

    if len(fid_original.shape) != 1: 
        # must be 1D array
        raise ValueError("Data must be a 1D array")
    
    if not "complex" in str(fid_original.dtype):
        raise TypeError("Data must be complex, not %s" % str(fid_original.dtype))
        
    dim0 = len(fid_original)
        
    fid = fid_original.copy()
    fid = fid/fid[0]
    
    # Determine initial shift fft and max peak in freq domain    
    dat = np.fft.fft(fid) / len(fid)
    
    b0shft = (np.abs(dat)).argmax()-search_range
    
    if b0shft > dim0/2: 
        b0shft -= dim0
    
    # Swapped sign (to -) on indx to correct measured shift
    # b0shft sign convention matches FITT with (-) sign
    
    initial_shift = -b0shft                    # initial FREQ shift
        
    #-------------------------------------------------------
    # Filter results and apply to the envelopes
    initial_roll = np.exp(complex(0,1) * 2.0 * np.pi * np.arange(dim0) * initial_shift / acqsw)    
    fid = fid * initial_roll    
    

    #----------------------------------------------------------------
    # Tested Oct 9, 2001 - Correlation method worked a lot
    # better at tweaking the data than did the 2nd derivative
    # minimizing method.
    #
    # Also, tested if a fixed width input to the shiftOpt2 works
    # as well as a variable width input, based on a data measure.
    # Sure enough, the results were quite similar and in some cases
    # a bit better with a fixed narrow width (long decay FID)
    #----------------------------------------------------------------

    #-------------------------------------------------------
    # Now tweak search across a smaller region
    
    delta = optimize_b0_shift_correlation(fid.copy(), -search_range*2, 0.0, acqsw, ideal_width)
    
    final_roll  = np.exp(1j * 2.0 * np.pi * np.arange(dim0) * delta / acqsw)
    final_shift = initial_shift + delta
    
    fid = fid * final_roll
    fid = fid / fid[0]
    
    return fid, final_shift



def shift_correlation_function(shift, info={}):
    '''
    Returns negative correlation because we are using an optimization call
    that minimizes things
    
    info should be a dictionary with keys of 
    "data"  - a numpy array 
    "acqsw" - float, acquisition sweep width
    "expo"  = a numpy array
    
    '''
    if not info: 
        return 0
    
    sw   = info['acqsw']
    data = info['data'].copy()
    
    npts = len(data)
    data = np.exp(1j*2.0*np.pi*np.arange(npts)*shift/sw) * data
    cc   = np.corrcoef(data, info['expo'])
#    print 'shift = ', shift, '   corr = ',-cc[1,0] 
    
    return -cc[1,0]




def optimize_b0_shift_correlation(data_original, minhz, maxhz, acqsw, width=5.0):
    '''
    This routine take a complex numpy array containing MRS time data (an FID
    data set) and calculates the phase roll (aka frequency shift) needed to
    best match the frequency shifted real part of the data to an ideal
    Gaussian decay envelope matched to the linewidth of the data. 
    
    data:  complex array, FID data to be shifted
    minhz: float, lower bound in Hz
    maxhz: float, upper bound in Hz
    acqsw: float, acquisition sweep width in Hz

    returns the float value for shifting the FID in Hz
    
    '''
    
    data = data_original.copy()
    
    npts = len(data)
    xx = np.arange(npts) / acqsw
    expo = util_generic_spectral.apodize(xx, width, 'Gaussian')   # from 4DFT
    info = {'data'  : (data/data[0]).real,
            'acqsw' : acqsw,
            'expo'  : expo }
    
    # Call parabolic interpolation, Brent's method
    # 1-d minimization routine
    
    shift_a = np.where(minhz < maxhz, minhz, maxhz)        # lower bound
    shift_c = np.where(minhz > maxhz, minhz, maxhz)        # upper bound
    shift_b = np.mean([shift_c,shift_a]) + 0.1             # "some" point in the middle   
    minshift, _ = minf.minf_parabolic_info(shift_a, 
                                           shift_b, 
                                           shift_c, 
                                           shift_correlation_function, 
                                           info)
        
    return minshift



def _test():
    import vespa.common.minf_parabolic_info as minf
#    import pylab

    xa = -0.1
    xb = -5
    xc = -20
    func = shift_correlation_function
    
    sw    = 3906.25
    dim0  = 4096
    width = 6
    range = 15
    global_shift = 3
    wat_shift = 0.0
    cho_shift = 95.4
    cr_shift  = 108.8
    naa_shift = 172.8
    wat_area  = 50.0
    cho_area  = 0.6
    cr_area   = 0.8
    naa_area  = 1.0
    
    xx = np.arange(dim0) / sw

    water = wat_area * np.exp(1j*2.0*np.pi*xx*wat_shift)
    naa   = naa_area * np.exp(1j*2.0*np.pi*xx*naa_shift)
    cho   = cho_area * np.exp(1j*2.0*np.pi*xx*cho_shift)
    cr    = cr_area  * np.exp(1j*2.0*np.pi*xx*cr_shift)
    
    roll = np.exp(1j*2.0*np.pi*xx*global_shift)
    
    data0 = water + naa + cr + cho

    lshape = util_generic_spectral.apodize(xx, width, 'Gaussian') # from 4DFT
    data1  = data0 * lshape 
    data   = data1 * roll
    
#    pylab.plot(data1, 'b')
#    pylab.plot(data, 'r')
#    pylab.plot(lshape, 'g')
#    pylab.show()
    
    info = { "data"  : data, 
             "acqsw" : sw, 
             "expo"  : lshape }

    minshift, iter = minf.minf_parabolic_info(xa, xb, xc, func, info)
                            
    print('Direct Minf optimize = ', minshift, '   iter = ', iter)

    fid, shift = optimize_b0_shift(data, sw, ideal_width=width,
                                             search_range=range)
    
    print('optimize_b0_shift = ', shift)
    

if __name__ == '__main__':
    _test()
