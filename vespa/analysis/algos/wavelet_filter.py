from __future__ import print_function

# Python modules

# 3rd party modules
import numpy as np

# Our modules
import vespa.common.wavelet_1d as wave1d

try:
    import pywt
    flag1 = 'pywavelets'
except:
    flag1 = 'local'



def wavelet_filter_pywt_original(data, thresh, scale, dyadmin):
    '''
    The wavelet filter uses the library PyWavelets (pywt).

    From the website (http://www.pybytes.com/pywavelets/):
    PyWavelets is a free Open Source wavelet transform software for Python
    programming language. It is written in Python, Pyrex/Cython and C for
    a mix of easy and powerful high-level interface and the best performance.

    The implementation below is meant to be used as a method for estimating the
    signals that compose a baseline with broad signals that change slowly. It is
    performed by using the discrete wavelet transform (DWT) on the data to be
    estimated. Then, some of the smaller scale dyad coefficients are set to zero,
    and the inverse DWT is performed.

    Args:

    data (array, float): dependent data to be filtered

    thresh (float):  value in points for the minimum scale wavelet to be
            left in the result based on the linewidth value passed in from
            the application

        scale (float): multiplier of thresh, to make result have even broader terms

        dyadmin (float): hard value in points of the smallest scale object allowed
            regardless of thresh and scale values.

    Result:

        res (array, float): the filtered waveform (ie. data)

        coeffs (array, float): the full set coefficients from the discrete
            wavelet transform

        filtered_coeffs (array, float): the set of coefficients left after
            removing small scale dyads

        base_level (int): the calculated value for dyad not to be removed

    Algorithm:

        We use wavedec because it returns a list of arrays that correspond to
        dyad coeffs. So to get rid of a whole dyad of coeffs, we just have to
        multiply the numpy array by zero. That is, gtting rid of one of the
        arrays in the coeffs list is the same as zeroing a whole 2^N category.

        To have some room to filter at different levels of broadness, I found
        that the pywt.dwt_max_level() method did not typically set the level to
        a very big dyad range since it returns a "useful" level of decomposition
        not the complete decomposition. Empirically, I found that for arrays of
        32 - 2048 points, if I leave 8 coeffs in the approximation coeffs array,
        and 8 more in the first details coefs array, then I have a workable
        number of decimation levels to work with.

        Thus the base_level calc below with the magic number 32. For 128 points,
        128 / 32 = 4 which gives us dyad arrays of 8, 16, 32, 64 to work with.
        and this scales to bigger data arrays too, with 128, 256 etc dyad levels.

        We use line width as a criterion for removing wavelet coefficients
        that have features on smaller scales. This results in a filtered signal
        that with broad components that change broadly.

        Here we calculate which of two values will determine the broadness of
        the reconstructed signal, 1) peak LW times a scaling multiplier, or
        2) a hard minimum scale value. We calc each and use the bigger one.

    Usage:

      res = wavelet_filter_pywt_original(data, thresh, scale, dyadmin)

    ''' 
    npts   = np.size(data)
    
    # Use the Coiflet 3 wavelet (18 data points), for historical reasons 
    wavelet = pywt.Wavelet('coif3')   
    
    # per - periodization - is like periodic-padding but gives the smallest 
    #       possible number of decomposition coefficients. Wavedec and waverec 
    #       methods must be performed with the same mode.
    _mode = 'per'

    # hard set level value empirically - see note above
    max_level  = int(np.log(len(data))/np.log(2))    
    base_level = int(len(data) / 32) 
    base_level = base_level if base_level < max_level else max_level
    
    coeffs  = pywt.wavedec(data, wavelet, _mode, level=base_level)

    # Calculate criterion for removing wavelet coefficients with features on 
    # smaller scales using either line width or a hard minimum scale value. 
    scale_mult = 1 if scale < 1 else scale
    lw_thresh  = np.where(thresh >  1, thresh, 1)  # this is in points
    lw_thresh *= scale_mult
    hard_min   = np.where(dyadmin > 1, dyadmin, 1) 
    width = np.where(hard_min > lw_thresh, hard_min, lw_thresh)
    
    # now we need to apply the calculated filtering width to the coeffs
    # and decide which ones to use in the reconstruction.
    filtered_coeffs = [coeffs[0]]
    for coeff in coeffs[1:]:
        width_of_each_dyad = npts/len(coeff)
        if width_of_each_dyad > width:
            filtered_coeffs.append(coeff.copy())
        else:
            filtered_coeffs.append(coeff.copy()*0)
    
    res = pywt.waverec(filtered_coeffs, wavelet, _mode)

    return res, coeffs, filtered_coeffs, base_level



def wavelet_filter_pywt_new(data, thresh, scale, dyadmin):
    '''
    This version of wavelet filter uses the library PyWavelets (pywt).

    From the website (http://www.pybytes.com/pywavelets/):
    PyWavelets is a free Open Source wavelet transform software for Python
    programming language. It is written in Python, Pyrex/Cython and C for
    a mix of easy and powerful high-level interface and the best performance.

    The implementation below is meant to be used as a method for estimating the
    signals that compose a baseline with broad signals that change slowly. It is
    performed by using the discrete wavelet transform (DWT) on the data to be
    estimated. Then, some of the smaller scale dyad coefficients are set to zero,
    and the inverse DWT is performed.

    Args:

        data (array,float):  dependent data to be filtered

        thresh (float):  value in points for the minimum scale wavelet to be
            left in the result based on the linewidth value passed in from
            the application

        scale (float): multiplier of thresh, to make result have even broader terms

        dyadmin (float): hard value in points of the smallest scale object allowed
            regardless of thresh and scale values.

    Result:

        res (array, float): the filtered waveform (ie. data)

        coeffs (array, float): the full set coefficients from the discrete
            wavelet transform

        filtered_coeffs (array, float): the set of coefficients left after
            removing small scale dyads

        base_level (int): the calculated value for dyad not to be removed

    Usage:

      res = wavelet_filter_pywt_new(data, thresh, scale, dyadmin)

    '''
    # Use the Coiflet 3 wavelet (18 data points), for historical reasons
    wavelet = 'coif3' 
    wavelet = pywt.Wavelet('coif3')
    _mode   = 'per'

    # hard set level value empirically - see note above    
    base_level = int((np.log(len(data)) / np.log(2)) ) 
    
    coeffs2 = pywt.wavedec(data, wavelet, _mode, base_level)

    coeffs = np.concatenate(coeffs2)
    

    # Calculate criterion for removing wavelet coefficients with features on 
    # smaller scales using either line width or a hard minimum scale value. 
    scale_mult = 1 if scale  < 1 else scale   # float, multiplier
    lw_thresh  = 1 if thresh < 1 else thresh  # in points
    lw_thresh *= scale_mult
    
    hard_thresh = 2 if dyadmin < 2 else dyadmin  # in points
    
    width = hard_thresh if hard_thresh > lw_thresh else lw_thresh
    width = np.floor(np.log(width)/np.log(2))
    
    # Convert base2 dyad widths into index into the DWT result
    # - if width=1, then index==base_level, which would zero out all coeffs
    #   in the DWT results
    # - added (-1) to this line to always leave at least 2 coeffs to be 
    #   transformed back
    # - empirically, this also keeps result a bit smoother, too
    
    indx = base_level - width  
    indx = int(2**indx) 
    
    # now we need to apply the calculated filtering width to the coeffs
    # and decide which ones to use in the reconstruction.
    filtered_coeffs = coeffs.copy()
    filtered_coeffs[indx:] = 0.0 
    
    # restack flattened coeffs so we can apply waverec() 
    tmp = [np.array(filtered_coeffs[0:1]),]
    for i in range(0,base_level):
        tmp.append(np.array(filtered_coeffs[2**(i):2**(i+1)]))
        
    res = pywt.waverec(tmp, wavelet, _mode)

    return res, coeffs, filtered_coeffs, base_level


def wavelet_filter_local(data, thresh, scale, dyadmin):
    '''
    This wavelet filter uses our native wavelets methods.

    The implementation below is meant to be used as a method for estimating the
    signals that compose a baseline with broad signals that change slowly. It
    is performed by using the discrete wavelet transform (DWT) on the data to
    be estimated. Then, some of the smaller scale dyad coefficients are set to
    zero, and the inverse DWT is performed.

    Args:

        data (array,float):  dependent data to be filtered

        thresh (float):  value in points for the minimum scale wavelet to be
            left in the result based on the linewidth value passed in from
            the application

        scale (float): multiplier of thresh, to make result have even broader terms

        dyadmin (float): hard value in points of the smallest scale object allowed
            regardless of thresh and scale values.

    Result:

        res (array, float): the filtered waveform (ie. data)

        coeffs (array, float): the full set coefficients from the discrete
            wavelet transform

        filtered_coeffs (array, float): the set of coefficients left after
            removing small scale dyads

        base_level (int): the calculated value for dyad not to be removed

    Usage:

      res = wavelet_filter_local(data, thresh, scale, dyadmin)

    ''' 
    # Use the Coiflet 3 wavelet (18 data points), for historical reasons 
    wavelet = 'coif3' 

    # hard set level value empirically - see note above    
    base_level = int((np.log(len(data)) / np.log(2)) ) 
    
    coeffs, widths, alerts = wave1d.dwt(wavelet, data, base_level)

    # Calculate criterion for removing wavelet coefficients with features on 
    # smaller scales using either line width or a hard minimum scale value. 
    scale_mult = 1 if scale  < 1 else scale   # float, multiplier
    lw_thresh  = 1 if thresh < 1 else thresh  # in points
    lw_thresh *= scale_mult
    
    hard_thresh = 2 if dyadmin < 2 else dyadmin  # in points
    
    width = hard_thresh if hard_thresh > lw_thresh else lw_thresh
    width = np.floor(np.log(width)/np.log(2))
    
    # Convert base2 dyad widths into index into the DWT result
    # - if width=1, then index==base_level, which would zero out all coeffs
    #   in the DWT results
    # - user needs to ensure that this does not happen
    # - if we subtract 1 from indx to prevent this from happening we are no
    #   longer working quite the same way as wavepy dwt/idwt calls
    
    indx = base_level - width  
    indx = int(2**indx)
    
    # now we need to apply the calculated filtering width to the coeffs
    # and decide which ones to use in the reconstruction.
    filtered_coeffs = coeffs.copy()
    filtered_coeffs[indx:] = 0.0 
    
    res = wave1d.idwt(wavelet, filtered_coeffs, widths, alerts)

    return res, coeffs, filtered_coeffs, base_level


def wavelet_filter(data, thresh, scale, dyadmin):
    # the extra return values below are good for error checking 
    # in the test code below. Here we just return result

    if flag1 == 'local':
        result, coefs, filts, levels = wavelet_filter_local(data, thresh, scale, dyadmin)

    else:
        result, coefs, filts, levels = wavelet_filter_pywt_new(data, thresh, scale, dyadmin)
        #result, coefs, filts, levels = wavelet_filter_pywt_original(data, thresh, scale, dyadmin)

    return result
    

#------------------------------------------------------------------------------
# Test Code 

data0 = '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'

def _test(data):
    """
    Ran this 2020-02-17 and got the results below. Seems like the different 
    versions of pywt work the same for all intents and purposes.
    
    NB. The difference between the old and new versions of wavelet_filter_pywt
        is how the hard limit for the 'level' parameter is set.  The new way
        is a bit more flexible for leaving enough coeffs in the 'approx' arrays
        as well as in the 'detail' arrays for us to zero as needed. 

    python 3 - pywt version 1.1.1
    
    max diff12 = 0.3857474631292188   mean coefs 1vs2 = -0.1602445608675241
    max diff13 = 0.3857474631281934   mean coefs 1vs3 = -0.1544178762518972
    max diff23 = 1.4902801215299633e-12   mean coefs 2vs3 = -0.12480753243166957
    
    python 2 - pywt version 0.3.0
    
    max diff12 = 0.38574746312822544   mean coefs 1vs2 = -0.16024456086773445
    max diff13 = 0.38574746312828   mean coefs 1vs3 = -0.154417876251769
    max diff23 = 7.038813976123492e-14   mean coefs 2vs3 = -0.12480753243188117

    """
    import pylab


    thresh  = 1      # should be in points
    scale   = 1      # multiplier for thresh
    dyadmin = 8      # in points
    
    result1, coef1, filt1, base_level1 = wavelet_filter_pywt_original( data64, thresh, scale, dyadmin)
    result2, coef2, filt2, base_level2 = wavelet_filter_pywt_new(data64, thresh, scale, dyadmin)
    result3, coef3, filt3, base_level3 = wavelet_filter_local(data64, thresh, scale, dyadmin)

    diff12 = result1 - result2
    diff13 = result1 - result3 
    diff23 = result2 - result3

    print("max diff12 = " + str(np.max(diff12))+'   mean coefs 1vs2 = '+str(np.mean(np.concatenate((coef1[0],coef2) ))))
    print("max diff13 = " + str(np.max(diff13))+'   mean coefs 1vs3 = '+str(np.mean(np.concatenate((coef1[0],coef3) ))))
    print("max diff23 = " + str(np.max(diff23))+'   mean coefs 2vs3 = '+str(np.mean(np.concatenate((coef2,coef3) ))))

    pylab.plot(data64)
    pylab.plot(result1, 'b', linewidth=2)
    pylab.plot(result2, 'r', linewidth=3)
    pylab.plot(result3, 'g', linewidth=2)
    pylab.show()
    
    bob = 10
    bob += 1
    

if __name__ == '__main__':

    import io, zlib, base64

    data512 = base64.b64decode(data0)
    data512 = zlib.decompress(data512)
    buf = io.BytesIO(data512)
    data512 = np.load(buf)
    data64 = data512[:64]

    _test(data64)
    
#    import cProfile
#    cProfile.run('_test2()')
 