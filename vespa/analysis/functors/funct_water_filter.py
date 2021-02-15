# Python imports


# 3rd party imports
import numpy as np
import scipy

# Vespa imports





WATFILT_MENU_ITEMS = ['None', 
                      'FIR - water filter',
                      'Hamming - water filter',
                      'SVD - water filter']

#----------------------------------------------------------
# FIR (Finite Impulse Response) and Hamming constants
    
FIR_LENGTH_MIN = 1
FIR_LENGTH_MAX = 99
FIR_LENGTH_DEFAULT = 11

FIR_HALF_WIDTH_MIN = 0
FIR_HALF_WIDTH_MAX = 500
FIR_HALF_WIDTH_DEFAULT = 50
FIR_HALF_WIDTH_STEP = 2

FIR_RIPPLE_MIN = 0
FIR_RIPPLE_MAX = 500
FIR_RIPPLE_DEFAULT = 40
FIR_RIPPLE_STEP = 2

FIR_EXTRAPOLATION_ALL = ['None', 'Linear', 'AR Model']
FIR_EXTRAPOLATION_DEFAULT = FIR_EXTRAPOLATION_ALL[0]

FIR_EXTRAPOLATION_POINTS_MIN = 2
FIR_EXTRAPOLATION_POINTS_MAX = 1000
FIR_EXTRAPOLATION_POINTS_DEFAULT = 10

HAM_LENGTH_MIN = 1
HAM_LENGTH_MAX = 99
HAM_LENGTH_DEFAULT = 11

HAM_EXTRAPOLATION_ALL = ['None', 'Linear', 'AR Model']
HAM_EXTRAPOLATION_DEFAULT = HAM_EXTRAPOLATION_ALL[0]

HAM_EXTRAPOLATION_POINTS_MIN = 2
HAM_EXTRAPOLATION_POINTS_MAX = 1000
HAM_EXTRAPOLATION_POINTS_DEFAULT = 10 

SVD_N_DATA_POINTS = 1024
SVD_N_SINGULAR_VALUES = 20
 
 

def water_filter_fir(chain):

    filter_length           = chain._block.set.fir_length
    half_width              = chain._block.set.fir_half_width
    ripple                  = chain._block.set.fir_ripple
    extrapolation_method    = chain._block.set.fir_extrapolation_method
    extrapolation_point_count = chain._block.set.fir_extrapolation_point_count

    water_time = 0

    # Get filter cutoff in Hz in terms of Nyquist frequency,
    # 1/2T, where T is the time between data samples
    # Always do as lowpass, to simplify in removal part
    cutoff = 2 * half_width / chain.sw

    # FIR ripple is in db, and is converted to the approximate width
    # of the transition region (normalized so that 1 corresonds to pi)
    # for use in kaiser FIR filter design.
    width = (ripple - 8) / (2.285 * filter_length)

    # Approximate the digital_filter function with firwin
    # Differences are as high as 5% at width=0, and ~0% for larger widths
    h2o_filter = firwin(filter_length, cutoff, width)

    # normalize
    h2o_filter = h2o_filter / sum(h2o_filter)

    # convolution function is always high-pass, so get water function
    water_time = np.convolve(chain.data, h2o_filter, 1)

    if extrapolation_method == 'Linear':
        k2 = (filter_length - 1) // 2
        x  = np.arange(k2)

        # Get slope of linear fit
        slope = np.polyfit(x, water_time[k2:k2+k2], 1)[0]

        for j in range(int(k2)):
            water_time[j] = water_time[k2] + ((k2 - j) * slope)

    elif extrapolation_method == 'AR Model':
        k = (filter_length - 1) // 2
        p = extrapolation_point_count

        rwdata = time_series_forecast(water_time[k:].real, p, k, backcast=True)
        iwdata = time_series_forecast(water_time[k:].imag, p, k, backcast=True)

        # will return NaN if constant function !
        if not np.isfinite(rwdata[0]):
            rwdata = np.zeros(k,float)+water_time[k].real
        if not np.isfinite(iwdata[0]):
            iwdata = np.zeros(k,float)+water_time[k].imag

        water_time[0:k] = rwdata + 1j * iwdata

    # Return the time data with the water estimate subtracted

    chain.data = chain.data - water_time     
  

def water_filter_hamming(chain):

    filter_length           = chain._block.set.ham_length
    extrapolation_method    = chain._block.set.ham_extrapolation_method
    extrapolation_point_count = chain._block.set.ham_extrapolation_point_count

    water_time = 0

    # get Hamming filter
    hamming_filter = np.hamming(filter_length + 1)[0:filter_length]

    # normalize
    hamming_filter = hamming_filter / sum(hamming_filter)

    # Convolution function is always high-pass, so get water function
    water_time = np.convolve(chain.data, hamming_filter, 1)

    if extrapolation_method == 'Linear':
        k2 = (filter_length - 1) // 2
        x  = np.arange(k2)

        # Get slope of linear fit
        slope = np.polyfit(x, water_time[k2:k2+k2], 1)[0]

        for j in range(int(k2)):
            water_time[j] = water_time[k2] + ((k2 - j) * slope)

    elif extrapolation_method == 'AR Model':
        k = (filter_length - 1) // 2
        p = extrapolation_point_count

        rwdata = time_series_forecast(water_time[k:].real, p, k, backcast=True)
        iwdata = time_series_forecast(water_time[k:].imag, p, k, backcast=True)

        # will return NaN if constant function !
        if not np.isfinite(rwdata[0]):
            rwdata = np.zeros(k,float)+water_time[k].real
            #rwdata = replicate(water_time[k].real, k)
        if not np.isfinite(iwdata[0]):
            iwdata = np.zeros(k,float)+water_time[k].imag
            #iwdata = replicate(water_time[k].imag, k)

        water_time[0:k] = rwdata + 1j * iwdata

    # Return the time data with the water estimate subtracted

    chain.data = chain.data - water_time      


def water_filter_hlsvd(chain):
    """
    The chain object contains both the array of FIDs that HLSVD generates
    and the array of summed FIDs that represent only the lines selected by
    checking the boxes in the list on the HLSVD tab. These are the two 
    choices that can be applied in this filter. 

    If the user has checked the Apply Threshold box on the Spectral tab's
    water filter, then the array of FIDs is parsed to create a water 
    filter from all lines whose frequency is >= the threshold value. The 
    use of this setting does NOT affect the check boxes on the HLSVD tab.

    If the Apply Threshold box is not checked, then the array of summed
    FIDs set from the HLSVD tab is used as-is. In this case the value in 
    the water threshold field is ignored.

    """
    # the HLSVD results were calculated aligned with the original raw time
    # data, no frequency shift applied. As we apply any frequency shift to
    # the raw data, we must also shift the HLSVD fids. However, if we use
    # the Spectral tab's cutoff to determine which HLSVD fids to remove, 
    # then we need to apply the threshold to HLSVD frequencies that have
    # had the frequency shift added to them. And in the end, the HLSVD fids
    # need to have a phase roll applied to line them up with the raw data. 

    chain.data = chain.data - chain.svd_fids_checked    



def do_water_filter_processing(chain):
    
    water_filter_method = chain._block.set.water_filter_method
    
    if water_filter_method == 'FIR - water filter':
        water_filter_fir(chain)
        
    elif water_filter_method == 'Hamming - water filter':
        water_filter_hamming(chain)
        
    elif water_filter_method == 'SVD - water filter':
        water_filter_hlsvd(chain)
        


#------------------------------------------------------------------------------
# Helper Functions

def time_series_forecast(x, p, nvalues, backcast=False, reflect=False):
    """
    This function computes future or past values of a stationary time-
    series (X) using a Pth order autoregressive model. The result is an
    nvalues-element numpy array whose type is identical to X.

    This function uses the last P elements [Xn-1, Xn-2, ... , Xn-p]
    of the time-series [x0, x1, ... , xn-1] to compute the forecast.
    More coefficients correspond to more past time-series data used
    to make the forecast.

    x   An n-element numpy array of floats containing time-series data

    p   A scalar that specifies the number of actual time-series values 
        to be used in the forecast. In general, a larger number of values 
        results in a more accurate result.

    nvalues    A scalar that specifies the number of future or past 
               values to be computed.

    backcast   If set then "backcasts" (backward-forecasts)are computed

    Based on ...
    The Analysis of Time Series, An Introduction (Fourth Edition)
    Chapman and Hall  ISBN 0-412-31820-2

    """
    nvalues = int(nvalues)

    if nvalues <= 0:
        raise ValueError("nvalues must be a scalar > 0")

    nx = len(x)

    if p<2 or p>(nx-1):
        raise ValueError("p must be a scalar in  [2, len(x)-1]")

    # reverse time-series for backcasting.
    if backcast:
        x = x[::-1]

    # last p elements of time-series.
    data = (x[nx-int(p):])[::-1]

    fcast = np.zeros(nvalues, float)

    # compute coeffs
    arcoeff = time_series_coef(x, int(p))

    for j in np.arange(nvalues):
        data = np.concatenate((np.array([(data * arcoeff).sum()]), data[0:int(p)-1]))
        fcast[j] = data[0]

    if backcast:
        return fcast[::-1]
    else:
        return fcast


def time_series_coef(x, p):
    """
    This function computes the coefficients used in a Pth order
    autoregressive time-series forecasting/backcasting model. The
    result is a P-element numpy array 

    Used to compute the coefficients of the Pth order autoregressive 
    model used in time-series forecasting/backcasting.
    arcoef = arcoef[0, 1, ... , p-1]

      x:    An n-element vector of type float or double containing time-
            series samples.

      p:    A scalar of type integer or long integer that specifies the
            number of coefficients to be computed.

    mse:    calculates the mean square error of the Pth order 
            autoregressive model

    Based on ...
    
    The Analysis of Time Series, An Introduction (Fourth Edition)
    Chapman and Hall
    ISBN 0-412-31820-2
    
    """

    nx = len(x)

    if p<2 or p>(nx-1):
        msg = "p must be a scalar in [2, len(x)-1]"
        return

    mse = (x*x).sum() / nx

    arcoef = np.zeros(p, float)
    str1 = np.concatenate((np.array([0.0]), x[0:nx], np.array([0.0])))
    str2 = np.concatenate((np.array([0.0]), x[1:],   np.array([0.0])))
    str3 = np.zeros(nx+1, float)

    for k in np.arange(p)+1:
        arcoef[k-1] = 2.0*(str1[1:nx-k+1] * str2[1:nx-k+1]).sum() / (str1[1:nx-k+1]**2 + str2[1:nx-k+1]**2).sum()

        mse = mse * (1.0 - arcoef[k-1]**2)

        if k>1:
            for i in np.arange(k-1)+1:
                arcoef[i-1] = str3[i] - (arcoef[k-1] * str3[k-i])

        # if k == p then skip the remaining calcs
        if k == p:
            break
            
        str3[1:k+1] = arcoef[0:k]
        for j in np.arange(nx-k-1)+1:
            str1[j] = str1[j]   - str3[k] * str2[j]
            str2[j] = str2[j+1] - str3[k] * str1[j+1]
    
    return arcoef
    

#---------------------------------------------------------------------------
# The sinc() function in SciPy 0.8.0 (and probably older versions) causes 
# errors if the array contains zeros. This was fixed in SciPy 0.9.0. Since 
# we'd like to continue to support older versions of scipy, we use local 
# copies of sinc() and firwin() (which relies on sinc()) if the installed 
# scipy is < 0.9.0.

if scipy.__version__ >= "0.9.0":
    # It's safe to use the scipy versions of these
    sinc = scipy.sinc
    import scipy.signal
    firwin = scipy.signal.firwin
else:
    # If an older scipy is installed, I use my local versions
    def firwin(N, cutoff, width=None, window='hamming'):
        """
        FIR Filter Design using windowed ideal filter method.

        Parameters
        ----------
        N      -- order of filter (number of taps)
        cutoff -- cutoff frequency of filter (normalized so that 1 corresponds to
                  Nyquist or pi radians / sample)

        width  -- if width is not None, then assume it is the approximate width of
                  the transition region (normalized so that 1 corresonds to pi)
                  for use in kaiser FIR filter design.
        window -- desired window to use. See get_window for a list
                  of windows and required parameters.

        Returns
        -------
        h      -- coefficients of length N fir filter.

        """

        from scipy.signal.signaltools import get_window
        if isinstance(width,float):
            A = 2.285*N*width + 8
            if (A < 21): beta = 0.0
            elif (A <= 50): beta = 0.5842*(A-21)**0.4 + 0.07886*(A-21)
            else: beta = 0.1102*(A-8.7)
            window=('kaiser',beta)

        win = get_window(window,N,fftbins=1)
        alpha = N//2
        m = np.arange(0,N)
        h = win*sinc(cutoff*(m-alpha))
        return h / np.sum(h,axis=0)

    def sinc(x):
        """
        Returns sin(pi*x)/(pi*x) at all points of array x.

        """
        w = np.pi * np.asarray(x)
        # w might contain 0, and so temporarily turn off warnings
        # while calculating sin(w)/w.
        old_settings = np.seterr(all='ignore')
        s = np.sin(w) / w
        np.seterr(**old_settings)
        return np.where(x==0, 1.0, s)        



