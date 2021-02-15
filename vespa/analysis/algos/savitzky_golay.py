# Python modules


# 3rd party modules
import numpy as np


def savitzky_golay(y, window_size, order, deriv=0):
    """
    Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techhniques.
    Parameters

    Obtained from: http://www.scipy.org/Cookbook/SavitzkyGolay
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = _savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError as msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = list(range(order+1))
    half_window = (window_size -1) // 2
    
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv]
    
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    
    return np.convolve( m, y, mode='valid')



def _test():
    
    import numpy as np
    import pylab
    import io, zlib, base64

    # bjs - I create a 256 point array, to make the timing trial a bit
    #       more real world, by appending "y" four times to itself

    y = 'eNqb7BfqGxDJyFDGUK2eklqcXKRupaBuk2ahrqOgnpZfVFKUmBefX5SSChJ3S8wpTgWKF2ckFqQC+RpmJjqaOgq1CuQDrqzbBlxNx1/tfzH/QZIr77H9L90PaxZ/Pm6flM419/qNn/YHbFJOtS84tl+34mL4hGdP7Gcf+xZ1ofrh/mbJ4A3Hbp6wnzpHuDpj/lt7kc1xPE9dPtu/9ND7r1j7w77PYc+6+A1f908/ZXb6QMwq+8COBucVQR/s5efu+hQ57/N+U96uZ4xrntmXnqk8tu/Ja/v6O3U/U1xf7q/h03E/4XDVvvFC1Q+jwFP2144wxV9ce3S/bXL8po8JV/bLLRWad+H1pf3qN+rnzt72wF7zLe8+g5uMB/7x3lrd337Bvl330YPyim37a9dtS6qP/LVfO4dth0L0JfttDSADL9u7LLbumLvyoz1/wsR4+yPP7T967LGZq31k/+6Lx2scFB7vV8+YfFhd9qn9iy07BL6eWbD/x9K6eyFLPtn/CH68dPaRD/tVtjyN6dJ+ZW947n6LU8u3/UZhs/+Hz31gf6vnskBoxfP955sC7h2/eMVeME2lZ8OnpfsZ/u6r//zo9X7plEkTN294YS8bI3kj2KHFvoqh+8+vhW/2n298yiF54sn+Ft5NN/3PNe6/bFDj4sP5z/5O8ktn8SvMB0Lv9RVEffi1/+PVvkcypS/2b44xs3q+/JB9CkuIRVPIGXulnIqqpTrMDlm/w8NmTflnf/B/2KZk5d/2VWar7cJvv7SP5igKOOTD5sD/YUnraa059uf+Xn287+ou+7jWKfNDqv/YAwAG+zSG'

    y = base64.b64decode(y)
    y = zlib.decompress(y)
    buf = io.BytesIO(y)
    y = np.load(buf)
    y = np.concatenate((y,y,y,y))       # create 512 pt array

    wsize = 31      # must be odd
    order = 4
    deriv = 0
    ntest = 1       # when doing 10 repeats it takes 0.018 sec total
    
    y = np.array(y)
            
    for i in list(range(ntest)):                                                                
        yout = savitzky_golay(y, wsize, order, deriv=deriv)                                                                            
    
    pylab.plot(y)
    pylab.plot(yout,  'b', linewidth=2)
    pylab.show()
    
    bob = 10
    bob += 1



if __name__ == '__main__':
    """
    Not as smooth as a LOWESS filtering, but 10-15 times faster.
    
    Need to do more playing with parameters and my real data to better dial
    in the results I'd like to get.
    
    """
    

#    from vespa.analysis import loess_1d     # from loess 2.0.11 on PyPI

    import cProfile
    
    cProfile.run('_test()')
#    _test()

