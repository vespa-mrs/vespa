# Python modules

# 3rd party modules
import numpy as np
from scipy import interpolate


def splinevar(y, nknots, kord, hz_per_pt=1.0):
    '''
    =========
    Arguments 
    =========
    **y:**  [array][float]

      dependent data to be approximated using splines

    **nknots:**  [integer]

      Literally, the number of spline knots to use in the fit to the data.
      But, since this only works well for splrep when knots are fixed, we
      use this as a surrogate to indicate the 's' smoothing value that 
      splrep can take to determine the tradeoff between closeness and 
      smoothness of the fit.

    **kord:**  [integer]

      object containing fitting parameters


    =========== 
    Description
    ===========   
    Returns the b-spline representation of the data in y using the 
    scipy.interpolate routines splrep and splev. This version allows
    the routine to determine the number and spacing of the knots based
    on a smoothing value set by an empircal guess based on the number
    of knots the user wanted to use in the first place.
    
    splrep finds the b-spline representation of a 1-D curve given the 
    set of data points (x[i], y[i]) it determine a smooth spline 
    approximation of degree k on the interval xb <= x <= xe. The 
    coefficients, c, and the knot points, t, are returned. Uses the 
    FORTRAN routine curfit from FITPACK.
    
    splev evaluates the value of the smoothing polynomial. This is a 
    wrapper around the FORTRAN routines splev and splder of FITPACK


    ======    
    Syntax
    ====== 
    ::

      res = splinevar(y, nknots, kord)

    '''
    ny = np.size(y)
    
    # the following weights were empirically optimized for 64 random points
    # for which we assumed a hetz per point of 1.0, we created a range of
    # smoothing parameters from 1-100 where smoother results are seen for
    # higher smoothing parameter values.
    #
    # We have to scale the weights to the actual hertz_per_points in the 
    # data, but only if that info is provided in the keyword.

    nk = np.where(nknots >  1, nknots, 1)
    nk = np.where(nknots < 100, nknots, 100)

    sarr = np.array( [ 16,  18,  20,  22,  24,  26,  28,  30,  33,  36, \
                       40,  44,  50,  54,  60,  64,  68,  72,  76,  80, \
                       84,  88,  92,  96, 100, 105, 110, 115, 120, 125, \
                      130, 135, 140, 145, 150, 155, 160, 165, 170, 175, \
                      180, 185, 190, 195, 200, 205, 210, 215, 220, 225, \
                      230, 235, 240, 245, 250, 255, 260, 265, 270, 275, \
                      280, 285, 290, 295, 300, 305, 310, 315, 320, 325, \
                      330, 335, 340, 345, 350, 355, 360, 365, 370, 375, \
                      380, 385, 390, 390, 400, 405, 410, 415, 420, 425, \
                      430, 435, 440, 445, 450, 460, 470, 480, 490, 500 ])

    # here we allow the user to provide info about the hertz per point for
    # each data point. This allows us to scale the smoothing factor to try
    # to keep smoothness equivalent across baselines with different sweep
    # widths and or data points. This is likely also important for zero
    # filled data too.
    sarr = sarr * hz_per_pt
    #sarr *= hz_per_pt
    
    # assumes a parameter from 1-100
    s = sarr[int(nk-1)]
    
    x  = np.arange(ny)
    wt = np.zeros(ny, float) + 1.0
    
    tck = interpolate.splrep(x, y, wt, k=kord, task=0, s=s)
    res = interpolate.splev(x, tck)
    
    return res


def splinefix(y, nknots, kord):
    '''
    =========
    Arguments 
    =========
    **y:**  [array][float]

      dependent data to be approximated using splines

    **nknots:**  [integer]

      The number of evenly spaced spline knots to use in the fit to 
      the data.

    **kord:**  [integer]

      object containing fitting parameters


    =========== 
    Description
    ===========   
    Returns the b-spline representation of the data in y using the 
    scipy.interpolate routines splrep and splev. This version tells 
    the routine the locations of the knots.
    
    splrep finds the b-spline representation of a 1-D curve given the 
    set of data points (x[i], y[i]) it determine a smooth spline 
    approximation of degree k on the interval xb <= x <= xe. The 
    coefficients, c, and the knot points, t, are returned. Uses the 
    FORTRAN routine curfit from FITPACK.
    
    splev evaluates the value of the smoothing polynomial. This is a 
    wrapper around the FORTRAN routines splev and splder of FITPACK

    ======    
    Syntax
    ====== 
    ::

      res = splinefix(y, nknots, kord)

    '''    
    ny = np.size(y)
    x  = np.arange(ny)
    wt = np.zeros(ny, float) + 1.0
    
    # Calculate an even distribution of the interior knots 
    knots  = (ny-1.0)*(np.arange(nknots)/(nknots-1.0))
    knots  = [round(k) for k in knots]
    
    tck = interpolate.splrep(x, y, wt, k=kord, t=knots[1:-1])
    res = interpolate.splev(x, tck)
    
    return res


def _test_splinevar():
    import pylab as pl
    
    kord = 3
    nknots = 10
    y = [-0.836854, -0.172280,  0.187117,   1.61544,   -0.176774,  0.653145, -0.546364, 0.194146,\
          0.925709,  1.20432,   1.53055,   -1.35556,    0.0514889, 1.02018,  -1.22616,  0.708497,\
          0.871673, -0.789721,  0.332079,   0.205603,  -0.169367, -0.318417, -0.295643, 0.522291,\
         -2.23105,   0.258274, -0.0877757, -1.64685,    0.286812,  0.299986,  1.10391,  0.742706,\
         -0.157581, -0.597687,  0.659809,  -0.0328137,  1.16512,  -1.04800,   0.817815,-1.40729,\
          0.519207, -0.733439,  0.325304,  -0.0428672, -0.871454,  0.771570,  0.00988832, -0.894773,\
         -0.649426, -0.00869429,1.87727,   -2.47856,   -1.68368,  -0.764296,  0.145749, 0.221329,\
          2.39680,   1.91128,   1.69614,    0.808025,   2.78748,   0.0275070, 0.104824, 1.78011]
    
    y = np.array(y)
    
    res1 = splinevar(y, 10, kord, hz_per_pt=0.5)
#    res2 = splinevar(y, 10, kord, hz_per_pt=1.0)
#    res3 = splinevar(y, 15, kord, hz_per_pt=1.0)
    pl.plot(y)
    pl.plot(res1, 'r')
#    pl.plot(res2, 'black')
#    pl.plot(res3, 'g')
    pl.show()



def _test_splinefix():
    import pylab as pl
    
    kord = 3
    nknots = 16
    y = [-0.836854, -0.172280,  0.187117,   1.61544,   -0.176774,  0.653145, -0.546364, 0.194146,
          0.925709,  1.20432,   1.53055,   -1.35556,    0.0514889, 1.02018,  -1.22616,  0.708497,
          0.871673, -0.789721,  0.332079,   0.205603,  -0.169367, -0.318417, -0.295643, 0.522291,
         -2.23105,   0.258274, -0.0877757, -1.64685,    0.286812,  0.299986,  1.10391,  0.742706,
         -0.157581, -0.597687,  0.659809,  -0.0328137,  1.16512,  -1.04800,   0.817815,-1.40729,
          0.519207, -0.733439,  0.325304,  -0.0428672, -0.871454,  0.771570,  0.00988832, -0.894773,
         -0.649426, -0.00869429,1.87727,   -2.47856,   -1.68368,  -0.764296,  0.145749, 0.221329,
          2.39680,   1.91128,   1.69614,    0.808025,   2.78748,   0.0275070, 0.104824, 1.78011]
    y = np.array(y)
    
    res = splinefix(y, nknots, kord)
    print(res)
    pl.plot(y)
    pl.plot(res, 'r')
    pl.show()



if __name__ == '__main__':
    _test_splinevar()
#    _test_splinefix()
