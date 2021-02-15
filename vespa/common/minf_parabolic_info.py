# Python modules


# 3rd party modules
import numpy as np

# Our modules


def minf_parabolic_info(xa, xb, xc, func, info, maxit=100, 
                                                tol=np.sqrt(1.e-7), 
                                                pn=0, 
                                                dirn=1):
    '''
    The Python version of this code was adapted from minf_parabolic.pro:
    http://idlastro.gsfc.nasa.gov/ftp/pro/math/minf_parabolic.pro

    The IDL version of this code is public domain:
    http://idlastro.gsfc.nasa.gov/idlfaq.html#A14

    This version is governed by the license found in the LICENSE file. 

    +-----------------------------------------------------------------------------
     NAME:
           MINF_PARABOLIC_INFO

            Same as MinF_Parabolic, but passes a structure to the optimization
            function, to allow additional info to be passed in.

     PURPOSE:

        Find a local minimum of a 1-D function up to specified tolerance.
        This routine assumes that the function has a minimum nearby.
        (recommend first calling minF_bracket, xa,xb,xc, to bracket minimum).
        Routine can also be applied to a scalar function of many variables,
        for such case the local minimum in a specified direction is found,
        This routine is called by minF_conj_grad, to locate minimum in the
        direction of the conjugate gradient of function of many variables.

     CALLING EXAMPLES:

        minf_parabolic_info, xa,xb,xc, xmin, fmin, FUNC="name", INFO=info    ;for 1-D func.
      or:
        minf_parabolic_info, xa,xb,xc, xmin, fmin, FUNC="name", INFO=info $
                          POINT=[0,1,1],   $
                          DIRECTION=[2,1,1]    ;for 3-D func.
     INPUTS:
        xa,xb,xc = scalars, 3 points which bracket location of minimum,
            that is, f(xb) < f(xa) and f(xb) < f(xc), so minimum exists.
            When working with function of N variables
            (xa,xb,xc) are then relative distances from POINT_NDIM,
            in the direction specified by keyword DIRECTION,
            with scale factor given by magnitude of DIRECTION.
     KEYWORDS:

        FUNC = function 
            Calling mechanism should be:  F = func( px, INFO=info )
            where:
                px = scalar or vector of independent variables, input.
                F = scalar value of function at px.
                INFO = structure

        INFO = structure containing additional info to be used in FUNC function


        POINT_NDIM = when working with function of N variables,
            use this keyword to specify the starting point in N-dim space.
            Default = 0, which assumes function is 1-D.
        DIRECTION = when working with function of N variables,
            use this keyword to specify the direction in N-dim space
            along which to bracket the local minimum, (default=1 for 1-D).
            (xa, xb, xc, x_min are then relative distances from POINT_NDIM)
        MAX_ITER = maximum allowed number iterations, default=100.
        TOLERANCE = desired accuracy of minimum location, default=sqrt(1.e-7).

     OUTPUTS:

        xmin = estimated location of minimum.
            When working with function of N variables,
            xmin is the relative distance from POINT_NDIM,
            in the direction specified by keyword DIRECTION,
            with scale factor given by magnitude of DIRECTION,
            so that min. Loc. Pmin = Point_Ndim + xmin * Direction.
        iterations = the number of iterations performed. When this value is
            > max_iter, the function was unable to find a satisfactory xmin.

     PROCEDURE:

        Brent's method to minimize a function by using parabolic interpolation,
        from Numerical Recipes (by Press, et al.), sec.10.2 (p.285).

     MODIFICATION HISTORY:
        Written, Frank Varosi NASA/GSFC 1992.
    '''    
    # WARNING:  This function does not have functional equivalence to the IDL 
    # version since IDL uses floating point precision and Python uses double.

    zeps = 1.e-7                    # machine epsilon, smallest addition.
    goldc = 1 - (np.sqrt(5)-1)/2    # complement of golden mean.

    xLo = np.where(xa < xc, xa, xc)
    xHi = np.where(xa > xc, xa, xc)
    
    xmin = xb
    fmin = func(pn + xmin * dirn, info)
    xv = xmin
    xw = xmin
    fv = fmin
    fw = fmin
    es = 0.0

    for iteration in range(maxit):
        goldstep = 1
        xm = (xLo + xHi)/2.0
        TOL1 = tol * abs(xmin) + zeps
        TOL2 = 2*TOL1

        if (abs(xmin - xm) <= (TOL2 - (xHi-xLo)/2.0)): 
            return xmin, iteration

        if (abs(es) > TOL1):
            r = (xmin-xw) * (fmin-fv)
            q = (xmin-xv) * (fmin-fw)
            p = (xmin-xv) * q + (xmin-xw) * r
            q = 2 * (q-r)
            if (q > 0): p = -p
            q = abs(q)
            etemp = es
            es = ds

            if (p > q*(xLo-xmin)) and (p < q*(xHi-xmin)) and (abs(p) < abs(q*etemp/2)):
                ds = p/q
                xu = xmin + ds
                if (xu-xLo < TOL2) or (xHi-xu < TOL2):
                    ds = TOL1 * (1-2*((xm-xmin) < 0))
                goldstep = 0

        if goldstep:
            if (xmin >= xm): 
                es = xLo-xmin
            else:
                es = xHi-xmin
            ds = goldc * es

        xu = xmin + (1-2*(ds < 0)) * np.where(abs(ds) > TOL1, abs(ds), TOL1)
        fu = func(pn + xu * dirn, info)

        if (fu <= fmin):
            if (xu >= xmin):
                xLo=xmin
            else:
                xHi=xmin
            xv = xw
            fv = fw
            xw = xmin
            fw = fmin
            xmin = xu
            fmin = fu
        else:
            if (xu < xmin): 
                xLo=xu 
            else: 
                xHi=xu
            if (fu <= fw) or (xw == xmin):
                xv = xw
                fv = fw
                xw = xu
                fw = fu
            elif (fu <= fv) or (xv == xmin) or (xv == xw):
                xv = xu
                fv = fu

    return xmin, (maxit + 1)


def _test():
    import vespa.analysis.algos.shiftopt
    import pylab
    import vespa.common.util.generic_spectral as util_generic_spectral

    xa = -11
    xb = -15
    xc = -20
    func = shiftopt.shift_correlation_function
    
    sw    = 3906.25
    dim0  = 4096
    width = 5.0
    shift = 15.2
    xx = np.arange(dim0) / sw

    lshape = util_generic_spectral.apodize(xx, width, 'Gaussian') # from 4DFT
    data0  = np.exp(1j*2.0*np.pi*xx*shift) 
    data   = data0 * lshape
    
    pylab.plot(data0, 'b')
    pylab.plot(data, 'r')
    pylab.plot(lshape, 'g')
    pylab.show()
    
    info = { "data"  : data, 
             "acqsw" : sw, 
             "expo"  : lshape }

    minshift, iter_ = minf_parabolic_info(xa, xb, xc, func, info)
                            
    print(minshift, '   ', iter_)

if __name__ == '__main__':
    _test()

"""
Minf unit test result:
minshift = 2.16607634704e-005
"""