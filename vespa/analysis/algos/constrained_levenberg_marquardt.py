# Python modules


# 3rd party modules
import numpy as np


def constrained_levenberg_marquardt(y, w, a, limits, function=None, 
                                                     itmax=50, 
                                                     tol=1.e-3):
    '''
    =========
    Arguments 
    =========
    **y:**  [array][float or complex]

      Array of dependent variable, can be real or complex

    **w:**  [array][float]

      Array of weights, the same length as y. For no weighting, w[i] = 1.0.
      For instrumental weighting, w[i] = 1.0/y[i], etc.

    **a:**  [array][float]

      Array with as many elements as the number of fitting parameter terms, 
      that contains the initial estimate for each parameter.

    **limits:**  [array 2D][float]

      Array with as many elements (x2) as the number of fitting parameter 
      terms that contain the min and max limits of the search space for each
      fitting term. For example, if there are 10 paramters, then the limits
      array would be 2 x 10, with the first 10 item array being the min limit
      value for each item and the second 10 item array being the max limit 
      value

    **function:**  [keyword][object][default=None]

      referece to the function to fit. The procedure must be written as 
      described under RESTRICTIONS, below.

    **itmax:**  [keyword][integer][default=50]

      maximum number of iterations. Default = 50

    **tol:**  [keyword][float][default=1.0e-3]

      the convergence tolerance. The routine returns when the relative 
      decrease in chi-squared is less than TOL in an iteration. Default 
      is 1.0e-3.

    =========
    Outputs 
    =========

    **yfit:**  [array][complex or float]

      Vector form of the optimization function for the final values in the
      parameter list.

    **a:**  [array][float]
    
      final vector of optimized parameters

    **sigmaa:**  [array][float]

      vector of standard deviations for the parameters

    **chi2:**  [float]

      The value of chi-squared on exit

    **wchi2:**  [float]

      The value of chi-squared on exit. This value can be NaN. 

    **badfit:**  [boolean]

      returns False if no math errors, True if there is a finite value
                     issue in matrix inversion or if the internal loop for parameter
               refinement reaches a max iteration of 100


    =========== 
    Description
    ===========   
    
    constrained_levenberg_marquardt is a non-linear least squares fit to a 
    function of an arbitrary number of parameters.  The function may be any 
    non-linear function.  If available, partial derivatives can be calculated 
    by the user function, else this routine will estimate partial derivatives
    with a forward difference approximation. 
    
    This routine was initially based on the IDL CURVEFIT procedure, which was
    based on "CURFIT", a least squares fit to a non-linear function, pages 
    237-239, Bevington, Data Reduction and Error Analysis for the Physical 
    Sciences. "This method is the Gradient-expansion algorithm which combines 
    the best features of the gradient search with the method of linearizing 
    the fitting function."
    
    Iterations are performed until the chi square changes by only TOL or until 
    ITMAX iterations have been performed. The initial guess of the parameter 
    values should be as close to the actual values as possible or the solution
    may not converge.
    
    Simple hard parameter contraints were added by Soher and Kiefer. Also 
    added a finite iteration limit to the curvature matrix inversion repeat 
    loop.  Also has a check for finite values in the repeat loop.  If either
    are twigged, the current A string is returned

    Fitting Function Restrictions -----------------------------------
           
    The function to be fit must be defined and passed into the procedure as
    a function reference in the 'function' keyword. This function, must accept
    values of 'a' (the fitted function's parameter values), and return F (the 
    function's value at X), and PDER (a 2D array of partial derivatives).
    For an example, see the gfunct() and _test() functions below.
           
    A call to 'function' is entered as:

    f, pder = funct(a)

    a    =  Vector of NTERMS function parameters, input.
    f    =  Vector of NPOINT values of function, y[i] = funct(x[i]), output.
    pder =  Array, (NPOINT, NTERMS), of partial derivatives of funct.
            PDER(I,J) = Derivative of function at ith point with
            respect to jth parameter.  

    Note. pder can be returned as None, in which case the forward difference
          approximation will be calculated by constrained_levenberg_marquardt()
          
    
    ======    
    Syntax
    ====== 
    ::
    
    yfit, a, sigmaa, chi2, wchi2, badfit = 
        constrained_levenberg_marquardt(y, w, a, limits, function, itmax, tol)
    


    '''
    if y.dtype in ('complex64', 'complex128'):
        complex_flag = True
    else:
        complex_flag = False
    
    if complex_flag:    
        y = np.concatenate([y.real, y.imag])
        w = np.concatenate([w, w])
    
    nterms = len(a)   # number of parameters 
    eps = np.sqrt(np.finfo(np.float).eps)

    # order limits [lowerlimit,upperlimit]
    limits = np.array([np.where(limits[0,:] < limits[1,:], limits[0,:], limits[1,:]), 
                       np.where(limits[0,:] > limits[1,:], limits[0,:], limits[1,:])])
 
    no_zero = np.where( (limits[0,:] != 0) + (limits[1,:] != 0) )[0] # '+' equals logical 'or'
    count_no_zero = np.size(no_zero)

    nfree = np.size(y)-nterms       # Degrees of freedom

    if nfree <= 0: 
        raise ValueError("constrained_levenberg_marquardt error - not enough degrees of freedom = %s " % str(nfree))

    flambda = 0.001                 # Initial lambda

    #Define the partial derivative array ---
    pder = np.zeros((np.size(y), nterms), float)

    inc = eps * np.abs(a) 
    badfit = False

    for iter in range(1, itmax+1):  # Iteration loop
        # Evaluate alpha and beta matrices ---
        # Evaluate function and estimate partial derivatives
        if (count_no_zero > 0):    # there are constraints
            # apply hard limits, keep in mind the increment for ---
            # the derivatives attempt constraint with current inc
            a[no_zero] = np.where( a[no_zero] > (limits[0,no_zero]+1.1*inc[no_zero]),
                                   a[no_zero],  (limits[0,no_zero]+1.1*inc[no_zero]) )
            a[no_zero] = np.where( a[no_zero] < (limits[1,no_zero]-1.1*inc[no_zero]),
                                   a[no_zero],  (limits[1,no_zero]-1.1*inc[no_zero]) )

            # if still out of constraints reset increment AND params ---
            indxb = np.where(a[no_zero] < limits[0,no_zero])[0]
            indxt = np.where(a[no_zero] > limits[1,no_zero])[0]

            if indxb:
                a[no_zero[indxb]]   = limits[0,no_zero[indxb]] * (1.0+eps)
                inc[no_zero[indxb]] = limits[0,no_zero[indxb]] * eps

            if indxt:
                a[no_zero[indxt]]   = limits[1,no_zero[indxt]] * (1.0-eps)
                inc[no_zero[indxt]] = limits[1,no_zero[indxt]] * eps

            yfit, pder1 = function(a)
            if complex_flag:    
                yfit = np.concatenate([yfit.real, yfit.imag])
                if pder1 is not None:
                    pder1 = np.concatenate([pder1.real, pder1.imag], axis=1)

            pder2 = pder1.copy()
            pder1 = None

            if pder1 is None:
                for term in range(nterms):
                    # Copy current parameters ---
                    p = a.copy()
    
                    # Increment size for forward difference derivative ---
                    p[term] = p[term] + inc[term]
                    yfit1, _ = function(p)
                    if complex_flag:    
                        yfit1 = np.concatenate([yfit1.real, yfit1.imag])
                    pder[:,term] = (yfit1-yfit)/inc[term]
    
                    # Try to set difference to 1e-5 for next iteration ---
                    inc[term] = np.size(y)*1e-5/np.sum(abs(pder[:,term]))
                pder1 = pder
            else:
                pder1 = pder1.T


        else:  # no constraints applied
            yfit, pder1 = function(a)
            if complex_flag:    
                yfit = np.concatenate([yfit.real, yfit.imag])
                if pder1 != None:
                    pder1 = np.concatenate([pder1.real, pder1.imag], axis=1)

            if pder1 == None:
                for term in range(nterms):
                    # Copy current parameters ---
                    p = a.copy()
    
                    # Increment size for forward difference derivative ---
                    p[term] = p[term] + inc[term]
                    yfit1, _ = function(p)
                    if complex_flag:    
                        yfit1 = np.concatenate([yfit1.real, yfit1.imag])
                    pder[:,term] = (yfit1-yfit)/inc[term]
                
                    # Try to set difference to 1e-5 for next iteration ---
                    inc[term] = np.size(y)*1e-5/np.sum(abs(pder1[:,term]))
                pder1 = pder
            else:
                pder1 = pder1.T

        beta   = np.dot((y-yfit)*w, pder1) 
        alpha  = np.dot(pder1.transpose(), np.outer(w, np.zeros(nterms,float)+1) * pder1)
        chisq1 = np.sum(w*(y-yfit)**2)/nfree # Present chi squared.

        bs = 0

        # Invert modified curvature matrix to find new parameters. ---
        while True:
            c = np.sqrt(np.outer(alpha.diagonal(), alpha.diagonal())) 

            arr = alpha/c
            arr = arr + flambda * arr * np.diag(np.zeros(len(arr))+1)
            try:
                arr = np.linalg.inv(arr)
            except np.linalg.LinAlgError:
                badfit = True
                chisqr = np.NAN
                break


            # Calculate new parameters ---
            b = a + np.dot(arr/c, np.transpose(beta))

            # Check for and set any constraints ---
            if count_no_zero > 0:
                b[no_zero] = np.where(b[no_zero] < limits[1,no_zero], b[no_zero], limits[1,no_zero])
                b[no_zero] = np.where(b[no_zero] > limits[0,no_zero], b[no_zero], limits[0,no_zero])

            yfit, pder1 = function(b)
            if complex_flag:    
                yfit = np.concatenate([yfit.real, yfit.imag])
            chisqr  = np.sum(w*(y-yfit)**2)/nfree     # New chisqr
            flambda = flambda*10.0                    # Assume fit got worse

            if not np.isfinite(chisqr):
                # Finite Error problem
                badfit = True

            bs = bs+1
            if bs > 100:
                # Kicker Loop Error - ccfit
                badfit = True

            if chisqr <= chisq1: break
            if badfit: break

        if badfit:
            break
        else:
            flambda = flambda/100.0      # Decrease flambda by factor of 10
            a = b                        # Save new parameter estimate.

        if (chisq1-chisqr)/chisq1 <= tol: break

    sigmaa = np.sqrt(arr.diagonal()/alpha.diagonal())   # Return sigma's
    wchi2  = chisqr                             # Return weighted chi squared
    chi2   = np.sum((y-yfit)**2)/nfree          # Return chi squared
    
    

    dim0 = len(yfit)/2
    yfit = yfit[0:int(dim0)] + 1j * yfit[int(dim0):int(2*dim0)]  # Convert to complex array

    return yfit, a, sigmaa, chi2, wchi2, badfit # Return result



def gfunct(a):  
    x  = np.arange(10)  
    bx = np.exp(a[1] * x)  
    f  = a[0] * bx + a[2]  
    
    # If the procedure is called with four parameters, calculate the  
    # partial derivatives.  
    pder = np.array([bx, a[0] * x * bx, np.ones(len(x))])  

    return f, pder



def _test():

    # Compute the fit to the function we have just defined. 
    # First, define the independent and dependent variables: 
    y = np.array([12.0, 11.0, 10.2, 9.4, 8.7, 8.1, 7.5, 6.9, 6.5, 6.1])  
      
    # Define a vector of weights.  
    w = 1.0/y  
      
    # Provide an initial guess of the function's parameters.  
    a = np.array([10.0,-0.1,2.0])  

    limits = [[ 8.0,-1.2,-2.0],
              [12.0, 1.2, 6.0]]
    limits = np.array(limits)
    
    # Compute the parameters.  
    yfit, a, sig, chis, wchis, badfit = constrained_levenberg_marquardt(y, w, a, limits, gfunct)
      
    # Print the parameters returned in A.  
    print('constrained_levenberg_marquardt : ', a)
    print('IDL CURVEFIT results            :  9.91120    -0.100883      2.07773')  


if __name__ == '__main__':
    _test()

