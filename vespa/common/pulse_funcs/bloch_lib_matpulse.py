# Python modules


# 3rd party modules
import numpy as np

# our modules
import vespa.common.pulse_funcs.constants as pf_constants


# Adopted from Karl Young (PyPulse.py) : 07/23/10 #
# Further adopted by David Todd for RFPulse : 01/24/11 (2011.01.24) #

GAMMA1H = pf_constants.GAMMA1H
PI      = np.pi
# epsilon to replace zeros in rotation angle calculation
reps = 0.00000001


GAMMA = 26753.0
TWOPI = 6.283185



def grad_field(mb1gl,extended,convert,offset,dwell,resol):
    """
    This function returns a gradient field, i.e. set
    of coordinates for calculating variation in tipping
    of magnetization - returned in 'spatial' (vn) and
    magnetic coordinates (gx).

    Matson, G.B., An Integrated Program For Amplitude-
    Modulated RF Pulse Generation And Remapping With
    Shaped Gradients, Mag. Res. Im., 12(8) (1994), 1205
    """
    # Calc.Nyquist
    nq = (1.0/(2.0*dwell))*1000.0
    
    # Array for spatial values
    vn = nq*np.arange(-resol/2,(resol/2))/float(resol/2)
    
    # For extended range
    if extended == True:
        vn = 4.0*vn
    
    # in mTesla
    gx = vn/GAMMA1H
    
    # units in radians/mT ([rad]*[KHz/mT]*[musec]/[KHz/MHz] = [rad/mT])
    constb = 2.0 * PI * GAMMA1H * dwell/1000.0
    
    # Convert to mm (mb1pmv == 1) for display
    if convert == 1:
        vn = (1000*vn)/(mb1gl*GAMMA1H)
        
    # Add offset
    gx += offset
    return(vn,gx)
   
    
def bloch_b1(b1,g1const,extended,convert,offset,dwell,resol,vn,gx):
    """
    This function returns a set of frequency (or space)
    dependent magnetizations calculated via the Bloch
    equations. 
    
    Input: 
    b1       - An RF pulse in the form of a complex array
    extended - Flag indicating whether we want extended range.    
    g1const  - Conversion factor from kHz to millimeters, as if applying constant gradient.
    convert  - Flag (0 or 1), indicating if conversion to mm
               should be done or not. (1 ==> convert to mm)
    offset   - Resonance Offset
    dwell    - Dwell time (in microseconds).
    resol    - Calculation Resolution (Number of freq locations over which bloch sim applied)
    vn       - spectral range for profile
    gx       - gradient array in mTesla off central resonance frequency

    Output:
    alph, bet - Cayley Klein parameters
    vn        

    Matson, G.B., An Integrated Program For Amplitude-
    Modulated RF Pulse Generation And Remapping With
    Shaped Gradients, Mag. Res. Im., 12(8) (1994), 1205
    
    Derived from: MBLB1CC.M
    """

    # Additional versions of the Bloch equations are available 
    # for more complex situations, Bloch with flow, and Bloch
    # with Gradients. See Jerry Matson's Matpulse code.

#    # Get gradient field to do calcs over
#    vn,gx = grad_field(g1const,extended,convert,offset,dwell,resol)
    gx2 = gx*gx

    # units in radians/mT ([rad]*[KHz/mT]*[musec]/[KHz/MHz] = [rad/mT])
    constb = 2*np.pi*pf_constants.GAMMA1H*dwell/1000.

    # Do Bloch calcs
    resol = len(gx)
    alph  = np.ones(resol,np.complex64)
    bet   = np.zeros(resol,np.complex64)

    for j in range(len(b1)):
        phib = -constb*(np.sqrt(gx2 + b1[j]*np.conj(b1[j])))

        # Avoid divide by 0
        phib = np.where(phib == 0.0,reps,phib)

        # and zero small b1's...
        if np.absolute(b1[j]) < reps:
            b1[j] = complex(0,0)

        tempdiv = constb / np.absolute(phib)
        nx = np.real(b1[j]) * tempdiv
        ny = np.imag(b1[j]) * tempdiv
        nz = gx * tempdiv


        phib2 = phib/2.0
        si = np.sin(phib2)
        alpha = np.cos(phib2) - 1.0j*nz*si
        beta = -1.0j*(nx + 1j*ny)*si

        alphn = alpha*alph - np.conj(beta)*bet
        betn  = beta*alph  + np.conj(alpha)*bet

        alph = alphn
        bet  = betn

    return(alph,bet)
    

def bloch_b2g2( b2, g3, f2, extended, offset, dwell, resol):
    """
    This function returns a set of frequency (or space)
    dependent magnetizations calculated via the Bloch
    equations. 
    
    Input: 
    b2       - An complex array with RF values in mT
    g3       - An array with gradient values in mT, same length as b2.
    f2       - An array with frequency values, same length as b2, may be zeros.
    extended - Flag indicating whether we want extended range.    
    offset   - Resonance Offset
    dwell - Dwell time (in microseconds).
    resol  - Calculation Resolution (Number of points used in Fourier transform calculation.?)

  NB. The logic below likely needs to be put into the block.py in common
      so depending on whether an RF waveform has a time varying gradient 
      or not we can pick which block algorithm to call.
      
  In MatPulse we have a module named, -----------------------------------------------
      MBLCALCC.M  - Script for calling Bloch equations
                  - Called by Calculate in mblbm
  
  In this script, if mblb1v  == 1 | mblb1gv == 1 then it is a b1/g1 calc
                  if mblb2gv == 1 | mblbgfv == 1)then it is a b2/g2 / f2 calc

  Next step in b2/g2/f2 is check if mblbgfv == 1  then needs f2 calc, done by

          mblg3 = mblg3c(mrgg2,mrlslew,dwell) ; % G3
          mblf2 = -mblf2mmv*mmggamma*mblg3/1000  ; % F2 values in kHz

  If mblbgfv != 1 then we calc g3 and f2 this way

          mblg3 = mblg3c(mrgg2,mrlslew,dwell) ;    % G3
          mblf2 = zeros(size(mblg3)) ;

        where mblg3c is calculated as (seems to be checking if gradient is allowed by slew rate
        
                function g3 = mblg3c(g2,slew,dwelln)
                
                % Function file for obtaining integral (ave) for g2 (=g3)
                % Called by mblcalcc
                % CORRECTED
                % Convention is that the gradient spends the first dwell 
                %     period at g2(1), and slews to subsequent values 
                % No checking for slew here
                % Set up g2l for g2(l) - g2(l-1) 
                
                n = length(g2) ;
                g2l = g2 ; 
                g2l(n) = [] ;
                g2l = [g2(1) g2l] ;
                delg2 = (g2 - g2l) ;
                
                % Note that slew in ms and dwell in us
                % Del area g3*dwell is delg2*(dwell/1000-abs(delg2/(2*slew)),
                % Leading to g3 = delg2(1-abs(delg3/(2*dwell/1000*slew))
                
                g3 = delg2.*(1-abs(delg2)/(2*dwelln/1000*slew)) ;
                g3 = g2l + g3 ;        

  Then set freq offset 
      if mblofsf == 1 ;           % Use offset
           mbltesla = mblofsv/(1000*mmggamma) ;
       else    mbltesla = 0 ;

  Finally, we can run the Bloc Equa for b2g2f2

    [alph,bet,vn] = mblbgfcc(b2,g3,f2,mblperv,offset)
    
    Bloch equa for b2 and g2 (with f2 and offset)
    Called by mblcalcc
    
    b2, g3, f2 [mT], ext. range flag, offset [mT]

    Specify globals ********  Moved to input variables
    global mmgresol dwell mmggamma
    
    """


    # Set up frequency axis & convert to mm ********************

    nq = (1.0/(2.0*dwell))*1000      # Nyquist freq in kHz (plot limits)
    
    # array for spatial values
    vn = nq*np.arange(-resol/2,(resol/2))/float(resol/2.0)

    # For extended range  
    if extended == True:         # Extended range
        vn = 4*vn  

    # in mTesla
    gf     = vn/pf_constants.GAMMA1H                            
    # units in radians/mT ([rad]*[KHz/mT]*[musec]/[KHz/MHz] = [rad/mT])
    constb = 2 * PI * pf_constants.GAMMA1H * dwell / 1000.0  

    # Convert to mm ****

    vn = (1000.0*vn)/(np.max(g3)*pf_constants.GAMMA1H)   # x in mm

    # Do actual Bloch equa calculations ***************************

    resol = len(gf) 
    alph  = np.ones(resol,np.complex64)
    bet   = np.zeros(resol,np.complex64)

    for j in range(len(b2)):

        gx = g3[j]*vn/1000 + f2[j]/pf_constants.GAMMA1H + offset   # Set gx in mT
        phib = -constb*(np.sqrt(gx*gx + b2[j]*np.conj(b2[j]))) 

        # Avoid divide by 0
        phib = np.where(phib == 0.0,reps,phib)

        # and zero small b1's...
        if np.absolute(b2[j]) < reps:   # For avoid divide by zero
            b2[j] = complex(0,0) 

        tempdiv = constb / np.absolute(phib)
        nx = np.real(b2[j]) * tempdiv         # real b1 along x
        ny = np.imag(b2[j]) * tempdiv         # imag b1 along y
        nz = gx * tempdiv 

        phib2 = phib/2.0
        si = np.sin(phib2) 

        alpha = np.cos(phib2) - 1.0j*nz*si 
        beta = -1.0j*(nx + 1j*ny)*si 

        alphn = alpha*alph - np.conj(beta)*bet 
        betn  = beta*alph  + np.conj(alpha)*bet 

        alph = alphn  
        bet = betn 

    return(alph,bet,vn)
        
        
def mag_recon(alpha,beta,minit):
    """
    This function returns a resultant magnetization
    vector given an initial magnetization vector
    (minit) and a 'pulse' in the form of Cayley-Klein
    parameters alpha and beta - note that these must match
    i.e. if there are n initial magnetization vectors, there must be
    n sets of alpha and beta parameters
    
    alpha, beta - output of Bloch equations in the Cayley-Klein formalism.
    
    minit       - initial magnetization. This depends on pulse type.
                  Excitation, Saturation, Inversion: Only a z component
                   Spin-Echo: Only a y component
    
    Matson, G.B., An Integrated Program For Amplitude-
    Modulated RF Pulse Generation And Remapping With
    Shaped Gradients, Mag. Res. Im., 12(8) (1994), 1205
    """

    # Do a little checking
    dimension = len(minit)
    if dimension != 3:
        raise ValueError("Initial magnetization vector needs to be 3D")

    # set up internal arrays for multiplication and output array
    mi = np.zeros([3],np.complex64)
    mo = np.zeros([3],np.complex64)
    mout = np.zeros([3],np.float32)
    
    # Check alpha, beta too

    mi[0] = minit[0] + 1j*minit[1]
    mi[1] = minit[0] - 1j*minit[1]
    mi[2] = minit[2] 
    
    alphac = np.conj(alpha)
    betac = np.conj(beta)
    tran = np.zeros([3,3],np.complex64)
    
    tran[0,0] = alphac*alphac
    tran[1,0] = -betac*betac
    tran[2,0] = -alphac*betac
    tran[0,1] = -beta*beta
    tran[1,1] = alpha*alpha
    tran[2,1] = -alpha*beta
    tran[0,2] = 2.0*alphac*beta
    tran[1,2] = 2.0*alpha*betac
    tran[2,2] = alpha*alphac-beta*betac

    mo = np.dot(tran,mi)

    # mx,my,mz
    mout[0] = np.real(mo[0])
    mout[1] = np.imag(mo[0])
    mout[2] = np.real(mo[2])

    return(mout)


def mag_recon_vector(_alpha, _beta, _minit):
    """
    This function returns a resultant magnetization
    vector given an initial magnetization vector
    (minit) and a 'pulse' in the form of Cayley-Klein
    parameters alpha and beta - note that these must match
    i.e. if there are n initial magnetization vectors, there must be
    n sets of alpha and beta parameters
    
    alpha, beta - output of Bloch equations in the Cayley-Klein formalism.
    
    minit       - initial magnetization. This depends on pulse type.
                  Excitation, Saturation, Inversion: Only a z component
                   Spin-Echo: Only a y component
    
    Matson, G.B., An Integrated Program For Amplitude-
    Modulated RF Pulse Generation And Remapping With
    Shaped Gradients, Mag. Res. Im., 12(8) (1994), 1205
    """
    dimension = len(_minit[0])
    if dimension != 3:
        raise ValueError("Initial magnetization vector needs to be 3D")
    
    nx = len(_alpha)

    _alphac = np.conj(_alpha)
    _betac  = np.conj(_beta)
    _tran   = np.zeros([3,3,nx],np.complex64)
    
    _tran[0,0,:] = _alphac*_alphac
    _tran[1,0,:] = -_betac*_betac
    _tran[2,0,:] = -_alphac*_betac
    _tran[0,1,:] = -_beta*_beta
    _tran[1,1,:] = _alpha*_alpha
    _tran[2,1,:] = -_alpha*_beta
    _tran[0,2,:] = 2.0*_alphac*_beta
    _tran[1,2,:] = 2.0*_alpha*_betac
    _tran[2,2,:] = _alpha*_alphac-_beta*_betac

    # set up internal arrays for multiplication and output array
    mi   = np.zeros([3],np.complex64)
    mo   = np.zeros([3],np.complex64)
    mout = np.zeros([3],np.float32)

    outz = np.zeros([nx,3], np.float64)
    for i in range(nx):
        
        minit = _minit[i]
        
        # Check alpha, beta too
    
        mi[0] = minit[0] + 1j*minit[1]
        mi[1] = minit[0] - 1j*minit[1]
        mi[2] = minit[2] 

        mo = np.dot(_tran[:,:,i],mi)
    
        # mx,my,mz
        mout[0] = np.real(mo[0])
        mout[1] = np.imag(mo[0])
        mout[2] = np.real(mo[2])

        outz[i,:] = mout 

    return(outz)




#------------------------------------------------------------------------------
# Testing

def hsinc(npts, ncycles, filter='hamming'):
    """
    Returns a sinc function of length npts, with ncycles sinc-cycles. This
    yields a time-bandwidth value of 4 * ncycles
    """    
    t = np.arange(npts) - (npts/2.0)
    t = t / (npts/2.0)
    val = 2*np.pi*ncycles*t + 0.00001
    res = np.sin(val) / val
    if filter == 'hamming':
        res = res * 4 * ncycles * (0.54 + 0.46*np.cos(np.pi*t)) / npts
        
    return res   


def bloch_twostep(b1, T, x, vn, gx):
    """
    Having this method here allows me to time the cumulative time it
    takes to get from b1 to frequency profile in cProfile call
    
    >python -m cProfile -s cumulative bloch_lib_matpulse.py
    
    without having to add up the times in bloch_b1 and mag_recon
    separately
    
    """
    nx = len(x)
    xax = vn
    
    [alpha,beta] = bloch_b1(b1,1.0, False,False,0.0,T*1000000,nx, vn, gx)

    minit1 = np.zeros([nx, 3], np.float64)          
    minit1[:,2] = 1.0       # for exc, inv, sat
#     outz = np.zeros([nx,3], np.float64)
#     for i in range(len(alpha)):
#         outz[i,:] = mag_recon( alpha[i], beta[i], minit1[i] ) 

#    minit2 = np.zeros([nx, 3], np.float64)   
#    minit2[:,1] = 1.0       # for spin-echo 
#    outxy = np.zeros([nx,3], np.float64)
#    # for i in range(len(alpha)):
#    #    outxy[i,:] = mag_recon( alpha[i], beta[i], minit2[i] )
#    outz = mag_recon_vector( alpha, beta, minit1 )
    
    outz = mag_recon_vector( alpha, beta, minit1 )
    
    return outz, xax  


def run_tests():
    """
    Blah blah blah blah blah blah blah blah blah blah blah blah blah blah blah 
    blah blah blah blah  

    """
    
    my_sinc = hsinc(120,6)
    b1 = np.concatenate((np.zeros(4), my_sinc, np.zeros(4)))
    b1 = 0.5*b1/np.max(b1);
    g  = np.concatenate((np.zeros(4), 1*np.ones(120), np.zeros(4)))

    b1 = np.array(b1)
    if not np.iscomplexobj(b1):
        b1r = b1.real.copy()
        b1i = b1r.copy() * 0
        b1 = b1r + 1j*b1i
    b1 = b1 / 10.0

 
    # Example - on resonance, 2000 spatial locations

    T = 0.00002
    x = np.arange(-25,25,0.025)
    f = np.arange(-1000.0,2000.0,200.0)
    t = np.arange(1,len(b1)+1)*T

    # Get gradient field to do calcs over
    khz_to_mm           = 1.0
    do_extended_flag    = False
    convert_to_mm       = False
    resonance_offset    = 0.0
    dwell               = T * 1000000.0
    calc_resolution     = 2000
    
    vn,gx = grad_field(khz_to_mm,
                       do_extended_flag,
                       convert_to_mm,
                       resonance_offset,
                       dwell,
                       calc_resolution)

    # Run Bloch Simulation
    for i in range(100):
        outz, xax = bloch_twostep(b1, T, x, vn, gx)

    
#     mx = np.copy(outz[:,0])
#     my = np.copy(outz[:,1])
#     mz = np.copy(outz[:,2])
#        
#     
#     # Display Results in matplotlib 
#         
#     from pylab import *
#       
#     mxy = mx + 1j*my
#        
#     subplot(3,1,1)
#     xlabel('Time [ms]')
#     plot(t*1000,b1)
#        
#     subplot(3,1,2)
#     plot(xax, abs(mxy), xax, mx, xax, my )
#     xlabel('Position [cm]')
#     ylabel('Magnetization |Mxy|')
#            
#     subplot(3,1,3)
#     plot(xax, mz)
#     xlabel('Position [cm]')
#     ylabel('Magnetization Mz')
#            
#     show()
   
    

if __name__ == "__main__":
#    run_tests()

    import cProfile
    cProfile.run('run_tests()')

    # to time this call user
    #>python -m cProfile -s cumulative bloch_lib_matpulse.py

