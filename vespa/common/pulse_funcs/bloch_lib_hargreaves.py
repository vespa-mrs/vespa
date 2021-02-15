# Python imports


import numpy as np


# Derived from Brian Hargreaves bloch.c code for his Matlab mex library
#
#   Bloch simulator.
#
#   [mx,my,mz] = bloch(b1,gr,tp,t1,t2,df,dp,mode,mx,my,mz)
#
#   Bloch simulation of rotations due to B1, gradient and
#   off-resonance, including relaxation effects.  At each time
#   point, the rotation matrix and decay matrix are calculated.
#   Simulation can simulate the steady-state if the sequence
#   is applied repeatedly, or the magnetization starting at m0.
#
#   INPUT:
#       b1 = (Mx1) RF pulse in G.  Can be complex.
#       gr = (Mx1,2,or 3) 1,2 or 3-dimensional gradient in G/cm.
#       tp = (Mx1) time duration of each b1 and gr point, in seconds,
#               or 1x1 time step if constant for all points
#               or monotonically INCREASING endtime of each
#               interval..
#       t1 = T1 relaxation time in seconds.
#       t2 = T2 relaxation time in seconds.
#       df = (Nx1) Array of off-resonance frequencies (Hz)
#       dp = (Px1,2,or 3) Array of spatial positions (cm).  
#           Width should match width of gr.
#       mode= Bitmask mode:
#           Bit 0:  0-Simulate from start or M0, 1-Steady State
#           Bit 1:  1-Record m at time points.  0-just end time.
#
#   (optional)
#       mx,my,mz (PxN) arrays of starting magnetization, where N
#           is the number of frequencies and P is the number
#           of spatial positions.
#
#   OUTPUT:
#       mx,my,mz = PxN arrays of the resulting magnetization
#               components at each position and frequency.
#
#   B. Hargreaves.  Nov 2003.
#
#  
#  To profile line by line: 
#  >python kernprofile.py -l -v bloch_lib_hargreaves.py > statsline1.txt
#
#  also have to uncomment the @profile
#

TWOPI = 6.283185

GAMMA = 26751.3     # 1H  = 267.513 10^6 rad/(sec.Tesla) = 267.513 x 100 rad/(sec.gauss) = 26751.3 rad/(sec.gauss)
                    # 13C = 6726.2
                    # 17O = -3626.4
                    # 19F = 25166.2
                    # 23Na = 7076.1
                    # 31P  = 10829.1
                    # 129Xe = -73.997
                    
def isiterable(p_object):
    try:
        it = iter(p_object)
    except TypeError: 
        return False
    return True


def times2intervals( endtimes ):
    """
    Function takes the given endtimes of intervals, and
    returns the interval lengths in an array, assuming that
    the first interval starts at 0.

    If the intervals are all greater than 0, then this
    returns True, otherwise it returns False.

    """
    allpos    = True
    lasttime  = 0.0
    intervals = []
    for endtime in endtimes:
        intervals.append(endtime-lasttime)
        lasttime = endtime
        if intervals[-1] <= 0:
            allpos = False
            
    return allpos, np.array(intervals)


#@profile            # for use with line_profile/kernprofile.py
def calcrotmat(nx, ny, nz, rmat):
    """ 
    Find the rotation matrix that rotates |n| radians about 
        the vector given by nx,ny,nz               

    From:  https://code.google.com/p/robotics-toolbox-python/source/browse/trunk/robotics-toolbox-python/robot/transform.py
    Approach: Uses Matrices from numpy
    Turns out to be 1.8 times slower than for-loop with original straight math

    Found on: http://mail.scipy.org/pipermail/numpy-discussion/2009-March/040806.html
    Approach: Uses numpy array manipulations and skew array.
    Turns out to be 2.5 times slower than than for-loop with original straight math

    This final version takes advantage of ufunc speed of cos, sin, *, etc acting
    on an array of numbers. This was ~ 3-4 times faster than for-loop with math.

    Using for loop to check for 0.0 in phi_array was 2x slower than using an 
    index replace method.
       
    """
    phi_array = np.sqrt(nx*nx+ny*ny+nz*nz)
 
    # First define Cayley-Klein parameters     
    hp = phi_array/2.0
    cp = np.cos(hp)
    
    # NB. this line occasionally causes a div by 0 error to be thrown, but we
    #     fix any issues this causes below on line 152, so leave be for now.
    sp = np.sin(hp)/phi_array   # /phi because n is unit length in defs. 
    
    ar = cp
    ai = -nz*sp
    br =  ny*sp
    bi = -nx*sp
 
    # Make auxiliary variables to speed this up    
    arar  = ar*ar  
    aiai  = ai*ai
    arai2 = 2*ar*ai
    brbr  = br*br
    bibi  = bi*bi
    brbi2 = 2*br*bi
    arbi2 = 2*ar*bi
    aibr2 = 2*ai*br
    arbr2 = 2*ar*br
    aibi2 = 2*ai*bi
  
    # Make rotation matrix.  

    rmat[:,0,0] =  arar-aiai-brbr+bibi
    rmat[:,0,1] =  arai2-brbi2
    rmat[:,0,2] =  arbr2+aibi2
    rmat[:,1,0] =  -arai2-brbi2
    rmat[:,1,1] =  arar-aiai+brbr-bibi
    rmat[:,1,2] =  arbi2-aibr2
    rmat[:,2,0] =  -arbr2+aibi2
    rmat[:,2,1] =  -aibr2-arbi2
    rmat[:,2,2] =  arar+aiai-brbr-bibi 

    indx = np.where(phi_array == 0.0)
    rmat[indx,:,:] = np.array( [[1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,1.0]] )
      
#     rmat = np.array([[arar-aiai-brbr+bibi,  arai2-brbi2,         arbr2+aibi2],
#                      [-arai2-brbi2,         arar-aiai+brbr-bibi, arbi2-aibr2],
#                      [-arbr2+aibi2,        -aibr2-arbi2,         arar+aiai-brbr-bibi]])
#                             
#     rmat = rmat.transpose([2,0,1])
# 
#     for i, phi in enumerate(phi_array):
#         if phi == 0.0:
#             rmat[i,:,:] = np.array( [[1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,1.0]] )

    return rmat


#@profile            # for use with line_profile/kernprofile.py
def blochsim( b1real, b1imag, xgrad, ygrad, zgrad, 
              tsteps, e1, e2, df, 
              dx, dy, dz, mx, my, mz, mode, gamma):
    """
    Go through time for one df and all dx,dy,dz.    
    
    Calling Format
    
    mxx,myy,mzz = blochsim(b1real, b1imag, xgrad, ygrad, zgrad, 
                           tsteps, e1, e2, dfreq, 
                           dxpos, dypos, dzpos, 
                           mx, my, mz, mode);

    """
    npos = len(dx)
    nf   = len(df)

    bvecs   = [np.zeros((3,)) for i in range(nf*npos)]    
    decmat  = np.zeros((3,3))       # Decay matrix for each time step.
    decvec  = np.zeros((3,))        # Recovery vector for each time step.

    amats    = [np.eye(3) for i in range(nf*npos)]  # A is the identity matrix.
    imats    = [np.eye(3) for i in range(nf*npos)]  # I is the identity matrix.
    mcurr0s  = [np.array([mx[j,i,0],my[j,i,0],mz[j,i,0]]) for j in range(nf) for i in range(npos)] # Set starting x,y,z magnetizations
    rotmats  = np.zeros([nf*npos,3,3],np.float32)
    
    for t in range(len(tsteps)):

        #  Rotation    

        df_array = np.repeat(df*TWOPI*tsteps[t], npos)

        rotz = -(xgrad[t] * dx + ygrad[t] * dy + zgrad[t] * dz) * tsteps[t]
        rotz = np.tile(rotz, nf) - df_array
        rotx = (- b1real[t] * gamma * tsteps[t])
        roty = (- b1imag[t] * gamma * tsteps[t]) # based on Hao Sun's UMich blochCim code

        if t == 83:
            bob = 10

        rotmats = calcrotmat(rotx, roty, rotz, rotmats)
        
        if (mode == 1):
            #
            # FIXME - bjs, may need to swap these out for the same code as in 
            #         the else statement below
            #
            arots = [np.dot(rotmat, amat) for rotmat, amat in zip(rotmats,amats)] 
            brots = [np.dot(rotmat, bvec) for rotmat, bvec in zip(rotmats,bvecs)]
        else:
            #mcurr1s = [np.dot(rotmat, mcurr0) for rotmat, mcurr0 in zip(rotmats,mcurr0s)]
            # this effectively does a dot product on lists of 3x3 and 3x1 data
            mcurr1s = np.einsum('...ij,...j->...i',rotmats,mcurr0s)

        #  Decay   
        if e1 is not None:
            decvec[2]   = 1-e1[t]
            decmat[2,2] = e1[t]
        if e2 is not None:
            decmat[0,0] = e2[t]
            decmat[1,1] = e2[t]
        
        if e1 is not None and e2 is not None:
            
            if (mode == 1):
                amats = [np.dot(decmat, arot) for arot in arots]
                bvecs = [np.dot(decmat, brot) for brot in brots]
                bvecs = [bvec+decvec for bvec in bvecs]
            else:
                mcurr0s = [np.dot(decmat, mcurr1) for mcurr1 in mcurr1s]
                mcurr0s = [mcurr0+decvec for mcurr0 in mcurr0s]
        else:
            mcurr0s = mcurr1s
                
    
        if mode == 2:      
            # Sample output at times. Only do this if transient!
            mcurr0 = np.array(mcurr0s)
            mcurr0.shape = nf, npos, 3
            mx[:,:,t] = mcurr0[:,:,0]
            my[:,:,t] = mcurr0[:,:,1]
            mz[:,:,t] = mcurr0[:,:,2]

    # If only recording the endpoint, either store the last
    # point, or calculate the steady-state endpoint. 

    if mode == 0:        
        # Indicates start at given m, or m0.
        mcurr0 = np.array(mcurr0s)
        mcurr0.shape = nf, npos, 3
        mx[:,:,0] = mcurr0[:,:,0]
        my[:,:,0] = mcurr0[:,:,1]
        mz[:,:,0] = mcurr0[:,:,2]
        
    elif mode == 1:         
        # Indicates to find steady-state magnetization 
        amats = [imat-amat for imat,amat in zip(imats,amats)]           # Now amat = (I-A) 
        imats = [np.linalg.inv(amat) for amat in amats]                 # Reuse imat as inv(I-A)   
        mvec  = [np.dot(imat,bvec) for imat,bvec in zip(imats,bvecs)]   # Now M = inv(I-A)*B       
        mvec  = np.array(mvec)
        mvec.shape = nf, npos, 3
        mx[:,:,0] = mvec[:,:,0]
        my[:,:,0] = mvec[:,:,1]
        mz[:,:,0] = mvec[:,:,2]
        
    return mx, my, mz



def blochsimfz( b1real, b1imag, xgrad, ygrad, zgrad,
                tsteps, t1, t2, dfreq, nfreq,
                dxpos, dypos, dzpos, npos,
                mx, my, mz, mode, gamma):
    """
    Comment?
    """
    if mode & 2:
        ntout = len(tsteps)
    else:
        ntout = 1

    mxout = mx.copy()
    myout = my.copy()
    mzout = mz.copy()

    # First calculate the e1 and e2 values at each time step. 
    if t1 is not None:
        e1 = np.exp(-tsteps/t1)
    else:
        e1 = None
    if t2 is not None:
        e2 = np.exp(-tsteps/t2)
    else:
        e2 = None

    gammadx = dxpos.copy()*gamma          # Convert to Hz/cm 
    gammady = dypos.copy()*gamma          # Convert to Hz/cm 
    gammadz = dzpos.copy()*gamma          # Convert to Hz/cm 

    #for ifreq in range(nfreq):

    if mode == 3:       # Steady state AND record all time points.  

        # First go through and find steady state, then
        # repeat as if transient starting at steady st.

        mxx, myy, mzz = blochsim(b1real, b1imag, xgrad, ygrad, zgrad, 
                                 tsteps, e1, e2, dfreq, 
                                 gammadx, gammady, gammadz, 
                                 mx, my, mz, 1, gamma);

        mxx, myy, mzz = blochsim(b1real, b1imag, xgrad, ygrad, zgrad, 
                                 tsteps, e1, e2, dfreq, 
                                 gammadx, gammady, gammadz, 
                                 mxx, myy, mzz, 2, gamma);


    else:

        mxx, myy, mzz = blochsim(b1real, b1imag, xgrad, ygrad, zgrad, 
                                 tsteps, e1, e2, dfreq, 
                                 gammadx, gammady, gammadz, 
                                 mx, my, mz, mode, gamma);
    mxout[:,:,:] = mxx
    myout[:,:,:] = myy
    mzout[:,:,:] = mzz

    return mxout, myout, mzout



def bloch(b1,gr,tp,t1,t2,df,dp,mode=0,mx=[],my=[],mz=[], do_mxy=False, gamma=GAMMA):
    """ 
    Calling format 
    
    [mx,my,mz] = bloch(b1,gr,tp,t1,t2,df,dp,mode=0,mx=[],my=[],mz=[])

    blochsimfz(b1r,b1i,gx,gy,gz,tp,t1,t2,df,nf,dx,dy,dz,npos,mx,my,mz,md);
     
    Bloch simulation of rotations due to B1, gradient and
    off-resonance, including relaxation effects.  At each time
    point, the rotation matrix and decay matrix are calculated.
    Simulation can simulate the steady-state if the sequence
    is applied repeatedly, or the magnetization starting at m0.
 
    INPUT:
        b1 = (Mx1) [Gauss] RF pulse.  Can be complex.
        gr = (Mx1,2,or 3) [G/cm] 1,2 or 3-dimensional gradient.
        tp = (Mx1) [sec] time duration of each b1 and gr point,
                or 1x1 time step if constant for all points
                or monotonically INCREASING endtime of each
                interval..
        t1 = T1 relaxation time [sec].
        t2 = T2 relaxation time [sec].
        df = (Nx1) Array of off-resonance frequencies [Hz]
        dp = (Px1,2,or 3) Array of spatial positions [cm].  
            Width should match width of gr.
        mode= Bitmask mode:
            Bit 0:  0-Simulate from start or M0, 1-Steady State
            Bit 1:  0-Just end time, 1-Record m at time points.
 
    (optional) - NB. bjs, swapped N anp P parameters here versus Matlab code
    
        mx,my,mz (NxP) arrays of starting magnetization, where N
            is the number of frequencies and P is the number
            of spatial positions.
 
    OUTPUT:
        mx,my,mz = NxP arrays of the resulting magnetization
                components at each position and frequency.
 
    B. Hargreaves.  Nov 2003.

    """
    ntime = len(b1)  # Number of Time, RF, and Grad points  

    # ====================== RF (B1) =========================
    # :  If complex, split up.  If real, allocate an imaginary part. 
    b1 = np.array(b1)
    if np.iscomplexobj(b1):
        b1r = b1.real.copy()
        b1i = b1.imag.copy()
    else:
        b1r = b1.real.copy()
        b1i = b1r.copy() * 0

    # ======================= Gradients ========================= 
    gr = np.array(gr)
    ngrad = gr.size                         # Number of Time, RF, and Grad points  
    gx = np.zeros((ntime,),'float')
    gy = np.zeros((ntime,),'float')
    gz = np.zeros((ntime,),'float')

    gx[0:ntime] = gr[0:ntime]               # X-gradient is first N points.  

    if ngrad >= 2*ntime:                    # Need to allocate Y-gradient. 
        gy[0:ntime] = gr[ntime:ntime*2]     # Assign from Nx3 input array. Assuming (at least) 2-Dimensional Gradient

    if (ngrad >= 3*ntime):                  # Need to allocate Z-gradient. 
        gz[0:ntime] = gr[ntime*2:ntime*3]   # Assign from Nx3 input array.  Assuming 3-Dimensional Gradient

    # Warning if Gradient length is not 1x, 2x, or 3x RF length. 
    if (ngrad != ntime) and (ngrad != 2*ntime) and (ngrad != 3*ntime):
        raise ValueError("Bloch: Gradient length differs from B1 length.")


    # === Time points ===== 
    #
    #  THREE Cases:
    #        1) Single value given -> this is the interval length for all.
    #        2) List of intervals given.
    #        3) Monotonically INCREASING list of end times given.
    #
    #    For all cases, the goal is for tp to have the intervals.
    #
    if isinstance(tp,float):                # === Case 1 === 
        tstep = tp
        tp = np.zeros((ntime,),'float') + tstep

    elif len(tp) == 1:                # === Case 1 === 
        tstep = tp
        tp = np.zeros((ntime,),'float') + tstep

    elif len(tp) != ntime:
        raise ValueError("Time-point length differs from B1 length")

    else:
        tp = np.array(tp)
        posflag, tp = times2intervals( tp )
        if posflag:
#            print "Times are monotonically increasing. "
            pass

    # === Relaxation Times ===== 
    if t1 is not None:
        t1 = float(t1)
    if t2 is not None:
        t2 = float(t2)

    # === Frequency Points ===== 
    if isiterable(df):
        df = np.array(df)
    else:
        df = np.array([df])
    nf = len(df)

    # === Position Points ===== 
    if isiterable(dp):
        dp = np.array(dp)
    else:
        dp = np.array([dp])
        
    if len(dp.shape) == 1:
        dp.shape = 1,dp.shape[0]
    nposN, nposM = dp.shape

    npos = nposM
    dx = np.zeros((npos,), 'float')
    dy = np.zeros((npos,), 'float')
    dz = np.zeros((npos,), 'float')

    if (nposN==3):                      # Assume 3 position dimensions given 
        dx[0:npos] = dp[0]
        dy[0:npos] = dp[1]
        dz[0:npos] = dp[2]

    elif (nposN==2):                    # Assume only 2 position dimensions given 
        dx[0:npos] = dp[0]
        dy[0:npos] = dp[1]

    else:              
        dx[0:npos] = dp[0]

    nfnpos = nf*npos;   # Just used to speed things up below.  

    # ===== Mode, defaults to 0 (simulate single endpoint, transient). ==== 
    md = int(mode)

    if (md & 2):
        ntout = ntime       # Include time points. 
    else:
        ntout = 1

    ntnfnpos = ntout*nfnpos;

#     if (md & 1)==0:
#         print "Simulation from Initial Condition."
#     else:
#         print "Simulation of Steady-State."
#     
#     if (md & 2)==0:
#         print "Simulation to Endpoint. "
#     else:
#         print "Simulation over Time."


    # ===== Allocate Output Magnetization vectors arrays.  
    mxin = np.zeros((nf, npos, ntout), 'float')
    myin = np.zeros((nf, npos, ntout), 'float')
    mzin = np.zeros((nf, npos, ntout), 'float')

    # ===== If Initial Magnetization is given... 
    if mx and my and mz and len(mx)==nfnpos and len(my)==nfnpos and len(mz)==nfnpos:
        # Set output magnetization to that passed.
        # If multiple time points, then just the first is set.               
        # print "Using Specified Initial Magnetization."
        for ipos in range(npos):
            for ifreq in range(nf):
                mxin[ifreq,ipos,0] = mx[ifreq,ipos]
                myin[ifreq,ipos,0] = my[ifreq,ipos]
                mzin[ifreq,ipos,0] = mz[ifreq,ipos]
    else:
        # if mx and my and mz:   # Magnetization given, but wrong size! 
        #     print "Initial magnetization passed, but not Npositions x Nfreq. "
        #     print " --> Using [0; 0; 1] for initial magnetization. "

        if do_mxy:
            for ipos in range(npos):
                for ifreq in range(nf):
                    mxin[ifreq,ipos,0] = 0
                    myin[ifreq,ipos,0] = 1
                    mzin[ifreq,ipos,0] = 0            
        else:
            for ipos in range(npos):
                for ifreq in range(nf):
                    mxin[ifreq,ipos,0] = 0
                    myin[ifreq,ipos,0] = 0
                    mzin[ifreq,ipos,0] = 1

    # ======= Do The Simulation! ====== 
    mxout, myout, mzout = blochsimfz(b1r,b1i,gx,gy,gz,tp,t1,t2,df,nf,dx,dy,dz,npos,mxin,myin,mzin,md,gamma)

    # ======= Reshape Output Matrices ====== 
    if (ntout > 1) and (nf > 1) and (npos > 1):
        outsize = nf, npos, ntout
    else:                   # Basically "squeeze" the matrix. 
        if nf <= 1:
            outsize = npos, ntout
        else:
            if npos <= 1:
                outsize = nf, ntout
            else:
                outsize = nf, npos, ntout

    mxout.shape = outsize
    myout.shape = outsize
    mzout.shape = outsize
    
    return mxout, myout, mzout



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

def run_tests():
    """
    Blah blah blah blah blah blah blah blah blah blah blah blah blah blah blah 
    blah blah blah blah  

    """
    
    # Set up B1 and gradient
    
    my_sinc = hsinc(120,6)
    b1 = np.concatenate((np.zeros(4), my_sinc, np.zeros(4)))
    b1 = 0.5*b1/np.max(b1)
    g  = np.concatenate((np.zeros(4), 0.25*np.ones(120), np.zeros(4)))

    # Example - on resonance, 2000 spatial locations
 
    T  = 0.00002
    x  = np.arange(-25,25,0.025)
    f  = np.array([0.0,])            # on resonance
    t  = np.arange(1,len(b1)+1)*T;

#     # Example - on/off resonance, 2000 spatial locations
# 
#     T = 0.00002
#     x = np.arange(-25,25,0.025)
#     f = np.arange(-1000.0,2000.0,200.0)
#     t = np.arange(1,len(b1)+1)*T;

    # Run Bloch Simulation
    #
    # Input units: Gauss, G/cm, sec, sec, sec, Hz, cm
    
    # mx, my, mz = bloch(b1,g,t,1.0,0.2,f,x,mode=0)     # this has T1/T2 decays
    # mx, my, mz = bloch(b1,g,t,None,None,f,x,mode=0)     # no decays
    # this is for testing time spent in functions with cProfile or kernprofile.py
    for i in range(1):
        mx, my, mz = bloch(b1,g,t,None,None,f,x,mode=0)     # no decays


    # Display Results in matplotlib - for on resonance freq only
    
    from pylab import subplot, xlabel, plot, ylabel, show
        
    mxy = mx + 1j*my
    ioff = 1
    subplot(3,1,1)
    xlabel('Time [ms]')
    plot(t*1000,b1)
    subplot(3,1,2)
    plot(x, abs(mxy[:,0]), x, real(mxy[:,0]), x, imag(mxy[:,0]) )
    xlabel('Position [cm]')
    ylabel('Magnetization |Mxy|')
    subplot(3,1,3)
    plot(x, mz[:,0])
    xlabel('Position [cm]')
    ylabel('Magnetization Mz')
    show()
 
#     # Display Results in matplotlib - for multiple frequencies
#  
#     from pylab import *
#      
#     mxy = mx + 1j*my
#     ioff = 1
#     subplot(3,1,1)
#     xlabel('Time [ms]')
#     plot(t*1000,b1)
#     subplot(3,1,2)
#     plot(x, abs(mxy[ioff,:]), x, real(mxy[ioff,:]), x, imag(mxy[ioff,:]) )
#     xlabel('Position [cm]')
#     ylabel('Magnetization |Mxy|')
#     subplot(3,1,3)
#     plot(x, mz[ioff,:])
#     xlabel('Position [cm]')
#     ylabel('Magnetization Mz')
#     show()    

if __name__ == "__main__":
    run_tests()

#     import cProfile
#     cProfile.run('run_tests()')
 
    # to time this call user
    #>python -m cProfile -s cumulative bloch_lib_hargreaves.py






