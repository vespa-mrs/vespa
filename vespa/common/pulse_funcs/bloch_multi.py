# Python modules

import multiprocessing
import math

# 3rd party modules
import numpy as np

# Our modules

from pylab import *

# GAMMA = 26753.0 - replaced with user defined value
TWOPI = 6.283185


def isiterable(p_object):
    try:
        it = iter(p_object)
    except TypeError: 
        return False
    return True


"""
Code for calculating bloch simulations. This code is independent of any GUI
and can be called from the command line. There's only one public function
which is bloch_multi().

This code iterates over each position and frequency offset in the simulation. 
For each, it returns the magnetization values (Mx,My,Mz).

It uses multiprocessing.

"""

# The naive algorithm for this code would be to process one simulation at
# a time. We found that experiments are processed maybe 10% faster if we 
# combine the bloch simulation's frequency offset dimension into chunks 
# prior to processing. 
#
# MAX_CHUNK_SIZE is the largest # of freq offsets that will be grouped together
# in one chunk. This is a hard limit.
#
# Choosing MAX_CHUNK_SIZE was the result of a lot of experimentation. 
# Obviously we can't test all possibilities but 4-10 gave good average 
# performance. 
# 
# Because the x,y,z position dimension is always fully included in each chunk
# the time saved by adding more than 1 freq offset line to the chunk is much
# less than expected. Thus for a 1000x200 nf x npos array, we went from 55 sec
# on 12 processors to 45 seconds on 12 processors when we set MAX_CHUNK_SIZE
# from 1 to 10 respectively.  At 50 or 100 per chunk, the total time approached
# 55 seconds again.  So, nothing too simple to figure out here, just empirical
# mysto-crud.
#
# To force this code to use the naive algorithm (each simulation = one chunk), 
# set MAX_CHUNK_SIZE = 1.
MAX_CHUNK_SIZE = 4


class _Chunk(object):
    """
    A chunk of frequency offset values to process. len(df) is 
    always < MAX_CHUNK_SIZE.
    
    b1real, b1imag, xgrad, ygrad, zgrad, tsteps, e1, e2, df, dx, dy, dz, mx, my, mz, mode, gamma
    
    """ 
    def __init__(self, b1real, b1imag, xgrad, ygrad, zgrad, tsteps, 
                       e1, e2, dx, dy, dz, mx, my, mz, mode, gamma):
        # Each chunk has a bunch of static info that is used in all 
        # calculations. And a small set of frequency offset values and their
        # indices for storage back into the results arrays on finish. 
        self.b1real = b1real
        self.b1imag = b1imag
        self.xgrad  = xgrad
        self.ygrad  = ygrad
        self.zgrad  = zgrad
        self.tsteps = tsteps
        self.e1     = e1
        self.e2     = e2
        self.df     = []
        self.ifreq  = []
        self.dx     = dx
        self.dy     = dy
        self.dz     = dz
        self.mx     = mx
        self.my     = my
        self.mz     = mz
        self.mode   = mode
        self.gamma  = gamma

    @property
    def nlines(self):
        return len(self.df)
        
    def __str__(self):
        # __str__ is useful for debugging
        lines = [ ]
        lines.append("---- chunk ----")
        lines.append("b1real: %d" % self.b1real.size)
        lines.append("b1imag: %d" % self.b1imag.size)
        lines.append("xgrad: %d" % self.xgrad.size)
        lines.append("ygrad: %d" % self.ygrad.size)
        lines.append("zgrad: %d" % self.zgrad.size)
        lines.append("tsteps: %d" % self.tsteps.size)
        lines.append("e1: %d" % self.e1.size)
        lines.append("e2: %d" % self.e2.size)
        lines.append("df: %d" % self.df.size)
        lines.append("dx: %d" % self.dx.size)
        lines.append("dy: %d" % self.dy.size)
        lines.append("dz: %d" % self.dz.size)
        lines.append("mx: %d" % self.mx.size)
        lines.append("my: %d" % self.my.size)
        lines.append("mz: %d" % self.mz.size)
        lines.append("mode: %d" % self.mode)
        lines.append("gamma: %f" % self.gamma)
        
        return "\n".join(lines)
        


def blochsim(b1real, b1imag, xgrad, ygrad, zgrad, 
             tsteps, ntime, e1, e2, df, 
             dx, dy, dz, mx, my, mz, mode, 
             cpu_count=1, gamma=26751.3 ):
    """
    Builds a set of Magnetization values for spatial locations (dx,dy,dz) and
    frequency offset values (df) given a B1 pulse and set of x,y,z gradient
    values. It runs a bloch simulation at each of nf x npos frequencies and 
    locations. It uses multiprocessing. 
    
    Returns a numpy array appropriate for the magnetization vector at each 
    location and frequency offset.

    """
    mxout = mx.copy()
    myout = my.copy()
    mzout = mz.copy()
    
    # PS - If you want to run this code without using multiprocessing (e.g.
    # in order to profile execution), use the 3 lines below in place of
    # the use of multiprocessing.Pool.
#     _initializer()
#     chunks = _build_chunks(b1real, b1imag, xgrad, ygrad, zgrad, tsteps, e1, e2, df, dx, dy, dz, mx, my, mz, mode, gamma )
#     results = [_process_chunk(chunk) for chunk in chunks]

    pool = multiprocessing.Pool(cpu_count, _initializer, [])
  
      
    # The 3rd param to imap_unordered() is a chunksize. These chunks are not
    # to be confused with the chunks returned by _build_chunks()! chunksize
    # just determines how many values will be grabbed from the iterator 
    # at once. Using a chunksize > 1 gives slightly better performance, but
    # only slightly. The advantage of using a chunksize == 1 is that 
    # _build_chunks() is called every time a worker needs new work, so we
    # can use it as a cheap callback/progress indicator.
    results = pool.imap_unordered(_process_chunk, 
                                  _build_chunks(b1real, b1imag, xgrad, ygrad, zgrad, tsteps, 
                                                e1, e2, df, dx, dy, dz, mx, my, mz, mode, gamma ), 
                                  1)
      
    pool.close()
    pool.join()
    
    # The lines from each bloch simulation are combined into one array that
    # has results from all nf x npos frequency offsets and spatial positions

    for result in results:

        ifreqs = result[0][0]
        mxs    = result[0][1]
        mys    = result[0][2]
        mzs    = result[0][3]

        for i, ifreq in enumerate(ifreqs):

            mxout[ifreq,:,:] = mxs[i,:,:]
            myout[ifreq,:,:] = mys[i,:,:]
            mzout[ifreq,:,:] = mzs[i,:,:]
    
    return mxout, myout, mzout


###############       Internal use only below this line     ###############

def _build_chunks(b1real, b1imag, xgrad, ygrad, zgrad, tsteps, 
                  e1, e2, df, dx, dy, dz, mx, my, mz, mode, gamma ):
    """
    A generator function. Given an experiment, iterates over the bloch 
    simulation's frequency offset dimension and returns a set of offsets
    chunked according to MAX_CHUNK_SIZE.
    
    See here for more info on generators:
    http://docs.python.org/tutorial/classes.html#generators
    """
    current = _Chunk(b1real, b1imag, xgrad, ygrad, zgrad, tsteps, 
                     e1, e2, dx, dy, dz, mx, my, mz, mode, gamma)
    
    nlines_processed = 0

    for ifreq, dfreq in enumerate(df):
        nlines = 1
        
        if current.nlines and ((current.nlines + nlines) > MAX_CHUNK_SIZE):
            # Chunk has enough in it, adding more would exceed the max.
            nlines_processed += current.nlines
            
            yield current
            current = _Chunk(b1real, b1imag, xgrad, ygrad, zgrad, tsteps, 
                             e1, e2, dx, dy, dz, mx, my, mz, mode, gamma)
        #else:
            # The current chunk is empty or there's still room in the current
            # chunk for the next collection of lines. 

        # Append the contents of this simulation to the current chunk.
        current.df.append(dfreq)
        current.ifreq.append(ifreq)
        

    # Return the final set of lines.
    yield current


def _initializer():
    # This function is subtle...it's called by each worker process, and is
    # passed the values of the global constants that I need in 
    # _process_chunk(). Under *nix, I can just declare them global and 
    # (thanks to the magic of fork()) the variables and their values will be 
    # copied to the worker processes'. Under Windows, this module is 
    # re-imported once for each worker, and as a result these globals are
    # recreated and re-initialized to 0 in each worker. This function sets 
    # them back to the values they need to be, and that's the only reason
    # it exists.
    pass

    
def _process_chunk(chunk):
    # This is what each worker executes
    
    b1real = chunk.b1real 
    b1imag = chunk.b1imag 
    xgrad  = chunk.xgrad  
    ygrad  = chunk.ygrad  
    zgrad  = chunk.zgrad  
    tsteps = chunk.tsteps 
    e1     = chunk.e1     
    e2     = chunk.e2     
    df     = np.array(chunk.df)
    ifreq  = chunk.ifreq  
    dx     = chunk.dx     
    dy     = chunk.dy     
    dz     = chunk.dz     
    mx     = chunk.mx     
    my     = chunk.my     
    mz     = chunk.mz     
    mode   = chunk.mode  
    gamma  = chunk.gamma

    mxout = np.zeros((len(df),mx.shape[1],mx.shape[2]), 'float')
    myout = np.zeros((len(df),mx.shape[1],mx.shape[2]), 'float')
    mzout = np.zeros((len(df),mx.shape[1],mx.shape[2]), 'float')

    npos = len(dx)
    nf   = len(df)

    bvecs   = [np.zeros((3,)) for i in range(nf*npos)]    
    decmat  = np.zeros((3,3))       # Decay matrix for each time step.
    decvec  = np.zeros((3,))        # Recovery vector for each time step.

    amats    = [np.eye(3) for i in range(nf*npos)]  # A is the identity matrix.
    imats    = [np.eye(3) for i in range(nf*npos)]  # I is the identity matrix.
    mcurr0s  = [np.array([mx[j,i,0],my[j,i,0],mz[j,i,0]]) for j in range(nf) for i in range(npos)] # Set starting x,y,z magnetizations

    for t in range(len(tsteps)):

        #  Rotation    

        df_array = np.repeat(df*TWOPI*tsteps[t], npos)

        rotz = -(xgrad[t] * dx + ygrad[t] * dy + zgrad[t] * dz) * tsteps[t]
        rotz = np.tile(rotz, nf) - df_array
        rotx = (- b1real[t] * gamma * tsteps[t])
        roty = (- b1imag[t] * gamma * tsteps[t]) # based on Hao Sun's UMich blochCim code

        rotmats = calcrotmat(rotx, roty, rotz)

        if (mode == 1):
            arots = [np.dot(rotmat, amat) for rotmat, amat in zip(rotmats,amats)] 
            brots = [np.dot(rotmat, bvec) for rotmat, bvec in zip(rotmats,bvecs)]
        else:
            mcurr1s = [np.dot(rotmat, mcurr0) for rotmat, mcurr0 in zip(rotmats,mcurr0s)]

        #  Decay   
        decvec[2]   = 1-e1[t]
        decmat[0,0] = e2[t]
        decmat[1,1] = e2[t]
        decmat[2,2] = e1[t]
        
        if (mode == 1):
            amats = [np.dot(decmat, arot) for arot in arots]
            bvecs = [np.dot(decmat, brot) for brot in brots]
            bvecs = [bvec+decvec for bvec in bvecs]
        else:
            mcurr0s = [np.dot(decmat, mcurr1) for mcurr1 in mcurr1s]
            mcurr0s = [mcurr0+decvec for mcurr0 in mcurr0s]
            

        if mode == 2:      
            # Sample output at times. Only do this if transient!
            mcurr0 = np.array(mcurr0s)
            mcurr0.shape = nf, npos, 3
            mxout[:,:,t] = mcurr0[:,:,0]
            myout[:,:,t] = mcurr0[:,:,1]
            mzout[:,:,t] = mcurr0[:,:,2]

    # If only recording the endpoint, either store the last
    # point, or calculate the steady-state endpoint. 

    if mode == 0:        
        # Indicates start at given m, or m0.
        mcurr0 = np.array(mcurr0s)
        mcurr0.shape = nf, npos, 3
        mxout[:,:,0] = mcurr0[:,:,0]
        myout[:,:,0] = mcurr0[:,:,1]
        mzout[:,:,0] = mcurr0[:,:,2]
        
    elif mode == 1:         
        # Indicates to find steady-state magnetization 
        amats = [imat-amat for imat,amat in zip(imats,amats)]           # Now amat = (I-A) 
        imats = [np.linalg.inv(amat) for amat in amats]                 # Reuse imat as inv(I-A)   
        mvec  = [np.dot(imat,bvec) for imat,bvec in zip(imats,bvecs)]   # Now M = inv(I-A)*B       
        mvec  = np.array(mvec)
        mvec.shape = nf, npos, 3
        mxout[:,:,0] = mvec[:,:,0]
        myout[:,:,0] = mvec[:,:,1]
        mzout[:,:,0] = mvec[:,:,2]
        


    # The results are a list of 2-tuples (index, Mx, My, Mz). index is an index 
    # into the frequency offset dimension of the magnetization array -- it's 
    # where these mx, my, mz values will reside in the overall results array.
    result = [ (ifreq, mxout, myout, mzout) ]
    
    return result




#==============================================================================



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


def calcrotmat(nx, ny, nz):
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
       
    """
    phi_array  = np.sqrt(nx*nx+ny*ny+nz*nz)

    rmat = []

    # First define Cayley-Klein parameters     
    hp = phi_array/2.0
    cp = np.cos(hp)
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
     
    rmat = np.array([[arar-aiai-brbr+bibi,  arai2-brbi2,         arbr2+aibi2],
                     [-arai2-brbi2,         arar-aiai+brbr-bibi, arbi2-aibr2],
                     [-arbr2+aibi2,        -aibr2-arbi2,         arar+aiai-brbr-bibi]])
                          
    rmat = rmat.transpose([2,0,1])

    for i, phi in enumerate(phi_array):

        if phi == 0.0:
            rmat[i,:,:] = np.array( [[1.0,0.0,0.0],
                                     [0.0,1.0,0.0],
                                     [0.0,0.0,1.0]] )
    return rmat




def blochsimfz(b1real, b1imag, xgrad, ygrad, zgrad,
               tsteps, ntime, t1, t2, dfreq, nfreq,
               dxpos, dypos, dzpos, npos,
               mx, my, mz, mode, cpu_count=0, gamma=26751.3):
    """
    Comment?
    """
    if mode & 2:
        ntout = ntime
    else:
        ntout = 1

    mxout = mx.copy()
    myout = my.copy()
    mzout = mz.copy()

    # First calculate the e1 and e2 values at each time step. 

    e1 = np.exp(-tsteps/t1)
    e2 = np.exp(-tsteps/t2)

    gammadx = dxpos.copy()*gamma          # Convert to Hz/cm 
    gammady = dypos.copy()*gamma          # Convert to Hz/cm 
    gammadz = dzpos.copy()*gamma          # Convert to Hz/cm 

    #for ifreq in range(nfreq):

    if mode == 3:       # Steady state AND record all time points.  

        # First go through and find steady state, then
        # repeat as if transient starting at steady st.

        mxx, myy, mzz = blochsim(b1real, b1imag, xgrad, ygrad, zgrad, 
                                 tsteps, ntime, e1, e2, dfreq, 
                                 gammadx, gammady, gammadz, 
                                 mx, my, mz, 1, cpu_count, gamma);

        mxx, myy, mzz = blochsim(b1real, b1imag, xgrad, ygrad, zgrad, 
                                 tsteps, ntime, e1, e2, dfreq, 
                                 gammadx, gammady, gammadz, 
                                 mxx, myy, mzz, 2, cpu_count, gamma);


    else:

        mxx, myy, mzz = blochsim(b1real, b1imag, xgrad, ygrad, zgrad, 
                                 tsteps, ntime, e1, e2, dfreq, 
                                 gammadx, gammady, gammadz, 
                                 mx, my, mz, mode, cpu_count, gamma);
    mxout[:,:,:] = mxx
    myout[:,:,:] = myy
    mzout[:,:,:] = mzz

    return mxout, myout, mzout



def bloch_multi(b1,gr,tp,t1,t2,df,dp,mode=0,mx=[],my=[],mz=[], cpu_count=0, gamma=26751.3):
    """ 
    Calling format 
    
    [mx,my,mz] = bloch(b1,gr,tp,t1,t2,df,dp,mode=0,mx=[],my=[],mz=[])

    blochsimfz(b1r,b1i,gx,gy,gz,tp,ntime,t1,t2,df,nf,dx,dy,dz,npos,mx,my,mz,md);
     
    Bloch simulation of rotations due to B1, gradient and
    off-resonance, including relaxation effects.  At each time
    point, the rotation matrix and decay matrix are calculated.
    Simulation can simulate the steady-state if the sequence
    is applied repeatedly, or the magnetization starting at m0.
 
    INPUT:
        b1 = (Mx1) RF pulse in G.  Can be complex.
        gr = (Mx1,2,or 3) 1,2 or 3-dimensional gradient in G/cm.
        tp = (Mx1) time duration of each b1 and gr point, in seconds,
                or 1x1 time step if constant for all points
                or monotonically INCREASING endtime of each
                interval..
        t1 = T1 relaxation time in seconds.
        t2 = T2 relaxation time in seconds.
        df = (Nx1) Array of off-resonance frequencies (Hz)
        dp = (Px1,2,or 3) Array of spatial positions (cm).  
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
    # cpu_count is the number of processing cores (virtual CPUs) available on 
    # this machine. We ask multiprocessing.Pool() to create CPU_COUNT workers.
    # cpu_count can be determined from a variety of sources. We accept any 
    # int > 0.

    if not cpu_count:
        # OK, the user didn't specify so we ask multiprocessing. Note that 
        # multiprocessing.cpu_count() isn't implemented on all platforms.
        # Where it's not implemented we default to 2 for no really strong reasons.
        try:
            cpu_count = multiprocessing.cpu_count()
        except NotImplementedError:
            cpu_count = 2




    
    print("---------------------------------------")
    print("3D-position + frequency Bloch Simulator")
    print("---------------------------------------")

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
        print("Gradient length differs from B1 length")



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
        print("Time-point length differs from B1 length")

    else:
        tp = np.array(tp)
        posflag, tp = times2intervals( tp )
        if posflag:
            print("Times are monotonically increasing. ")

    # === Relaxation Times ===== 
    t1 = float(t1)
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

    if (md & 1)==0:
        print("Simulation from Initial Condition.")
    else:
        print("Simulation of Steady-State.")

    if (md & 2)==0:
        print("Simulation to Endpoint. ")
    else:
        print("Simulation over Time.")


    # ===== Allocate Output Magnetization vectors arrays.  
    mxin = np.zeros((nf, npos, ntout), 'float')
    myin = np.zeros((nf, npos, ntout), 'float')
    mzin = np.zeros((nf, npos, ntout), 'float')

    # ===== If Initial Magnetization is given... 
    if mx and my and mz and len(mx)==nfnpos and len(my)==nfnpos and len(mz)==nfnpos:
        # Set output magnetization to that passed.
        # If multiple time points, then just the first is set.               
        print("Using Specified Initial Magnetization.")
        for ipos in range(npos):
            for ifreq in range(nf):
                mxin[ifreq,ipos,0] = mx[ifreq,ipos]
                myin[ifreq,ipos,0] = my[ifreq,ipos]
                mzin[ifreq,ipos,0] = mz[ifreq,ipos]
    else:
        if mx and my and mz:   # Magnetization given, but wrong size! 
            print("Initial magnetization passed, but not Npositions x Nfreq. ")
            print(" --> Using [0; 0; 1] for initial magnetization. ")

        for ipos in range(npos):
            for ifreq in range(nf):
                mxin[ifreq,ipos,0] = 0
                myin[ifreq,ipos,0] = 0
                mzin[ifreq,ipos,0] = 1

    # ======= Do The Simulation! ====== 
    
    print("Calling blochsimfz_par() function.")
    mxout, myout, mzout = blochsimfz(b1r, b1i, x, gy, gz,
                                     tp, ntime, t1, t2,
                                     df, nf, dx, dy, dz,
                                     npos, mxin, myin, mzin,
                                     md, cpu_count, gamma=gamma)

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


def main():
    """
    Blah blah blah blah blah blah blah blah blah blah blah blah blah blah blah 
    blah blah blah blah  

    """
    gamma = 26751.3   # 1H Hz/gauss 
    
    my_sinc = hsinc(120,6)
 
    T  = 0.00002
    b1 = np.concatenate((np.zeros(4), my_sinc, np.zeros(4)))
    g  = np.concatenate((np.zeros(4), 1*np.ones(120), np.zeros(4)))
    b1 = 0.5*b1/np.max(b1);
    x  = np.arange(-5,5,0.05)
    f  = np.arange(-1000.0,2000.0,200.0)
    t  = np.arange(1,len(b1)+1)*T;
 
    mx, my, mz = bloch_multi(b1,g,t,1,.2,f,x,mode=0, cpu_count=0, gamma=gamma)
 
    mxy = mx + 1j*my
    ioff = int(len(f)/2)-1
    subplot(3,1,1)
    xlabel('Time [ms]')
    plot(t*1000,b1)
    subplot(3,1,2)
    plot(x, abs(mxy[ioff,:]), x, real(mxy[ioff,:]), x, imag(mxy[ioff,:]) )
    xlabel('Position [cm]')
    ylabel('Magnetization |Mxy|')
    subplot(3,1,3)
    plot(x, mz[ioff,:])
    xlabel('Position [cm]')
    ylabel('Magnetization Mz')
    show()

    bob = 10
    bob = 2*bob
    

if __name__ == "__main__":
    main()
      
    #import cProfile
    #cProfile.run('main()')






