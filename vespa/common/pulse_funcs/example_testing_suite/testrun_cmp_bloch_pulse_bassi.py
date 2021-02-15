# Python imports


from pylab import *

import numpy as np

import vespa.common.bloch_call as bloch_call

PI = np.pi
    

def bassi_pulse_design(npts, dwell, alpha, beta, b0, f0, vref=11.7, b1ref=11.7, kappa=2.0, x0=0.0, deltax=0.02):
    r"""

    Derived from Warnking paper MRM 52:1190-1199 (2004)

    Note. this code was rewritten to be all in one function (ie. no separate
    aBt or bBt functs so as to make port to C/C++ easier

    vref [uT] = reference voltage for scaling, default for Siemens 1ms hard 180 pulse = 11.7uT

    Pulse Shape Parameters and Variables and Their Units and Definitions
    ----------------------------------------------------------------------------------------------
    Variable        Units       Description                             Definition
    ----------------------------------------------------------------------------------------------
    a, a0, aB       rad         Amplitude parameter                     Eqs. [1] and [2]
    alpha           deg         Effective flip angle                    alpha = acos(1-2Pe)
    AM              rad/s       RF amplitude modulation function        Eq. [1]
    b, b0, bB       rad         Bandwidth parameter                     Eqs. [1] and [2], b = pi*mu
    beta            rad         Truncation factor of driving function   g(+/-Tp/2)=+/-2*beta
    c               1           Normalized off-resonance                Eqs. [3] and [5]
    f0              1           Maximal relative bandwidth increase     Eq. [15]
    FM              rad/s       RF frequency modulation function        Eq. [2]
    g               rad         HS driving function                     Eqs. [1] and [2]
    gdot            rad/s       Derivative of the driving function
    gamma           rad/s*T     Gyromagnetic ratio                      2.675222exp8 rad/s*T
    G               T/m         Slice-select gradient                   Eq. [18]
    kappa           1           BASSI amplitude scaling parameter       Eqs. [14] and [15]
    L               1           relative RF energy, as function of Pc   Eq. 9
    mu              1           Bandwidth parameter, notation from (14) mu = beta/pi
    Pe              1           Population inversion                    Pe=(1-cos(alpha))/2
    Tp              s           Pulse duration
    x0              m           Center of inversion / saturation band
    deltax          m           Width of inversion / saturation slab

    The population Pe inversion corresponds to an effective flip angle (i.e., the 
    flip angle of a plane rotation that would produce the same longitudinal 
    magnetization), of alpha = arccos(1-2Pe).


    Pulse parameters, as well as peak B1, peak gradient amplitude (Gmax) and width of the frequency 
    sweep (BW), for the pulses compared in the simulations and phantom experiments
    -----------------------------------------------------------------------------------------------
    Pulse               Tp      alpha       beta    b0      f0      Peak B1     Gmax       BW
                       (ms)     (deg)       (rad)                    (uT)       mT/m)     (kHz)
    -----------------------------------------------------------------------------------------------
    HS                  8.1     90          5.3     157     1       23.0        4.9         20.8
    FOCI                8.1     90          5.3     157     3.2     23.0        15.7        66.7
    BASSI(kappa=2)      8.1     90          5.3     168     3       23.1        15.7        66.7
    VERSE-HS            8.1     90          5.3     503             22.8        19.0        81.0

    HS                  10.24   170         5.3     37.7    1       22.9        0.9         4.0
    FOCI                10.24   175         6.3     22.8    22      23.0        14.7        62.7
    BASSI(kappa=2)      10.24   175         6.3     22.8    22      23.0        14.7        62.7
    BASSI(kappa=1.6)    10.24   175         6.3     22.8    22      23.0        14.7        62.7
    BASSI(GOIA)         10.24   175         6.3     22.8    22      23.0        14.7        62.7
    VERSE-HS            10.24   175         5.3     251             23.0        23.6        100.5


    """
    gamma = 267522200.       # 2.675222 x 10^8 rad/s*T

    Tp = npts * dwell
    t  = np.arange(npts) * dwell
    t  = t - (Tp/2.0)

    Pe = 0.5*(1-np.cos(alpha*PI/180.0))      # from comment above
    
    # Equation 9
    L = (1.0/(2*PI))*np.log(1.0/(np.sqrt(1-Pe)))

    # Equation 15,  bB(t) = bBt
    bBt = np.power(np.power((np.power(np.cosh(2*beta*t/Tp),kappa)*(b0-(PI*L)) + (PI*L)),-2) + np.power(f0*b0,-2), -0.5)
    bB0 = np.power(np.power((np.power(np.cosh(2*beta*0/Tp),kappa)*(b0-(PI*L)) + (PI*L)),-2) + np.power(f0*b0,-2), -0.5)
    
    # Equation 14
    aBt = np.power((bBt-PI*L)/(b0-PI*L),1.0/kappa) * np.sqrt(b0*b0 - 4*np.power(np.arccosh(np.cosh(b0/2)*np.sqrt(1-Pe)),2))
    aB0 = np.power((bB0-PI*L)/(b0-PI*L),1.0/kappa) * np.sqrt(b0*b0 - 4*np.power(np.arccosh(np.cosh(b0/2)*np.sqrt(1-Pe)),2))
    
    # Equation 16
    #
    # Acc'd to paper, if we plug in vref and b1ref output is in Volts
    # For now, we set vref=b1ref thus removing that term from the equation
    # and it returns in T (I think)
    
    amt = (vref/b1ref)*(2*aB0*beta/(PI*gamma*Tp)) * (aBt.copy()/aB0)*(1.0/np.cosh(2*beta*t/Tp)) # for sech()
        
    # Equation 18
    #
    # Acc'd to paper this returns in Tesla/meter
    
    g = bBt.copy()*4*beta /(PI*gamma*deltax*Tp)
    
    # Equation 17
    # 
    # first we create the function, then integrate under it for each point
    fmt  = ((2*bBt.copy()*beta)/(PI*Tp))*((-2*x0/deltax) + np.tanh(2*beta*t/Tp))
    phit = np.cumsum(fmt*dwell)
    phit = phit % (PI*2)
    
    return amt, fmt, g, phit



def bassi_pulse_design_c(npts, dwell, alpha, beta, b0, f0, vref=11.7, b1ref=11.7, kappa=2.0, x0=0.0, deltax=0.02):
    r"""

    Note. this code was rewritten to be all in one function (ie. no separate
    aBt or bBt functs so as to make port to C/C++ easier

    """
    gamma = 267522200.       # 2.675222 x 10^8 rad/s*T

#     if Tp > 0.008192:
#         dwell = 1.0/500000.0
#     else:
#         dwell = 1.0/1000000.0
#     
#     npts = int(np.round(Tp/dwell))

    Tp = npts * dwell

    ts   = np.ndarray((npts,), 'float')
    aBt  = np.ndarray((npts,), 'float')
    bBt  = np.ndarray((npts,), 'float')
    amt  = np.ndarray((npts,), 'float')
    fmt  = np.ndarray((npts,), 'float')
    g    = np.ndarray((npts,), 'float')
    phit = np.ndarray((npts,), 'float')
    aB0  = 0.0
    bB0  = 0.0
    
    for i in range(npts):
        ts[i] = i*dwell - (Tp/2.0)

    Pe = 0.5*(1-np.cos(alpha*PI/180.0))      # from comment above
    
    # Equation 9
    L = (1.0/(2*PI))*np.log(1.0/(np.sqrt(1-Pe)))

    # from Eqns 15 and 16 respectively for t=0
    bB0 = np.power(np.power((np.power(np.cosh(2*beta*0/Tp),kappa)*(b0-(PI*L)) + (PI*L)),-2) + np.power(f0*b0,-2), -0.5)
    aB0 = np.power((bB0-PI*L)/(b0-PI*L),1.0/kappa) * np.sqrt(b0*b0 - 4*np.power(np.arccosh(np.cosh(b0/2)*np.sqrt(1-Pe)),2))

    print(bB0, aB0)

    # Eqn 15 (bB(t) = bBt), Eqn 14, Eqn 16, Eqn 18, Eqn 17
    #
    # Equation 16 - returns in T (I think)
    #  - Acc'd to paper, if we plug in vref and b1ref output is in Volts. 
    #    For now we set vref=b1ref thus removing that term from the equation 
    # Equation 18 - Acc'd to paper this returns in Tesla/meter
    # Equation 17 - First we create the function, then integrate under it for each point
    last_sum = 0.0
    for i in range(npts):
        bBt[i]  = np.power(np.power((np.power(np.cosh(2*beta*ts[i]/Tp),kappa)*(b0-(PI*L)) + (PI*L)),-2) + np.power(f0*b0,-2), -0.5)
        aBt[i]  = np.power((bBt[i]-PI*L)/(b0-PI*L),1.0/kappa) * np.sqrt(b0*b0 - 4*np.power(b_acosh(np.cosh(b0/2)*np.sqrt(1-Pe)),2))
        amt[i]  = (vref/b1ref)*(2*aB0*beta/(PI*gamma*Tp)) * (aBt[i]/aB0)*(1.0/np.cosh(2*beta*ts[i]/Tp)) # for sech()
        g[i]    = bBt[i]*4*beta /(PI*gamma*deltax*Tp)
        fmt[i]  = ((2*bBt[i]*beta)/(PI*Tp))*((-2*x0/deltax) + np.tanh(2*beta*ts[i]/Tp))
        phit[i] = (last_sum + fmt[i]*dwell) 
        last_sum = phit[i]
        phit[i] = WrapTwoPI(phit[i])
        

    return amt, fmt, g, phit

def b_acosh(x):
    return np.log(x + np.sqrt(np.power(x,2) - 1))


_PI     = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348;
_TWO_PI = 6.2831853071795864769252867665590057683943387987502116419498891846156328125724179972560696;

# Floating-point modulo
# The result (the remainder) has same sign as the divisor.
# Similar to matlab's mod(); Not similar to fmod() -   Mod(-3,4)= 1   fmod(-3,4)= -3
def Mod(x,y):

    if 0.0 == y:
        return x

    m = x - y * np.floor(x/y)

    # handle boundary cases resulted from floating-point cut off:
    if y > 0:                   # modulo range: [0..y)
        if m >= y:              # Mod(-1e-16             , 360.    ): m= 360.
            return 0
        if m <0:
            if (y+m) == y:
                return 0        # just in case...
            else:
                return y+m      # Mod(106.81415022205296 , _TWO_PI ): m= -1.421e-14 
    else:                       # modulo range: (y..0]
        if m <= y:              # Mod(1e-16              , -360.   ): m= -360.
            return 0
        if m > 0:
            if (y+m) == y:
                return 0    # just in case...
            else:
                return y+m  # Mod(-106.81415022205296, -_TWO_PI): m= 1.421e-14 
    return m


# wrap [rad] angle to [-PI..PI)
def WrapPosNegPI(fAng):
    return Mod(fAng + _PI, _TWO_PI) - _PI

# wrap [rad] angle to [0..TWO_PI)
def WrapTwoPI(fAng):
    return Mod(fAng, _TWO_PI)

# wrap [deg] angle to [-180..180)
def WrapPosNeg180(fAng):
    return Mod(fAng + 180.0, 360.0) - 180.0

# wrap [deg] angle to [0..360)
def Wrap360(fAng):
    return Mod(fAng ,360.0)


#------------------------------------------------------------------------------
# Testing starts here

def main():
    """
    Pulse parameters, as well as peak B1, peak gradient amplitude (Gmax) and width of the frequency 
    sweep (BW), for the pulses compared in the simulations and phantom experiments
    -----------------------------------------------------------------------------------------------
    Pulse               Tp      alpha       beta    b0      f0      Peak B1     Gmax       BW
                       (ms)     (deg)       (rad)                    (uT)       mT/m)     (kHz)
    -----------------------------------------------------------------------------------------------
    HS                  8.1     90          5.3     157     1       23.0        4.9         20.8
    FOCI                8.1     90          5.3     157     3.2     23.0        15.7        66.7
    BASSI(kappa=2)      8.1     90          5.3     168     3       23.1        15.7        66.7
    VERSE-HS            8.1     90          5.3     503             22.8        19.0        81.0

    HS                  10.24   170         5.3     37.7    1       22.9        0.9         4.0
    FOCI                10.24   175         6.3     22.8    22      23.0        14.7        62.7
    BASSI(kappa=2)      10.24   175         6.3     22.8    22      23.0        14.7        62.7
    BASSI(kappa=1.6)    10.24   175         6.3     22.8    22      23.0        14.7        62.7
    BASSI(GOIA)         10.24   175         6.3     22.8    22      23.0        14.7        62.7
    VERSE-HS            10.24   175         5.3     251             23.0        23.6        100.5

    WET Multi-pulse settings, assumes 8192 pts, 1us dwell, B1max=23ish, Gmax=13
    Flip [deg]      b0
      90           140
      83.6         140
      99.7         140
      74.7         140
     160.0          36
              
      

    """
    import sys
    sys.path.insert(0, 'C:\\Users\\bsoher\\code\\repository_svn\\orphans\\soher\\rf_pulse_design\\hargreaves_bloch')
#    import bloch    # NB this is Hargreaves Bloch
    from . import bloch_multi

    gamma = 26751.3  # Hz/gauss

#     npts   = 8100
#     dwell  = 0.000001
#     Tp     = 0.0081    # npts * dwell
#     alpha  = 90.0
#     beta   = 5.3
#     b0     = 168.0
#     f0     = 3.0
#     deltax = 0.1
#     x0     = 0.0
#     kappa  = 2.0

#  This gives a 2ms 90 pulse with 17kHz bw, 10 mT/m gradient
#
#     npts   = 512       # 8100 with dwell 0.000001
#     dwell  = 0.000004
#     Tp     = npts * dwell
#     alpha  = 90.0
#     beta   = 5.3
#     b0     = 11.0  #168.0
#     f0     = 3.0
#     deltax = 0.04
#     x0     = 0.0
#     kappa  = 2.0

#  This gives a 4.5ms 178 pulse with 6.7kHz bw, 15uT rf and 3.8 mT/m gradient
#
    npts   = 1000       
    dwell  = 0.0000045
    Tp     = npts * dwell
    alpha  = 178.0
    beta   = 3.5
    b0     = 9.7
    f0     = 5.0
    deltax = 0.04
    x0     = 0.0
    kappa  = 2.0


#  This gives a 2ms 90 pulse with 24kHz bw, 14 mT/m gradient
#
#     npts   = 1024       # 8100 with dwell 0.000001
#     dwell  = 0.000002
#     Tp     = npts * dwell
#     alpha  = 90.0
#     beta   = 5.3
#     b0     = 15.0  #168.0
#     f0     = 3.0
#     deltax = 0.04
#     x0     = 0.0
#     kappa  = 2.0


    amt,  fmt,  g,  phit  = bassi_pulse_design(npts, dwell, alpha, beta, b0, f0, kappa=2.0, deltax=deltax, x0=x0)
#    amti, fmti, gi, phiti = bassi_pulse_design_c(npts, dwell, alpha, beta, b0, f0, kappa=2.0, deltax=deltax, x0=x0)
 
#    print "Inline - funct = "+str(sum(amti-amt))+" "+str(sum(fmti-fmt))+" "+str(sum(gi-g))+" "+str(sum(phiti-phit))
 #   print "g[0] = "+str(g[0]*1000)+"  max amt = "+str(max(amt))

#     subplot(411)
#     plot(amt*1000000)                   # plot in uT
#     subplot(412)
#     plot(fmt / (2.0 * PI * 1000))       # plot in kHz
#     subplot(413)
#     plot(g * 1000.0)                    # plot in mT/m
#     subplot(414)
#     plot(real(amt*1000000 * np.exp(1j*phit)))                  # plot in rad
#     subplot(414)
#     plot(imag(amt*1000000 * np.exp(1j*phit)))                  # plot in rad

    show()

    b1c = amt*10000.0 * np.exp(1j*phit)        # convert T to G
    gr  = g * 1000 * 0.1                       # convert to mT/m (G/cm)
    xlft = -0.5*deltax*100 + x0*100 - 2
    xrgt =  0.5*deltax*100 + x0*100 + 2
    x = np.arange(xlft,xrgt,float(xrgt-xlft)/800.0)
    T = dwell  #0.000001
    f = np.array([0.0])
    t = np.arange(1,len(b1c)+1)*T;

    bw = gr[0] * 4358.0 * (deltax*100)    # G/cm  Hz/G  cm 

    mx,my,mz = bloch_multi.bloch_multi(b1c,gr,t,1,.2,f,x,mode=0, gamma=gamma);
    
    mxy = mx + 1j*my

    print("Bandwidth = "+str(bw))
    subplot(311)
    plot(t,b1c.real*100,t,b1c.imag*100,t,np.abs(b1c*100)) # plot in uT
    subplot(312)
    plot(g * 1000.0)                    # plot in mT/m
    subplot(313)
    plot(x*10,mz)                       # scale cm to mm
#     subplot(414)
#     plot(x*10,mxy, x*10,np.abs(mxy))    # scale cm to mm
    ylim([-1.1,1.1]) 
    show()

    #--------------------------------------------------
    # Matpulse engine
    #--------------------------------------------------

#     calc_resolution = 5000
#     # dwell time in microseconds 
#     dwell_time = dwell * 1000000
# 
#     rfwave = b1c / 10   # convert G -> mT
#     rfgrad = gr * 10    # convert G -> mT
# 
#     bout = bloch_call.calc_all_profiles_b2g2(rfwave, rfgrad, dwell_time, calc_resolution)
# 
# 
#     mz_std    = bout[0]
#     mxy_std   = bout[1]
#     xaxis_std = bout[2]
#     mz_ext    = bout[3]
#     mxy_ext   = bout[4]
#     xaxis_ext = bout[5]     
#      
#     bob = 10
#     bob += 1
#      
#     mxy = np.zeros(len(mz_std), np.complex128) 
#     mxy = np.zeros(len(mz_std), np.complex128) 
#     for i in range(len(mz_std)):
#         mxy[i] = complex(mz_std[i][0], mz_std[i][1])
#     
#     subplot(4,1,1)
#     plot(t,rfwave.real*1000,t,rfwave.imag*1000,t,np.abs(rfwave*1000)) # plot in uT
#     xlabel('Time [ms]')
#     subplot(4,1,2)
#     plot(t,rfgrad)                                # plot in mT/m
#     xlabel('Time [ms]')
#     subplot(4,1,3)
#     plot(xaxis_std, abs(mxy), xaxis_std, real(mxy))     # x-axis in mm
#     xlabel('Position [mm]')
#     ylabel('Magnetization |Mxy|')
#     subplot(4,1,4)
#     plot(xaxis_std, mz_std[:,2])                        # x-axis in mm
#     xlabel('Position [mm]')
#     ylabel('Magnetization Mz')
#     show()

    bob = 10
    bob = 2*bob
    

if __name__ == "__main__":
    main()
      
    #import cProfile
    #cProfile.run('main()')
