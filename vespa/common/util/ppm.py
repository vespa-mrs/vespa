""" 
=========== 
Description
===========   

Functions for converting from/to PPM, points and hertz, including: 

 ppm2pts()
 ppm2hz()
 pts2ppm()
 pts2hz()
 hz2ppm()
 hz2pts()

===================
Arguments (general)
===================
**val:**  [array/scalar][float]

  The values (in ppm/points/hertz) that are to be converted. They can be
  a numpy array or a scalar value of floats.  

**data:**  [object]

  Object that contains data information used to convert the values. These 
  may be provided as actual attributes of the object or as properties that
  return the values. The only restriction at the moment is that they are 
  available as calls at the top level of the object.
  
  The information required may include some or all of the following: 
  
  data.frequency     - [float] central frequency, in MHz (eg. 123.8 at 3 Tesla)
  data.sw            - [float] sweep width, in Hz (eg. 1536.0)
  data.raw_dims      - [float] data shape of raw data array. This typically has
                          the spectral points then any spatial dimensions in an
                          array or tuple (eg. [1024,1,1,1]) 
  data.spectral_dims - [float] data shape of the data array after spectral 
                          processing. Data may have zero filling applied. This 
                          shape also has the spectral points then any spatial 
                          dimensions in an array or tuple (eg. [2048,1,1,1])
  data.resppm        - [float] resonance ppm value, the ppm value of the 
                          central frequency value.
  data.raw_hpp       - [float] hertz per point of raw data. 
  data.spectral_hpp  - [float] hertz per point of spectral data. This may differ 
                          from raw_hpp if the data has been zero filled to a 
                          larger spectral dimension.

===================
Keywords (general)
===================
**acq:**  [float][default=False]

  Flag to do calculations based on data acquisition values rather than data
  values that may be changed by spectral processing. This is mainly whether
  the raw_dims or spectral_dims are used.   

**rel:**  [float][default=False]

  Flag to do calculatios based on relative change versus absolute change from
  a reference. Typical references are the central frequency ppm value (which
  is 4.7 ppm for 1H data). 
  
======    
Syntax
====== 
::

  val = ppm2pts( val, data )
  val = ppm2hz( val, data, acq=True)

  and similarly ...


"""

# Python modules


# 3rd party modules
import numpy as np


def ppm2pts(val, data, acq=False, rel=False):
    """ 
    Returns the number of points away from 0.0 ppm based on an assumed ppm 
    value for the center data point. 
    
    If rel=True, assumes that center point is 0.0 ppm and calculates the 
    relative points away represented by the ppm value.
    
    """
    if acq:
        dim0   = data.raw_dims[0]
        ctrppm = data.resppm
        hpp    = data.raw_hpp
    else:
        dim0   = data.spectral_dims[0]
        ctrppm = data.resppm
        hpp    = data.spectral_hpp
    
    if rel:
        pts = data.frequency * val / hpp
    else:
        pts = (dim0 / 2) - (data.frequency * (val - ctrppm) / hpp)
        pts = np.where(pts > 0, pts, 0)
    
    return pts


def ppm2hz(val, data, acq=False, rel=False):
    """ 
    Returns the absolute number of hz away from 0.0 ppm based on an assumed ppm 
    value for the center data point.

    If rel=True, assumes that center point is 0.0 ppm and calculates the 
    relative hertz away represented by the ppm value.
    
    """
    if rel:
        ppm = pts2hz(ppm2pts(val,data),data)
    else:
        hpp = data.raw_hpp if acq else data.spectral_hpp
        ppm = ppm2pts(val, data, rel=rel) * hpp
    
    return ppm


def pts2ppm(val, data, acq=False, rel=False):
    """ 
    Returns the number of ppm from the points based on an assumed ppm value for
    the center data point.

    If rel=True, assumes that center point is 0.0 ppm and calculates the 
    relative ppm away represented by the points value.
    
    """
    if acq:
        dim0   = data.raw_dims[0]
        ctrppm = data.resppm
        hpp    = data.raw_hpp
    else:
        dim0   = data.spectral_dims[0]
        ctrppm = data.resppm
        hpp    = data.spectral_hpp
    
    if rel:
        ppm = val * hpp / data.frequency
    else:
        ppm = ( ((dim0 / 2) - val) * (hpp / data.frequency) ) + ctrppm
    
    return ppm


def pts2hz(val, data, acq=False, rel=False):
    """ 
    Returns the number of hertz away from 0.0 ppm from the points based on an 
    assumed ppm value for the center point.

    If rel=True, assumes that center point is 0.0 ppm and calculates the 
    relative hz away represented by the points value.
    
    """
    hpp = data.raw_hpp if acq else data.spectral_hpp
    
    if rel:
        hz = val * hpp
    else:
        hz = (ppm2pts(0.0, data) - val) * hpp
    
    return hz


def hz2ppm(val, data, acq=False, rel=False):
    """ 
    Returns the number of ppm from hertz based on an assumed ppm value for the 
    center point.
    
    If rel=True, it is assumed that the hertz value is relative to 0.0 ppm 
    equals 0.0 hertz. Thus we convert the hz value to points, take the distance 
    in points from the 0.0 ppm point and convert that to ppm
    
    """
    if rel:
        val = pts2ppm(hz2pts(val,data),data)
    else:
        hpp = data.raw_hpp if acq else data.spectral_hpp
        
        val = pts2ppm(val/hpp, data)
    
    return val


def hz2pts(val, data, acq=False, rel=False):
    """ 
    Returns the number of points away from 0.0 hertz (0.0 ppm) based on an 
    assumed ppm value for the center point.
    
    If rel=True, it is assumed that the hertz value is relative to 0.0 ppm 
    equals 0.0 hertz. Thus we convert the hz value to points, take the distance 
    in points from the 0.0 ppm point and convert that to points.
    
    """
    hpp = data.raw_hpp if acq else data.spectral_hpp
    
    if rel:
        pts = val / hpp
    else:
        pts = ppm2pts(0.0, data) - (val / hpp)
    
    return pts

