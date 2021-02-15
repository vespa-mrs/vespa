# Python modules



# 3rd party modules
import numpy as np


# Our modules
import vespa.common.constants as common_constants
import vespa.common.util.misc as util_misc
from vespa.common.constants import Deflate


class DataBasis(object):
    """ 
    This is a fundamental object that represents the data being
    manipulated and displayed in the visualization tab of Simulation. 
    
    It consists of the MRS spectral paramters and a numpy array that is used
    to contain various types of spectral data.
                
    """
    # This is the version of this object's XML output format. 
    XML_VERSION = "1.0.0"
    
    # Some default attributes follow. Creating them as static properties of
    # the class ensures that every time someone accesses one of these 
    # constants, they get a fresh (unique) copy of the list.
    
    DEFAULT_DIMS = \
        util_misc.StaticProperty(lambda: [common_constants.DEFAULT_SPECTRAL_POINTS, 1])
    
    def __init__(self, attributes=None):
        # FIXME PS - docstring is wrong
        """     
        Attributes
        -------------------------------------
        
        raw         object, typically an mrs_data_raw containing the raw 
                    k-space data for this Analysis data set.
        
        data        numpy.ndarray, complex, the result of applying all of the 
                    various transforms (svd, tdfd, baseline) to the time data
                    
        dims        list, int, dimension of the data array in case it differs 
                    from the time data.
                    
        sw          float, sweep width of data array in case it differs from 
                    raw data sweep width.
        """

        self.sw         = common_constants.DEFAULT_SWEEP_WIDTH
        self.resppm     = common_constants.DEFAULT_PROTON_CENTER_PPM
        self.frequency  = common_constants.DEFAULT_PROTON_CENTER_FREQUENCY

        self.data       = np.ndarray(self.DEFAULT_DIMS[::-1], dtype='complex64')

        if attributes is not None:
            self.inflate(attributes)



    @property
    def hpp(self):
        """Current hertz per point. It's read only."""
        return (self.sw / self.dims[0]) if self.dims[0] else 0.0

    @property
    def dims(self):
        """Current shape of the data array. It's read only."""
        return list(self.data.shape)[::-1]

    @property
    def raw_hpp(self):
        """Current hertz per point. It's read only."""
        return (self.sw / self.dims[0]) if self.dims[0] else 0.0

    @property
    def raw_dims(self):
        """Current shape of the data array. It's read only."""
        return list(self.data.shape)[::-1]

    @property
    def spectral_hpp(self):
        """Current hertz per point. It's read only."""
        return (self.sw / self.dims[0]) if self.dims[0] else 0.0

    @property
    def spectral_dims(self):
        """Current shape of the data array. It's read only."""
        return list(self.data.shape)[::-1]

    def ppm2pts(self, val, acq=False, rel=False):
        """
        Returns the point index along spectrum for given ppm value.
        - Assumes center point <--> resppm for rel False
        - Assumes center point <--> 0.0 ppm for rel True

        """
        dim0 = self.raw_dims[0] if acq else self.spectral_dims[0]
        hpp = self.raw_hpp if acq else self.spectral_hpp
        pts = self.frequency*val/hpp if rel else (dim0/2) - (self.frequency*(val-self.resppm)/hpp)
        pts = np.where(pts > 0, pts, 0)
        return pts

    def ppm2hz(self, val, acq=False, rel=False):
        """
        Returns the absolute number of hz away from 0.0 ppm based on an assumed ppm
        value for the center data point.

        If rel=True, assumes that center point is 0.0 ppm and calculates the
        relative hertz away represented by the ppm value.

        """
        hpp = self.raw_hpp if acq else self.spectral_hpp
        ppm = self.pts2hz(self.ppm2pts(val)) if rel else self.ppm2pts(val, rel=rel) * hpp
        return ppm

    def pts2ppm(self, val, acq=False, rel=False):
        """
        Returns the ppm value of the given point index along spectrum.
        - Assumes center point <--> resppm for rel False
        - Assumes center point <--> 0.0 ppm for rel True

        """
        dim0 = self.raw_dims[0] if acq else self.spectral_dims[0]
        hpp = self.raw_hpp if acq else self.spectral_hpp
        ppm = val*hpp/self.frequency if rel else (((dim0/2)-val)*(hpp/self.frequency))+self.resppm
        return ppm

    def pts2hz(self, val, acq=False, rel=False):
        """
        Returns the number of hertz away from 0.0 ppm from the points based on an
        assumed ppm value for the center point.

        If rel=True, assumes that center point is 0.0 ppm and calculates the
        relative hz away represented by the points value.

        """
        hpp = self.raw_hpp if acq else self.spectral_hpp
        hz = val * hpp if rel else (self.ppm2pts(0.0) - val) * hpp
        return hz

    def hz2ppm(self, val, acq=False, rel=False):
        """
        Returns the number of ppm from hertz based on an assumed ppm value for the
        center point.

        If rel=True, it is assumed that the hertz value is relative to 0.0 ppm
        equals 0.0 hertz. Thus we convert the hz value to points, take the distance
        in points from the 0.0 ppm point and convert that to ppm

        """
        hpp = self.raw_hpp if acq else self.spectral_hpp
        val = self.pts2ppm(self.hz2pts(val)) if rel else self.pts2ppm(val / hpp)
        return val

    def hz2pts(self, val, acq=False, rel=False):
        """
        Returns the number of points away from 0.0 hertz (0.0 ppm) based on an
        assumed ppm value for the center point.

        If rel=True, it is assumed that the hertz value is relative to 0.0 ppm
        equals 0.0 hertz. Thus we convert the hz value to points, take the distance
        in points from the 0.0 ppm point and convert that to points.

        """
        hpp = self.raw_hpp if acq else self.spectral_hpp
        pts = val / hpp if rel else self.ppm2pts(0.0) - (val / hpp)
        return pts


    def __str__(self):
        return self.__unicode__()


    def __unicode__(self):
        lines = [ ]
        lines.append("Data Dimensions: %s" % self.dims)
        lines.append("Spectral Sweep Width (Hz): %f" % self.sw)
        lines.append("Spectrometer Frequency (MHz): %f" % self.frequency)
        lines.append("Data type: %s" % str(self.data.dtype))

        # __unicode__() must return a Unicode object. In practice the code
        # above always generates Unicode, but we ensure it here.
        return '\n'.join(lines)
        
    
    def deflate(self, flavor=Deflate.ETREE):
        if flavor == Deflate.ETREE:
            raise NotImplementedError
            
        elif flavor == Deflate.DICTIONARY:
            return self.__dict__.copy()
    
        
    def inflate(self, source):
        if hasattr(source, "makeelement"):
            # Quacks like an ElementTree.Element
            raise NotImplementedError

        elif hasattr(source, "keys"):
            # Quacks like a dict
            for key in list(source.keys()):
                if hasattr(self, key):
                    setattr(self, key, source[key])
