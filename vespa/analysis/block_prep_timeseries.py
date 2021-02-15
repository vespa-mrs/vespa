# Python modules


# 3rd party modules
import numpy as np
from xml.etree.cElementTree import Element

# Our modules
import vespa.analysis.block_prep_identity as block_prep_identity
import vespa.analysis.chain_prep_timeseries as chain_prep_timeseries
import vespa.analysis.block as block

import vespa.common.util.xml_ as util_xml
import vespa.common.mrs_data_raw as mrs_data_raw

from vespa.common.constants import Deflate




class _Settings(object):
    """
    Settings object contains the parameter inputs used for processing in the 
    Chain object in this Block. Having a separate object helps to delineate 
    inputs/outputs and to simplify load/save of preset values.

    This object can also save/recall these values to/from an XML node.

    """
    XML_VERSION = "1.0.0"
    
    def __init__(self, attributes=None):
        """
        Based on 'block_prep_fidsum' class which deals with "raw" types of 
        data objects that need to do a bit of massaging of the data as it 
        comes in (e.g. align and sum individual FIDs for an SVS data set).
        
        Timeseries specific attributes
        -----------------------------------
        fids_to_average        int, number of raw fids to sum together to
                                 form the 2nd dimension in the data object.
                                 Note. this changes the data dimensions and
                                 will cause all downstream results to reset.
                                 
        Attributes in common with Fidsum class
        -------------------------------------------
        coil_combine_method    string, name of algorithm to combine data from coils
        fid_left_shift         int, applied before auto phase/shift algorithms are 
                                 calculated, not applied to final data calculation
        gaussian_apodization   float, applied before auto phase/shift algorithms are 
                                 calculated, not applied to final data calculation
        apply_peak_shift       bool, flag on whether to perform/apply correction
        reference_peak_center  float, input for b0 shift algorithm
        peak_search_width      float, used in b0 shift algorithm
        apply_phase0           bool, flag on whether to perform/apply correction
        phase0_range_start     float, input for phase correction
        phase0_range_end       float, input for phase correction
        global_phase0          float, degrees, phase 0 amount added to all FIDs.
                                 This value is applied prior to the phase 0 
                                 correction algorithm and also as part of the 
                                 final data calculation. Note. each FID can also 
                                 have a local delta ph0 added as well. 
        global_phase1          float, degrees, phase 1 amount added to all FIDs
                                 This value is applied to final data calculation


        
        """
        self.coil_combine_method    = 'Siemens'
        self.fids_to_average        = 1
        self.fid_left_shift         = 0
        self.gaussian_apodization   = 2.0
        self.apply_peak_shift       = True
        self.reference_peak_center  = 2.01
        self.peak_search_width      = 0.2
        self.apply_phase0           = True
        self.phase0_range_start     = 3.5
        self.phase0_range_end       = 0.5
        self.global_phase0          = 0.0
        self.global_phase1          = 0.0
        
        if attributes is not None:
            self.inflate(attributes)


    def __str__(self):
        lines = []
        lines.append("------- {0} Object -------".format(self.__class__.__name__))
        lines.append("coil_combine_method       : " + str(self.coil_combine_method))
        lines.append("fids_to_average           : " + str(self.fids_to_average))
        lines.append("fid_left_shift            : " + str(self.fid_left_shift))
        lines.append("gaussian_apodization      : " + str(self.gaussian_apodization))
        lines.append("apply_peak_shift          : " + str(self.apply_peak_shift))
        lines.append("reference_peak_center     : " + str(self.reference_peak_center))
        lines.append("peak_search_width         : " + str(self.peak_search_width))
        lines.append("apply_phase0              : " + str(self.apply_phase0))
        lines.append("phase0_range_start        : " + str(self.phase0_range_start))
        lines.append("phase0_range_end          : " + str(self.phase0_range_end))
        lines.append("global_phase0             : " + str(self.global_phase0))
        lines.append("global_phase1             : " + str(self.global_phase1))

        # __unicode__() must return a Unicode object. In practice the code
        # above always generates Unicode, but we ensure it here.
        return '\n'.join(lines)


    def deflate(self, flavor=Deflate.ETREE):
        if flavor == Deflate.ETREE:
            e = Element("settings", {"version" : self.XML_VERSION})

            util_xml.TextSubElement(e, "coil_combine_method",   self.coil_combine_method)
            util_xml.TextSubElement(e, "fids_to_average",       self.fids_to_average)
            util_xml.TextSubElement(e, "fid_left_shift",        self.fid_left_shift)
            util_xml.TextSubElement(e, "gaussian_apodization",  self.gaussian_apodization)
            util_xml.TextSubElement(e, "apply_peak_shift",      self.apply_peak_shift)
            util_xml.TextSubElement(e, "reference_peak_center", self.reference_peak_center)
            util_xml.TextSubElement(e, "peak_search_width",     self.peak_search_width)
            util_xml.TextSubElement(e, "apply_phase0",          self.apply_phase0)
            util_xml.TextSubElement(e, "phase0_range_start",    self.phase0_range_start)
            util_xml.TextSubElement(e, "phase0_range_end",      self.phase0_range_end)
            util_xml.TextSubElement(e, "global_phase0",         self.global_phase0)
            util_xml.TextSubElement(e, "global_phase1",         self.global_phase1)

            return e
            
        elif flavor == Deflate.DICTIONARY:
            return self.__dict__.copy()


    def inflate(self, source):
        if hasattr(source, "makeelement"):
            # Quacks like an ElementTree.Element

            for name in ("reference_peak_center", "gaussian_apodization", 
                         "peak_search_width", "global_phase0", "global_phase1", 
                         'phase0_range_start', 'phase0_range_end'):
                item = source.findtext(name)
                if item is not None:
                    setattr(self, name, float(item))

            for name in ("fids_to_average", "fid_left_shift",):
                item = source.findtext(name)
                if item is not None:
                    setattr(self, name, int(item))
                
            for name in ("apply_peak_shift", "apply_phase0", ):
                item = source.findtext(name)
                if item is not None:
                    setattr(self, name, util_xml.BOOLEANS[item])

            for name in ("coil_combine_method",):
                item = source.findtext(name)
                if item is not None:
                    setattr(self, name, item)

        elif hasattr(source, "keys"):
            # Quacks like a dict
            for key in list(source.keys()):
                if hasattr(self, key):
                    setattr(self, key, source[key])



class BlockPrepTimeseries(block_prep_identity.BlockPrepIdentity):
    """ 
    Building block to hold the state of a step in an MRS processing chain.
    Includes the functionality to save/recall this object to/from an XML node.

    Contains inputs/results for preprocessing of the raw data from the previous
    block ('raw') in the dataset.blocks list. This step modifies coil/average
    data into a single summed FID array for one dataset.
    
    """
    XML_VERSION = "1.0.0"

    
    def __init__(self, attributes=None):
        """
        Block objects have a self.set attribute that contains a _Settings object
        that contains the input attributes for the processing done in the Chain.
        This simplifies using a Block object as a preset. Results from the Chain
        are stored in this object at the level of self.set

        Base class sets references to:  self.id, self.data, self.chain and self.behave_as_preset
        
        """
        super().__init__(attributes)
        
        # processing parameters
        self.set = _Settings()

        # results storage
        self.measure_time = None          # store here in case we average FIDs, filled by chain!
        self.frequency_shift = None
        self.phase_0 = None
        self.data = None
                                 
        if attributes is not None:
            self.inflate(attributes)

        self.chain = None


    ##### Standard Methods and Properties #####################################

    def __str__(self):
        lines = []
        lines.append("------- {0} Object -------".format(self.__class__.__name__))
        lines.append("\n")
        lines += _Settings.__str__(self).split('\n')
        lines.append("\n")
        lines.append("------- Main Object -------")
        lines.append("Data shape                    : %s" % str(self.dims))

        return '\n'.join(lines)


    def create_chain(self, dataset):
        self.chain = chain_prep_timeseries.ChainPrepTimeseries(dataset, self)


    def set_dims(self, dataset):
        """
        Given a Dataset object, this is an opportunity for this block object 
        to ensure that its dims match those of the parent dataset. 
        """
        block.Block.set_dims(self, dataset)
        
        raw = dataset.blocks['raw']

        # local reference to input data
        data = dataset.get_source_data('prep')

        # this is the calculated proper size for self.data
        raw_dims  = list(data.shape)   # so we can compare to self.dims LIST

        # if we average FIDs, the dimensions change here        
        raw_dims[-2] = int(raw_dims[-2]/self.set.fids_to_average)

        if self.dims is None:
            self._reset_dimensional_data(dataset)
        elif (self.dims)[::-1] != raw_dims:                  #FIXME bjs - need reverse here, ARRRRGH, why?
            self._reset_dimensional_data(dataset)

        # calculate measure_time array based on whether we average or not
        measure_time = list(raw.measure_time)
        nfids = len(measure_time)
        navgs = self.set.fids_to_average
        measure_time = measure_time[0::navgs]
        if (nfids % navgs) != 0:
            del measure_time[-1]
        self.measure_time = np.array(measure_time)
        

    def _reset_dimensional_data(self, dataset):
        """
        Resets (to zero) and resizes dimensionally-dependent data
        
        """
        # local reference to input data
        raw = dataset.get_source_data('prep')

        nfids = raw.shape[-2]
        
        nfids = int(nfids/self.set.fids_to_average)
        
        data_shape = list(raw.shape)
        data_shape[-2] = nfids

        self.frequency_shift = np.zeros([nfids])
        self.phase_0         = np.zeros([nfids])
        self.measure_time    = np.arange(nfids)

        self.data = np.zeros(data_shape, dtype=raw.dtype)
        if self.chain is not None:
            self.chain.reset_results_arrays()
        
    
    def concatenate(self, new):
        raise NotImplementedError


    def deflate(self, flavor=Deflate.ETREE):
        if flavor == Deflate.ETREE:

            e = Element("block_prep_timeseries", {"id" : self.id,
                                                  "version" : self.XML_VERSION})

            util_xml.TextSubElement(e, "behave_as_preset", self.behave_as_preset)
            
            # Now I deflate the attribs that are specific to this class
            e.append(self.set.deflate())
            
            if not self.behave_as_preset:

                e.append(util_xml.numpy_array_to_element(self.measure_time,'measure_time'))
                e.append(util_xml.numpy_array_to_element(self.frequency_shift,'frequency_shift'))
                e.append(util_xml.numpy_array_to_element(self.phase_0,'phase_0'))
                e.append(util_xml.numpy_array_to_element(self.data, 'data'))

            return e

        elif flavor == Deflate.DICTIONARY:
            return self.__dict__.copy()


    def inflate(self, source):
        if hasattr(source, "makeelement"):

            val = source.findtext("behave_as_preset")   # default is False
            if val is not None:
                self.behave_as_preset = util_xml.BOOLEANS[val]

            # Quacks like an ElementTree.Element
            self.set = _Settings(source.find("settings"))

            if not self.behave_as_preset:
            
                # Now I inflate the attribs that are specific to this class
                temp = source.find("frequency_shift")
                self.frequency_shift = util_xml.element_to_numpy_array(temp)
                temp = source.find("phase_0")
                self.phase_0 = util_xml.element_to_numpy_array(temp)
                temp = source.find("data")
                self.data = util_xml.element_to_numpy_array(temp)

                temp = source.find("measure_time")
                if temp is not None:
                    self.measure_time = util_xml.element_to_numpy_array(temp)
                else:
                    self.measure_time = np.arange(len(self.frequency_shift))


        elif hasattr(source, "keys"):
            # Quacks like a dict
            for key in list(source.keys()):
                if key == "set":
                    setattr(self, key, source[key])


    ##### Private Methods #####################################

    
