# Python modules


# 3rd party modules
import numpy as np
from xml.etree.cElementTree import Element

# Our modules
import vespa.analysis.block_prep_identity as block_prep_identity
import vespa.analysis.chain_prep_wbnaa as chain_prep_wbnaa
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
        The BlockPrepWbnaa class is a specialized version of BlockPrepFidsum 
        deals with "raw" data from Oded Gonen's Whole Brain NAA sequence.
         
        This version allows us to do a bit of massaging of the data that is
        slightly more specific that the generic Fidsum version of this class. 
        (e.g. align and sum individual FIDs for an SVS data set).
        
        Attributes 
        -------------------------------------------
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
        global_phase1          float, degrees, phase 1 amount added to all FIDs
                                 This value is applied to final data calculation
        
        """
        self.apply_data_exclusion   = False
        self.fid_left_shift         = 0
        self.gaussian_apodization   = 2.0
        
        self.apply_peak_shift       = True
        self.reference_peak_center  = 2.01
        self.peak_search_width      = 0.2
        self.fid_left_shift_b0      = 56
        
        self.apply_phase0           = True
        self.phase0_range_start     = 2.2
        self.phase0_range_end       = 1.8
        self.fid_left_shift_phase0  = 56
        self.ref_spectrum_source    = 'singlet_centered_in_range'
        self.ref_peak_line_width    = 18
        self.constant_phase0_offset = 70     # degrees
        self.global_phase1          = 0.0
        
        if attributes is not None:
            self.inflate(attributes)


    def __str__(self):
        lines = []
        lines.append("------- {0} Object -------".format(self.__class__.__name__))
        lines.append("apply_data_exclusion      : " + str(self.apply_data_exclusion))
        lines.append("fid_left_shift            : " + str(self.fid_left_shift))
        lines.append("gaussian_apodization      : " + str(self.gaussian_apodization))
        lines.append("apply_peak_shift          : " + str(self.apply_peak_shift))
        lines.append("reference_peak_center     : " + str(self.reference_peak_center))
        lines.append("peak_search_width         : " + str(self.peak_search_width))
        lines.append("fid_left_shift_b0         : " + str(self.fid_left_shift_b0))
        lines.append("apply_phase0              : " + str(self.apply_phase0))
        lines.append("phase0_range_start        : " + str(self.phase0_range_start))
        lines.append("phase0_range_end          : " + str(self.phase0_range_end))
        lines.append("fid_left_shift_phase0     : " + str(self.fid_left_shift_phase0))
        lines.append("ref_spectrum_source       : " + str(self.ref_spectrum_source))
        lines.append("ref_peak_line_width       : " + str(self.ref_peak_line_width))
        lines.append("constant_phase0_offset    : " + str(self.constant_phase0_offset))
        lines.append("global_phase1             : " + str(self.global_phase1))
        return '\n'.join(lines)


    def deflate(self, flavor=Deflate.ETREE):
        if flavor == Deflate.ETREE:
            e = Element("settings", {"version" : self.XML_VERSION})

            util_xml.TextSubElement(e, "apply_data_exclusion",  self.apply_data_exclusion)
            util_xml.TextSubElement(e, "fid_left_shift",        self.fid_left_shift)
            util_xml.TextSubElement(e, "gaussian_apodization",  self.gaussian_apodization)
            util_xml.TextSubElement(e, "apply_peak_shift",      self.apply_peak_shift)
            util_xml.TextSubElement(e, "reference_peak_center", self.reference_peak_center)
            util_xml.TextSubElement(e, "peak_search_width",     self.peak_search_width)
            util_xml.TextSubElement(e, "fid_left_shift_b0",     self.fid_left_shift_b0)
            util_xml.TextSubElement(e, "apply_phase0",          self.apply_phase0)
            util_xml.TextSubElement(e, "phase0_range_start",    self.phase0_range_start)
            util_xml.TextSubElement(e, "phase0_range_end",      self.phase0_range_end)
            util_xml.TextSubElement(e, "fid_left_shift_phase0", self.fid_left_shift_phase0)
            util_xml.TextSubElement(e, "ref_spectrum_source",   self.ref_spectrum_source)
            util_xml.TextSubElement(e, "ref_peak_line_width",   self.ref_peak_line_width)
            util_xml.TextSubElement(e, "constant_phase0_offset",self.constant_phase0_offset)
            util_xml.TextSubElement(e, "global_phase1",         self.global_phase1)

            return e
            
        elif flavor == Deflate.DICTIONARY:
            return self.__dict__.copy()


    def inflate(self, source):
        if hasattr(source, "makeelement"):
            # Quacks like an ElementTree.Element

            for name in ("reference_peak_center", "gaussian_apodization", 
                         "peak_search_width", "global_phase1", 
                         'phase0_range_start', 'phase0_range_end'):
                item = source.findtext(name)
                if item is not None:
                    setattr(self, name, float(item))

            for name in ("fid_left_shift", "fid_left_shift_b0",
                         "fid_left_shift_phase0", "ref_peak_line_width",
                         "constant_phase0_offset"):
                item = source.findtext(name)
                if item is not None:
                    setattr(self, name, int(item))

            for name in ("apply_peak_shift", "apply_phase0", "apply_data_exclusion", ):
                item = source.findtext(name)
                if item is not None:
                    setattr(self, name, util_xml.BOOLEANS[item])

            for name in ("ref_spectrum_source",):
                item = source.findtext(name)
                if item is not None:
                    setattr(self, name, item)


        elif hasattr(source, "keys"):
            # Quacks like a dict
            for key in list(source.keys()):
                if hasattr(self, key):
                    setattr(self, key, source[key])




class BlockPrepWbnaa(block_prep_identity.BlockPrepIdentity):
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
        self.exclude_indices = None
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
        self.chain = chain_prep_wbnaa.ChainPrepWbnaa(dataset, self)


    def set_dims(self, dataset):
        """
        Given a Dataset object, this is an opportunity for this block object 
        to ensure that its dims match those of the parent dataset. 
        """
        block.Block.set_dims(self, dataset)

        # local reference to input data
        raw = dataset.get_source_data('prep')

        # this is the calculated proper size for self.data
        fidsum_dims = [raw.shape[-1],1,1,1]

        if not self.dims or self.dims != fidsum_dims: 
            self._reset_dimensional_data(dataset)


    def _reset_dimensional_data(self, dataset):
        """
        Resets (to zero) and resizes dimensionally-dependent data
        
        """
        # local reference to input data
        raw = dataset.get_source_data('prep')

        n_fids = raw.shape[-2]

        self.exclude_indices = []
        self.frequency_shift = np.zeros([n_fids])
        self.phase_0         = np.zeros([n_fids])

        self.data = np.zeros((1,1,1,raw.shape[-1]), dtype=raw.dtype)
        if self.chain is not None:
            self.chain.reset_results_arrays()
        
    
    def concatenate(self, new):
        raise NotImplementedError


    def toggle_exclude_index(self, dataset, index):
        """
        Given an index, add it to the exclude data list. If it is already
        in the list, then remove it. If it is outside the list dimensions, 
        ignore it. If index is None, reset the list to empty
        
        """
        if index is None:
            self.exclude_indices = []
            return
        
        # local reference to input data
        raw = dataset.get_source_data('prep')
        nfids = raw.shape[-2]
        
        if index < 0 or index >= nfids: return
        
        if index in self.exclude_indices:
            self.exclude_indices.remove(index)
        else:
            self.exclude_indices.append(index)
            self.exclude_indices.sort()


    def deflate(self, flavor=Deflate.ETREE):
        if flavor == Deflate.ETREE:

            e = Element("block_prep_wbnaa", {"id" : self.id, 
                                             "version" : self.XML_VERSION})

            util_xml.TextSubElement(e, "behave_as_preset", self.behave_as_preset)
            
            # Now I deflate the attribs that are specific to this class
            e.append(self.set.deflate())
            
            if not self.behave_as_preset:

                # These attributes are all python lists.
                for attribute in ("exclude_indices",  ):
                    for value in getattr(self, attribute):
                        util_xml.TextSubElement(e, attribute, value)

                # these attributes are numpy arrays
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
            val_set = source.find("settings")
            self.set = _Settings(val_set)

            if not self.behave_as_preset:

                self.exclude_indices = [int(val.text) for val in source.getiterator("exclude_indices")]
            
                # Now I inflate the attribs that are specific to this class
                temp = source.find("frequency_shift")
                self.frequency_shift = util_xml.element_to_numpy_array(temp)
                temp = source.find("phase_0")
                self.phase_0 = util_xml.element_to_numpy_array(temp)
                temp = source.find("data")
                self.data = util_xml.element_to_numpy_array(temp)


        elif hasattr(source, "keys"):
            # Quacks like a dict
            for key in list(source.keys()):
                if key == "set":
                    setattr(self, key, source[key])


    ##### Private Methods #####################################

    
    
    
  
  
