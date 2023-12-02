# Python modules


# 3rd party modules
import numpy as np
from xml.etree.cElementTree import Element

# Our modules
import vespa.analysis.block_prep_identity as block_prep_identity
import vespa.analysis.chain_prep_fidsum as chain_prep_fidsum
import vespa.analysis.block as block

import vespa.common.util.xml_ as util_xml

from vespa.analysis.mrs_user_prior import UserPrior
from vespa.common.constants import Deflate



class _Settings(object):
    """
    Settings object contains the parameter inputs used for processing in the 
    Chain object in this Block. Having a separate object helps to delineate 
    inputs/outputs and to simplify load/save of preset values.

    This object can also save/recall these values to/from an XML node.

    Version 1.0.0 - only had Vespa-Analysis correction algorithm in it
            1.1.0 - added Suspect and some FID-A algorithms for correction and
                      data exclusion. Renamed most attributes and added lots.

    """
    XML_VERSION = "1.1.0"

    
    def __init__(self, attributes=None):
        """
        The BlockPrepFidsum class  deals with "raw" types of data objects 
        that need to do a bit of massaging of the data as it comes in 
        (e.g. align and sum individual FIDs for an SVS data set).
        
        Attributes 
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
        global_phase1          float, degrees, phase 1 amount added to all FIDs
                                 This value is applied to final data calculation

        
        """
        self.coil_combine_method                = 'Siemens'
        self.coil_combine_external_filename     = ''
        self.coil_combine_external_dataset_id   = ''
        self.coil_combine_external_dataset      = None
        self.coil_combine_noise_whitening       = False

        self.apply_data_exclusion               = False
        self.exclude_method                     = 'Manual'
        self.exclusion_input_adjust             = False
        self.fida_bad_threshold                 = 4.0       # from fid-a run_pressproc_auto.m (and GE auto)
        self.fida_bad_domain                    = 't'       # from fid-a
        self.fida_n_worst                       = 3         # guess by bjs

        self.correction_method                  = 'Manual'
        self.auto_correct                       = False
        self.correction_input_adjust            = False
        self.vespa_reference_peak_center        = 2.01      # was reference_peak_center
        self.vespa_peak_search_width            = 0.2       # was peak_search_width
        self.vespa_phase0_range_start           = 3.5       # was phase0_range_start
        self.vespa_phase0_range_end             = 0.5       # was phase0_range_end
        self.vespa_target_method                = 'Average all'

        self.vespa_preprocess_prior             = UserPrior()          

        self.suspect_initial_guess_freq         = 0.1       # in Hz
        self.suspect_initial_guess_phase        = 0.1       # in deg here, but suspect wants rads
        self.suspect_optimization_range_start   = 3.3       # empirical
        self.suspect_optimization_range_end     = 2.7       # empirical
        self.suspect_target_method              = 'Avg first 4'

        self.rats_initial_guess_freq            = 0.1       # in Hz
        self.rats_initial_guess_phase           = 0.1       # in deg here, but suspect wants rads
        self.rats_optimization_range_start      = 3.3       # empirical
        self.rats_optimization_range_end        = 2.7       # empirical
        self.rats_baseline_order                = 2         # from suspect
        self.rats_target_method                 = 'Avg first 4'

        self.global_left_shift                  = 0         # was fid_left_shift
        self.global_phase0                      = 0.0       # was phase0_global_const
        self.global_phase1                      = 0.0
        self.global_gaussian_apodization        = 0.0       # was gaussian_apodization
        self.chop_data                          = True
        self.zero_phase1                        = True
        self.apply_peak_shift                   = True
        self.apply_phase0                       = True

        if attributes is not None:
            self.inflate(attributes)


    def __str__(self):
        lines = []
        lines.append("------- {0} Object -------".format(self.__class__.__name__))
        lines.append("coil_combine_method               : " + str(self.coil_combine_method))
        lines.append("coil_combine_external_filename    : " + str(self.coil_combine_external_filename))
        lines.append("coil_combine_external_dataset_id  : " + str(self.coil_combine_external_dataset_id))
        lines.append("coil_combine_noise_whitening      : " + str(self.coil_combine_noise_whitening))
        lines.append("apply_data_exclusion              : " + str(self.apply_data_exclusion))
        lines.append("exclude_method                    : " + str(self.exclude_method))
        lines.append("exclusion_input_adjust            : " + str(self.exclusion_input_adjust))
        lines.append("fida_bad_threshold                : " + str(self.fida_bad_threshold))
        lines.append("fida_bad_domain                   : " + str(self.fida_bad_domain))
        lines.append("fida_n_worst                      : " + str(self.fida_n_worst))
        lines.append("correction_method                 : " + str(self.fida_bad_threshold))
        lines.append("auto_correct                      : " + str(self.auto_correct))
        lines.append("correction_input_adjust           : " + str(self.correction_input_adjust))
        lines.append("vespa_reference_peak_center       : " + str(self.vespa_reference_peak_center))
        lines.append("vespa_peak_search_width           : " + str(self.vespa_peak_search_width))
        lines.append("vespa_phase0_range_start          : " + str(self.vespa_phase0_range_start))
        lines.append("vespa_phase0_range_end            : " + str(self.vespa_phase0_range_end))
        lines.append("vespa_target_method               : " + str(self.vespa_target_method))
        lines.append("suspect_initial_guess_freq        : " + str(self.suspect_initial_guess_freq))
        lines.append("suspect_initial_guess_phase       : " + str(self.suspect_initial_guess_phase))
        lines.append("suspect_optimization_range_start  : " + str(self.suspect_optimization_range_start))
        lines.append("suspect_optimization_range_end    : " + str(self.suspect_optimization_range_end))
        lines.append("suspect_target_method             : " + str(self.suspect_target_method))
        lines.append("rats_initial_guess_freq           : " + str(self.rats_initial_guess_freq))
        lines.append("rats_initial_guess_phase          : " + str(self.rats_initial_guess_phase))
        lines.append("rats_optimization_range_start     : " + str(self.rats_optimization_range_start))
        lines.append("rats_optimization_range_end       : " + str(self.rats_optimization_range_end))
        lines.append("rats_baseline_order               : " + str(self.rats_baseline_order))
        lines.append("rats_target_method                : " + str(self.rats_target_method))
        lines.append("global_phase0                     : " + str(self.global_phase0))
        lines.append("global_phase1                     : " + str(self.global_phase1))
        lines.append("global_left_shift                 : " + str(self.global_left_shift))
        lines.append("chop_data                         : " + str(self.chop_data))
        lines.append("zero_phase1                       : " + str(self.zero_phase1))
        lines.append("apply_peak_shift                  : " + str(self.apply_peak_shift))
        lines.append("apply_phase0                      : " + str(self.apply_phase0))

        lines.append(str(self.vespa_preprocess_prior))

        return '\n'.join(lines)


    def deflate(self, flavor=Deflate.ETREE):
        if flavor == Deflate.ETREE:
            e = Element("settings", {"version" : self.XML_VERSION})

            util_xml.TextSubElement(e, "coil_combine_method",               self.coil_combine_method)
            util_xml.TextSubElement(e, "coil_combine_external_filename",    self.coil_combine_external_filename)
            # In the next line, we *have* to save the uuid values from the
            # actual object rather than from the attribute above, in
            # order for the associated dataset uuid to reflect the new id
            # that is given in the top level dataset. Associated datasets are
            # given new temporary uuid values so that if the main dataset is
            # saved and immediately loaded back in, we do not get collisions
            # between the newly opened datasets and already existing ones.
            if self.coil_combine_external_dataset is not None:
                util_xml.TextSubElement(e, "coil_combine_external_dataset_id", self.coil_combine_external_dataset.id)
            util_xml.TextSubElement(e, "coil_combine_noise_whitening",      self.coil_combine_noise_whitening)
            
            util_xml.TextSubElement(e, "apply_data_exclusion",              self.apply_data_exclusion)
            util_xml.TextSubElement(e, "exclude_method",                    self.exclude_method)
            util_xml.TextSubElement(e, "exclusion_input_adjust",            self.exclusion_input_adjust)
            util_xml.TextSubElement(e, "fida_bad_threshold",                self.fida_bad_threshold)
            util_xml.TextSubElement(e, "fida_bad_domain",                   self.fida_bad_domain)
            util_xml.TextSubElement(e, "fida_n_worst",                      self.fida_n_worst)

            util_xml.TextSubElement(e, "correction_method",                 self.correction_method)
            util_xml.TextSubElement(e, "auto_correct",                      self.auto_correct)
            util_xml.TextSubElement(e, "correction_input_adjust",           self.correction_input_adjust)
            util_xml.TextSubElement(e, "vespa_reference_peak_center",       self.vespa_reference_peak_center)
            util_xml.TextSubElement(e, "vespa_peak_search_width",           self.vespa_peak_search_width)
            util_xml.TextSubElement(e, "vespa_phase0_range_start",          self.vespa_phase0_range_start)
            util_xml.TextSubElement(e, "vespa_phase0_range_end",            self.vespa_phase0_range_end)
            util_xml.TextSubElement(e, "vespa_target_method",               self.vespa_target_method)
            util_xml.TextSubElement(e, "suspect_initial_guess_freq",        self.suspect_initial_guess_freq)
            util_xml.TextSubElement(e, "suspect_initial_guess_phase",       self.suspect_initial_guess_phase)
            util_xml.TextSubElement(e, "suspect_optimization_range_start",  self.suspect_optimization_range_start)
            util_xml.TextSubElement(e, "suspect_optimization_range_end",    self.suspect_optimization_range_end)
            util_xml.TextSubElement(e, "suspect_target_method",             self.suspect_target_method)
            util_xml.TextSubElement(e, "rats_initial_guess_freq",           self.rats_initial_guess_freq)
            util_xml.TextSubElement(e, "rats_initial_guess_phase",          self.rats_initial_guess_phase)
            util_xml.TextSubElement(e, "rats_optimization_range_start",     self.rats_optimization_range_start)
            util_xml.TextSubElement(e, "rats_optimization_range_end",       self.rats_optimization_range_end)
            util_xml.TextSubElement(e, "rats_baseline_order",               self.rats_baseline_order)
            util_xml.TextSubElement(e, "rats_target_method",                self.rats_target_method)

            util_xml.TextSubElement(e, "global_left_shift",                 self.global_left_shift)
            util_xml.TextSubElement(e, "global_phase0",                     self.global_phase0)
            util_xml.TextSubElement(e, "global_phase1",                     self.global_phase1)
            util_xml.TextSubElement(e, "global_gaussian_apodization",       self.global_gaussian_apodization)
            util_xml.TextSubElement(e, "chop_data",                         self.chop_data)
            util_xml.TextSubElement(e, "zero_phase1",                       self.zero_phase1)
            util_xml.TextSubElement(e, "apply_peak_shift",                  self.apply_peak_shift)
            util_xml.TextSubElement(e, "apply_phase0",                      self.apply_phase0)

            e.append(self.vespa_preprocess_prior.deflate())

            return e
            
        elif flavor == Deflate.DICTIONARY:
            return self.__dict__.copy()


    def inflate(self, source):
        if hasattr(source, "makeelement"):
            # Quacks like an ElementTree.Element

            for name in ("fida_bad_threshold",
                         "vespa_reference_peak_center",
                         "vespa_peak_search_width",
                         "vespa_phase0_range_start",
                         "vespa_phase0_range_end",
                         "suspect_initial_guess_freq",
                         "suspect_initial_guess_phase",
                         "suspect_optimization_range_start",
                         "suspect_optimization_range_end",
                         "rats_initial_guess_freq",
                         "rats_initial_guess_phase",
                         "rats_optimization_range_start",
                         "rats_optimization_range_end",
                         "global_gaussian_apodization",
                         "global_phase0",
                         "global_phase1",):
                item = source.findtext(name)
                if item is not None:
                    setattr(self, name, float(item))

            for name in ("fida_n_worst",
                         "input_left_shift",
                         "rats_baseline_order",
                         "global_left_shift",):
                item = source.findtext(name)
                if item is not None:
                    setattr(self, name, int(item))
                
            for name in ("coil_combine_noise_whitening", 
                         "apply_data_exclusion",
                         "exclusion_input_adjust",
                         "correction_input_adjust",
                         "auto_correct",
                         "chop_data",
                         "zero_phase1",
                         "apply_peak_shift", 
                         "apply_phase0",):
                item = source.findtext(name)
                if item is not None:
                    setattr(self, name, util_xml.BOOLEANS[item])

            for name in ("coil_combine_method",
                         "coil_combine_external_filename",
                         "coil_combine_external_dataset_id",
                         "exclude_method", 
                         "correction_method",
                         "vespa_target_method",
                         "suspect_target_method",
                         "rats_target_method",):
                item = source.findtext(name)
                if item is not None:
                    setattr(self, name, item)

            # Version 1.0.0 -> 1.1.0 check for old tag names and reassign
            for old in (("gaussian_apodization",  "global_gaussian_apodization", "float"),
                        ("reference_peak_center", "vespa_reference_peak_center", "float"),
                        ("peak_search_width",     "vespa_peak_search_width", "float"),
                        ("phase0_range_start",    "vespa_phase0_range_end", "float"),
                        ("phase0_range_end",      "vespa_phase0_range_start", "float"),
                        ("phase0_global_const",   "global_phase0", "float"),
                        ("fid_left_shift",        "global_left_shift", "int"),):
                old_tag, new_tag, dtype = old
                item = source.findtext(old_tag)
                if item is not None:
                    if dtype == "float":
                        setattr(self, new_tag, float(item))
                    elif dtype == "int":
                        setattr(self, new_tag, int(item))

            self.vespa_preprocess_prior.inflate(source.find("user_prior"))


        elif hasattr(source, "keys"):
            # Quacks like a dict
            for key in list(source.keys()):
                if hasattr(self, key):
                    setattr(self, key, source[key])



class BlockPrepFidsum(block_prep_identity.BlockPrepIdentity):
    """ 
    Building block to hold the state of a step in an MRS processing chain.
    Includes the functionality to save/recall this object to/from an XML node.

    Contains inputs/results for preprocessing of the raw data from the previous
    block ('raw') in the dataset.blocks list. This step modifies coil/average
    data into a single summed FID array for one dataset.

    Version 1.0.0 - only had Vespa-Analysis correction algorithm in it
            1.1.0 - See Settings for algorithm changes. At this level the major
                      change is that the results for coil_combine_weights, and 
                      coil_combine_phases were moved from Settings to here.
    """
    XML_VERSION = "1.1.0"

    def __init__(self, attributes=None):
        """
        Block objects have a self.set attribute that contains a _Settings object
        that contains the input attributes for the processing done in the Chain.
        This simplifies using a Block object as a preset. Results from the Chain
        are stored in this object at the level of self.set

        Base class sets references to:  self.id, self.data, self.chain and self.behave_as_preset

        """
        super().__init__(attributes)
        
        # input parameters for Chain processing
        self.set = _Settings()

        # results storage
        self.coil_combine_weights   = None
        self.coil_combine_phases    = None
        self.exclude_indices        = []
        self.frequency_shift        = None
        self.phase_0                = None
        self.data                   = None
                                 
        if attributes is not None:
            self.inflate(attributes)



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
        self.chain = chain_prep_fidsum.ChainPrepFidsum(dataset, self)


    def get_associated_datasets(self, is_main_dataset=True):
        """
        Returns a list of datasets associated with this object. 
        
        The 'is_main_dataset' flag allows the method to know if it is the top
        level dataset gathering associated datasets, or some dataset that is
        only associated with the top dataset. This is used to stop circular
        logic conditions where one or more datasets refer to each other.
        """
        # Call base class first
        datasets = block_prep_identity.BlockPrepIdentity.get_associated_datasets(self, is_main_dataset)
        
        if self.set.coil_combine_external_dataset:
            if self.set.coil_combine_method == 'External Dataset':
                datasets += self.set.coil_combine_external_dataset.get_associated_datasets(is_main_dataset=False)                
                datasets += [self.set.coil_combine_external_dataset]
        else:
            return []

        return datasets


    def set_associated_datasets(self, datasets):                          
        for dataset in datasets:
            if dataset.id == self.set.coil_combine_external_dataset_id:
                self.set.coil_combine_external_dataset = dataset
    

    def set_dims(self, dataset):
        """
        Given a Dataset object, this is an opportunity for this block object 
        to ensure that its dims match those of the parent dataset. 

        """
        block.Block.set_dims(self, dataset)
        raw = dataset.get_source_data('prep')       # local reference to input data

        if not self.dims or self.shape != [1,1,1,raw.shape[-1]]:
            self._reset_dimensional_data(dataset)


    def _reset_dimensional_data(self, dataset):
        """
        Resets (to zero) and resizes dimensionally-dependent data
        
        """
        # local reference to input data
        raw = dataset.get_source_data('prep')

        ncoil = raw.shape[1]
        nfids = raw.shape[2]

        self.coil_combine_weights = np.zeros([ncoil], dtype=np.float32)
        self.coil_combine_phases  = np.zeros([ncoil], dtype=np.float32)
        self.exclude_indices = []
        self.frequency_shift = np.zeros([nfids], dtype=np.float32)
        self.phase_0         = np.zeros([nfids], dtype=np.float32)

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

        raw = dataset.get_source_data('prep')   # local reference to input data
        nfids = raw.shape[-2]
        
        if index < 0 or index >= nfids: return
        
        if index in self.exclude_indices:
            self.exclude_indices.remove(index)
        else:
            self.exclude_indices.append(index)
            self.exclude_indices.sort()
        

    def deflate(self, flavor=Deflate.ETREE):
        if flavor == Deflate.ETREE:

            e = Element("block_prep_fidsum",{"id" : self.id,
                                             "version" : self.XML_VERSION})

            util_xml.TextSubElement(e, "behave_as_preset", self.behave_as_preset)

            if self.behave_as_preset:
                # presets are only input parameters, set results temporarily to None
                save_filename  = self.set.coil_combine_external_filename
                save_datasetid = self.set.coil_combine_external_dataset_id
                save_dataset   = self.set.coil_combine_external_dataset

                self.set.coil_combine_external_filename     = ''
                self.set.coil_combine_external_dataset_id   = ''
                self.set.coil_combine_external_dataset      = None

            # Now I deflate the attribs that are specific to this class
            e.append(self.set.deflate())

            if self.behave_as_preset:
                self.set.coil_combine_external_filename   = save_filename
                self.set.coil_combine_external_dataset_id = save_datasetid
                self.set.coil_combine_external_dataset    = save_dataset

            if not self.behave_as_preset:

                # These attributes are all python lists.
                for attribute in ("exclude_indices",  ):
                    for value in getattr(self, attribute):
                        util_xml.TextSubElement(e, attribute, value)

                # these attributes are numpy arrays
                if self.coil_combine_weights is not None:
                    e.append(util_xml.numpy_array_to_element(self.coil_combine_weights, 'coil_combine_weights'))
                if self.coil_combine_phases is not None:
                    e.append(util_xml.numpy_array_to_element(self.coil_combine_phases, 'coil_combine_phases'))
                if self.frequency_shift is not None:
                    e.append(util_xml.numpy_array_to_element(self.frequency_shift,'frequency_shift'))
                if self.phase_0 is not None:
                    e.append(util_xml.numpy_array_to_element(self.phase_0,'phase_0'))
                if self.data is not None:
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

                val = source.find("coil_combine_weights")
                if val is not None:
                    self.coil_combine_weights = util_xml.element_to_numpy_array(val)
                else:
                    # backward compatibility for Version 1.0.0
                    val = val_set.find("coil_combine_weights")
                    if val is not None:
                        self.coil_combine_weights = util_xml.element_to_numpy_array(val)

                val = source.find("coil_combine_phases")
                if val is not None:
                    self.coil_combine_phases = util_xml.element_to_numpy_array(val)
                else:
                    # backward compatibility for Version 1.0.0
                    val = val_set.find("coil_combine_phases")
                    if val is not None:
                        self.coil_combine_phases = util_xml.element_to_numpy_array(val)

                # lists
                self.exclude_indices = [val.text for val in source.getiterator("exclude_indices")]

                # Now I inflate the attribs that are specific to this class
                temp = source.find("frequency_shift")
                if temp is not None:
                    self.frequency_shift = util_xml.element_to_numpy_array(temp)
                temp = source.find("phase_0")
                if temp is not None:
                    self.phase_0 = util_xml.element_to_numpy_array(temp)
                temp = source.find("data")
                if temp is not None:
                    self.data = util_xml.element_to_numpy_array(temp)

        elif hasattr(source, "keys"):
            # Quacks like a dict
            for key in list(source.keys()):
                if key == "set":
                    setattr(self, key, source[key])


    ##### Private Methods #####################################

    
