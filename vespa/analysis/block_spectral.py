# Python modules
from functools import reduce

# 3rd party modules
import numpy as np
from xml.etree.cElementTree import Element

# Our modules
import vespa.analysis.block_spectral_identity as block_spectral_identity
import vespa.analysis.chain_spectral as chain_spectral
import vespa.analysis.svd_output as svd_output_module
import vespa.analysis.functors.funct_water_filter as funct_water
import vespa.common.util.xml_ as util_xml
import vespa.common.util.generic_spectral as util_spectral
from vespa.common.constants import Deflate



def _group_svd_outputs(frequencies, damping_factors, amplitudes, phases, in_model):
    """
    Given five 3D arrays with the same shape, returns a numpy object array
    with the same shape. Each array element is an SvdOutput instance with the
    frequencies, damping factors, etc. copied from the appropriate element of
    the corresponding array parameter.

    """
    shape = frequencies.shape

    # This line multiplies the elements of shape together; e.g. (3, 4, 5) = 60.
    size = reduce(lambda x, y: x * y, shape)

    # Start with a flat array; it's easier.
    svd_outputs = np.empty(size, dtype=object)

    # Here's some clever code for quickly and concisely iterating over a 3D
    # list or numpy array. Note that izip is from itertools in the Python
    # standard lib. I swiped this solution from here:
    # http://stackoverflow.com/questions/189087/how-can-i-in-python-iterate-over-multiple-2d-lists-at-once-cleanly
    i = 0
    for f1, d1, a1, p1, i1 in zip(frequencies, damping_factors, amplitudes, phases, in_model):
        for f2, d2, a2, p2, i2 in zip(f1, d1, a1, p1, i1):
            for f3, d3, a3, p3, i3 in zip(f2, d2, a2, p2, i2):
                svd_outputs[i] = svd_output_module.SvdOutput(f3, d3, a3, p3, i3)
                i += 1

    # Now give the array its proper shape.
    svd_outputs.shape = shape

    return svd_outputs


def _split_svd_outputs(svd_outputs):
    """
    Given the 3D array of SVD output from BlockSpectral, returns five
    3D arrays of the same shape representing the SVD output components
    (frequencies, damping factors, amplitudes, phases, and 'in model' flags).

    """
    shape = svd_outputs.shape

    frequencies = np.zeros(shape, dtype=object)
    damping_factors = np.zeros(shape, dtype=object)
    amplitudes = np.zeros(shape, dtype=object)
    phases = np.zeros(shape, dtype=object)
    in_model = np.zeros(shape, dtype=object)

    for position, svd_output in np.ndenumerate(svd_outputs):
        frequencies[position] = svd_output.frequencies
        damping_factors[position] = svd_output.damping_factors
        amplitudes[position] = svd_output.amplitudes
        phases[position] = svd_output.phases
        in_model[position] = svd_output.in_model

    return frequencies, damping_factors, amplitudes, phases, in_model





class _Settings(object):
    """
    Settings object contains the parameter inputs used for processing in the 
    Chain object in this Block. Having a separate object helps to delineate 
    inputs/outputs and to simplify load/save of preset values.

    This object can also save/recall these values to/from an XML node.

    Version 1.0.0 - had 'kiss off' attrib naming and separate water/ecc filter objects
            1.2.0 - removed separate water/ecc filter objects, new attibute names
            1.3.0 - removed 'kiss off' attib names, deprecated water/ecc filter conversion

    """
    XML_VERSION = "1.3.0"

    def __init__(self, attributes=None):
        """
        Most of these values appear in a GUI when the app first starts.

        SPECTRAL Processing Parameters
        -----------------------------------------

        flip - Reverse the X-axis
        chop - Spectral chop data
        fft  - Apply the FFT?

        zero_fill_multiplier    1 <= zero_fill_multiplier <= 32, must be an
                                  integral power of 2. Valid values are
                                  1, 2, 4, 8, 16 and 32.
        apodization             pulldown menu, algorithm for filter
        apodization_width       float, broadening width in Hz
                                  0 <= apodization_width <= 100; in Hertz
        amplitude               relative scaling, 0 <= amplitude <= 1e12
        phase_1_pivot           in ppm, -1000 <= pivot <= 1000
        dc_offset               offset in vertical values,
                                  -1e5 <= dc_offset <= 1e5
        left_shift_correct      Boolean -- correct spectrum for missing points
                                  at beginning of FID
        left_shift_value        Number of points to drop from the start of
                                  the FID before FFT
                                  0 <= left_shift_value <= acqdim
        svd_apply_threshold     Boolean -- use the svd_threshold value to set
                                  which results peaks check boxes are set ON
        svd_threshold           int, [Hz] offset from on-resonance (4.7 ppm
                                  for 1H) below which we include lines to remove
                                  from the data in an SVD filter operation
        
        svd_last_n_singular_values
        svd_last_n_data_points  ints, these save the last values selected by 
                                the user on the GUI sliders. They serve to 
                                allow the user to change the default value 
                                for the slider when we save a Preset file.
                                A de novo object will have funct_water_filter
                                default values. But with each slider GUI 
                                change, these attributes are updated. When 
                                a Preset file is read in, these attributes
                                are set to last user value, and THEN the 
                                default values in the SVD results arrays are 
                                set using these values. 
        """
        #------------------------------------------------------------
        # Spectral Processing Variables
        #------------------------------------------------------------

        # General Parameters
        self.flip                    = False
        self.chop                    = True
        self.fft                     = True
        self.zero_fill_multiplier    = 1
        self.apodization             = 'gaussian'
        self.apodization_width       = 1.0           # Hz
        self.amplitude               = 1.0
        self.phase_1_pivot           = 2.01          # PPM
        self.dc_offset               = 0.0
        self.left_shift_value        = 0             # was kiss_off_point_count v1.2.0
        self.left_shift_correct      = False         # was kiss_off_correction  v1.2.0
        # SVD Filter parameters
        self.svd_apply_threshold     = True
        self.svd_threshold           = 30            # value
        self.svd_threshold_unit      = 'Hz'          # units for svd_threshold, 'Hz' or 'PPM'
        self.svd_exclude_lipid       = False
        self.svd_exclude_lipid_start = 1.9
        self.svd_exclude_lipid_end   = 0.0     

        # ecc filter attributes
        self.ecc_method              = 'None'
        self.ecc_filename            = ''
        self.ecc_dataset_id          = ''
        self.ecc_dataset             = None
        self.ecc_raw                 = None

        # water filter attributes
        self.water_filter_method            = 'None'
        self.fir_length                     = funct_water.FIR_LENGTH_DEFAULT
        self.fir_half_width                 = funct_water.FIR_RIPPLE_DEFAULT
        self.fir_ripple                     = funct_water.FIR_HALF_WIDTH_DEFAULT
        self.fir_extrapolation_method       = funct_water.FIR_EXTRAPOLATION_DEFAULT
        self.fir_extrapolation_point_count  = funct_water.FIR_EXTRAPOLATION_POINTS_DEFAULT
        self.ham_length                     = funct_water.HAM_LENGTH_DEFAULT
        self.ham_extrapolation_method       = funct_water.HAM_EXTRAPOLATION_DEFAULT
        self.ham_extrapolation_point_count  = funct_water.HAM_EXTRAPOLATION_POINTS_DEFAULT
        self.svd_last_n_data_points         = funct_water.SVD_N_DATA_POINTS
        self.svd_last_n_singular_values     = funct_water.SVD_N_SINGULAR_VALUES


        if attributes is not None:
            self.inflate(attributes)


    ##### Standard Methods and Properties #####################################

    def __str__(self):
        return self.__unicode__()


    def __unicode__(self):
        lines = [ ]
        lines.append("--- Block Spectral Settings ---")
        lines.append("flip                    : " + str(self.flip))
        lines.append("chop                    : " + str(self.chop))
        lines.append("fft                     : " + str(self.fft))

        lines.append("zero_fill_multiplier    : " + str(self.zero_fill_multiplier))
        lines.append("apodization             : " + str(self.apodization))
        lines.append("apodization_width       : " + str(self.apodization_width))
        lines.append("amplitude               : " + str(self.amplitude))
        lines.append("phase_1_pivot           : " + str(self.phase_1_pivot))
        lines.append("dc_offset               : " + str(self.dc_offset))
        lines.append("left_shift_value        : " + str(self.left_shift_value))
        lines.append("left_shift_correct      : " + str(self.left_shift_correct))
        lines.append("svd_apply_threshold     : " + str(self.svd_apply_threshold))
        lines.append("svd_threshold           : " + str(self.svd_threshold))
        lines.append("svd_threshold_unit      : " + str(self.svd_threshold_unit))
        lines.append("svd_exclude_lipid       : " + str(self.svd_exclude_lipid))
        lines.append("svd_exclude_lipid_start : " + str(self.svd_exclude_lipid_start))
        lines.append("svd_exclude_lipid_end   : " + str(self.svd_exclude_lipid_end))

        lines.append("ecc_method              : " + str(self.ecc_method))
        lines.append("ecc_filename            : " + str(self.ecc_filename))
        lines.append("ecc_dataset_id          : " + str(self.ecc_dataset_id))
        if self.ecc_raw is not None:
            lines.append("ecc_raw (length)        : " + str(len(self.ecc_raw)))
        else:
            lines.append("ecc_raw (length)        : None")

        lines.append("water_filter_method           : " + str(self.water_filter_method))
        lines.append("fir_length                    : " + str(self.fir_length))
        lines.append("fir_half_width                : " + str(self.fir_half_width))
        lines.append("fir_ripple                    : " + str(self.fir_ripple))
        lines.append("fir_extrapolation_method      : " + str(self.fir_extrapolation_method))
        lines.append("fir_extrapolation_point_count : " + str(self.fir_extrapolation_point_count))
        lines.append("ham_length                    : " + str(self.ham_length))
        lines.append("ham_extrapolation_method      : " + str(self.ham_extrapolation_method))
        lines.append("ham_extrapolation_point_count : " + str(self.ham_extrapolation_point_count))
        lines.append("svd_last_n_data_points        : " + str(self.svd_last_n_data_points))
        lines.append("svd_last_n_singular_values    : " + str(self.svd_last_n_singular_values))

        # __unicode__() must return a Unicode object. In practice the code
        # above always generates Unicode, but we ensure it here.
        return '\n'.join(lines)


    def deflate(self, flavor=Deflate.ETREE):
        if flavor == Deflate.ETREE:
            e = Element("settings", {"version" : self.XML_VERSION})

            util_xml.TextSubElement(e, "flip",                    self.flip)
            util_xml.TextSubElement(e, "chop",                    self.chop)
            util_xml.TextSubElement(e, "fft",                     self.fft)
            util_xml.TextSubElement(e, "zero_fill_multiplier",    self.zero_fill_multiplier)
            util_xml.TextSubElement(e, "apodization",             self.apodization)
            util_xml.TextSubElement(e, "apodization_width",       self.apodization_width)
            util_xml.TextSubElement(e, "amplitude",               self.amplitude)
            util_xml.TextSubElement(e, "phase_1_pivot",           self.phase_1_pivot)
            util_xml.TextSubElement(e, "dc_offset",               self.dc_offset)
            util_xml.TextSubElement(e, "left_shift_value",        self.left_shift_value)
            util_xml.TextSubElement(e, "left_shift_correct",      self.left_shift_correct)
            util_xml.TextSubElement(e, "svd_apply_threshold",     self.svd_apply_threshold)
            util_xml.TextSubElement(e, "svd_threshold",           self.svd_threshold)
            util_xml.TextSubElement(e, "svd_threshold_unit",      self.svd_threshold_unit)
            util_xml.TextSubElement(e, "svd_exclude_lipid",       self.svd_exclude_lipid)
            util_xml.TextSubElement(e, "svd_exclude_lipid_start", self.svd_exclude_lipid_start)
            util_xml.TextSubElement(e, "svd_exclude_lipid_end",   self.svd_exclude_lipid_end)

            util_xml.TextSubElement(e, "ecc_method",               self.ecc_method)
            util_xml.TextSubElement(e, "ecc_filename",             self.ecc_filename)
            
            # We *have* to save uuid from actual object versus the attribute above, in order for 
            # the associated dataset uuid to reflect its new_id. Associated datasets are given 
            # new temporary uuid to avoid collisions if, say, the main dataset is saved and then
            # immediately loaded back in.
            if self.ecc_dataset is not None:
                util_xml.TextSubElement(e, "ecc_dataset_id",        self.ecc_dataset.id)
            if self.ecc_raw is not None:
                e.append(util_xml.numpy_array_to_element(self.ecc_raw, 'ecc_raw'))

            util_xml.TextSubElement(e, "water_filter_method",           self.water_filter_method)
            util_xml.TextSubElement(e, "fir_length",                    self.fir_length)
            util_xml.TextSubElement(e, "fir_half_width",                self.fir_half_width)
            util_xml.TextSubElement(e, "fir_ripple",                    self.fir_ripple)
            util_xml.TextSubElement(e, "fir_extrapolation_method",      self.fir_extrapolation_method)
            util_xml.TextSubElement(e, "fir_extrapolation_point_count", self.fir_extrapolation_point_count)
            util_xml.TextSubElement(e, "ham_length",                    self.ham_length)
            util_xml.TextSubElement(e, "ham_extrapolation_method",      self.ham_extrapolation_method)
            util_xml.TextSubElement(e, "ham_extrapolation_point_count", self.ham_extrapolation_point_count)
            util_xml.TextSubElement(e, "svd_last_n_data_points",        self.svd_last_n_data_points)
            util_xml.TextSubElement(e, "svd_last_n_singular_values",    self.svd_last_n_singular_values)

            return e

        elif flavor == Deflate.DICTIONARY:
            return self.__dict__.copy()


    def inflate(self, source):
        if hasattr(source, "makeelement"):
            # Quacks like an ElementTree.Element

            # in some cases below, we need to check if a value returns None
            # because in mrs_dataset versions prior to 1.1.0, HLSVD was its own
            # tab not a sub-tab of spectral, thus some of these attributes were
            # not present in the settings

            for item in ("flip",
                         "chop",
                         "fft",
                         "left_shift_correct",
                         "svd_apply_threshold",
                         "svd_exclude_lipid" ):
                val = source.findtext(item)
                if val:
                    setattr(self, item, util_xml.BOOLEANS[val])

            for item in ("apodization",):
                setattr(self, item, source.findtext(item))      # this code allows 'None' value

            for item in ("apodization_width",
                         "amplitude",
                         "phase_1_pivot",
                         "dc_offset",
                         "svd_threshold",
                         "fir_half_width",
                         "fir_ripple",
                         "svd_exclude_lipid_start",
                         "svd_exclude_lipid_end",
                         ):
                val = source.findtext(item)
                if val:
                    setattr(self, item, float(val))

            for item in ("zero_fill_multiplier",
                         "left_shift_value",
                         "fir_length",
                         "fir_extrapolation_point_count",
                         "ham_length",
                         "ham_extrapolation_point_count",
                         "svd_last_n_data_points",
                         "svd_last_n_singular_values" ):
                val = source.findtext(item)
                if val:
                    setattr(self, item, int(val))

            val = source.find("ecc_raw")
            if val is not None:
                self.ecc_raw = util_xml.element_to_numpy_array(val)

            for item in ("ecc_method",
                         "ecc_filename",
                         "ecc_dataset_id",
                         "water_filter_method",
                         "fir_extrapolation_method",
                         "ham_extrapolation_method",
                         "svd_threshold_unit"):
                val = source.findtext(item)
                if val is not None:
                    setattr(self, item, val)

            # As of VERSION 1.2.0, ECC and Water filters were no longer their own
            # object. As of VERSION 1.3.0 we no longer find/convert these objects.


            # Version 1.2.0 -> 1.3.0 check for old tag names and reassign
            for old in (("kiss_off_point_count", "left_shift_value",   "int"),
                        ("kiss_off_correction",  "left_shift_correct", "bool"),):
                old_tag, new_tag, dtype = old
                item = source.findtext(old_tag)
                if item is not None:
                    if dtype == "float":
                        setattr(self, new_tag, float(item))
                    elif dtype == "int":
                        setattr(self, new_tag, int(item))
                    elif dtype == 'bool':
                        setattr(self, new_tag, util_xml.BOOLEANS[item])



        elif hasattr(source, "keys"):
            # Quacks like a dict
            for key in list(source.keys()):
                if hasattr(self, key):
                    setattr(self, key, source[key])




class BlockSpectral(block_spectral_identity.BlockSpectralIdentity):
    """
    This is a building block object that can be used to create a list of
    processing blocks.

    This object represents the settings and results involved in processing
    data from the Time Domain to the Frequency Domain (TDFD) for the spectral
    dimension.

    In here we also package all the functionality needed to save and recall
    these values to/from an XML node.

    """
    # The XML_VERSION enables us to change the XML output format in the future
    XML_VERSION = "1.3.0"

    def __init__(self, attributes=None):
        """

        Parameters:

            id (str): A permanent, unique identifying string for this object.
                Typically serves as a "source_id" for some other object. It
                is part of the provenance for this processing functor chain


        SPECTRAL Processing Parameters
        -----------------------------------------

        set                 Settings object

        Initialization of _phase_0, _phase_1 and _frequency_shift is
        handled in the call to _reset_dimensional_data().

        _phase_0            Phase 0 for all voxels, 0 <= phase 0 <= 360, in degrees
        _phase_1            Phase 1 for all voxels, -10000.0 <= phase 1 <= 10000.0
        _frequency_shift    Frequency shift, -10000 <= frequency shift <= 10000
        data                Complex data after spectral processing

        As of Vespa >= 0.6.0, the SVD filter is part of the spectral block.

        # Singular value decomposition filtering (SVD)
        # - original algorithm from Dutch Hankel-Lanczos HLSVD
        # - more recently from Laudadio/Sima/Munk HLSVD-PRO

        See hlsvdpro/hlsvd.py for limits.

        _data_point_count            (replaces IDL variable ndp)
        _signal_singular_value_count (replaces IDL variable nssv)

        _water_filter_do_fit         Determines whether to do hlsvd when
                                      water removal is applied

        _svd_outputs                 3D array of results from HLSVD algorithm

        """
        super().__init__(attributes)

        # processing parameters
        self.set = _Settings()

        # svd algorithm settings
        # When populated, each of these is a 3D numpy array in the shape of
        # the last 3 spectral dims. For single voxel data, this shape is
        # always (1, 1, 1). The value at each array index ([x, y, z]) is
        # a scalar. (_water_filter_do_fit contains booleans, the others
        # contain ints.)
        self._data_point_count            = None
        self._signal_singular_value_count = None
        self._water_filter_do_fit         = None

        # flags for widget value locks
        self.frequency_shift_lock   = False
        self.phase_lock             = False
        self.phase_1_lock_at_zero   = True
        self.left_shift_correct     = False

        # dimensional data - list of lists
        self._phase_0               = None
        self._phase_1               = None
        self._frequency_shift       = None

        # results storage
        self.data = None

        # _svd_outputs is a 3D numpy object array shaped like the spectral
        # dims. Each array element is an instance of the SvdOutput class.
        self._svd_outputs = None

        if attributes is not None:
            self.inflate(attributes)

            # Prior to version 0.6.0, the HLSVD tab was separate and not
            # always present. As a result, the HLSVD-related values below
            # don't always get set during inflate. If that's the case, we
            # set them to reasonable defaults here.
            # We have to test explicitly for None because if it is not None
            # it is a numpy array, and numpy arrays don't like being cast
            # to bool.
            if self._svd_outputs is None and self.data is not None:
                self._set_default_svd_inputs()
                self._set_default_svd_outputs()

        self.chain = None


    ##### Standard Methods and Properties #####################################

    @property
    def dims(self):
        """Data dimensions in a list, read only."""
        return list(self.data.shape[::-1]) if self.data is not None else None


    def __str__(self):
        lines = []
        lines.append("------- {0} Object -------".format(self.__class__.__name__))
        lines.append("\n")
        lines += _Settings.__str__(self).split('\n')
        lines.append("\n")
        lines.append("------- Main Object -------")
        lines.append("No additional data")
        return lines


    def create_chain(self, dataset):
        self.chain = chain_spectral.ChainSpectral(dataset, self)


    def set_dims(self, dataset):
        """
        Given a Dataset object, this is an opportunity for this block object
        to ensure that its dims match those of the parent dataset.
        """
        if not self.dims or self.dims != dataset.spectral_dims:
            self._reset_dimensional_data(dataset)
            if self.chain is not None:
                self.chain.reset_results_arrays()


    def get_phase_0(self, xyz):
        """ Returns 0th order phase for the voxel at the xyz tuple """
        x, y, z = xyz
        return self._phase_0[x, y, z]

    def set_phase_0(self, phase_0, xyz):
        """ Sets 0th order phase for the voxel at the xyz tuple """
        x, y, z = xyz
        if self.phase_lock:
            self._phase_0 += (phase_0 - self._phase_0[x, y, z])
            self._phase_0 = (self._phase_0+180)%360-180
            #self._phase_0 %= 360
        else:
            phase_0 = (phase_0+180)%360-180
            self._phase_0[x, y, z] = phase_0


    def get_phase_1(self, xyz):
        """ Returns 1st order phase for the voxel at the xyz tuple """
        x, y, z = xyz
        return self._phase_1[x, y, z]

    def set_phase_1(self, phase_1, xyz):
        """ Sets 1st order phase for the voxel at the xyz tuple """
        x, y, z = xyz
#        phase_1  = util_spectral.clip(phase_1, -1e5, 1e5)
        phase_1  = np.clip(phase_1, -1e5, 1e5)

        if self.phase_lock:
            if self.phase_1_lock_at_zero:
                self._phase_1 *= 0.0
            else:
                self._phase_1 += (phase_1 - self._phase_1[x, y, z])
        else:
            if self.phase_1_lock_at_zero:
                self._phase_1[x, y, z] = 0.0
            else:
                self._phase_1[x, y, z] = phase_1


    def get_frequency_shift(self, xyz):
        """ Returns frequency_shift for the voxel at the xyz tuple """
        x, y, z = xyz
        return self._frequency_shift[x, y, z]

    def set_frequency_shift(self, frequency_shift, xyz):
        """ Sets frequency_shift for the voxel at the xyz tuple """
        x, y, z = xyz
        if self.frequency_shift_lock:
            self._frequency_shift += (frequency_shift - self._frequency_shift[x, y, z])
        else:
            self._frequency_shift[x, y, z] = frequency_shift


    def get_data_point_count(self, xyz):
        """ Returns data_point_count for the voxel at the xyz tuple """
        x, y, z = xyz
        return self._data_point_count[x, y, z]

    def set_data_point_count(self, data_point_count, xyz):
        """ Sets data_point_count for the voxel at the xyz tuple """
        x, y, z = xyz
        self._data_point_count[x, y, z] = data_point_count


    def get_signal_singular_value_count(self, xyz):
        """ Returns singular_value_count for the voxel at the xyz tuple """
        x, y, z = xyz
        return self._signal_singular_value_count[x, y, z]

    def set_signal_singular_value_count(self, singular_value_count, xyz):
        """ Sets singular_value_count for the voxel at the xyz tuple """
        x, y, z = xyz
        self._signal_singular_value_count[x, y, z] = singular_value_count


    def get_do_fit(self, xyz):
        """ Returns do_fit for the voxel at the xyz tuple """
        x, y, z = xyz
        return self._water_filter_do_fit[x, y, z]

    def set_do_fit(self, do_fit, xyz):
        """ Sets do_fit for the voxel at the xyz tuple """
        x, y, z = xyz
        self._water_filter_do_fit[x, y, z] = do_fit


    def get_svd_output(self, xyz):
        """Returns SVD output for the voxel at the xyz tuple. See the SvdOutput
        class for details."""
        x, y, z = xyz
        return self._svd_outputs[x, y, z]


    def set_svd_output(self, svd_output, xyz):
        """Given an instance of SvdOutput, saves that instance for the voxel at
        the xyz tuple."""
        x, y, z = xyz
        self._svd_outputs[x, y, z] = svd_output


    def get_associated_datasets(self, is_main_dataset=True):
        """
        Returns a list of datasets associated with this object

        'is_main_dataset' signals that this is the top level dataset gathering 
        associated datasets, and is used to stop circular references

        """
        datasets = block_spectral_identity.BlockSpectralIdentity.get_associated_datasets(self, is_main_dataset)

        if self.set.ecc_dataset:
            if self.set.ecc_method != 'None':
                datasets += self.set.ecc_dataset.get_associated_datasets(is_main_dataset=False)
                datasets += [self.set.ecc_dataset]
        else:
            return []

        return datasets


    def set_associated_datasets(self, datasets):
        """
        When we open a VIFF format file, main._import_file() calls this method
        to parse/store any datasets associated with this one as described below.
        
        """
        for dataset in datasets:
            if dataset.id == self.set.ecc_dataset_id:
                self.set.ecc_dataset = dataset



    def deflate(self, flavor=Deflate.ETREE):
        
        if flavor == Deflate.ETREE:
            e = Element("block_spectral", {"id" : self.id, 
                                           "version" : self.XML_VERSION})

            util_xml.TextSubElement(e, "behave_as_preset", self.behave_as_preset)

            if self.behave_as_preset:
                # presets are for generalized parameters only, We set the ECC
                # filter file reference temporarily to None and then back.
                save_ecc_filename               = self.set.ecc_filename
                save_ecc_datasetid              = self.set.ecc_dataset_id
                save_ecc_dataset                = self.set.ecc_dataset
                save_ecc_raw                    = self.set.ecc_raw
                self.set.ecc_filename           = ''
                self.set.ecc_dataset_id         = ''
                self.set.ecc_dataset            = None
                self.set.ecc_raw                = None

            e.append(self.set.deflate())

            if self.behave_as_preset:
                self.set.ecc_filename           = save_ecc_filename
                self.set.ecc_dataset_id         = save_ecc_datasetid
                self.set.ecc_dataset            = save_ecc_dataset
                self.set.ecc_raw                = save_ecc_raw

            util_xml.TextSubElement(e, "frequency_shift_lock", self.frequency_shift_lock)
            util_xml.TextSubElement(e, "phase_lock",           self.phase_lock)
            util_xml.TextSubElement(e, "phase_1_lock_at_zero", self.phase_1_lock_at_zero)

            e.append(util_xml.numpy_array_to_element(self._phase_0,'phase_0'))
            e.append(util_xml.numpy_array_to_element(self._phase_1,'phase_1'))
            e.append(util_xml.numpy_array_to_element(self._frequency_shift,'frequency_shift'))

            if not self.behave_as_preset:

                for dim in self.dims:
                    util_xml.TextSubElement(e, "dim", dim)

                e.append(util_xml.numpy_array_to_element(self.data, 'data'))

                e.append(util_xml.numpy_array_to_element(self._data_point_count, 'data_point_count'))
                e.append(util_xml.numpy_array_to_element(self._signal_singular_value_count, 'signal_singular_value_count'))
                e.append(util_xml.numpy_array_to_element(self._water_filter_do_fit, 'water_filter_do_fit'))

                # We write SVD output as correlated 3D arrays to save some
                # space in the XML (as opposed to deflating each element in the
                # 3D self._svd_outputs array).

                # rbt = Really Big Tuple
                rbt = _split_svd_outputs(self._svd_outputs)
                frequencies, damping_factors, amplitudes, phases, in_model = rbt

                e.append(util_xml.array3d_to_element(frequencies, 'frequencies'))
                e.append(util_xml.array3d_to_element(damping_factors, 'damping_factors'))
                e.append(util_xml.array3d_to_element(amplitudes, 'amplitudes'))
                e.append(util_xml.array3d_to_element(phases, 'phases'))
                e.append(util_xml.array3d_to_element(in_model, 'index'))

            return e

        elif flavor == Deflate.DICTIONARY:
            return self.__dict__.copy()


    def inflate(self, source):
        if hasattr(source, "makeelement"):
            # Quacks like an ElementTree.Element

            self.id = source.get("id")

            val = source.findtext("behave_as_preset")   # default is False
            if val is not None:
                self.behave_as_preset = util_xml.BOOLEANS[val]

            # Look for settings under the old name as well as the standard name.
            self.set = util_xml.find_settings(source, "block_spectral_settings")
            self.set = _Settings(self.set)

            self.frequency_shift_lock = util_xml.BOOLEANS[source.findtext("frequency_shift_lock")]
            self.phase_lock           = util_xml.BOOLEANS[source.findtext("phase_lock")]
            self.phase_1_lock_at_zero = util_xml.BOOLEANS[source.findtext("phase_1_lock_at_zero")]

            temp = source.find("phase_0")
            self._phase_0 = util_xml.element_to_numpy_array(temp)
            temp = source.find("phase_1")
            self._phase_1 = util_xml.element_to_numpy_array(temp)
            temp = source.find("frequency_shift")
            self._frequency_shift = util_xml.element_to_numpy_array(temp)

            if not self.behave_as_preset:

                temp = source.find("data")
                self.data = util_xml.element_to_numpy_array(temp)

                # in the code below, we need to check if a value returns None
                # because in mrs_dataset versions prior to 1.1.0, HLSVD was its own
                # tab not a sub-tab of spectral, thus some of these results were
                # not present in the block_spectral

                temp = source.find("data_point_count")
                if temp is not None:
                    self._data_point_count = util_xml.element_to_numpy_array(temp)
                temp = source.find("signal_singular_value_count")
                if temp is not None:
                    self._signal_singular_value_count = util_xml.element_to_numpy_array(temp)
                temp = source.find("water_filter_do_fit")
                if temp is not None:
                    self._water_filter_do_fit = util_xml.element_to_numpy_array(temp)

                # As of Vespa >= 0.6.0, there's always an SVD subtab, so SVD
                # output info is always written to the XML. However, prior to
                # that the SVD tab was optional. That means that we're not
                # guaranteed to find SVD nodes in the XML.
                frequencies = source.find("frequencies")
                if frequencies is not None:
                    # If frequencies is present, so are all the others.
                    frequencies = util_xml.element_to_array3d(frequencies)

                    damping_factors = source.find("damping_factors")
                    damping_factors = util_xml.element_to_array3d(damping_factors)

                    amplitudes = source.find("amplitudes")
                    amplitudes = util_xml.element_to_array3d(amplitudes)

                    phases = source.find("phases")
                    phases = util_xml.element_to_array3d(phases)

                    # The concept "in model" was originally called "index".
                    in_model = source.find("index")
                    in_model = util_xml.element_to_array3d(in_model)

                    self._svd_outputs = _group_svd_outputs(frequencies,
                                                           damping_factors,
                                                           amplitudes, phases,
                                                           in_model)
                #else:
                    # If there's no SVD data, we leave it up to the caller
                    # to intialize self._svd_outputs.


        elif hasattr(source, "keys"):
            # Quacks like a dict
            for key in list(source.keys()):
                if hasattr(self, key):
                    setattr(self, key, source[key])


######################    Private methods

    def _set_default_svd_inputs(self):
        """Populates self._data_point_count and the other SVD inputs with a 3D
        array of zeros of the appropriate shape. The shape relies on the dims,
        so this method will fail if called when dims are not yet set.
        """
        dims = self.dims

        n_pts  = min(dims[0], self.set.svd_last_n_data_points)
        n_sing =              self.set.svd_last_n_singular_values

        shape = dims[1:]

        self._data_point_count            = n_pts  + np.zeros(shape, int )
        self._signal_singular_value_count = n_sing + np.zeros(shape, int )
        self._water_filter_do_fit         = np.empty(shape, bool)
        self._water_filter_do_fit.fill(True)


    def _set_default_svd_outputs(self):
        """Populates self._svd_outputs with a 3D array of empty SvdOutput
        objects of the appropriate shape. The shape relies on the dims, so
        this method will fail if called when dims are not yet set.
        """
        shape = self.dims[1:]

        # This line multiplies the elements of shape together; e.g. (3, 4, 5) = 60.
        size = reduce(lambda x, y: x * y, shape)
        self._svd_outputs = np.array([svd_output_module.create_zero() for i in range(size)])
        self._svd_outputs.shape = shape


    def _reset_dimensional_data(self, dataset):
        """Resets (to zero) and resizes dimensionally-dependent data"""

        dims = dataset.spectral_dims

        if self.data is None and not self.behave_as_preset:
            # there are no results to maintain
            self.data             = np.zeros(tuple(dims[::-1]), dtype='complex64')
            self._phase_0         = np.zeros(dims[1:])
            self._phase_1         = np.zeros(dims[1:])
            self._frequency_shift = np.zeros(dims[1:])

            self._set_default_svd_inputs()
            self._set_default_svd_outputs()
        else:
            # In the case of a Preset file, there is no data but we still want
            # the phase/b0 values we saved, so we create an empty data array
            # here as a placeholder. That's what this flag is for.

            # data has to match full dimensionality
            if self.dims != dims:
                self.data = np.zeros(tuple(dims[::-1]), dtype='complex64')

            # maintain results if only dims[0] has changed
            if self.dims[1:] != list(self._phase_0.shape):
                self._phase_0         = np.zeros(dims[1:])
                self._phase_1         = np.zeros(dims[1:])
                self._frequency_shift = np.zeros(dims[1:])

            self._set_default_svd_inputs()
            self._set_default_svd_outputs()





#--------------------------------------------------------------------
# test code

def _test():

    import vespa.common.util.time_ as util_time

    test = BlockSpectral([128,1,1,1])

    class_name = test.__class__.__name__
    filename = "_test_output_"+class_name+".xml"
    element = test.deflate()
    root = ElementTree.Element("_test_"+class_name, { "version" : "1.0.0" })
    util_xml.TextSubElement(root, "timestamp", util_time.now().isoformat())
    root.append(element)
    tree = ElementTree.ElementTree(root)
    tree.write(filename, "utf-8")

    tom = 10


if __name__ == '__main__':
    _test()


