# Python modules


# 3rd party modules
import numpy as np
from xml.etree.cElementTree import Element

# Our modules
import vespa.analysis.constants as constants
import vespa.analysis.block_fit_identity as block_fit_identity
import vespa.analysis.chain_fit_giso as chain_fit_giso
import vespa.common.mrs_prior as mrs_prior
import vespa.common.util.xml_ as util_xml

from vespa.common.constants import Deflate
from vespa.analysis.constants import FitLineshapeModel
from vespa.analysis.constants import FitMacromoleculeMethod

# FYI GISO = Grouped Individual Singlets Optimization
#   or maybe Generalized Individual Singlet Optimization





class _Settings(object):
    """ 
    Settings object contains the parameter inputs used for processing in the 
    Chain object in this Block. Having a separate object helps to delineate 
    inputs/outputs and to simplify load/save of preset values.

    This object can also save/recall these values to/from an XML node.
    
    """
    XML_VERSION = "1.0.0"
    
    def __init__(self,  attributes=None):

        #----------------------------------------
        # GISO algorithm input settings

        self.prior_ppm_start                    = 0.0
        self.prior_ppm_end                      = 6.1
        self.prior_mask_source                  = constants.FitPriorMaskSource.ALL_ON
        self.prior_list                         = []
        self.prior_area_scale                   = []
        self.prior_peak_ppm                     = []
        self.prior_search_ppm                   = []
        self.prior_db_ppm                       = []
        self.prior_fix_t2                       = []
        self.prior_search_ph0                   = []
        self.prior_xrange                       = [0,1]     # min/max x-range
        self.prior_yrange                       = [0,1]     # min/max y-range
        self.prior_zrange                       = [0,1]     # min/max z-range
        self.prior_ignore_mask                  = False      
        
        self.lineshape_model                    = FitLineshapeModel.GAUSS 

        self.initial_b0_shift_method            = constants.FitInitialB0ShiftMethod.MANUAL
        self.initial_b0_value                   = 0.0       # initial b0 shift in hz
        self.initial_baseline_method            = constants.FitInitialBaselineMethod.NONE
        self.initial_baseline_lowess_width        = 35.0    # Hz
        self.initial_baseline_lowess_delta        = 3.2     # 
        self.initial_baseline_lowess_ignore_width = 15.0    # Hz
        
        self.initial_cr_cho_separation          = True
        self.initial_peak_search_abs            = False
        self.initial_small_peak_areas           = constants.FitInitialSmallPeakAreas.NAA_RATIO
        self.initial_small_peak_freqs           = constants.FitInitialSmallPeakFreqs.REF_PEAK
        self.initial_linewidth_method           = constants.FitInitialLinewidthMethod.MANUAL 
        self.initial_linewidth_value            = 8.0       # initial linewidth in hz
        self.initial_linewidth_range_start      = 1.9       # NB this is *LW* optimization range
        self.initial_linewidth_range_end        = 4.1       #  start/end in ppm
        self.initial_linewidth_fudge            = 1.0       # LW fudge value (linear scaler)
        self.initial_phase_method               = constants.FitInitialPhaseMethod.MANUAL  
        self.initial_phase0_value               = 0.0       # initial phase0 in deg
        self.initial_phase1_value               = 0.0       # initial phase1 in deg
        self.initial_apply_ko_filter            = False
        self.initial_ko_linewidth_minimum       = 20
        self.initial_ko_points                  = 0
        self.initial_phase1_fid_constant        = 0.0       # for ultra-short FID data phase1 unwrap
        self.initial_lac_method                 = 1
        self.initial_peak_negative_flag         = 0         # ie. for Lac at 72ms, not implemented in GUI yet    

        self.baseline_method                    = constants.FitBaselineMethod.DEFAULT
        self.baseline_smoothing_flag            = False
        self.baseline_skip_last_smooth          = True
        self.baseline_smoothing_width           = 20.0      # in Hz
        self.baseline_smoothing_delta           = 3.2       # lowess smoothing pt skip delta during wavelet baseline calc
        self.baseline_underestimate_method      = constants.FitBaselineUnderestimateMethod.DEFAULT
        self.baseline_underestimate             = 0.0       # in %
        self.baseline_underestimate_last        = 0.0       # in %
        self.baseline_underestimate_steps       = 3         # in %

        self.baseline_spline_nknots             = 10
        self.baseline_spline_spacing            = 10        # in points
        self.baseline_spline_order              = 3
        self.baseline_wavelet_scale             = 4
        self.baseline_wavelet_min_dyad          = 5.0   # float, in Hz min threshold
        self.baseline_wavelet_shrinkage_flag    = 0     # flg - do wavelet shrinkage 0/1
        self.baseline_wavelet_invariance_lag    = 2     # wavelet invariance lag FACTOR (+/- points)
        self.baseline_wavelet_invariance_recombine = 0  # recombination method for invariance lags

        self.macromol_model                           = FitMacromoleculeMethod.DEFAULT
        self.macromol_single_basis_dataset            = None
        self.macromol_single_basis_dataset_id         = '' 
        self.macromol_single_basis_dataset_fname       = '' 
        self.macromol_single_basis_dataset_initval     = FitMacromoleculeMethodInitVal.DEFAULT
        self.macromol_single_basis_dataset_start_area = 1.0
        self.macromol_single_basis_dataset_start_fudge = 1.0        # multiplier, times start_area
        self.macromol_single_basis_dataset_ppm_start   = 1.5        # ppm
        self.macromol_single_basis_dataset_ppm_end     = 1.8        # ppm
        self.macromol_single_basis_dataset_limit_max   = 1.5        # multiplier times start value
        self.macromol_single_basis_dataset_limit_min   = 0.000001   # multiplier times start value, non-negative

        self.optimize_method                    = constants.FitOptimizeMethod.CONSTRAINED_LEVENBERG_MARQUARDT
        self.optimize_scaling_flag              = False
        self.optimize_stop_tolerance            = 0.005         # fit tolerance
        self.optimize_max_iterations            = 100           # max iterations
        self.optimize_global_iterations         = 10             # global iterations for optimization metab/base
        self.optimize_weights_method            = constants.FitOptimizeWeightsMethod.LOCAL_WEIGHTING
        self.optimize_weights_scale_factor      = 1000.0        # multiplier for weighted optimization
        self.optimize_weights_width_factor      = 6.0           # multiplier for width of weighted regions
        self.optimize_weights_water_flag        = True          # Off/On zero water region
        self.optimize_weights_water_start       = 4.1           # water suppression affected region
        self.optimize_weights_water_end         = 5.3           #  start/end in ppm
        self.optimize_weights_lipid_flag        = False         # Off/On zero lipid region
        self.optimize_weights_lipid_start       = -0.5          # region with lipid - lower weight here
        self.optimize_weights_lipid_end         = 1.1           #  start/end in ppm
        self.optimize_weights_small_peak_factor = 1.0           # scale multiple for small peak regions

        self.optimize_bounds_area_max           = 150.0                 
        self.optimize_bounds_area_min           = 50.0                 
        self.optimize_bounds_area_max_small     = 150.0                 
        self.optimize_bounds_area_min_small     = 50.0                 
        self.optimize_enable_bounds_area_small  = False         
        self.optimize_bounds_range_ppm          = 8.0           # Hz    - can be global or only for 'large' singlet peaks
        self.optimize_bounds_range_ppm_small    = 2.0           # Hz    - alternative range for 'small' peaks
        self.optimize_enable_bounds_ppm_small   = False         
        self.optimize_bounds_range_phase0       = 45.0          # deg
        self.optimize_bounds_range_phase1       = 2000.0        # deg
        self.optimize_bounds_max_linewidth      = 1.0
        self.optimize_bounds_min_linewidth      = 0.04
        self.optimize_bounds_mmol_range_area    = 50.0          # % from start values
        self.optimize_bounds_mmol_range_ppm     = 8.0           # Hz
        self.optimize_constrain_ppm_naa_naag    = False
        self.optimize_constrain_ppm_cr_pcr      = False
        self.optimize_constrain_ppm_gpc_pcho    = False
        self.optimize_constrain_ppm_cr2_pcr2    = False
        self.optimize_constrain_ppm_glu_gln     = False
        self.optimize_constrain_ppm_tau_glc     = False
        
        self.optimize_limits_range_area         = 50.0          # %
        self.optimize_limits_range_ppm          = 8.0           # Hz
        self.optimize_limits_range_phase0       = 45.0          # deg
        self.optimize_limits_range_phase1       = 2000.0        # deg
        self.optimize_limits_max_linewidth      = 1.0
        self.optimize_limits_min_linewidth      = 0.04

        self.confidence_intervals_flag          = False
        self.confidence_alpha                   = 0.85
        self.confidence_area_flag               = False
        self.confidence_ppm_flag                = False
        self.confidence_linewidth_flag          = False
        self.confidence_phase_flag              = False

        self.cramer_rao_flag                    = True      # flag to do cramer-rao bound calculations
        self.cramer_rao_ppm_start               = -8.0      # ppm range for noise measurement to calc
        self.cramer_rao_ppm_end                 = -6.0      #  the sigma squared variance.

        # GISO algorithm input objects
        self.prior              = mrs_prior.Prior()

        if attributes is not None:
            self.inflate(attributes)

    ##### Standard Methods and Properties #####################################

    def __str__(self):
        lines = []
        lines.append("------- {0} Object -------".format(self.__class__.__name__))
        lines.append("Not currently written")
        return '\n'.join(lines)


    def deflate(self, flavor=Deflate.ETREE):
        if flavor == Deflate.ETREE:
            e = Element("settings", {"version" : self.XML_VERSION})
                                            
            # These attributes are all lists.
            for attribute in ("prior_list", "prior_area_scale", 
                              "prior_peak_ppm", "prior_search_ppm",
                              "prior_db_ppm", "prior_fix_t2",
                              "prior_search_ph0",
                              "prior_xrange", "prior_yrange",
                              "prior_zrange", ):
                for value in getattr(self, attribute):
                    util_xml.TextSubElement(e, attribute, value)
                    
            # These atttributes are all scalars and map directly to 
            # XML elements of the same name.
            for attribute in ("prior_ppm_start", "prior_ppm_end", 
                    "prior_mask_source", "prior_ignore_mask",
                    "lineshape_model", "initial_b0_shift_method",
                    "initial_b0_value", "initial_baseline_method",
                    "initial_baseline_lowess_width",
                    "initial_baseline_lowess_delta",
                    "initial_baseline_lowess_ignore_width",
                    "initial_cr_cho_separation",
                    "initial_peak_search_abs", "initial_small_peak_areas",
                    "initial_small_peak_freqs", "initial_linewidth_method",
                    "initial_linewidth_value",
                    "initial_linewidth_range_start",
                    "initial_linewidth_range_end", "initial_linewidth_fudge",
                    "initial_phase_method",
                    "initial_phase0_value", "initial_phase1_value",
                    "initial_apply_ko_filter",
                    "initial_ko_linewidth_minimum", "initial_ko_points",
                    "initial_phase1_fid_constant",
                    "initial_lac_method",
                    "initial_peak_negative_flag",
                    "baseline_method",
                    "baseline_smoothing_flag", "baseline_skip_last_smooth",
                    "baseline_smoothing_width", "baseline_smoothing_delta",
                    "baseline_underestimate", "baseline_spline_nknots",
                    "baseline_spline_spacing", "baseline_spline_order",
                    "baseline_wavelet_scale", "baseline_wavelet_min_dyad",
                    "baseline_wavelet_shrinkage_flag",
                    "baseline_wavelet_invariance_lag",
                    "baseline_wavelet_invariance_recombine",
                    "macromol_model",
                    "macromol_single_basis_dataset_start_area",
                    "optimize_method", "optimize_scaling_flag",
                    "optimize_stop_tolerance", "optimize_max_iterations",
                    "optimize_limits_range_area", "optimize_limits_range_ppm",
                    "optimize_limits_range_phase0",
                    "optimize_limits_range_phase1",
                    "optimize_limits_max_linewidth",
                    "optimize_limits_min_linewidth",
                    "optimize_global_iterations", "optimize_weights_method",
                    "optimize_weights_scale_factor",
                    "optimize_weights_width_factor",
                    "optimize_weights_water_flag",
                    "optimize_weights_water_start",
                    "optimize_weights_water_end",
                    "optimize_weights_lipid_flag",
                    "optimize_weights_lipid_start",
                    "optimize_weights_lipid_end",
                    "optimize_weights_small_peak_factor",
                    "confidence_intervals_flag", "confidence_alpha",
                    "confidence_area_flag", "confidence_ppm_flag",
                    "confidence_linewidth_flag", "confidence_phase_flag",
                    "cramer_rao_flag", "cramer_rao_ppm_start",
                    "cramer_rao_ppm_end",
                    ):
                util_xml.TextSubElement(e, attribute, getattr(self, attribute))

            if self.macromol_single_basis_dataset:
                # In the next line, we *have* to save the uuid values from the 
                # actual object rather than from the attribute above, in  
                # order for the associated dataset uuid to reflect the new id
                # that is given in the top level dataset. Associated datasets are
                # given new temporary uuid values so that if the main dataset is 
                # saved and immediately loaded back in, we do not get collisions
                # between the newly opened datasets and already existing ones.
                
                util_xml.TextSubElement(e, "macromol_single_basis_dataset_id", self.macromol_single_basis_dataset.id)


            e.append(self.prior.deflate())
            
            return e
            
        elif flavor == Deflate.DICTIONARY:
            return self.__dict__.copy()


    def inflate(self, source):
        if hasattr(source, "makeelement"):
            # Quacks like an ElementTree.Element
            
            # We inflate in attributes grouped by type since there's so
            # doggone many attrs on this class.
            
            # Booleans
            for attribute in ("prior_ignore_mask", 
                              "initial_cr_cho_separation",
                              "initial_peak_search_abs",
                              "initial_apply_ko_filter", 
                              "baseline_smoothing_flag",
                              "baseline_skip_last_smooth",
                              "optimize_scaling_flag",                              
                              "optimize_weights_lipid_flag",
                              "optimize_weights_water_flag",
                              "confidence_intervals_flag",
                              "confidence_area_flag", 
                              "confidence_ppm_flag",
                              "confidence_linewidth_flag",
                              "confidence_phase_flag", 
                              "cramer_rao_flag",
                             ):
                item = source.findtext(attribute)
                if item is not None:
                    setattr(self, attribute, util_xml.BOOLEANS[item])
                        
                        
            # floats
            for attribute in ("prior_ppm_start", 
                              "prior_ppm_end", 
                              "initial_b0_value", 
                              "initial_baseline_lowess_width",
                              "initial_baseline_lowess_delta",
                              "initial_baseline_lowess_ignore_width",
                              "initial_linewidth_value",
                              "initial_linewidth_range_start",
                              "initial_linewidth_range_end", 
                              "initial_linewidth_fudge",
                              "initial_phase0_value", 
                              "initial_phase1_value", 
                              "initial_ko_linewidth_minimum", 
                              "initial_phase1_fid_constant",
                              "initial_lac_method",
                              "baseline_smoothing_width", 
                              "baseline_smoothing_delta",
                              "baseline_underestimate", 
                              "baseline_wavelet_min_dyad",
                              "macromol_single_basis_dataset_start_area",
                              "optimize_stop_tolerance", 
                              "optimize_limits_range_area",
                              "optimize_limits_range_ppm",
                              "optimize_limits_range_phase0",
                              "optimize_limits_range_phase1",
                              "optimize_limits_max_linewidth",
                              "optimize_limits_min_linewidth",
                              "optimize_weights_scale_factor",
                              "optimize_weights_width_factor",
                              "optimize_weights_water_start",
                              "optimize_weights_water_end",
                              "optimize_weights_lipid_start",
                              "optimize_weights_lipid_end",
                              "optimize_weights_small_peak_factor", 
                              "confidence_alpha",
                              "cramer_rao_ppm_start", 
                              "cramer_rao_ppm_end",
                             ):
                item = source.findtext(attribute)
                if item is not None:
                    setattr(self, attribute, float(item))

            # ints
            for attribute in ("initial_ko_points",
                              "initial_peak_negative_flag",
                              "baseline_spline_nknots", 
                              "baseline_spline_spacing",
                              "baseline_spline_order",
                              "baseline_wavelet_scale",
                              "baseline_wavelet_shrinkage_flag",
                              "baseline_wavelet_invariance_lag",
                              "baseline_wavelet_invariance_recombine",
                              "optimize_max_iterations",
                              "optimize_global_iterations",
                             ):
                item = source.findtext(attribute)
                if item is not None:
                    setattr(self, attribute, int(float(item)))

            # No translation required for these attrs
            for attribute in ("prior_mask_source", 
                              "lineshape_model", 
                              "initial_b0_shift_method",
                              "initial_baseline_method",
                              "initial_small_peak_areas",
                              "initial_small_peak_freqs",
                              "initial_linewidth_method",
                              "initial_phase_method", 
                              "baseline_method",
                              "macromol_model",
                              "macromol_single_basis_dataset_id",
                              "optimize_method", 
                              "optimize_weights_method",
                             ):
                item = source.findtext(attribute)
                if item is not None:
                    setattr(self, attribute, item)

            # lists
            self.prior_list       = [val.text for val in source.findall("prior_list")]
            self.prior_area_scale = [float(val.text) for val in source.findall("prior_area_scale")]
            self.prior_peak_ppm   = [float(val.text) for val in source.findall("prior_peak_ppm")]
            self.prior_search_ppm = [float(val.text) for val in source.findall("prior_search_ppm")]
            self.prior_db_ppm     = [util_xml.BOOLEANS[val.text] for val in source.findall("prior_db_ppm")]
            self.prior_fix_t2     = [float(val.text) for val in source.findall("prior_fix_t2")]
            self.prior_search_ph0 = [float(val.text) for val in source.findall("prior_search_ph0")]
            
            nmet = len(self.prior_list)
            if not self.prior_area_scale:
                self.prior_area_scale = [1.0 for i in range(nmet)]
            if not self.prior_search_ppm:
                self.prior_search_ppm = [0.1 for i in range(nmet)]
            if not self.prior_db_ppm:
                self.prior_db_ppm = [False for i in range(nmet)]
            if not self.prior_fix_t2:
                self.prior_fix_t2 = [1000.0 for i in range(nmet)]
            if not self.prior_search_ph0:
                self.prior_search_ph0 = [0.0 for i in range(nmet)]

            self.prior_xrange     = [int(val.text) for val in source.findall("prior_xrange")]
            self.prior_yrange     = [int(val.text) for val in source.findall("prior_yrange")]
            self.prior_zrange     = [int(val.text) for val in source.findall("prior_zrange")]
            
            # subobjects
            self.prior = mrs_prior.Prior(source.find("prior"))
            

        elif hasattr(source, "keys"):
            # Quacks like a dict
            for key in list(source.keys()):
                if hasattr(self, key):
                    setattr(self, key, source[key])




class BlockFitGiso(block_fit_identity.BlockFitIdentity):
    """ 
    Building block to hold the state of a step in an MRS processing chain.
    Includes the functionality to save/recall this object to/from an XML node.

    Contains inputs/results for converting the fitted metabolite areas to 
    concentration values.
    
    """
    XML_VERSION = "1.0.0"
    
    def __init__(self, attributes=None):
        """

        General Parameters
        -----------------------------------------

        id          A permanent, unique identifying string for this 
                    object. Typically serves as a "source_id" for 
                    some other object. It is part of the provenance
                    for this processing functor chain
                   
        source_id   The unique identifier used to find the input data
                    for this object. It may refer to one whole object
                    that has only one result, OR it could refer to s
                    single results inside an object that has multiple
                    results.

        """        
        super().__init__(attributes)

        #----------------------------------------
        # GISO algorithm input settings
        self.set = _Settings()
        
        #----------------------------------------
        # GISO result objects
        self.fit_results        = None      # fit parameters from optimization
        self.fit_stats          = None      # stats from the optimization for each voxel
        self.fit_baseline       = None      # baseline result (may be vector or parameterized)
        self.confidence         = None      # confidence limits results storage
        self.cramer_rao         = None      # cramer-rao bounds results storage
        self.initial_values     = None      # initial values for metab optimization   
        self.prior_mask         = None

        if attributes is not None:
            self.inflate(attributes)

        self.chain = None


    ##### Standard Methods and Properties #####################################

    @property
    def data(self):
        self.fit_results


    @property
    def dims(self):
        """Dimensions for fitting results in a list, read only."""
        return list(self.fit_results.shape) if self.fit_results is not None else None


    @property
    def nparam(self):
        nparam = len(self.set.prior_list)*3+2 + self.nmmol
        return nparam  
        
    @property
    def nmmol(self):
        nmmol = 0
        if self.set.macromol_model == FitMacromoleculeMethod.SINGLE_BASIS_DATASET:
            nmmol += 2
        return nmmol
        

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
        self.chain = chain_fit_giso.ChainFitGiso(dataset, self)


    def set_dims(self, dataset):
        """
        Given a Dataset object, this is an opportunity for this block object 
        to ensure that its results dims match those of the parent dataset.
        Because the fitting results are parameters, they do not have to match
        in the spectral dimensions, just in the spatial dimensions. But, we 
        have to calculate the parameter dimension size from current settings
        to compare to the exisiting array dimensions.
        """
        if self.dims is None or self.fit_baseline is None:
            # need this here for initialization
            self._reset_dimensional_data(dataset)
            if self.chain is not None:
                self.chain.reset_results_arrays()
        else:
    
            # Prior depends on spectral data settings
            self.set.prior.calculate_full_basis_set(None, None, dataset)
#             self.set.prior.calculate_full_basis_set(self.set.prior_ppm_start,
#                                                     self.set.prior_ppm_end, 
#                                                     dataset)
    
            # if the spectral, spatial or fit parameter dimensions have changed,
            # then set a flag to reset results arrays
            spectral_dims  = dataset.spectral_dims
            result_dims    = list(spectral_dims)
            result_dims[0] = self.nparam
            spectral_flag  = self.fit_baseline.shape[0] != spectral_dims[0]
            result_flag    = self.dims != result_dims
    
            if spectral_flag or result_flag:
                self._reset_dimensional_data(dataset)
                if self.chain is not None:
                    self.chain.reset_results_arrays()


    def check_parameter_dimensions(self, dataset):
        """
        Checks the "nparam" dimension in the results to see if the number of 
        parameters in the model has changed. Only resets results if this
        dimension has changed.
        
        """
        if self.fit_results.shape[0] != self.nparam:
            self._reset_dimensional_data(dataset)


    def _reset_dimensional_data(self, dataset):
        """
        Resets (to zero) and resizes dimensionally-dependent data
        
        fit_results     - fit parameters from optimization
        fit_stats       - optimization chisqr, wtchisqr and finitemath values
        fit_baseline    - baseline result (may be vector or parameterized)
        confidence      - confidence limits results storage
        cramer_rao      - cramer-rao bounds results storage
        initial_values  - initial values for metab optimization
        prior_mask      - spatial (x,y,z) mask to indicate what voxels to fit
        
        Note. B0 shift is accounted for in Spectral object/tabs
        
        """
        dims = dataset.spectral_dims
        
        nparam = self.nparam        
        
        if self.fit_results is None:
        
            self.fit_results    = np.zeros((nparam, dims[1], dims[2], dims[3]))      
            self.fit_stats      = np.zeros((3, dims[1], dims[2], dims[3]))      
            self.fit_baseline   = np.zeros(tuple(dataset.spectral_dims), dtype='complex64')
            self.confidence     = np.zeros((nparam, dims[1], dims[2], dims[3]))      
            self.cramer_rao     = np.zeros((nparam, dims[1], dims[2], dims[3]))      
            self.initial_values = np.zeros((nparam, dims[1], dims[2], dims[3]))      
            self.prior_mask     = np.ones((dims[1], dims[2], dims[3])) 
        
        else:
            param_dims = list(dims)
            param_dims[0] = nparam

            # maintain results if no dimension has changed
            if self.fit_baseline.shape[::-1] != dims or \
                self.fit_results.shape[::-1] != param_dims:
                
                self.fit_baseline   = np.zeros(tuple(dims), dtype='complex64')        
        
                self.fit_results    = np.zeros((nparam, dims[1], dims[2], dims[3]))      
                self.fit_stats      = np.zeros((3, dims[1], dims[2], dims[3]))      
                self.fit_baseline   = np.zeros(tuple(dataset.spectral_dims), dtype='complex64')
                self.confidence     = np.zeros((nparam, dims[1], dims[2], dims[3]))      
                self.cramer_rao     = np.zeros((nparam, dims[1], dims[2], dims[3]))      
                self.initial_values = np.zeros((nparam, dims[1], dims[2], dims[3]))      
            
                self.prior_mask     = np.ones((dims[1], dims[2], dims[3])) 


    def get_fitted_water_area(self):
        
        flag_water = False
        for indx, item in enumerate(self.set.prior_list):
            if 'water' in item.lower() or 'h2o' in item.lower():
                flag_water = True
                indx_water = indx
                break
        
        if not flag_water:
            return None
        else:
            dims = self.fit_results.shape
            res  = self.fit_results[indx_water,:,:,:]
            res  = res.copy() 
            if dims[1] * dims[2] * dims[3] == 1:
                res.shape = 1,1,1,1
            return res  # return water areas


    def get_associated_datasets(self, is_main_dataset=True):
        """
        Returns a list of datasets associated with this object

        'is_main_dataset' signals that this is the top level dataset gathering 
        associated datasets, and is used to stop circular references

        """
        datasets = block_fit_identity.BlockFitIdentity.get_associated_datasets(self, is_main_dataset)

        if is_main_dataset:
            # at this point, only a macromolecule model may need an associated dataset
            if self.set.macromol_single_basis_dataset:
                datasets += self.set.macromol_single_basis_dataset.get_associated_datasets(is_main_dataset=False)
                datasets += [self.set.macromol_single_basis_dataset]

        return datasets


    def set_associated_datasets(self, datasets): 
        """
        When we open a VIFF format file, main._import_file() calls this method
        to parse/store any datasets associated with this one as described below.
        
        """
        for dataset in datasets:
            if dataset.id == self.set.macromol_single_basis_dataset_id:
                self.set.macromol_single_basis_dataset = dataset


    def deflate(self, flavor=Deflate.ETREE):
        if flavor == Deflate.ETREE:
            e = Element("block_fit_giso", { "id" : self.id,
                                            "version" : self.XML_VERSION})
            
            util_xml.TextSubElement(e, "behave_as_preset", self.behave_as_preset)
            
            e.append(self.set.deflate())

            if not self.behave_as_preset:
            
                for dim in self.dims:
                    util_xml.TextSubElement(e, "dim", dim)

                if self.fit_results is not None:
                    e.append(util_xml.numpy_array_to_element(self.fit_results,'fit_results'))
                if self.fit_stats is not None:
                    e.append(util_xml.numpy_array_to_element(self.fit_stats,'fit_stats'))
                if self.fit_baseline is not None:
                    e.append(util_xml.numpy_array_to_element(self.fit_baseline,'fit_baseline'))
                if self.confidence is not None:
                    e.append(util_xml.numpy_array_to_element(self.confidence,'confidence'))
                if self.cramer_rao is not None:
                    e.append(util_xml.numpy_array_to_element(self.cramer_rao,'cramer_rao'))
                if self.initial_values is not None:
                    e.append(util_xml.numpy_array_to_element(self.initial_values,'initial_values'))
                if self.prior_mask is not None:
                    e.append(util_xml.numpy_array_to_element(self.prior_mask,'prior_mask'))

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
            self.set = util_xml.find_settings(source, "block_giso_settings")
            self.set = _Settings(self.set)

            if not self.behave_as_preset:

                # Explicit tests for None necessary in the code below. 
                temp = source.find("fit_results")
                if temp is not None:
                    self.fit_results = util_xml.element_to_numpy_array(temp)

                temp = source.find("fit_stats")
                if temp is not None:
                    self.fit_stats = util_xml.element_to_numpy_array(temp)

                temp = source.find("fit_baseline")
                if temp is not None:
                    self.fit_baseline = util_xml.element_to_numpy_array(temp)

                temp = source.find("confidence")
                if temp is not None:
                    self.confidence = util_xml.element_to_numpy_array(temp)

                temp = source.find("cramer_rao")
                if temp is not None:
                    self.cramer_rao = util_xml.element_to_numpy_array(temp)

                temp = source.find("initial_values")
                if temp is not None:
                    self.initial_values = util_xml.element_to_numpy_array(temp)

                temp = source.find("prior_mask")
                if temp is not None:
                    self.prior_mask = util_xml.element_to_numpy_array(temp)


        elif hasattr(source, "keys"):
            # Quacks like a dict
            for key in list(source.keys()):
                if hasattr(self, key):
                    setattr(self, key, source[key])


    def results_as_html(self, voxel, lw=0.0, lwmin=0.0, lwmax=0.0, 
                        data_source="", image=None):
        """
        Given a voxel, linewidth params, and a data source (often a filename), 
        returns HTML-formatted results for that voxel. The HTML is appropriate 
        for the wx.Html control (which understand limited HTML) as well as for
        writing to a file.

        If the image param is populated, it should be a tuple of 
        (mime_type, image_data). The former should be a string like "image/png".
        The latter should be base64-encoded image data.
        """
        # First we assemble the data we need.
        nmet = len(self.set.prior_list)
        
        names = self.set.prior_list
        res   = self.fit_results[:,voxel[0],voxel[1],voxel[2]]
        crao  = self.cramer_rao[:,voxel[0],voxel[1],voxel[2]]
        conf  = self.confidence[:,voxel[0],voxel[1],voxel[2]]
        stats = self.fit_stats[:,voxel[0],voxel[1],voxel[2]]
        # both cramer-rao and confidence intervals may be off/on
        if len(crao) != len(res):
            crao = res * 0
        if len(conf) != len(res):
            conf = res * 0

        table1 = [['Area Results', 'Area', ' CrRao[%]', ' CnfInt[%]']]
        for i, item in enumerate(names):
            table1.append([item, res[i], crao[i], conf[i]])

        if self.set.macromol_model == FitMacromoleculeMethod.SINGLE_BASIS_DATASET:
            table1.append(['MMol', res[nmet*3+2], 0.0, 0.0])
            
        table1 = _pretty_space_table(table1, places=4)

        table2 = [['PPM Results', 'PPM', ' CrRao[ppm]', ' CnfInt[ppm]']]
        for i,item in enumerate(names):
            table2.append([item, res[i+nmet], crao[i+nmet], conf[i+nmet]])

        if self.set.macromol_model == FitMacromoleculeMethod.SINGLE_BASIS_DATASET:
            table2.append(['MMol', res[nmet*3+3], 0.0, 0.0])

        table2 = _pretty_space_table(table2, places=4)

        table2a = [['Tb Results', 'Tb[sec]', 'LW[Hz]', 'CrRao[sec]', ' CnfInt[sec]']]
        for i,item in enumerate(names):
            table2a.append([item, res[i+nmet*2], lw[i], crao[i+nmet*2], conf[i+nmet*2]])

        if self.set.macromol_model == FitMacromoleculeMethod.SINGLE_BASIS_DATASET:
            table2a.append(['MMol', res[nmet*3+4], 0.0, 0.0])

        table2a = _pretty_space_table(table2a, places=4)

        table3 =     [['Global Results', 'Value', ' CrRao[delta]', ' CnfInt[%]']]
        table3.append(['Phase0',     res[nmet*3+0], crao[nmet*3+0], conf[nmet*3+0] ])
        table3.append(['Phase1',     res[nmet*3+1], crao[nmet*3+1], conf[nmet*3+1] ])
        table3 = _pretty_space_table(table3, places=5)
        
        table4 = [['Calculation Results', ' Value', '  Max LW', '  Min LW']]
        table4.append(['ChiSquare', stats[0], ' ', ' '])
        table4.append(['Weighted ChiSquare', stats[1], ' ', ' '])
        matherr = str(stats[2] != 0)
        table4.append(['Math Finite Error', matherr, ' ', ' '])
        table4 = _pretty_space_table(table4, places=5)

        # Now that the data is assembled, we HTML-ify it.
        html = ElementTree.Element("html")
        head = ElementTree.SubElement(html, "head")
        style = util_xml.TextSubElement(head, "style", _CSS)
        style.set("type", "text/css")

        body = ElementTree.SubElement(html, "body")
        
        util_xml.TextSubElement(body, "h2", "Analysis Fit-Giso Results")

        e_div = ElementTree.SubElement(body, "div")

        if data_source:
            e_tt = util_xml.TextSubElement(e_div, "tt", "Data Source: ")
            util_xml.TextSubElement(e_tt, "small", data_source)
            ElementTree.SubElement(e_div, "br")
            
        voxel = [x + 1 for x in voxel]
        util_xml.TextSubElement(e_div, "tt", 'Voxel: (%d,%d,%d)' % tuple(voxel))

        if image:
            # If there's image data, we assume that this will be written to 
            # a file for display in a proper browser, so we can use slightly
            # fancier HTML.
            mime_type, image_data = image

            e_div = ElementTree.SubElement(body, "div", 
                                           { "id" : "image",
                                             "style" : "float: right; width: 50%",
                                           }
                                          )

            e_div.append(ElementTree.Comment(_IE_INCAPABLE_MSG))

            # In order to keep the HTML + image in one file, we use the
            # little-known "data" scheme.
            # ref: http://en.wikipedia.org/wiki/Data_URI_scheme
            src = "data:%s;base64,%s" % (mime_type, image_data)

            ElementTree.SubElement(e_div, "img", 
                                           {"style" : "width: 90%",
                                            "src" : src
                                            })

        e_div = ElementTree.SubElement(body, "div", {"id" : "table"})

        tables = (table1, table2, table2a, table3, table4)

        for table in tables:
            title = table[0]
            e_pre = ElementTree.SubElement(e_div, "pre")
            e_u   = ElementTree.SubElement(e_pre, "u")
            util_xml.TextSubElement(e_u, "b", title)
            e     = util_xml.TextSubElement(e_div, "pre", '\n'.join(table[1:]))
            
        # Keep in mind that HTML is whitespace sensitive, and if you call 
        # util_xml.indent() on the HTML, it will change the formatting.

        return ElementTree.tostring(html)


    def results_as_csv(self, voxel, lw=0.0, lwmin=0.0, lwmax=0.0, source="", dsetname=""):
        """
        Given a voxel, linewidth params, and a data source (often a filename), 
        returns CSV-formatted (comma separated variables)string containing both
        the voxel fitting results and header string descriptions for each 
        column.

        """
        hdr = []
        val = []

        hdr.append('Filename')   
        source = source.replace(",","_")     # some users have commas in filenames 
        val.append(source)

        hdr.append('Dataset Name')        
        dsetname = dsetname.replace(",","_") # some users have commas in filenames    
        val.append(dsetname)
        
        hdr.append('Voxel')
        val.append(str(voxel[0])+' '+str(voxel[1])+' '+str(voxel[2]))
        
        nmet = len(self.set.prior_list)
        
        names = self.set.prior_list
        res   = self.fit_results[:,voxel[0],voxel[1],voxel[2]]
        crao  = self.cramer_rao[ :,voxel[0],voxel[1],voxel[2]]
        conf  = self.confidence[ :,voxel[0],voxel[1],voxel[2]]
        stats = self.fit_stats[  :,voxel[0],voxel[1],voxel[2]]
        # both cramer-rao and confidence intervals may be off/on
        if len(crao) != len(res):
            crao = res * 0
        if len(conf) != len(res):
            conf = res * 0

        for i, item in enumerate(names):
            hdr.append('Area '+item)        
            hdr.append('CrRao[%]')        
            hdr.append('CnfInt[%]')        
            val.append(str(res[i]))
            val.append(str(crao[i]))
            val.append(str(conf[i]))

        for i,item in enumerate(names):
            hdr.append('PPM '+item)        
            hdr.append('CrRao[%]')        
            hdr.append('CnfInt[%]')        
            val.append(str(res[i+nmet]))
            val.append(str(crao[i+nmet]))
            val.append(str(conf[i+nmet]))

        for i,item in enumerate(names):
            hdr.append('Tb[sec]'+item)        
            hdr.append('LW[Hz]'+item)        
            hdr.append('CrRao[%]')        
            hdr.append('CnfInt[%]')        
            val.append(str(res[i+nmet*2]))
            val.append(str(lw[i]))
            val.append(str(crao[i+nmet*2]))
            val.append(str(conf[i+nmet*2]))

        hdr.append('Phase0 ')        
        hdr.append('CrRao[%]')        
        hdr.append('CnfInt[%]')        
        val.append(str(res[nmet*3+0]))
        val.append(str(crao[nmet*3+0]))
        val.append(str(conf[nmet*3+0]))

        hdr.append('Phase1 ')        
        hdr.append('CrRao[%]')        
        hdr.append('CnfInt[%]')        
        val.append(str(res[nmet*3+1]))
        val.append(str(crao[nmet*3+1]))
        val.append(str(conf[nmet*3+1]))

        hdr.append('ChiSquare ')        
        val.append(str(stats[0]))

        hdr.append('WtChiSquare ')        
        val.append(str(stats[1]))

        hdr.append('Math Finite Error ')        
        matherr = str(stats[2] != 0)
        val.append(str(matherr))

        return val, hdr
    
    
    def results_as_csv_all_voxels_truncated(self, voxels, datafile=' ', measure_times=[]):
        """
        Given a voxel, linewidth params, and a data source (often a filename), 
        returns CSV-formatted (comma separated variables)string containing both
        the voxel fitting results and header string descriptions for each 
        column.

        """

        source = datafile.replace(",","_")     # some users have commas in filenames 

        nmet = len(self.set.prior_list)
        names = self.set.prior_list
        
        hdr = []
        val = []

        hdr.append('Datafile')   
        hdr.append('Voxel')
        hdr.append('Time [sec]')

        for i, item in enumerate(names):
            hdr.append('Area '+item)        
            hdr.append('PPM '+item)        
            hdr.append('Tb[sec]'+item)        
        hdr.append('Phase0 ')        
        hdr.append('Phase1 ')        


        for i,voxel in enumerate(voxels):

            res   = self.fit_results[:,voxel[0],voxel[1],voxel[2]]
            itime = measure_times[i]
            
            val.append(datafile)
            val.append(str(voxel[0])+' '+str(voxel[1])+' '+str(voxel[2]))
            val.append(itime)
            
            for i, item in enumerate(names):
                val.append(str(res[i+nmet*0]))
                val.append(str(res[i+nmet*1]))
                val.append(str(res[i+nmet*2]))

            val.append(str(res[nmet*3+0]))
            val.append(str(res[nmet*3+1]))
    

        return val, hdr    
    
    
    def results_as_csv_all_voxels_areas(self, voxels, datafile=' ', measure_times=[]):
        """
        Given a voxel, linewidth params, and a data source (often a filename), 
        returns CSV-formatted (comma separated variables)string containing both
        the voxel fitting results and header string descriptions for each 
        column.

        """

        source = datafile.replace(",","_")     # some users have commas in filenames 

        nmet = len(self.set.prior_list)
        names = self.set.prior_list
        
        lines = []

        times = []
        times.append(source)
        times.append('_Time [min] ')
        for measure_time in measure_times:
            times.append(str((measure_time - measure_times[0])/60.0))       # in minutes
        lines.append(times)

        
        for j, name in enumerate(names):
            val = []
            val.append(source)
            val.append(name+' Area')
            for i,voxel in enumerate(voxels):
                val.append( str(self.fit_results[j+nmet*0,voxel[0],voxel[1],voxel[2]]) )
            lines.append(val)

        return lines    



    def results_in_table(self, voxel, lw=0.0, lwmin=0.0, lwmax=0.0, nozeros=False, noppm=False, fixphase=False):
        """
        Given a voxel, linewidth params, and a data source (often a filename), 
        returns formatted results for that voxel. 
        
        """
        # First we assemble the data we need.
        nmet  = len(self.set.prior_list)
        
        names = self.set.prior_list
        res   = self.fit_results[:,voxel[0],voxel[1],voxel[2]]
        crao  = self.cramer_rao[:,voxel[0],voxel[1],voxel[2]]
        conf  = self.confidence[:,voxel[0],voxel[1],voxel[2]]
        stats = self.fit_stats[:,voxel[0],voxel[1],voxel[2]]
        
        do_crao = False if max(crao) == 0.0 and nozeros else True
        do_conf = False if max(conf) == 0.0 and nozeros else True
        
        # both cramer-rao and confidence intervals may be off/on
        if len(crao) != len(res):
            crao = res * 0
        if len(conf) != len(res):
            conf = res * 0

        hdr1  = ['Metab Results', 'Area']
        hdr2  = ['PPM Results', 'PPM']
        hdr2a = ['Tb Results', 'Tb [sec]', 'LW [Hz]']
        hdr3  = ['Global Results', 'Value']
        hdr4  = ['Calc Results', ' Value']
        if do_crao: 
            hdr1.append('CrRao[%]')
            hdr2.append('CrRao[delta]')
            hdr2a.append('CrRao[delta]')
            hdr3.append('CrRao[delta]')
            hdr4.append(' ')
        if do_conf: 
            hdr1.append('CnfInt[%]')
            hdr2.append('CnfInt[delta]')
            hdr2a.append('CnfInt[delta]')
            hdr3.append('CnfInt[delta]')
            hdr4.append(' ')

        nsect = [0,]
        table1 = [hdr1,]
        
        for i, item in enumerate(names):
            tmp = [item, res[i]]
            if do_crao: tmp.append(crao[i])
            if do_conf: tmp.append(conf[i])
            table1.append(tmp)

        if self.set.macromol_model == FitMacromoleculeMethod.SINGLE_BASIS_DATASET:
            tmp = ['MMol', res[nmet*3+2]]
            if do_crao: tmp.append(0.0)
            if do_conf: tmp.append(0.0)
            table1.append(tmp)
            
        nsect.append(len(table1))

        if not noppm:
            table1.append(hdr2)
            for i,item in enumerate(names):
                tmp = [item, res[i+nmet]]
                if do_crao: tmp.append(crao[i+nmet])
                if do_conf: tmp.append(conf[i+nmet])
                table1.append(tmp)
     
            if self.set.macromol_model == FitMacromoleculeMethod.SINGLE_BASIS_DATASET:
                tmp = ['MMol', res[nmet*3+3]]
                if do_crao: tmp.append(0.0)
                if do_conf: tmp.append(0.0)
                table1.append(tmp)

            table1.append(hdr2a)
            for i,item in enumerate(names):
                tmp = [item, res[i+nmet*2], lw[i]]
                if do_crao: tmp.append(crao[i+nmet*2])
                if do_conf: tmp.append(conf[i+nmet*2])
                table1.append(tmp)
     
            if self.set.macromol_model == FitMacromoleculeMethod.SINGLE_BASIS_DATASET:
                tmp = ['MMol', res[nmet*3+4]]
                if do_crao: tmp.append(0.0)
                if do_conf: tmp.append(0.0)
                table1.append(tmp)

        nsect.append(len(table1))
        table1.append(hdr3)
        
        if not fixphase:
            tmp = ['Phase0', res[nmet*3+0]]
            if do_crao: tmp.append(crao[nmet*3+0])
            if do_conf: tmp.append(conf[nmet*3+0])
            table1.append(tmp)
    
            tmp = ['Phase1', res[nmet*3+1]]
            if do_crao: tmp.append(crao[nmet*3+1])
            if do_conf: tmp.append(conf[nmet*3+1])
            table1.append(tmp)
        else:
            tmp = ['Phase0', '0.0 ('+str(np.round(res[nmet*3+0],1))+')']
            if do_crao: tmp.append(crao[nmet*3+0])
            if do_conf: tmp.append(conf[nmet*3+0])
            table1.append(tmp)
    
            tmp = ['Phase1', '0.0 ('+str(np.round(res[nmet*3+1],1))+')']
            if do_crao: tmp.append(crao[nmet*3+1])
            if do_conf: tmp.append(conf[nmet*3+1])
            table1.append(tmp)
            
        nsect.append(len(table1))
        table1.append(hdr4)

        tmp = ['ChiSquare', stats[0]]
        if do_crao: tmp.append(' ')
        if do_conf: tmp.append(' ')
        table1.append(tmp)
        
        tmp = ['Weighted ChiSquare', stats[1]]
        if do_crao: tmp.append(' ')
        if do_conf: tmp.append(' ')
        table1.append(tmp)

        matherr = str(stats[2] != 0)
        tmp = ['Math Finite Error', matherr]
        if do_crao: tmp.append(' ')
        if do_conf: tmp.append(' ')
        table1.append(tmp)
        
        table1 = _pretty_space_table_inplace(table1, places=5)

        return table1, nsect
    
 
