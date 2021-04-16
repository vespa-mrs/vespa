# Python modules

# 3rd party modules
import numpy as np

# Our modules
import vespa.analysis.functors.funct_fidsum_all as funct_fidsum_all
from vespa.analysis.chain_base import Chain



class ChainPrepFidsum(Chain):
    """ 
    Building block object used to create a processing chain for MRS data.
    
    Combines multi-coil/average FID data into a single averaged FID data set.
    
    """
    def __init__(self, dataset, block):
        """
        Chain objects organize Algo (algorithm) calls by setting up access to
        input data and parameters, and creating standard output values for View.

        Base class sets convenience references to:  self._block and self._dataset

        self.data is always initialized as []

        """
        super().__init__(dataset, block)

        # results specific to this Chain object
        self.coil_combine_weights = None
        self.coil_combine_phases = None
        self.exclude_indices = []
        self.frequency_shift = None
        self.phase_0         = None

        self.time_summed   = None
        self.time_adjusted = None
        self.freq_summed   = None
        self.freq_current  = None
        self.freq_adjusted = None
        self.freq_raw      = None

        self.reset_results_arrays()
        
 
    def reset_results_arrays(self):
        """
        A separate method so it can be called outside __init__. Should
        create/set enough results to keep View happy if run() fails.

        """
        _, ncoils, navgs, npts = self._dataset.raw_shape
        sum_shape = [1,1,1,npts]
        adj_shape = [1,1,navgs, npts]

        # check if results arrays exist
        if self.time_summed   is None: self.time_summed   = np.zeros(sum_shape, complex)
        if self.freq_summed   is None: self.freq_summed   = np.zeros(sum_shape, complex)
        if self.freq_current  is None: self.freq_current  = np.zeros(sum_shape, complex)
        if self.time_adjusted is None: self.time_adjusted = np.zeros(adj_shape, complex)
        if self.freq_adjusted is None: self.freq_adjusted = np.zeros(adj_shape, complex)

        # check if results arrays are right sized
        if self.time_summed.shape   != sum_shape: self.time_summed   = np.zeros(sum_shape, complex)
        if self.freq_summed.shape   != sum_shape: self.freq_summed   = np.zeros(sum_shape, complex)
        if self.freq_current.shape  != sum_shape: self.freq_current  = np.zeros(sum_shape, complex)
        if self.time_adjusted.shape != adj_shape: self.time_adjusted = np.zeros(adj_shape, complex)
        if self.freq_adjusted.shape != adj_shape: self.freq_adjusted = np.zeros(adj_shape, complex)

        if self.freq_raw is not None:
            if self.freq_raw.shape != adj_shape:
                self.freq_raw = None        # set to None signals a recalc in chain

 
    def run(self, voxels, entry='all', freq_raw=False):
        """
        Run is typically called every time a processing setting is changed
        in the parent (block) object. Run processes a single voxel at a time.

        This object maintains previous run() results values until next run().
        This allows the View to update without having to re-run the pipeline.

        The 'entry' keyword adds flexibility to Block-Chain-View relationship.

        NB. we are processing all raw data coils/average etc. into the final
            self.data attribute, so despite having multiple raw FIDs, we are
            really only processing one voxel, so no for loop
        
        """
        self.set = self._block.set                          # dereference the Block with this
        self.raw = self._dataset.get_source_data('prep')    # local reference to input data

        self.do_freq_raw = freq_raw     # flag to redo fft on pre-correction data set
        
        # update coil combination data from associated dataset - if any
        method = self.set.coil_combine_method
        if method in ['External Dataset','External Dataset with Offset']:
            external_dataset = self.set.coil_combine_external_dataset
            if external_dataset is not None:
                #
                # If nfid dims don't match, we take first set of weights from
                # external dataset and repeat them to match current data
                # - this solves issue of using water data to combine metab data
                #
                ext_block = external_dataset.blocks["prep"]
                wt = ext_block.coil_combine_weights.copy()
                ph = ext_block.coil_combine_phases.copy()
                if wt.shape[0] != self.raw.shape[2]:        # compare nfid for both
                    wt = np.tile(np.squeeze(wt[0,:]), self.raw.shape[2])
                    ph = np.tile(np.squeeze(ph[0,:]), self.raw.shape[2])
                    wt.shape = (self.raw.shape[2], self.raw.shape[1])
                    ph.shape = (self.raw.shape[2], self.raw.shape[1])
                self._block.coil_combine_weights = wt
                self._block.coil_combine_phases  = ph
            else:
                if self.freq_raw is None:
                    self.freq_raw = self.freq_current.copy() * 0
#                print(' error - chain_pref_fidsum, no external dataset coil_combine values')
                plot_results = { 'freq_current'  : self.freq_current.copy(),
                                 'freq_summed'   : self.freq_current.copy() * 0,
                                 'freq_adjusted' : self.freq_adjusted.copy() * 0,
                                 'freq_raw'      : self.freq_raw.copy() * 0}
                return plot_results

        # use current values to initialize chain results
        self.weights         = self._block.coil_combine_weights
        self.phases          = self._block.coil_combine_phases
        self.exclude_indices = self._block.exclude_indices
        self.frequency_shift = self._block.frequency_shift
        self.phase_0         = self._block.phase_0

        self.voxel = voxels[0]

        # select the chain processing functor based on the entry point
        if entry == 'all':
            funct_fidsum_all.do_processing_all(self)
        elif entry == 'correct_adjust':
            funct_fidsum_all.do_processing_correct_adjust(self)
        elif entry == 'adjust':
            funct_fidsum_all.do_processing_adjust(self)
        else:
            print('oooops!')

        # These are algorithm Results (and intermediary results)
        self._block.coil_combine_weights = self.coil_combine_weights
        self._block.coil_combine_phases  = self.coil_combine_phases
        self._block.exclude_indices      = self.exclude_indices
        self._block.frequency_shift      = self.frequency_shift
        self._block.phase_0              = self.phase_0
        self._block.data[0,0,0,:]        = self.time_summed.copy()

        # These are GUI results use for PlotPanel displays
        plot_results = { 'freq_current'  : self.freq_current.copy(),
                         'freq_summed'   : self.freq_summed.copy(),
                         'freq_adjusted' : self.freq_adjusted.copy(),
                         'freq_raw'      : self.freq_raw.copy() }
                        
        return plot_results