# Python modules

# 3rd party modules
import numpy as np

# Our modules
import vespa.analysis.functors.funct_spectral_all as funct_spectral_all
from vespa.analysis.chain_base import Chain



class ChainSpectral(Chain):
    """
    Building block object used to create a processing chain for MRS data.

    Processes coil-combined, averaged FID data into frequency domain. Applies
    (optionally) ECC, water removal, phase, dc, apodization, left-shift.

    """

    def __init__(self, dataset, block):
        """
        Chain objects organize Algo (algorithm) calls by setting up access to
        input data and parameters, and creating standard output values for View.

        Base class sets convenience references to:  self._block and self._dataset

        self.data is always initialized as []

        """
        super().__init__(dataset, block)

        self.raw_dims             = self._dataset.raw_dims
        self.raw_dim0             = self._dataset.raw_dims[0]
        self.raw_hpp              = self._dataset.raw_hpp
        
        # processing functor - provides entry points for chain
        self.functor_all = funct_spectral_all.do_processing_all

        self.reset_results_arrays()



    def reset_results_arrays(self):
        """
        A separate method so it can be called outside __init__. Should
        create/set enough results to keep View happy if run() fails.

        """
        spectral_dim0 = self._dataset.spectral_dims[0]
        if len(self.data) != spectral_dim0:
            self.pre_roll        = np.zeros(self.raw_dim0, complex)
            self.kodata          = np.zeros(self.raw_dim0, complex)
            self.freq            = np.zeros(spectral_dim0, complex)
            self.data            = np.zeros(spectral_dim0, complex)

            self.time_fids       = np.zeros((20,self.raw_dim0), complex)
            self.sum_time_fids   = np.zeros(self.raw_dim0,      complex)

            self.svd_data              = np.zeros(spectral_dim0,      complex)
            self.svd_peaks_checked     = np.zeros((20,spectral_dim0), complex)
            self.svd_fids_all          = np.zeros((20,self.raw_dim0), complex)
            self.svd_peaks_checked_sum = np.zeros(spectral_dim0,      complex)
                        


    def run(self, voxels, entry='all'):
        """
        Run is typically called every time a processing setting is changed
        in the parent (block) object. Run processes a single voxel at a time.

        This object maintains previous run() results values until next run().
        This allows the View to update without having to re-run the pipeline.

        The 'entry' keyword adds flexibility to Block-Chain-View relationship.

        """

        # Get 'global' parameters, that DO NOT change with voxel, from Dataset
        #  - these processing/data parameters have to be updated at run time 
        self.spectral_dims        = self._dataset.spectral_dims
        self.spectral_dim0        = self._dataset.spectral_dims[0]
        self.spectral_hpp         = self._dataset.spectral_hpp
        self.zero_fill_multiplier = self._dataset.zero_fill_multiplier
        self.phase_1_pivot        = self._dataset.phase_1_pivot
        
        for voxel in voxels:
            # local copy of input data
            self.data = self._dataset.get_source_data('spectral')
            self.data = self.data[voxel[2],voxel[1],voxel[0],:]
            self.data = self.data.copy()

            # copy 'global' parameters, that DO change with voxel, from Dataset
            self.frequency_shift = self._dataset.get_frequency_shift(voxel)
            self.phase0          = self._dataset.get_phase_0(voxel)
            self.phase1          = self._dataset.get_phase_1(voxel)
            
            # copy block parameters, that DO change with voxel, from Block
            svd_output = self._block.get_svd_output(voxel)

            self.ndp             = self._block.get_data_point_count(voxel)
            self.nssv            = self._block.get_signal_singular_value_count(voxel)
            self.do_fit          = self._block.get_do_fit(voxel)
            self.svd_output      = svd_output
            self.voxel           = voxel

            # select the chain processing functor based on the entry point
            if entry == 'all':
                self.functor_all(self)
            else:
                print('oooops!')

            # save data and parameter results into the Block results arrays
            self._block.data[voxel[2],voxel[1],voxel[0],:] = self.freq.copy()
            self._block.set_svd_output(self.svd_output, voxel)
            self._block.set_do_fit(self.do_fit, voxel)

        # Return values specific to calling Tab that contains this Block.Chain
        # Used to update its self.view (plot_panel_spectrum object).

        plot_results = { 'svd_data'               : self.svd_data.copy(),
                         'svd_peaks_checked'      : self.svd_peaks_checked.copy(),
                         'svd_peaks_checked_sum'  : self.svd_peaks_checked_sum.copy(),
                         'svd_fids_checked_sum'   : self.svd_fids_checked.copy(),
                         'freq'                   : self.freq.copy()   }
                        
        return plot_results
