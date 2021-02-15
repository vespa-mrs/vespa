# Python modules

# 3rd party modules
import numpy as np

# Our modules
import vespa.analysis.functors.funct_timeseries_all as funct_timeseries_all
from vespa.analysis.chain_base import Chain



class ChainPrepTimeseries(Chain):
    """ 
    Building block object used to create a processing chain for MRS data where
    we have multiple FID data sets that describe data that change through time.

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

        # these are results for display in the Tab
        self.time_all = None
        self.freq_all = None

        self.reset_results_arrays()
        
 
    def reset_results_arrays(self):
        """
        A separate method so it can be called outside __init__. Should
        create/set enough results to keep View happy if run() fails.
        
        """
        raw_dims = self._dataset.raw_dims
        raw_dims = list(raw_dims[::-1])
        
        if len(self.data) != raw_dims:
            self.freq_all = np.zeros(raw_dims, complex)
            self.time_all = np.zeros(raw_dims, complex)


    def run(self, voxels, entry='all', do_calculate=False):
        """
        Run is typically called every time a processing setting is changed
        in the parent (block) object. Run processes a single voxel at a time.

        This object maintains previous run() results values until next run().
        This allows the View to update without having to re-run the pipeline.

        The 'entry' keyword adds flexibility to Block-Chain-View relationship.
        
        """
        self.voxel = voxels[0]

        self.calculate_flag = do_calculate

        # Note. in this case we are processing all raw data into the data 
        #   attribute, so despite having multiple raw FIDs, we are really 
        #   only processing one voxel, so no for loop


        # local reference to input data
        self.raw = self._dataset.get_source_data('prep')
        
    
        # select the chain processing functor based on the entry point
        if entry == 'all':
            funct_timeseries_all.do_processing_all(self)
        else:
            print('oooops!')

        # save data and parameter results into the Block results arrays
        self._block.data[0,0,:,:] = self.time_all.copy()
 
        # Return values specific to calling Tab that contains this Block.Chain
        # Used to update its self.view (plot_panel_spectrum object).

        plot_results = { 'freq_current' : self.freq_all[0,0,self.voxel,:].copy(),
                         'freq_all'     : self.freq_all[0,0,self.voxel,:].copy()      }
                        
        return plot_results