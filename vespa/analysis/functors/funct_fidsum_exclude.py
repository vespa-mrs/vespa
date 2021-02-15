# Python imports

# 3rd party imports
import numpy as np

# Vespa imports
from vespa.analysis.algos.fida_rmbadaverages    import op_rmbadaverages
from vespa.analysis.algos.fida_rmNworstaverages import op_rmNworstaverages

EXCLUDE_MENU_ITEMS = ['Manual',
                      'Remove Bad Averages (fid-a)',
                      'Remove N Worst Averages (fid-a)',
                      ]



def exclude_remove_bad_averages_fida(data, chain):
    """ 
    Remove corrupted FIDs from data set - uses op_rmbadaverages from FID-A
    
    Inputs:
        data (ndarray, complex): shape=(1, 1, nfids, npts) kspace data for acquired FIDs
        chain (object): control values for processing
    Outputs:
        exclude_index (ndarray, int): Indices of FIDs to remove, can be empty array

    """
    ds      = chain._dataset
    set     = chain.set
    nsd     = set.fida_bad_threshold     # Setting the number of standard deviations
    domain  = set.fida_bad_domain        # 't' or 'f'

    _, ncoil, nfid, npts = data.shape

    exclude_indices = []
    metrics  = []
    max_iter = 8

    for iter in range(max_iter):

        data_in         = np.delete(np.squeeze(data.copy()), exclude_indices, axis=0)   # [nfid, npts] here
        current_indices = np.delete(np.arange(nfid, dtype=int), exclude_indices)

        bad_indices, metric = op_rmbadaverages(data_in, ds.sw, nsd=nsd, domain=domain)

        if len(bad_indices)>0:
            metrics.append(metric)
            if exclude_indices == []:
                exclude_indices = current_indices[bad_indices]
            else:
                exclude_indices = np.unique(np.r_[exclude_indices, current_indices[bad_indices]])
            if len(exclude_indices) > (nfid-3):
                # no good fids is an error, do a reset
                exclude_indices = []
                break
        else:
            break


    return exclude_indices
  

def exclude_remove_n_worst_averages_fida(data, chain):
    """ 
    Remove corrupted FIDs from data set - uses op_rmNworstaverages from FID-A
    
    Inputs:
        data (ndarray, complex): shape=(1, 1, nfids, npts) kspace data for acquired FIDs
        chain (object): control values for processing
    Outputs:
        exclude_index (ndarray, int): Indices of FIDs to remove, can be empty array

    """
    ds  = chain._dataset
    n   = chain.set.fida_n_worst

    exclude_indices, metric = op_rmNworstaverages(np.squeeze(data.copy()), n, ds.sw)   # [nfid, npts] here

    return exclude_indices

