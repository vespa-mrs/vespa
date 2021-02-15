# Python modules
from abc import ABC, abstractmethod


class Chain(ABC):
    """
    An abstract base class for Chain objects. It can't be instantiated, but all
    chains inherit from it and must have the abstract methods shown below.

    Each Block object has a chain object reference, the set of Chain objects
    perform the MRS worflow for a Dataset.

    """
    @abstractmethod
    def __init__(self, dataset, block):
        """ all subclasses must include this method """
        self._dataset = dataset
        self._block   = block
        self.data     = []

        # Set local values for data acquisiton parameters.
        #  - these do not change over time, so we can set them here
        
        self.sw        = dataset.sw
        self.frequency = dataset.frequency
        self.resppm    = dataset.resppm
        self.echopeak  = dataset.echopeak
        self.is_fid    = dataset.is_fid
        self.seqte     = dataset.seqte
        self.seqtr     = dataset.seqtr
        self.nucleus   = dataset.nucleus


    @abstractmethod
    def run(self, voxels, entry='all'):
        """ all subclasses must include this method """
        pass


    def reset_results_arrays(self):
        """ reminder that subclasses may want to override this method """
        pass
