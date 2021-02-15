# Python modules

# 3rd party modules

# Our modules
from vespa.analysis.chain_base import Chain


class ChainSpectralIdentity(Chain):
    """
    This is a no-op Chain object that is a placeholder in a Dataset whose
    Blocks have not been fully instantiated/inflated yet.

    """
    def __init__(self, dataset, block):
        """ Prepare the base class. """
        super().__init__(dataset, block)

        
    def run(self, voxels, entry='all'):
        pass
