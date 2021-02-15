# Python modules

# 3rd party modules

# Our modules
from vespa.analysis.chain_base import Chain


class ChainRaw(Chain):
    """ 
    This is a building block object that can be used to create a processing
    pipeline for MRS data processing.
    
    """
    def __init__(self, dataset, block):
        """
        A BlockRaw object only holds loaded data, there is no associated data
        processing. But each Block has to have a Chain object.
        
        """
        super().__init__(dataset, block)


    def run(self, voxels, entry='all'):
        """ This method is required """
        pass