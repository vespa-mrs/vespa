# Python modules

# 3rd party modules

# Our modules
import vespa.analysis.block_prep_fidsum as block_prep_fidsum
from vespa.common.constants import Deflate




class BlockPrepEditFidsum(block_prep_fidsum.BlockPrepFidsum):
    """ 
    Building block to hold the state of a step in an MRS processing chain.
    Includes the functionality to save/recall this object to/from an XML node.

    Contains inputs/results for preprocessing of the raw data from the previous
    block ('raw') in the dataset.blocks list. This step modifies coil/average
    data into a single summed FID array for one dataset.

    """
    # The XML_VERSION enables us to change the XML output format in the future 
    XML_VERSION = "1.0.0"

    
    def __init__(self, attributes=None):
        """ Set up standard functionality of the base class """
        super.__init__(attributes)


    ##### Standard Methods and Properties #####################################

    def __str__(self):

        lines = super().__str__().split('\n')
        lines[0] = "------- {0} Object -------".format(self.__class__.__name__)

        return '\n'.join(lines)


    def deflate(self, flavor=Deflate.ETREE):

        e = super().deflate(self, flavor)

        if flavor == Deflate.ETREE:
            # Alter the tag name & XML version info   
            e.tag = "block_prep_edit_fidsum"
            e.set("version", self.XML_VERSION)
            return e

        elif flavor == Deflate.DICTIONARY:
            return e



    ##### Private Methods #####################################

    
