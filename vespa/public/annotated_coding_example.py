# Python modules
                     # bjs - this is in all our Modules

# 3rd party modules                         # bjs - these three comments are all in our modules even if empty

# Our modules
import vespa.common.util.misc as util_misc
from vespa.common.constants import Deflate


"""
This file was created as an example of how we structure and comment code in the
Vespa project, and in general in the Soher lab.  It came from the Vespa project,
analysis/src/block.py specifically.  
"""


class Block(object):                                    # bjs - CamelCase for classes
    """
    This is an base class for blocks. It's purpose is to define the minimal 
    interface that all blocks must support. It also provides an implementation      # bjs - docstring when class is finalized
    for some methods. 

    Subclasses must override create_chain(), deflate() and inflate().

    Subclasses may want to override XML_VERSION, get_associated_datasets(),
    set_associated_datasets(), set_dims() and is_identity.

    """
    
    # The XML_VERSION enables us to change the XML output format in the future
    XML_VERSION = "1.0.0"

    def __init__(self, attributes=None):                        # bjs - we add ID attributes to most classes
        self.id = util_misc.uuid()
        self.set = None
        self.chain = None
        
        self.behave_as_preset = False


    @property
    def is_identity(self):                              # bjs - properties come right after 'init' method
        """Returns True if this class is an identity class. Identity classes
        are named for the mathematical concept of an identity transform; it
        doesn't change the data it touches.
        """
        # Subclasses should feel free to override this, although this code will 
        # return the correct answer if all classes follow our current naming
        # convention.
        return ("identity" in str(type(self)).lower())


    def create_chain(self, dataset):                            # bjs - methods are in lower case with underscore
        raise NotImplementedError


    def get_associated_datasets(self, is_main_dataset=True): 
        """
        Returns a list of datasets associated with this object. 
        
        The 'is_main_dataset' flag allows the method to know if it is the top
        level dataset gathering associated datasets, or some dataset that is
        only associated with the top dataset. This is used to stop circular
        logic conditions where one or more datasets refer to each other.
        """
        return []


    def set_associated_datasets(self, datasets):                          
        pass
        

    def set_dims(self, dataset):
        """
        Given a Dataset object, this is an opportunity for the block object 
        to ensure that its dims match those of the parent dataset. The block
        is typically only interested in dataset.raw_dims or 
        dataset.spectral_dims.

        The default implementation does nothing.

        This function is typically called after some dimension change on a 
        block (probably via the GUI).
        """
        pass


    def deflate(self, flavor=Deflate.ETREE):                        # bjs - most classes have inflate/deflate methods
        raise NotImplementedError


    def inflate(self, source):
        raise NotImplementedError

