# Python modules
from abc import ABC, abstractmethod

# 3rd party modules

# Our modules
import vespa.common.util.misc as util_misc
from vespa.common.constants import Deflate


class Block(ABC):
    """
    An abstract base class for Block objects. It can't be instantiated, but all
    blocks inherit from it.

    Subclasses MUST override:  create_chain(), deflate() and inflate()

    Subclasses MAY override: XML_VERSION, get_associated_datasets(),
        set_associated_datasets(), set_dims() and is_identity.

    """
    
    # The XML_VERSION enables us to change the XML output format in the future
    XML_VERSION = "1.0.0"

    @abstractmethod
    def __init__(self, attributes=None):
        self.id = util_misc.uuid()
        self.behave_as_preset = False
        self.chain = None
        self.data = None


    @property
    def is_identity(self):
        """
        Returns True if this class is an identity class. Named for mathematical
        concept of an identity transform; it doesn't change the data it touches.
        Override this if you want, but this code will return correct answer if
        naming is followed.

        """
        return ("identity" in str(type(self)).lower())


    @property
    def dims(self):
        """Data dimensions in a list, read only."""
        return list(self.data.shape[::-1]) if self.data is not None else None

    @property
    def shape(self):
        """Data shape in a list, read only."""
        return list(self.data.shape) if self.data is not None else None


    @abstractmethod
    def create_chain(self, dataset):
        pass


    def get_associated_datasets(self, is_main_dataset=True):
        """
        Returns a list of datasets associated with this object. 
        
        The 'is_main_dataset' flag signals method that it is the top level dataset
        gathering associated datasets. It is used to stop circular logic where one
        or more datasets refer to each other.
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


    @abstractmethod
    def deflate(self, flavor=Deflate.ETREE):
        pass


    @abstractmethod
    def inflate(self, source):
        pass

