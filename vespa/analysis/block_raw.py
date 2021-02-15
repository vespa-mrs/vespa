# Python modules

# 3rd party modules
from xml.etree.cElementTree import Element

# Our modules
import vespa.analysis.chain_raw as chain_raw
import vespa.analysis.block as block
import vespa.common.mrs_data_raw as mrs_data_raw
import vespa.common.util.xml_ as util_xml
from vespa.common.constants import Deflate



class _Settings(object):
    """
    Settings object contains the parameter inputs used for processing in the 
    Chain object in this Block. Having a separate object helps to delineate 
    inputs/outputs and to simplify load/save of preset values.

    This object can also save/recall these values to/from an XML node.

    """   
    XML_VERSION = "1.0.0"

    def __init__(self, attributes=None):
        """ Currently there are no input parameters set in this object. """
        pass


    def deflate(self, flavor=Deflate.ETREE):
        if flavor == Deflate.ETREE:
            e = Element("settings", {"version" : self.XML_VERSION})
            return e
            
        elif flavor == Deflate.DICTIONARY:
            return self.__dict__.copy()


    def inflate(self, source):
        if hasattr(source, "makeelement"):
            # Quacks like an ElementTree.Element
            pass

        elif hasattr(source, "keys"):
            # Quacks like a dict
            for key in list(source.keys()):
                if hasattr(self, key):
                    setattr(self, key, source[key])



class BlockRaw(block.Block, mrs_data_raw.DataRaw):
    """ 
    Building block to hold the state of a step in an MRS processing chain.
    Includes the functionality to save/recall this object to/from an XML node.

    Raw Blocks hold data loaded from file. They don't have 'inputs' for a
    Chain object. They do have the attributes inheirited from DataRaw. 
    
    For some subclasses, one or more DataRaws objects can be held in a Block
    object, such as the On/Off/Add/Diff objects of an Edited data file.

    """
    XML_VERSION = "1.0.0"
    
    def __init__(self, attributes=None):
        block.Block.__init__(self, attributes)
        mrs_data_raw.DataRaw.__init__(self, attributes)  

        # processing parameters
        self.set = _Settings(attributes)


    ##### Standard Methods and Properties #####################################

    @property
    def dims(self):
        """Data dimensions in a list, read only."""
        return list(self.data.shape[::-1]) if self.data is not None else None


    def __str__(self):
        lines = mrs_data_raw.DataRaw.__str__(self).split('\n')
        # Replace the heading line
        lines[0] = "------- {0} Object -------".format(self.__class__.__name__)
        lines.append("No printable data ")
        return '\n'.join(lines)


    def create_chain(self, dataset):
        self.chain = chain_raw.ChainRaw(dataset, self)


    def concatenate(self, new):
        # This is a method from DataRaw that's not supported here. 
        raise NotImplementedError


    def deflate(self, flavor=Deflate.ETREE):
        if flavor == Deflate.ETREE:
            
            # Call base class - then update for subclass
            e = mrs_data_raw.DataRaw.deflate(self, flavor)
            e.tag = "block_raw"
            e.set("version", self.XML_VERSION)
            
            # Now I deflate the attribs that are specific to this class
            e.append(self.set.deflate())
            
            return e

        elif flavor == Deflate.DICTIONARY:
            raise NotImplementedError


    def inflate(self, source):

        # Make my base class do its inflate work
        mrs_data_raw.DataRaw.inflate(self, source)

        # Now I inflate the attribs that are specific to this class
        if hasattr(source, "makeelement"):
            # Quacks like an ElementTree.Element
            self.set = util_xml.find_settings(source, "block_raw_settings")
            self.set = _Settings(self.set)

        elif hasattr(source, "keys"):
            # Quacks like a dict
            for key in list(source.keys()):
                if key == "set":
                    setattr(self, key, source[key])


        


