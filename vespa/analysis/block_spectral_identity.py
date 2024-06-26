# Python modules

# 3rd party modules
from xml.etree.cElementTree import Element

# Our modules
import vespa.analysis.chain_spectral_identity as chain_spectral_identity
import vespa.analysis.block as block
from vespa.common.constants import Deflate



class _Settings(object):

    XML_VERSION = "1.0.0"

    def __init__(self, attributes=None):
        """
        This is a no-op Block.Settings object that is a placeholder in a
        Dataset whose Blocks have not been fully instantiated/inflated yet.

        """
        pass


    @property
    def zero_fill_multiplier(self):
        return 1

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



class BlockSpectralIdentity(block.Block):
    """
    This is a no-op Block object that is a placeholder in a Dataset whose
    Blocks have not been fully instantiated/inflated yet.

    """
    XML_VERSION = "1.0.0"

    def __init__(self, attributes=None):

        super().__init__(attributes)
        self.set = _Settings(attributes)
        self.data = None


    ##### Standard Methods and Properties #####################################

    def __str__(self):
        lines = []
        lines.append("------- {0} Object -------".format(self.__class__.__name__))
        lines.append("Nothing here")
        return '\n'.join(lines)


    def create_chain(self, dataset):
        self.chain = chain_spectral_identity.ChainSpectralIdentity(dataset, self)


    def deflate(self, flavor=Deflate.ETREE):
        if flavor == Deflate.ETREE:
            e = Element("block_spectral_identity", { "id":self.id,
                                                     "version":self.XML_VERSION})
            e.append(self.set.deflate())
            return e

        elif flavor == Deflate.DICTIONARY:
            raise NotImplementedError


    def inflate(self, source):
        if hasattr(source, "makeelement"):
            # Quacks like an ElementTree.Element
            self.set = _Settings(source.find("settings"))

        elif hasattr(source, "keys"):
            # Quacks like a dict
            for key in list(source.keys()):
                if key == "set":
                    setattr(self, key, source[key])

