# Python modules


# 3rd party modules
import xml.etree.cElementTree as ElementTree

# Our modules
from vespa.common.mrs_data_raw import DataRawFidsum
from vespa.common.constants import Deflate


class DataRawWbnaa(DataRawFidsum):
    """
    A subclass of DataRawFidsum that differentiate between summed and regular, 
    unsummed raw data.
    
    """
    XML_VERSION = "1.0.0"

    def __init__(self, attributes=None):
        super().__init__(attributes)


    def deflate(self, flavor=Deflate.ETREE):
        if flavor == Deflate.ETREE:
            # Make my base class do its deflate work
            e = super().deflate(self, flavor)

            # Alter the tag name & XML version info   
            e.tag = "raw_wbnaa"
            e.set("version", self.XML_VERSION)

            # If this class had any custom attributes, I'd add them here.

            return e

        elif flavor == Deflate.DICTIONARY:
            return self.__dict__.copy()
