# Python modules


# 3rd party modules
import xml.etree.cElementTree as ElementTree

# Our modules
from vespa.common.mrs_data_raw import DataRaw
from vespa.common.constants import Deflate

class DataRawTimeseries(DataRaw):
    """
    This subclass differentiates between multiple SVS data loaded 'into the
    screen' that are time-dependent in the second dimension and those that are
    not time-dependent.

    """
    XML_VERSION = "1.0.0"

    def deflate(self, flavor=Deflate.ETREE):

        if flavor == Deflate.ETREE:
            e = super().deflate(flavor)
            e.tag = "raw_timeseries"
            e.set("version", self.XML_VERSION)

            return e

        elif flavor == Deflate.DICTIONARY:
            return self.__dict__.copy()
