# Python modules
import xml.etree.cElementTree as ElementTree

# 3rd party modules

# Our modules
import vespa.common.constants as constants
import vespa.common.util.xml_ as util_xml
from vespa.common.constants import Deflate


class MasterParameters(object):
    # The XML_VERSION enables us to change the XML output format in the future 
    XML_VERSION = "1.0.0"

    def __init__(self, attributes=None):

        self.calc_resolution  = 5000 # num_points for Bloch and rfpulse profile calcs.
        self.pulse_bandwidth_convention = constants.PulseConvention.HALF_HEIGHT

        if attributes is not None:
            self.inflate(attributes)
            
            
    def as_html(self):
        """Returns this object rendered as HTML in an ElementTree.Element. The
        HTML is suitable for display in a wx.html control."""
        table = ElementTree.Element("table", {"border" : "1px"})
        tbody = ElementTree.SubElement(table, "tbody")

        lines = ( ("Calculation resolution", self.calc_resolution),
                      ("Pulse bandwidth convention", 
                        self.pulse_bandwidth_convention["display"])
                )

        for line in lines:
            description, data = line
            tr = ElementTree.SubElement(tbody, "tr")
            util_xml.TextSubElement(tr, "td", description)
            util_xml.TextSubElement(tr, "td", data, { "align" : "right" })

        return table


    def deflate(self, flavor=Deflate.ETREE):
        if flavor == Deflate.ETREE:
            # mpe = master parameters element
            mpe = ElementTree.Element("master_parameters",
                                      {"version" : self.XML_VERSION})

            util_xml.TextSubElement(mpe, "calc_resolution",
                                    str(self.calc_resolution))
            util_xml.TextSubElement(mpe, "pulse_bandwidth_convention",
                                    self.pulse_bandwidth_convention["db"])

            return mpe
        elif flavor == Deflate.DICTIONARY:
            return self.__dict__.copy()


    def inflate(self, source):

        if hasattr(source, "makeelement"):
            # Quacks like an ElementTree.Element
            self.calc_resolution = int(source.findtext("calc_resolution"))
            self.pulse_bandwidth_convention = source.findtext("pulse_bandwidth_convention")
            self.pulse_bandwidth_convention = constants.PulseConvention.get_type_for_value(self.pulse_bandwidth_convention, "db")
        elif hasattr(source, "keys"):
            # Quacks like a dict
            for key in list(source.keys()):
                if hasattr(self, key):
                    setattr(self, key, source[key])
                    
                    
                    