# Python modules

import xml.etree.cElementTree as ElementTree

# 3rd party modules

# Our modules
from vespa.common.constants import Deflate
import vespa.common.util.xml_ as util_xml


class Spin(object):
    def __init__(self, attributes=None):
        self.id = 0
        self.isotope = ""
        self.chemical_shift = None
        
        if attributes is not None:
            self.inflate(attributes)
            
       
    def __str__(self):
        return self.__unicode__()


    def __unicode__(self):
        lines = [ ]
        id_ = self.id if self.id else 0
        lines.append("--- Spin %d ---" % id_) 
        lines.append("Isotope: %s" % self.isotope) 
        shift = ("%f" % self.chemical_shift) if self.chemical_shift else "[None]"
        lines.append("Shift: %s" % shift) 
        
        # __unicode__() must return a Unicode object. In practice the code
        # above always generates Unicode, but we ensure it here.
        return '\n'.join(lines)


    def deflate(self, flavor=Deflate.ETREE):
        if flavor == Deflate.ETREE:
            spin_element = ElementTree.Element("spin")

            util_xml.TextSubElement(spin_element, "isotope", self.isotope)            
            util_xml.TextSubElement(spin_element, "chemical_shift",
                                    self.chemical_shift)
            
            return spin_element
            
        elif flavor == Deflate.DICTIONARY:
            return self.__dict__.copy()

            
    def inflate(self, source):
        if hasattr(source, "makeelement"):
            # Quacks like an ElementTree.Element
            isotope = source.findtext("isotope")
            if isotope:
                self.isotope = isotope
            chemical_shift = source.findtext("chemical_shift")
            if chemical_shift:
                self.chemical_shift = float(chemical_shift)
        
        elif hasattr(source, "keys"):
            # Quacks like a dict
            for key in list(source.keys()):
                if hasattr(self, key):
                    setattr(self, key, source[key])
               

