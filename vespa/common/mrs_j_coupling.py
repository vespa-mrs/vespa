# Python modules

import xml.etree.cElementTree as ElementTree

# 3rd party modules


# Our modules
from vespa.common.constants import Deflate
import vespa.common.util.xml_ as util_xml


class JCoupling(object):
    def __init__(self, attributes=None):
        self.id = 0
        self.value = None
        self.spin1 = None
        self.spin2 = None

        if attributes is not None:
            self.inflate(attributes)


    def __str__(self):
        return self.__unicode__()


    def __unicode__(self):
        lines = [ ]
        lines.append("--- J Coupling %d ---" % self.id)
        value = "[empty]" if self.value is None else ("%f" % self.value) 
        lines.append("value: %s" % value) 
        spin = str(self.spin1) if self.spin1 else "[empty]"
        lines.append("spin 1: %s" % spin) 
        spin = str(self.spin2) if self.spin2 else "[empty]"
        lines.append("spin 2: %s" % spin) 
        
        # __unicode__() must return a Unicode object. In practice the code
        # above always generates Unicode, but we ensure it here.
        return '\n'.join(lines)


    def deflate(self, flavor=Deflate.ETREE):
        if flavor == Deflate.ETREE:
            return util_xml.TextElement("j_coupling", self.value)
            
        elif flavor == Deflate.DICTIONARY:
            d = self.__dict__.copy()
            # I deflate my child objects (spins) too. This is different from
            # deflating to XML where the j_coupling elements don't reference
            # their spin subelements. Why the difference? Preserving the
            # relationship as I do here makes life easier for whoever is
            # handling the deflated object, so it's preferable. However, if I 
            # did that when writing to XML it would make the XML very verbose.
            if d["spin1"]:
                d["spin1"] = d["spin1"].deflate(flavor)
            if d["spin2"]:
                d["spin2"] = d["spin2"].deflate(flavor)

            return d


    def inflate(self, source):
        if hasattr(source, "makeelement"):
            # Quacks like an ElementTree.Element
            if source.text:
                self.value = float(source.text)
        
        elif hasattr(source, "keys"):
            # Quacks like a dict
            for key in list(source.keys()):
                if hasattr(self, key):
                    setattr(self, key, source[key])
        
        
