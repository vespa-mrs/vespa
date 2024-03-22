# Python modules
import xml.etree.cElementTree as ElementTree

# 3rd party modules

# Our modules
import vespa.common.util.xml_ as util_xml
from vespa.common.constants import Deflate



class PriorMetabolite(object):
    """
    Describes a metabolite using lists of ppm, area and phase values. The Nth
    element of each list corresponds to the Nth element in the other lists.

    """
    XML_VERSION = "1.0.0"

    def __init__(self, attributes=None):

        self.spins     = 0
        self.dims      = ['',0,0,0]     # re. Experiment format ['name', loop1, loop2, loop3]
        self.group     = []             # currently unused but we're treating it as text
        self.ppms      = []
        self.areas     = []
        self.phases    = []
    
        if attributes is not None:
            self.inflate(attributes)


    ##### Standard Methods and Properties #####################################

    @property
    def name(self):
        return self.dims[0]
        

    def __str__(self):

        lines = []
        lines.append("------- {0} Object -------".format(self.__class__.__name__))
        lines.append("- (%d lines)" % len(self.ppms))

        lines.append("Name: %s" % self.name)
        lines.append("Spins: %s" % str(self.spins))
        
        if self.ppms:
            for i in range(len(self.ppms)):
                line  = "Line " +str(i)+" : "
                line += " PPM="  +str(self.ppms[i])
                line += " Area=" +str(self.areas[i])
                line += " Phase="+str(self.phases[i])
                line += " Loop1="+str(self.dims[1])
                line += " Loop2="+str(self.dims[2])
                line += " Loop3="+str(self.dims[3])
                lines.append(line)
        
        return '\n'.join(lines)


    def deflate(self, flavor=Deflate.ETREE):
        if flavor == Deflate.ETREE:
            e = ElementTree.Element("prior_metabolite", {"version" : self.XML_VERSION})
            
            util_xml.TextSubElement(e, "name", self.name)
            for dim in self.dims[1:]:
                util_xml.TextSubElement(e, "dim", dim)
            util_xml.TextSubElement(e, "spins", self.spins)

            # Currently not using groups, so we expand here for the zip() below
            groups = [''] * len(self.ppms) if len(self.group) != len(self.ppms) else self.group
                 
            for group, ppm, area, phase in zip(groups, self.ppms, self.areas, self.phases):
                line_element = ElementTree.SubElement(e, "line")
                util_xml.TextSubElement(line_element, "group", group)
                util_xml.TextSubElement(line_element, "ppm",   ppm)
                util_xml.TextSubElement(line_element, "area",  area)
                util_xml.TextSubElement(line_element, "phase", phase)

            return e
            
        elif flavor == Deflate.DICTIONARY:
            d = { }
            
            for attribute in ("spins", "group", "ppms", "areas", "phases", ):
                d[attribute] = getattr(self, attribute)
            
            d["name"] = self.name
            d["dims"] = self.dims[1:]
            
            return d

            
    def inflate(self, source):
        if hasattr(source, "makeelement"):
            # Quacks like an ElementTree.Element
            dims = [int(round(float(dim.text))) for dim in source.findall("dim")]
            self.dims   = [source.findtext("name")] + dims
            self.spins  = int(source.findtext("spins"))
            self.group  = [group.text        for group in source.iter("group")]
            self.ppms   = [float(ppm.text)   for ppm   in source.iter("ppm")]
            self.areas  = [float(area.text)  for area  in source.iter("area")]
            self.phases = [float(phase.text) for phase in source.iter("phase")]
        
        elif hasattr(source, "keys"):
            # Quacks like a dict
            for key in list(source.keys()):
                # the "name" key gets a little special handling
                if hasattr(self, key) and (key != "name"):
                    setattr(self, key, source[key])
            
            self.dims = [source["name"]] + source["dims"]


    def sim_display(self, dims):
        lines = [ ]
        for i in range(len(self.ppms)):
            line = "%s\t%s\t%s\t%s\t" % tuple(dims)
            line += (str(i)+'\t'+str(self.ppms[i])+'\t'+str(self.areas[i])+'\t'+str(self.phases[i]))
            lines.append(line)
        
        return '\n'.join(lines)                


    def get_rows(self):
        lines = []
        for ppm, area, phase in zip(self.ppms,self.areas,self.phases):
            lines.append([ppm,area,phase])
        
        return lines                


#--------------------------------------------------------------------
# test code

def _test():

    bob = PriorMetabolite()
    
    tom = 10  


if __name__ == '__main__':
    _test()