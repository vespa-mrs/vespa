# Python modules

# 3rd party modules
import xml.etree.cElementTree as ElementTree

# Our modules
from vespa.common.mrs_generic_basis import GenericBasis
import vespa.common.util.xml_ as util_xml

from vespa.common.constants import Deflate





class Macromolecule(GenericBasis):
    # FIXME PS - need a docstring


    # The XML_VERSION enables us to change the XML output format in the future 
    XML_VERSION = "1.0.0"

    def __init__(self, attributes=None):
        super().__init__(attributes)



    ##### Standard Methods and Properties #####################################

    def deflate(self, flavor=Deflate.ETREE):
        if flavor == Deflate.ETREE:
            # Make my base class do its deflate work
            e = super().deflate(flavor)

            # Alter the tag name & XML version info   
            e.tag = "macromolecule"
            e.set("version", self.XML_VERSION)
            
            return e
        elif flavor == Deflate.DICTIONARY:
            # My base class does all the work
            return super().deflate(flavor)
        

    def inflate(self, source):
        # Make my base class do its inflate work
        super().inflate(source)


    ##### Object Specific Methods and Properties #####################################

    def default_prior(self):
    
        # set up a default macromolecule model
        
        values_ppm        = [2.346, 2.543, 2.142, 1.638, 1.357, 0.90, 3.81]
        values_area       = [0.2,   0.1,   0.15,  0.3,   0.4,   0.6,  2.0]
        values_phase      = [0,     0,     0,     0,     0,     0,    0.0]
        values_linewidth  = [20.0,  20,    30,    35,    35,    45,   120]
        limits_ppm        = [0.1,   0.1,   0.02,  0.1,   0.1,   0.1,  0.1]
        limits_area       = [1.0,   1,     0.4,   1,     1,     1,    1]
        limits_phase      = [1.0,   1,     1,     1,     1,     1,    1]
        limits_linewidth  = [10,    10,    10,    10,    10,    10,   40]
    
        # create an element tree with prior info in it
    
        e = ElementTree.Element("macromolecule", {"version" : self.XML_VERSION})

        util_xml.TextSubElement(e, "source",    'default')
        util_xml.TextSubElement(e, "source_id", 'default')

        for ppm,     area,     phase,     lw, \
            ppm_lim, area_lim, phase_lim, lw_lim in zip(values_ppm, 
                                                        values_area, 
                                                        values_phase, 
                                                        values_linewidth, 
                                                        limits_ppm, 
                                                        limits_area, 
                                                        limits_phase, 
                                                        limits_linewidth):
            line_element = ElementTree.SubElement(e, "line")
            util_xml.TextSubElement(line_element, "ppm",       ppm)
            util_xml.TextSubElement(line_element, "area",      area)
            util_xml.TextSubElement(line_element, "phase",     phase)
            util_xml.TextSubElement(line_element, "lw",        lw)
            util_xml.TextSubElement(line_element, "ppm_lim",   ppm_lim)
            util_xml.TextSubElement(line_element, "area_lim",  area_lim)
            util_xml.TextSubElement(line_element, "phase_lim", phase_lim)
            util_xml.TextSubElement(line_element, "lw_lim",    lw_lim)

        return e
                
    
    


#--------------------------------------------------------------------
# test code

def _test():

    import vespa.common.util.time_ as util_time
    
    test = Macromolecule()

    class_name = test.__class__.__name__
    filename = "_test_output_"+class_name+".xml"
    element = test.deflate()
    root = ElementTree.Element("_test_"+class_name, { "version" : "1.0.0" })
    util_xml.TextSubElement(root, "timestamp", util_time.now().isoformat())
    root.append(element)
    tree = ElementTree.ElementTree(root)
    tree.write(filename, "utf-8")
    
    tom = 10


if __name__ == '__main__':
    _test()    
    
