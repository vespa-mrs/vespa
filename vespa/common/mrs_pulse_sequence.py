# Python modules

import xml.etree.cElementTree as ElementTree
import copy

# 3rd party modules


# Our modules
import vespa.common.constants as constants
import vespa.common.util.xml_ as util_xml
import vespa.common.util.time_ as util_time
import vespa.common.rfp_pulse_design as rfp_pulse_design
from vespa.common.constants import Deflate


class PulseSequence(object):
    """
    PulseSequence() - an object describing a simulated pulse
    sequence for use in the Simulate application.
    """
    # The XML_VERSION enables us to change the XML output format in the future 
    XML_VERSION = "1.0.0"

    def __init__(self, attributes=None):
        self.id = None
        self.name = ""
        self.is_public = False
        self.created = None
        self.creator = ""
        self.comment = ""
        self.loop_labels = [ ]
        self.sequence_code = ""
        self.binning_code = constants.DEFAULT_BINNING_CODE
        self.user_static_parameters = [ ]
        self.pulse_projects = [ ]

        if attributes is not None:
            self.inflate(attributes)

        if self.comment is None:
            self.comment = ""
        
        if self.created is None:
            self.created = util_time.now()
            
                
    __doc = """True if the pulse seq is frozen, False otherwise. Note that 
    this is only accurate in the context of the pulse seq management dialog 
    and may be inaccurate elsewhere or even raise an error. Use only within the 
    pulse seq management dialog and its children!
    """
    def __get_is_frozen(self):
        # self.experiment_names exists only within the pulse seq management
        # dialogs.
        return self.is_public or bool(self.experiment_names)
    is_frozen = property(__get_is_frozen, doc=__doc)


    def __str__(self):
        return self.__unicode__()


    def __unicode__(self):
        lines = [ ]
        id_ = self.id if self.id else "[No Id]"
        lines.append("--- Pulse Sequence %s ---" % id_)
        lines.append("Name: %s" % self.name)
        lines.append("Public: %s" % ("True" if self.is_public else "False"))
        s = self.created.isoformat() if self.created else "[empty]"
        lines.append("Created: %s" % s)
        lines.append("Creator: %s" % self.creator)
        lines.append("Comment (abbr.): %s" % self.comment[:40])
        parameters = [ "[%s (%s): %s]" % (parameter.name, parameter.type, parameter.default) \
                        for parameter in self.user_static_parameters]
        lines.append("User static params: %s" % ", ".join(parameters))
        pulse_designs = [ ]
        for pulse_design in self.pulse_projects:
            pulse_designs.append("%s (%s)" % (pulse_design.name, \
                                               pulse_design.id))
        lines.append("Pulse Designs: %s" % ", ".join(pulse_designs))

        for i, label in enumerate(self.loop_labels):
            if label:
                lines.append("Loop %d: %s" % (i, label))
            
        code = "None"    
        if self.sequence_code:
            code = ("%d characters" % len(self.sequence_code))
        lines.append("Sequence code: %s" % code)
        
        code = "None"    
        if self.binning_code:
            code = ("%d characters" % len(self.binning_code))
        lines.append("Binning code: %s" % code)

        # __unicode__() must return a Unicode object. In practice the code
        # above always generates Unicode, but we ensure it here.
        return '\n'.join(lines)


    def clone(self):
        """Creates & returns a new pulse_sequence that looks just like this 
           one, but with a different id.
        """
        pulse_sequence = copy.deepcopy(self)
        pulse_sequence.id = None
        pulse_sequence.created = util_time.now()
        pulse_sequence.is_public = False
        
        return pulse_sequence


    def deflate(self, flavor=Deflate.ETREE):
        if flavor == Deflate.ETREE:
            # pse = pulse sequence element
            pse = ElementTree.Element("pulse_sequence", {"id" : self.id,
                                                         "version" : self.XML_VERSION})
                                            
            for attribute in ("name", "created", "creator", "comment", ):
                util_xml.TextSubElement(pse, attribute, getattr(self, attribute))
            
            for label in self.loop_labels:
                if label:
                    util_xml.TextSubElement(pse, "loop_label", label)
            
            for parameter in self.user_static_parameters:
                pse.append(parameter.deflate(flavor))

            for pulse_design in self.pulse_projects:
                pse.append(pulse_design.deflate(flavor))

            if self.sequence_code:
                util_xml.TextSubElement(pse, "sequence_code", self.sequence_code)

            if self.binning_code:
                util_xml.TextSubElement(pse, "binning_code", self.binning_code)

            return pse
            
        elif flavor == Deflate.DICTIONARY:
            d = self.__dict__.copy()
            
            d["user_static_parameters"] = [param.deflate(flavor) for param 
                                            in self.user_static_parameters]

            d["pulse_projects"] = [item.deflate(flavor) for item in self.pulse_projects]
            
            return d


    def inflate(self, source):
        if hasattr(source, "makeelement"):
            # Quacks like an ElementTree.Element
            self.id = source.get("id")
            self.name = source.findtext("name")
            self.creator = source.findtext("creator")
            self.created = source.findtext("created")
            self.created = util_time.datetime_from_iso(self.created)
            self.comment = source.findtext("comment")
            self.sequence_code = source.findtext("sequence_code")
            self.binning_code = source.findtext("binning_code")

            self.loop_labels = [label.text for label
                                           in source.findall("loop_label")
                                           if label.text]
            
            self.user_static_parameters = \
                                [UserStaticParameter(parameter) 
                                    for parameter 
                                    in source.findall("user_static_parameter")]

            tmp = [rfp_pulse_design.PulseDesign(item) for item in source.findall("pulse_design")]
            self.pulse_projects = tmp


        elif hasattr(source, "keys"):
            # Quacks like a dict
            for key in list(source.keys()):
                if hasattr(self, key):
                    setattr(self, key, source[key])

            self.user_static_parameters = \
                                [UserStaticParameter(parameter) for parameter 
                                    in source.get("user_static_parameters", [ ])]

            self.pulse_projects = \
                                [rfp_pulse_design.PulseDesign(item) for item 
                                    in source.get("pulse_projects", [ ])]
            

                
class UserStaticParameter(object):
    def __init__(self, attributes=None):
        self.id = None
        self.name = ""
        self.type = ""
        self.default = ""

        if attributes is not None:
            self.inflate(attributes)


    def __str__(self):
        return self.__unicode__()


    def __unicode__(self):
        lines = [ ]
        id_ = str(self.id) if self.id else "[No Id]" 
        lines.append("--- Pulse Seq User Static Param %s ---" % id_)
        lines.append("Name: %s" % self.name)
        lines.append("Type: %s" % self.type)
        lines.append("Default: %s" % self.default)
    
        # __unicode__() must return a Unicode object. In practice the code
        # above always generates Unicode, but we ensure it here.
        return '\n'.join(lines)


    def deflate(self, flavor=Deflate.ETREE):
        if flavor == Deflate.ETREE:
            # p_element = parameter element
            p_element = ElementTree.Element("user_static_parameter")
                                            
            util_xml.TextSubElement(p_element, "name", self.name)
            util_xml.TextSubElement(p_element, "type", self.type)
            util_xml.TextSubElement(p_element, "default", self.default)

            return p_element
            
        elif flavor == Deflate.DICTIONARY:
            return self.__dict__.copy()


    def inflate(self, source):
        if hasattr(source, "makeelement"):
            # Quacks like an ElementTree.Element
            self.id = source.get("id")
            self.name = source.findtext("name")
            self.type = source.findtext("type")
            self.default = source.findtext("default")

        elif hasattr(source, "keys"):
            # Quacks like a dict
            for key in list(source.keys()):
                if hasattr(self, key):
                    setattr(self, key, source[key])
            
