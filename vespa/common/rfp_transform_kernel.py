# Python modules

import copy
import xml.etree.cElementTree as ElementTree

# Our Modules
import vespa.common.util.misc as util_misc
import vespa.common.util.xml_ as util_xml
import vespa.common.util.time_ as util_time

from vespa.common.constants import Deflate

# _NBSP is a Unicode Non-Breaking SPace. Used in HTML output.
_NBSP = "\xA0"

# Default code for a pulse create algorithm. Creates a complex array of ones.
DEFAULT_ALGORITHM1 = """import numpy as np

def run(trans_desc):
    
    npts = pulse_desc.time_steps

    rf = np.zeros(npts, dtype=np.COMPLEX64)
    gr = np.zeros(npts, dtype=np.FLOAT32)

    return rf, gr
"""

DEFAULT_ALGORITHM = """import numpy as np
from vespa.common.transform_run_exception import TransformRunException

def run(trans_desc):
    
    param = trans_desc.parameters
    
    npts      = param["time_steps"]  # int
    duration  = param["duration"]    # float, msec
    
    if npts < 5:
        error_msg =  "Cannot create pulse waveform with less than five points"
        raise TransformRunException(error_msg, 1)
    
    dwell_time = 0.001 * float(duration)/float(npts)
    
    ncycles = 6
    
    t = (np.arange(npts) - (npts/2.0)) / (npts/2.0)
    y = 2 * np.pi * ncycles * t + 0.00001   # 0.0001 to avoid div by zero
    b1 = np.sin(y) / y
    
    # apply a hamming filter
    b1 = b1 * 4 * ncycles * (0.54 + 0.46*np.cos(np.pi*t)) / npts
    
    b1x = np.arange(npts) * dwell_time
    g1x = np.arange(npts) * dwell_time
    
    # scale empirically to 90 degree tip
    b1 = 0.1266667 * b1/np.max(b1)
    gr = np.ones(npts)

    b1 = np.array(b1)
    if not np.iscomplexobj(b1):
        b1r = b1.real.copy()
        b1i = b1r.copy() * 0
        b1 = b1r + 1j*b1i

    print "max b1 = "+ str(np.max(b1))

    #b1  = b1.tolist()
    #gr  = gr.tolist()
    #b1x = xaxisb.tolist()
    #g1x = xaxisg.tolist()

    return b1, b1x, gr, g1x, None
"""



class TransformKernel(object):
    """
    TransformKernel: Contains information to determine a Create or Modify type
    of Transform. This includes algorithm and gui setup values.
    
    """
    # The XML_VERSION enables us to change the XML output format in the future 
    XML_VERSION = "1.0.0"

    def __init__(self, attributes=None):
        
        # Identifying information
        # - all kernels have these entries
        self.id         = ''
        self.type       = 'Create Transform'
        self.name       = ''
        self.menu_label = ''
        self.is_public  = False
        self.created    = None
        self.creator    = ""
        self.comment    = ""

        # The kernel attributes below are needed in order to display and 
        # control required and suggested parameters in the Kernel Editor.
        # Values listed below are DEFAULT values and should be in viable
        # ranges for your algorithm. These are not the final values actually
        # displayed in the main program. Those are stored at the Transform
        # object level.


        # suggested global parameters - can be hidden
        self.hide_file1      = True     # use to browse for one or two file 
        self.hide_file2      = True     # names to use in your algorithm
        self.file1_label     = ''
        self.file2_label     = ''
        self.hide_time_steps = False
        self.hide_duration   = False    # creation
        self.hide_tip_angle  = False
        self.hide_bandwidth  = False    # creation
        self.time_steps      = 128      # int, number of points
        self.duration        = 1.0      # float, total time in msec
        self.tip_angle       = 90.0     # float, in degrees
        self.bandwidth       = 1.0      # float, in kHz
        self.deprecated = False
        
        # When a TransformKernel is validated, all required, suggested AND
        # user defined parameters from the GUI are stored in the 
        # transform_kernel_controls list below.
        #
        # This may include different default values than those listed above as 
        # the user may change them in the Kernel editor. These default values 
        # take precedence when the TabTransform instantiates a finished kernel.
        #
        # Only the "fileX" parameters (if included) are treated differently in
        # the TabTransform since they are a pain to deal with at run time.
        self.transform_kernel_controls = []

        # we set up the default controls as defined above 

        globs = [[self.hide_time_steps,  'Long', 'Time Steps [int]', self.time_steps, 'time_steps'],
                 [self.hide_duration,  'Double', 'Duration [msec]',  self.duration,   'duration'],
                 [self.hide_tip_angle, 'Double', 'Tip Angle [deg]',  self.tip_angle,  'tip_angle'],
                 [self.hide_bandwidth, 'Double', 'Bandwidth [kHz]',  self.bandwidth,  'bandwidth']]
        
        for row in globs:
            if not row[0]:
                control = TransformKernelControl({'type':row[1], 'name':row[2], 'default':row[3], 'variable':row[4]})
                self.transform_kernel_controls.append(control)
        
        self.algorithm_code = DEFAULT_ALGORITHM

        # referrers is a (possibly empty) list of 2-tuples of (id, name).
        # This contains all of the pulse sequences that refer to this 
        # pulse project.
        self.referrers = [ ]

        if attributes is not None:
            self.inflate(attributes)

        if self.comment is None:
            self.comment = ""
        
        if self.created is None:
            self.created = util_time.now()


    __doc = """True if the transform is frozen, False otherwise. Note that 
    this is only accurate in the context of the transform editor dialog 
    and may be inaccurate elsewhere or even raise an error. Use only within the 
    transform editor dialog and its children!
    """
    def __get_is_frozen(self):
        # see mrs_pulse_sequence for more code/info
        return self.is_public 
    is_frozen = property(__get_is_frozen, doc=__doc)
    

    def __str__(self):
        return self.__unicode__()

    def __unicode__(self):
        lines = [ ]
        id_ = self.id if self.id else "[No Id]"
        lines.append("--- Transform Kernel %s ---" % id_)
        lines.append("Type           : %s" % self.type)
        lines.append("Name           : %s" % self.name)
        lines.append("Menu Label: %s" % self.menu_label)
        lines.append("Public         : %s" % ("True" if self.is_public else "False"))
        s = self.created.isoformat() if self.created else "[empty]"
        lines.append("Created        : %s" % s)
        lines.append("Creator        : %s" % self.creator)
        lines.append("Comment (abbr.): %s" % self.comment[:40])
        
        lines.append("Hide File1     : %s" % str(self.hide_file1))
        lines.append("Hide File2     : %s" % str(self.hide_file2))
        lines.append("Hide Time Step : %s" % str(self.hide_time_steps))
        lines.append("Hide Duration  : %s" % str(self.hide_Duration))
        lines.append("Hide Tip Angle : %s" % str(self.hide_tip_angle))
        lines.append("Hide Bandwidth : %s" % str(self.hide_bandwidth))

        lines.append("File1 Label    : %s" % str(self.file1_label))
        lines.append("File1 Label    : %s" % str(self.file2_label))
        lines.append("Time Step      : %s" % str(self.time_steps))
        lines.append("Duration       : %s" % str(self.duration))
        lines.append("Tip Angle      : %s" % str(self.tip_angle))
        lines.append("Bandwidth      : %s" % str(self.bandwidth))
        lines.append("Deprecated     : %s" % str(self.deprecated))
        
        controls = [ "[%s (%s): %s  var_name=%s]" % (control.name, control.type, control.default, control.variable) \
                        for control in self.transform_kernel_controls]
        lines.append("User Parameters: %s" % ", ".join(controls))
            
        code = "None"    
        if self.algorithm_code:
            code = ("%d characters" % len(self.algorithm_code))
        lines.append("Algorithm Code : %s" % code)

        # __unicode__() must return a Unicode object. In practice the code
        # above always generates Unicode, but we ensure it here.
        return '\n'.join(lines)


    def as_html(self):
        """
        Returns this kernel rendered as a string of HTML suitable for
        display in a wx.html control.
        """
        # Note that this is being built for the wx.html control which has
        # no CSS support so our formatting options here are pretty limited.
        html = ElementTree.Element("html")
        body = ElementTree.SubElement(html, "body")
        
        h1 = util_xml.TextSubElement(body, "h1", "Transform Kernel ")
        
        p = util_xml.TextSubElement(body, "p", "Name: ")
        util_xml.TextSubElement(p, "i", self.name+"  ") # Name is italicized
        util_xml.TextSubElement(p, "p", "   Menu Label: ")
        util_xml.TextSubElement(p, "i", self.menu_label)
        
        if self.id:
            p = util_xml.TextSubElement(body, "p", "UUID: ")
            util_xml.TextSubElement(p, "tt", self.id)
            # After UUID we add an "empty" paragraph to provide some space
            # between the UUID and the next element. Since completely empty
            # paragraphs tend to get ignored by HTML renderers, we add a 
            # <p> that contains a non-breaking space.
            util_xml.TextSubElement(body, "p", _NBSP)

        p = util_xml.TextSubElement(body, "p", "Global Parameters ---- ")

        # Basic project attrs go into a table.
        table = ElementTree.SubElement(body, "table", {"border" : "1px"})
        tbody = ElementTree.SubElement(table, "tbody")

        lines = [ ]
            
        lines.append( ("Type", "%s" % self.type))
        lines.append( ("Public",  "%s" % ("True" if self.is_public else "False")))
        if self.created:
            timestamp = self.created.strftime(util_time.DISPLAY_TIMESTAMP_FORMAT)
            lines.append( ("Created", "%s" % timestamp))
        lines.append( ("Creator", "%s" % self.creator))
        lines.append( ("Hide File1", "%s" % self.hide_file1))
        lines.append( ("Hide File2", "%s" % self.hide_file2))
        lines.append( ("Hide Time Steps", "%s" % self.hide_time_steps))
        lines.append( ("Hide Duration",  "%s" % self.hide_duration))
        lines.append( ("Hide Tip Angle", "%s" % self.hide_tip_angle))
        lines.append( ("Hide Bandwidth", "%s" % self.hide_bandwidth))
        lines.append( ("File1 Label", "%s" % self.file1_label))
        lines.append( ("File2 Label", "%s" % self.file2_label))
        
        for line in lines:
            description, data = line
            tr = ElementTree.SubElement(tbody, "tr")
            util_xml.TextSubElement(tr, "td", description)
            util_xml.TextSubElement(tr, "td", data, { "align" : "right" })

        p = util_xml.TextSubElement(body, "p", "User Defined Parameters ---- ")

        # User parameters go into a table.
        table2 = ElementTree.SubElement(body, "table", {"border" : "1px"})
        tbody2 = ElementTree.SubElement(table2, "tbody")

        lines = [ ]
        lines.append(("Name", "Type", "Default", "Variable"))
        for control in self.transform_kernel_controls:
            lines.append((control.name, control.type, control.default, control.variable))    
        
        for line in lines:
            name, type_, default, variable = line
            tr = ElementTree.SubElement(tbody2, "tr")
            util_xml.TextSubElement(tr, "td", name)
            util_xml.TextSubElement(tr, "td", type_)
            util_xml.TextSubElement(tr, "td", default, { "align" : "right" })
            util_xml.TextSubElement(tr, "td", variable)
            
        if self.comment:
            util_xml.TextSubElement(body, "h2", "Comment")
            util_xml.TextSubElement(body, "p", self.comment)
            
        # The two commented-out lines below are handy for debugging, should 
        # you need that. =)
        # util_xml.indent(html)
        # print  ElementTree.tostring(html)

        return ElementTree.tostring(html)
    

    def clone(self):
        """
        Creates & returns a new transform that looks just like this 
        one, but with a different id.
        """
        transform = copy.deepcopy(self)
        transform.id = util_misc.uuid()
        transform.created = util_time.now()
        transform.is_public = False
        
        return transform
        

    def deflate(self, flavor=Deflate.ETREE):
        if flavor == Deflate.ETREE:
            
            e = ElementTree.Element("transform_kernel", 
                                    {"id" : self.id, "version" : self.XML_VERSION})
                                            
            for attribute in ("type", 
                              "name", 
                              "menu_label", 
                              "created", 
                              "creator", 
                              "comment", ):
                util_xml.TextSubElement(e, attribute, getattr(self, attribute))

            # These atttributes are all scalars and map directly to 
            # XML elements of the same name.
            for attribute in ("hide_file1", 
                              "hide_file2",
                              "hide_time_steps",
                              "hide_duration",
                              "hide_tip_angle",
                              "hide_bandwidth",
                              "file1_label",
                              "file2_label",
                              "tip_angle", 
                              "time_steps", 
                              "duration", 
                              "bandwidth",
                              "deprecated"):
                util_xml.TextSubElement(e, attribute, getattr(self, attribute))

            for control in self.transform_kernel_controls:
                e.append(control.deflate(flavor))

            if self.algorithm_code:
                util_xml.TextSubElement(e, "algorithm_code", self.algorithm_code)

            return e
            
        elif flavor == Deflate.DICTIONARY:
            d = self.__dict__.copy()
            
            d["transform_kernel_controls"] = [control.deflate(flavor) for control in self.transform_kernel_controls]
            return d


    def inflate(self, source):
        if hasattr(source, "makeelement"):
            # Quacks like an ElementTree.Element
            self.id = source.get("id")
            self.type = source.findtext("type")
            self.name = source.findtext("name")
            self.menu_label = source.findtext("menu_label")
            self.creator = source.findtext("creator")
            self.created = source.findtext("created")
            self.created = util_time.datetime_from_iso(self.created)
            self.comment = source.findtext("comment")
            self.algorithm_code = source.findtext("algorithm_code")

            # Booleans
            for attribute in ("hide_file1", 
                              "hide_file2",
                              "hide_time_steps",
                              "hide_duration",
                              "hide_tip_angle",
                              "hide_bandwidth",
                              "deprecated"
                             ):
                item = source.findtext(attribute)
                if item is not None:
                    setattr(self, attribute, util_xml.BOOLEANS[item])

            # floats
            for attribute in ("tip_angle", 
                              "duration", 
                              "bandwidth",
                             ):
                item = source.findtext(attribute)
                if item is not None:
                    setattr(self, attribute, float(item))

            # ints
            for attribute in ("time_steps", ):
                item = source.findtext(attribute)
                if item is not None:
                    setattr(self, attribute, int(float(item)))

            # No translation required for these attrs
            for attribute in ("file1_label",
                              "file2_label",
                             ):
                item = source.findtext(attribute)
                if item is not None:
                    setattr(self, attribute, item)
                    
            val = [TransformKernelControl(control) for control in source.findall("transform_kernel_control")]
            self.transform_kernel_controls = val
            
        elif hasattr(source, "keys"):
            # Quacks like a dict
            for key in list(source.keys()):
                if hasattr(self, key):
                    setattr(self, key, source[key])

            val = [TransformKernelControl(control) for control in source.get("transform_kernel_control", [])]
            self.transform_kernel_controls = val


class TransformKernelControl(object):

    # The XML_VERSION enables us to change the XML output format in the future 
    XML_VERSION = "1.0.0"
    
    def __init__(self, attributes=None):
        self.name = ""
        self.type = ""
        self.default = ""
        self.variable = ""

        if attributes is not None:
            self.inflate(attributes)


    def __str__(self):
        return self.__unicode__()


    def __unicode__(self):
        lines = [ ]
        lines.append("--- Transform Kernel Control ---" )
        lines.append("Name: %s" % self.name)
        lines.append("Type: %s" % self.type)
        lines.append("Default: %s" % self.default)
        lines.append("Variable: %s" % self.variable)
    
        # __unicode__() must return a Unicode object. In practice the code
        # above always generates Unicode, but we ensure it here.
        return '\n'.join(lines)


    def deflate(self, flavor=Deflate.ETREE):
        if flavor == Deflate.ETREE:
            # e = parameter element
            e = ElementTree.Element("transform_kernel_control", {"version" : self.XML_VERSION})
                                            
            util_xml.TextSubElement(e, "name", self.name)
            util_xml.TextSubElement(e, "type", self.type)
            util_xml.TextSubElement(e, "default", self.default)
            util_xml.TextSubElement(e, "variable", self.variable)

            return e
            
        elif flavor == Deflate.DICTIONARY:
            return self.__dict__.copy()


    def inflate(self, source):
        if hasattr(source, "makeelement"):
            # Quacks like an ElementTree.Element
            self.name = source.findtext("name")
            self.type = source.findtext("type")
            self.default = source.findtext("default")
            self.variable = source.findtext("variable")

        elif hasattr(source, "keys"):
            # Quacks like a dict
            for key in list(source.keys()):
                if hasattr(self, key):
                    setattr(self, key, source[key])

