# Python modules

import xml.etree.cElementTree as ElementTree
import copy

# 3rd party modules
import numpy as np

# Our modules
import vespa.common.constants as constants
import vespa.common.rfp_transform as rfp_transform
import vespa.common.rfp_machine_specs as rfp_machine_specs
import vespa.common.util.time_ as util_time
import vespa.common.util.xml_ as util_xml
import vespa.public.minimalist_pulse as minimalist_pulse

from vespa.common.rfp_machine_specs import MachineSpecs
from vespa.common.constants import Deflate

#1234567890123456789012345678901234567890123456789012345678901234567890123456789

# _NBSP is a Unicode Non-Breaking SPace. Used in HTML output.
_NBSP = "\xA0"

class PulseDesign(object):
    """
    A container to hold the final pulse, frequency profile (and any gradients),
    as well as the provenance of the pulse, including all the prior
    transforms (and parameters) and their results.

    Attributes:
        calc_resolution (int): default 5000, num_points for rf pulse Bloch
            profile calculations.

        pulse_bandwidth_type (str): default 'half_height', this is location on
            rf profile at which the bandwidth is meaured. For now, 'half-height'
            is the only setting. In future, may add 'minimum' or 'maximum'.

        gyromagnetic_nuclei (str): default '1H'

        bloch_range_value (float): default 4.0, +/- width of x-axis displayed
            in the rf pulse plots. This can be in either kHz or cm units
            depending on the value in bloch_range_units attribute.

        bloch_range_units (str): default 'cm', unit of the value in the
            bloch_range_value attribute. Used together with the calc_resolution
            to calculate the point locations for the Bloch calculation.

        bloch_offset_value (float): default 0.0, offset frequency value used in
            the Bloch calculations.

        machine_specs (object): useful info about the MR scanner specifications
            for which this rf pulse is designed. Info includes things like B0
            field, max gradient strength, max gradient slew, max B1 field, etc.

        transforms (list): a list of rfp_transform objects that describe the
            steps used to create this rf pulse Design.

        refferers (list): a (possibly empty) list of 2-tuples of (id, name).
            This contains all of the Vespa-Simulation 'pulse sequences' that
            refer to this pulse design.


    """
    # The XML_VERSION enables us to change the XML output format in the future 
    XML_VERSION = "1.0.0"

    def __init__(self, attributes=None):
        self.name = ''
        self.id = ''                # uuid
        self.is_public = False
        self.creator = ''
        self.created = util_time.now()
        self.comment = ''  # User's note(s) for describing this pulse.

        self.calc_resolution      = constants.DEFAULT_CALC_RESOLUTION
        self.pulse_bandwidth_type = constants.DEFAULT_PULSE_BANDWIDTH_TYPE
        self.gyromagnetic_nuclei  = constants.DEFAULT_GYROMAGNETIC_NUCLEI
        self.bloch_range_value    = constants.DEFAULT_BLOCH_RANGE_VALUE
        self.bloch_range_units    = constants.DEFAULT_BLOCH_RANGE_UNITS
        self.bloch_offset_value   = constants.DEFAULT_BLOCH_OFFSET_VALUE
 
        self.machine_specs = MachineSpecs()
        self.transforms = []
        self.referrers = []

        if attributes is not None:
            self.inflate(attributes)

        if self.comment is None:
            self.comment = ""


    @property
    def create_transform(self):
        """Returns the create transform if one exists, None otherwise."""
        return self.transforms[0] if len(self.transforms) else None


    @property
    def has_results(self):
        return any([bool(transform.result) for transform in self.transforms])

    @property
    def is_frozen(self):
        """
        A pulse design is frozen when it's public or when one or more
        pulse sequences refers to it.
        """
        return bool(self.referrers) or self.is_public

    def __str__(self):
        return self.__unicode__()

    def __unicode__(self):
        lines = [ ]
        id_ = self.id if self.id else ""
        lines.append("--- Pulse Design %s ---" % id_)
        lines.append("Name: %s" % self.name)
        lines.append("Public: %s" % ("True" if self.is_public else "False"))
        s = self.created.isoformat() if self.created else "[empty]"
        lines.append("Created: %s" % s)
        lines.append("Comment (abbr.): %s" % self.comment[:40])
        lines.append("Creator: %s" % self.creator)        
        
        # Present a summary of the machine specs
        machine_type = self.machine_specs.machine_type
        if machine_type in constants.MachineType.ALL:
            machine_type = machine_type["display"]
        machine_type += (" %.3fT" % self.machine_specs.field_strength)
        
        lines.append("Global - Calc Resolution : %s" % str(self.calc_resolution)) 
        lines.append("Global - Pulse Bandwidth Type : %s" % self.pulse_bandwidth_type) 
        lines.append("Global - Gyromagnetic Nuclei : %s" % self.gyromagnetic_nuclei) 

        lines.append("Machine Settings: %s" % machine_type)
        lines.append("%d Transforms: (not listed)" % len(self.transforms))

        # __unicode__() must return a Unicode object.
        return '\n'.join(lines)
        

    def as_html(self):
        """Returns this design rendered as a string of HTML suitable for
        display in a wx.html control."""
        # Note that this is being built for the wx.html control which has
        # no CSS support so our formatting options here are pretty limited.
        html = ElementTree.Element("html")
        body = ElementTree.SubElement(html, "body")
        
        h1 = util_xml.TextSubElement(body, "h1", "Pulse Design ")
        # Name is italicized
        util_xml.TextSubElement(h1, "i", self.name)
        
        if self.id:
            p = util_xml.TextSubElement(body, "p", "UUID: ")
            util_xml.TextSubElement(p, "tt", self.id)
            # After UUID we add an "empty" paragraph to provide some space
            # between the UUID and the next element. Since completely empty
            # paragraphs tend to get ignored by HTML renderers, we add a 
            # <p> that contains a non-breaking space.
            util_xml.TextSubElement(body, "p", _NBSP)

        # Basic design attrs go into a table.
        table = ElementTree.SubElement(body, "table", {"border" : "1px"})
        tbody = ElementTree.SubElement(table, "tbody")

        lines = [ ]
            
        lines.append( ("Public",  "%s" % ("True" if self.is_public else "False")))
        if self.created:
            timestamp = self.created.strftime(util_time.DISPLAY_TIMESTAMP_FORMAT)
            lines.append( ("Created", "%s" % timestamp))
        lines.append( ("Creator", "%s" % self.creator))
        lines.append( ("Calc Resolution", "%s" % str(self.calc_resolution)))
        lines.append( ("Pulse Bandwidth Type", "%s" % str(self.pulse_bandwidth_type)))
        lines.append( ("Gyromagnetic Nuclei", "%s" % str(self.gyromagnetic_nuclei)))
        
        for line in lines:
            description, data = line
            tr = ElementTree.SubElement(tbody, "tr")
            util_xml.TextSubElement(tr, "td", description)
            util_xml.TextSubElement(tr, "td", data, { "align" : "right" })
            
        if self.comment:
            util_xml.TextSubElement(body, "h2", "Comment")
            util_xml.TextSubElement(body, "p", self.comment)

        # The machine settings object renders itself
        util_xml.TextSubElement(body, "h2", "Machine Settings")
        body.append(self.machine_specs.as_html())

        # For now, only transform titles are listed. Maybe later we'll
        # want to add for info about each transform, and if we do that
        # we should probably let the transform objects render themselves.
        util_xml.TextSubElement(body, "h2", "Transforms")
        ol = ElementTree.SubElement(body, "ol")
        for transform in self.transforms:
            txt = transform.type + " - " + transform.transform_kernel.name
            util_xml.TextSubElement(ol, "li", txt)
            
        # The two commented-out lines below are handy for debugging, should 
        # you need that. =)
        # util_xml.indent(html)
        # print  ElementTree.tostring(html)

        return ElementTree.tostring(html)


    def clone(self):
        """
        Creates & returns a new pulse_design that looks just like this one, 
        but with a different value for 'created'. NB. The id (uuid) should be
        modified externally.

        """
        pulse_design = copy.deepcopy(self)
        pulse_design.id = None
        pulse_design.created = util_time.now()
        pulse_design.is_public = False
        
        return pulse_design


    def delete_transform(self, index):
        '''
        Deletes the transform at the given position in the list.
        The index has to be >= 0 

        '''
        length = len(self.transforms)
        if index >= length or index < 0:
            return
        self.transforms.pop(index)
        return


    def get_pulse(self, index=-1):
        """
        Given an index into the list of transforms in this design,
        returns a MinimalistPulse representation of that transform. 
        If the design has no transforms, or the indicated 
        transform has no result, None is returned.
        
        The index defaults to -1, which refers to the last transform.
        It's used directly in the list of transforms, so you can pass
        e.g. -2 to get the next-to-last transform, etc. Passing an
        index >= the number of transforms raises an IndexError.

        The waveform in the object that's returned is a copy of the list 
        inside the pulse design. Changes to the list returned from here 
        won't affect the pulse design.

        """
        pulse = None
        
        if self.transforms:
            transform = self.transforms[index]
            
            if transform.result:
                pulse = minimalist_pulse.MinimalistPulse()
                pulse.name = self.name
                pulse.dwell_time  = transform.result.dwell_time
                pulse.rf_waveform = transform.result.rf_waveform[:]
                pulse.rf_xaxis    = transform.result.rf_xaxis[:]
                if transform.result.gradient is not None:
                    pulse.gradient   = transform.result.gradient[:]
                if transform.result.grad_xaxis is not None:
                    pulse.grad_xaxis = transform.result.grad_xaxis[:]

        return pulse


    def get_bandwidth(self, index=-1):
        """
        Given an index into the list of transforms in this design,
        returns the bandwidth of that transform. If there is no 
        bandwidth parameter in this transform, it searches backwards through
        the transform list until it finds one. Otherwise it returns None.
        
        If the design has no transforms, None is returned.
        
        The index defaults to -1, which refers to the last transform.
        It's used directly in the list of transforms, so you can pass
        e.g. -2 to get the next-to-last transform, etc. Passing an
        index >= the number of transforms raises an IndexError.

        """
        val = None

        ntrans = len(self.transforms)
        if abs(index) > ntrans:
            return val
        
        for indx in reversed(list(range(ntrans + index + 1))):
            val = self.transforms[indx].get_parameter('bandwidth')
            if val != '':
                break
        
        return val
                
                
    def get_tip_angle(self, index=-1):
        """ See get_bandwidth() method description """
        val = None

        ntrans = len(self.transforms)
        if abs(index) > ntrans:
            return val
        
        for indx in reversed(list(range(ntrans + index + 1))):
            val = self.transforms[indx].get_parameter('tip_angle')
            if val != '':
                break
        
        return val
            

    def get_previous_result(self, transform):
        idx = self.transforms.index(transform)
        if idx > 0:
            return self.transforms[idx-1].result
        else:
            return None


    def deflate(self, flavor=Deflate.ETREE):
        if flavor == Deflate.ETREE:
            # e = pulse design element
            e = ElementTree.Element("pulse_design", { "id" : self.id, "version" : self.XML_VERSION})

            for attribute in ("name", "created", "creator", "comment", ):
                util_xml.TextSubElement(e, attribute, getattr(self, attribute))

            util_xml.TextSubElement(e, "calc_resolution", str(self.calc_resolution))
            util_xml.TextSubElement(e, "pulse_bandwidth_type", str(self.pulse_bandwidth_type))
            util_xml.TextSubElement(e, "gyromagnetic_nuclei", str(self.gyromagnetic_nuclei))

            e.append(self.machine_specs.deflate(flavor))

            for transform in self.transforms:
                e.append(transform.deflate(flavor))

            return e

        elif flavor == constants.Deflate.DICTIONARY:
            d = self.__dict__.copy()
            
            d["machine_specs"] = self.machine_specs.deflate(flavor)
            d["transforms"] = [transform.deflate(flavor) for transform in self.transforms]
            
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

            self.calc_resolution = int(source.findtext("calc_resolution"))
            self.pulse_bandwidth_type = source.findtext("pulse_bandwidth_type")
            
            val = source.findtext("gyromagnetic_nuclei")
            if val is not None:
                self.gyromagnetic_nuclei = val
            
            # We need an explicit test for None here; casting to bool doesn't
            # work as expected. See "caution" at the end of this section:
            # http://docs.python.org/release/2.6.6/library/xml.etree.elementtree.html#the-element-interface
            machine_specs_element = source.find("machine_specs")
            if machine_specs_element is not None:
                self.machine_specs = MachineSpecs(machine_specs_element)

            items = source.findall("transform")
            if items:
                self.transforms = [rfp_transform.Transform(item) for item in items]
                
        elif hasattr(source, "keys"):
            # Quacks like a dict
            for key in list(source.keys()):
                if hasattr(self, key):
                    setattr(self, key, source[key])


    def update_results(self):
        """
        Forces an update for the results on each transform.
        
        NB. This may be deprecated since we need display values that are only
            available at the tab_transform object level.
        """

        for transform in self.transforms:
            if transform.result:
                outputs = {}
                outputs['calc_resolution']     = self.calc_resolution
                outputs['gyromagnetic_nuclei'] = self.gyromagnetic_nuclei
                outputs['bloch_range_value']   = 25.0  #self.bloch_range_value
                outputs['bloch_range_units']   = 'cm'  #self.bloch_range_units
                outputs['update_profiles']     = True
                transform.result.update_profiles(outputs)

    

# def _convert_project_to_design(db, rfpulse_project):
#
#     # some constants
#     kernel_id    = '3cb7a9e2-bb23-4934-ae17-96c01b8d93c7'     # the Import From RFPulse transform
#     type_complex = constants.DataTypes.any_type_to_internal(complex)
#     type_float   = constants.DataTypes.any_type_to_internal(float)
#
#     # start conversion
#     min_pulse = rfpulse_project.get_pulse()           # returns MinimalistPulse of last tab
#
#     if min_pulse is not None:
#         flag_add_transform = True
#
#         # encode numpy arrays to strings
#         rf_waveform    = util_xml.encode_numeric_list(min_pulse.rf_waveform, type_complex)
#         rf_xaxis       = util_xml.encode_numeric_list(min_pulse.rf_xaxis,    type_float)
#         gradient       = 'None'
#         grad_xaxis     = 'None'
#         if min_pulse.gradient is not None and min_pulse.gradient != []:
#             gradient   = util_xml.encode_numeric_list(min_pulse.gradient,    type_float)
#         if min_pulse.grad_xaxis is not None and min_pulse.grad_xaxis != []:
#             grad_xaxis = util_xml.encode_numeric_list(min_pulse.grad_xaxis,  type_float)
#     else:
#         flag_add_transform = False
#
#     # create a new PulseDesign and copy basic info
#     pulse_design               = PulseDesign()
#     pulse_design.id            = rfpulse_project.id
#     pulse_design.name          = rfpulse_project.name
#     pulse_design.created       = rfpulse_project.created
#
#     default_template           = db.fetch_default_machine_specs_template()
#     pulse_design.machine_specs = rfp_machine_specs.specs_from_template(default_template)
#     pulse_design.name          = db.find_unique_name(pulse_design, "rfpulse_convert")
#
#     if flag_add_transform:
#
#         # add an ImportFromRFPulse transform kernel
#         transform = rfp_transform.Transform()
#         transform.transform_kernel = db.fetch_transform_kernel(kernel_id)
#         transform.reset_parameters()
#         transform.reset_result()
#         pulse_design.transforms.append(transform)
#
#         # copy results from PulseProject to PulseDesign
#         transform.parameters[0].value = rf_waveform           # rf_waveform_encoded
#         transform.parameters[1].value = rf_xaxis              # rf_xaxis_encoded
#         transform.parameters[2].value = gradient              # gradient_encoded
#         transform.parameters[3].value = grad_xaxis            # grad_xaxis_encoded
#
#         transform.result.rf_waveform    = np.array(min_pulse.rf_waveform)
#         transform.result.rf_xaxis       = np.array(min_pulse.rf_xaxis)
#         if min_pulse.gradient is not None and min_pulse.gradient != []:
#             transform.result.gradient   = np.array(min_pulse.gradient)
#         if min_pulse.grad_xaxis is not None and min_pulse.grad_xaxis != []:
#             transform.result.grad_xaxis = np.array(min_pulse.grad_xaxis)
#
#     # Append a comment marking the converted object
#     comment = "Converted RFPulse Project %s (%s) to Pulse Design on %s\n" % \
#                 (rfpulse_project.id,
#                  rfpulse_project.name,
#                  util_time.now(util_time.DISPLAY_TIMESTAMP_FORMAT))
#     if rfpulse_project.comment:
#         comment += "\n" + rfpulse_project.comment
#     pulse_design.comment = comment
#
#     return pulse_design

