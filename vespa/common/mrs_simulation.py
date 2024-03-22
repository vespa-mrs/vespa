# Python modules

import xml.etree.cElementTree as ElementTree

# 3rd party modules
import numpy as np

# Our modules
from vespa.common.constants import Deflate
import vespa.common.mrs_metabolite as mrs_metabolite
import vespa.common.constants as constants
import vespa.common.util.xml_ as util_xml
import vespa.common.util.time_ as util_time


class Simulation(object):
    """Represents one point in an experiment's results space."""

    # The XML_VERSION for this object is in the Experiment class.

    def __init__(self, attributes=None):
        # Started & completed hold timestamps. If they're None, the 
        # Simulation hasn't been run.
        self.started = None
        self.completed = None
        self.metabolite = None
        self.dims = [ ]
        self.ppms = np.array([ ], float)
        self.areas = np.array([ ], float)
        self.phases = np.array([ ], float)
        
        if attributes:
            self.inflate(attributes)


    def __lt__(self, other):
        """Python docs say, "The sort routines are guaranteed to use __lt__() 
        when making comparisons between two objects. So, it is easy to add a 
        standard sort order to a class by defining an __lt__() method."
        ref: http://wiki.python.org/moin/HowTo/Sorting/
        
        In the Simulation class, we typically sort by the dims attribute. 
        However, the first value in dims is a metabolite object, and we only
        want to sort on the name of the metabolite. Subsequently we want to 
        sort on however many other values there are in dims. thus we get the
        code below that compare two lists composed of the metabolite name and 
        other dims values using the less than operator.
        """
        # The metab ordering needs to be kept in sync with that of the
        # database function fetch_metabolites(). They're currently both
        # case-insensitive.
        name = self.metabolite.name if self.metabolite else ""
        me = [name.lower()] + self.dims
        me.reverse()
        
        name = other.metabolite.name if other.metabolite else ""
        other = [name.lower()] + other.dims
        other.reverse()
        
        return me < other


    def __str__(self):
        return self.__unicode__()


    def __unicode__(self):
        lines = [ ]
        lines.append("--- Simulation ---")
        lines.append("Started: %s" % self.started)
        lines.append("Completed: %s" % self.completed)

        if self.metabolite:        
            lines.append("metabolite: %s (%s)" % (self.metabolite.name,
                                                  self.metabolite.id))
        else:
            lines.append("metabolite: None")
            
        lines.append("dims: %s" % str(self.dims))
        
        lines.append("PPMs: " + str(self.ppms))
        lines.append("Areas: " + str(self.areas))
        lines.append("Phases: " + str(self.phases))
        
        # __unicode__() must return a Unicode object. In practice the code
        # above always generates Unicode, but we ensure it here.
        return '\n'.join(lines)


    def add_simulation(self, new):
        """ appends all lines into this Simulation instance"""
        if isinstance(new,Simulation):
            self.ppms = np.concatenate((self.ppms,new.ppms))
            self.areas = np.concatenate((self.areas,new.areas))
            self.phases = np.concatenate((self.phases,new.phases))

    def subtract_simulation(self, new):
        """ 
        Appends all lines into this Simulation with 180 degree phase change 
        
        """
        if isinstance(new,Simulation):
            self.ppms = np.concatenate((self.ppms,new.ppms))
            self.areas = np.concatenate((self.areas,new.areas))
            self.phases = np.concatenate((self.phases,new.phases+180.0))


    def deflate(self, flavor=Deflate.ETREE):
        if flavor == Deflate.ETREE:
            simulation_element = ElementTree.Element("simulation")

            if self.started:
                util_xml.TextSubElement(simulation_element, "started", 
                                        self.started)
            #else:
                # I don't even bother to create an element.
            
            if self.completed:
                util_xml.TextSubElement(simulation_element, "completed", 
                                        self.completed)
            #else:
                # I don't even bother to create an element.
            
            if self.metabolite:
                util_xml.TextSubElement(simulation_element, "metabolite_id", 
                                        self.metabolite.id)
            else:
                # Create an empty metab element
                ElementTree.SubElement(simulation_element, "metabolite_id")
                                    
            for dim in self.dims:
                util_xml.TextSubElement(simulation_element, "dim", dim)
                
            if self.ppms.size:
                ppms = util_xml.numpy_array_to_element(self.ppms, "ppms")
                simulation_element.append(ppms)
            else:
                # Create an empty element
                ElementTree.SubElement(simulation_element, "ppms")
                
            if self.areas.size:
                areas = util_xml.numpy_array_to_element(self.areas, "areas")
                simulation_element.append(areas)
            else:
                # Create an empty element
                ElementTree.SubElement(simulation_element, "areas")
            
            if self.phases.size:
                phases = util_xml.numpy_array_to_element(self.phases, "phases")
                simulation_element.append(phases)
            else:
                # Create an empty element
                ElementTree.SubElement(simulation_element, "phases")

            return simulation_element
            
        elif flavor == Deflate.DICTIONARY:
            d = self.__dict__.copy()
            
            return d

            
    def inflate(self, source, xml_version="", metabolites={ }):
        """Populates the object from source.
        
        When inflated from XML (an ElementTree), a simulation is a little
        different than other objects that have an inflate() method. In XML,
        simulations contain references to metabolites (via the id) rather
        than a description of the metabolite itself. That's because including
        the full metabolite definition would make for a huge XML file 
        containing a lot of not-very-useful repetition.
        
        To make it easier to build Simulation instances, inflate() accepts
        a dictionary of metabolite objects keyed by id. If the metabolite 
        to which this simulation refers is in the dictionary, the 
        Simulation object resolves the reference and populates its
        .metabolite attribute.
        
        In addition, when inflating from XML, the xml_version parameter must
        be populated. When inflating from a dict, it's ignored.
        """
        if hasattr(source, "makeelement"):
            # Quacks like an ElementTree.Element
            
            assert(xml_version)
            
            self.dims = [ ]
            
            started_element = source.find("started")
            # Be careful here; bool() of an ElementTree element returns
            # False if the element has no kids, even if it contains text.
            if started_element is not None:
                self.started = util_time.datetime_from_iso(started_element.text)
            completed_element = source.find("completed")
            if completed_element is not None:
                self.completed = util_time.datetime_from_iso(completed_element.text)

            # I don't expect to find a metabolite element here (see function's
            # docstring), but I check for one anyway.
            metabolite_element = source.find("metabolite")
            if metabolite_element is not None:
                self.metabolite = mrs_metabolite.Metabolite(metabolite_element)
            else:
                # This is the more likely case.
                metabolite_id = source.findtext("metabolite_id")
                self.metabolite = metabolites.get(metabolite_id, None)
            
            self.dims = [ ]
            for dim_element in source.findall("dim"):
                self.dims.append(float(dim_element.text))
                
            # Ensure that self.dims is the correct length.
            correct_dims_length = constants.RESULTS_SPACE_DIMENSIONS - 1 
            if len(self.dims) < correct_dims_length:
                self.dims += [0.0] * (correct_dims_length - len(self.dims))
                
            if xml_version == "1.0.0":
                # In XML version 1.0.0, we wrote out each spectrum and the
                # lines therein as separate elements.
                # s_e = spectrum_element
                s_e = source.find("spectrum")
                if s_e is not None:
                    self.ppms = [float(ppm.text) for ppm 
                                                 in s_e.findall("ppm")]
                    self.areas = [float(area.text) for area 
                                                   in s_e.findall("area")]
                    self.phases = [float(phase.text) for phase 
                                                     in s_e.findall("phase")]
            else:
                # In XML versions > 1.0.0, we dispense with the spectrum
                # element and write ppms, areas and phases as base64 encoded
                # binary arrays.
                ppms = source.find("ppms")
                if ppms.text:
                    self.ppms = util_xml.element_to_numpy_array(ppms)
                areas = source.find("areas")
                if areas.text:
                    self.areas = util_xml.element_to_numpy_array(areas)
                phases = source.find("phases")
                if phases.text:
                    self.phases = util_xml.element_to_numpy_array(phases)

        elif hasattr(source, "keys"):
            # Quacks like a dict
            for key in list(source.keys()):
                if hasattr(self, key):
                    setattr(self, key, source[key])

        # Before returning to the caller, I ensure that my attrs are of the
        # correct type. The code below might look tricky, but it's not. 
        # For each of ppm, areas and phases it does two things --
        #    1. If it's None, set it to [ ]
        #    2. It it's a tuple or a list, turn it into a numpy array.
        # A value of None is unusual here but possible (in the case where a
        # spectrum contains no lines).
        for name in ("ppms", "areas", "phases"):
            value = getattr(self, name)
            if value is None:
                value = [ ]
                setattr(self, name, value)                
                
            if isinstance(value, (list, tuple)):
                setattr(self, name, np.array(value))


    def summary(self):
        lines = [ ]
        dims = tuple([self.metabolite.name] + self.dims)
        spectrum = list(zip(self.ppms, self.areas, self.phases))
        for i, (ppm, area, phase) in enumerate(spectrum):
            line = "%s\t%s\t%s\t%s\t" % dims
            line += "%d\t%s\t%s\t%s" % (i, ppm, area, phase)
            lines.append(line)
        
        return '\n'.join(lines)                
