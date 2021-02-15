# Python modules

import xml.etree.cElementTree as ElementTree
import copy

# 3rd party modules


# Our modules
from vespa.common.constants import Deflate
import vespa.common.mrs_metabolite as mrs_metabolite
import vespa.common.mrs_pulse_sequence as mrs_pulse_sequence
import vespa.common.mrs_simulation as mrs_simulation
import vespa.common.constants as constants
import vespa.common.util.time_ as util_time
import vespa.common.util.xml_ as util_xml


def deconstruct_loop(dim):
    """A convenience function. Given a loop dimension, returns a tuple
    of (start, step, length). The opposite of expand_loop()."""
    length = len(dim)
    
    start = dim[0] if length else None
    
    step = (dim[1] - dim[0]) if (length > 1) else 0.0

    return start, step, length
    

def expand_loop(start, step, length):
    """A convenience function. Given a loop's start, step size, and length, 
    returns the loop values.  The opposite of deconstruct_loop()."""
    return [start + (step * i) for i in range(length)]
    
DEFAULT_LOOP = expand_loop(0.0, 0.0, 1)     # start/step need to be float or we get an int in our dims array

class Experiment(object):

    """This is the primary object in the Vespa-Simulation application. 
    It contains a description of the pulse sequence, sequence parameters
    and metabolites that were simulated, as well as the results from 
    each simulation. Contained in here is the description of the object
    and the ability to clone, inflate and deflate the object. 
    """
    # See http://scion.duhs.duke.edu/vespa/project/wiki/ViffVersionHistory
    # for an explanation of the different XML versions this code supports.
    XML_VERSION = "1.1.0"

    SIMULATIONS_XML_VERSION = "1.1.0"

    def __init__(self, attributes=None):
    
        """In here the default values for a Experiment object are set.
        Most of these values appear in the Simulate tab when an Experiment 
        is instatiated and loaded into an Experiment Tab
        
        
        Experiment Parameters
        -----------------------------------------
        
        id           - string, uuid for the experiment simulation
        name         - string, unique short description 
        comment      - string, longer winded description
        is_public    - boolean, flag indicating whether or not the experiment
                       has been exported.
        created      - DateTime object indicating when the object was created
        investigator - string, person running simulation
        b0           - float, MR main field strength in MHz
        isotope      - string, description of isotopes to be observed in 
                       the simulation, e.g. "1H" for protons
        
        peak_search_ppm_low   - float, binning search range
        peak_search_ppm_high  - float, binning search range
        blend_tolerance_ppm   - float, frequency bin width in [ppm]
        blend_tolerance_phase - float, phase bin width in [degrees phase]
        
        pulse_sequence - object, description of simulated pulse sequence
        user_static_parameters - list of the values of user static parameters 
                     (represented as strings) for the pulse sequence
                     associated with this experiment. If the pulse seq doesn't
                     define any parameters, this list is empty.
        
        simulations    - list, of Simulation objects containing results 
                         for all metabolites and timings run 
                         
        dims - An N element list of lists where N == RESULTS_SPACE_DIMENSIONS.
               In a completed experiment, dims is a subset or view of 
               simulations. 
               The first element of dims is a list of metabolites and the
               remaining elements are lists of loop timings.
        """

        self.id = 0
        self.name = ""
        self.is_public = False
        self.created = None
        self.comment = ""
        self.investigator = ""
        self.b0 = constants.DEFAULT_B0
        self.isotope = constants.DEFAULT_ISOTOPE
        self.peak_search_ppm_low = 0.0
        self.peak_search_ppm_high = 10.0
        self.blend_tolerance_ppm = 0.0015
        self.blend_tolerance_phase = 50.0
        self.pulse_sequence = None        
        self.user_static_parameters = [ ]
        self.metabolites = [ ]        
        self.dims = [ [ ] for i in range(constants.RESULTS_SPACE_DIMENSIONS - 1)]

        # Each element in simulations is a Simulation object.
        self.simulations = [ ]
        
        if attributes is not None:
            self.inflate(attributes)
            
        if self.comment is None:
            self.comment = ""

        if self.created is None:
            self.created = util_time.now()
            

    def __str__(self):
        return self.__unicode__()


    def __unicode__(self):
        lines = [ ]
        lines.append("--- Experiment %s ---" % self.id) 
        lines.append("Name: %s" % self.name) 
        lines.append("Public: %s" % ("True" if self.is_public else "False"))
        s = self.created.isoformat() if self.created else "[empty]"
        lines.append("Created: %s" % s)        
        lines.append("Comment (abbr.): %s" % self.comment[:40]) 
        lines.append("PI: %s" % self.investigator)
        s = ("%f" % self.b0) if (self.b0 is not None) else "[empty]"
        lines.append("b0: %s" % s)
        
        if self.peak_search_ppm_low is None:
            s = "[empty]"
        else:
            s = "%f" % self.peak_search_ppm_low
        s += " / "
        if self.peak_search_ppm_high is None:
            s += "[empty]"
        else:
            s += "%f" % self.peak_search_ppm_high

        lines.append("Peak Search PPM low/high: %s" % s) 
        
        if self.blend_tolerance_ppm is None:
            s = "[empty]"
        else:
            s = "%f" % self.blend_tolerance_ppm
        s += " / "
        if self.blend_tolerance_phase is None:
            s += "[empty]"
        else:
            s += "%f" % self.blend_tolerance_phase

        lines.append("Blend tol. PPM/phase: %s" % s) 
        
        if self.pulse_sequence:
            s = "%s (%s)" % (self.pulse_sequence.id, self.pulse_sequence.name)
        else:
            s = "[empty]"
        lines.append("Pulse seq.: %s" % s)
        
        s = ", ".join(self.user_static_parameters)
        lines.append("%d User static parameters: %s" % \
                                        (len(self.user_static_parameters), s))

        metabolites = [metabolite.name for metabolite in self.metabolites]
        s = ", ".join(metabolites)
        lines.append("%d Metabolites: %s" % (len(metabolites), s))

        lines.append("%d Simulations: (not shown)" % len(self.simulations))
        
        # __unicode__() must return a Unicode object. In practice the code
        # above always generates Unicode, but we ensure it here.
        return '\n'.join(lines)

        
    @property
    def all_indices(self):
        """Generator that returns a tuple for each unique numeric index in
        the results space e.g. (1.0, 2.0, .04)
        """
        for dim3 in self.dims[2]:
            for dim2 in self.dims[1]:
                for dim1 in self.dims[0]:
                    yield (dim1, dim2, dim3)

    @property
    def n_indices(self):
        """ Returns the number of for loop calls based on elements in dims """
        return len(self.dims[2]) * len(self.dims[1]) * len(self.dims[0])


    def clone(self):
        """Creates & returns a new experiment that looks just like this one, 
        but with a different id.
        """
        experiment = copy.deepcopy(self)
        experiment.id = None
        experiment.created = util_time.now()
        experiment.is_public = False
        
        return experiment


    def deflate(self, flavor=Deflate.ETREE):
        if flavor == Deflate.ETREE:
            # e_e = experiment_element
            e_e = ElementTree.Element("experiment", { "id" : self.id,
                                                      "version" : self.XML_VERSION})
            for attribute in ("name", "created", "investigator", "comment", ):
                util_xml.TextSubElement(e_e, attribute, getattr(self, attribute))

            # b0 gets a more explicit name in the export
            util_xml.TextSubElement(e_e, "main_field_frequency", self.b0)

            for attribute in ("isotope",
                              "peak_search_ppm_low", "peak_search_ppm_high",
                              "blend_tolerance_ppm", "blend_tolerance_phase", ):
                util_xml.TextSubElement(e_e, attribute, getattr(self, attribute))

            for parameter in self.user_static_parameters:
                e = ElementTree.SubElement(e_e, "user_static_parameter_value")
                if parameter:
                    e.text = parameter
                
            e_e.append(self.pulse_sequence.deflate(flavor))
            
            for metabolite in self.metabolites:
                e_e.append(metabolite.deflate(flavor))
                
            if self.simulations:
                # s_e = Simulations Element
                attributes = { "version" : self.SIMULATIONS_XML_VERSION }
                s_e = ElementTree.SubElement(e_e, "simulations", attributes)
                
                for simulation in self.simulations:
                    s_e.append(simulation.deflate(flavor))
            else:
                # Simulations are what describe the dims in XML, so we have
                # to make fake ones.
                # FIXME PS - this is ugly.
                for dim3 in experiment.dims[2]:
                    for dim2 in experiment.dims[1]:
                        for dim1 in experiment.dims[0]:
                            for metabolite in metabolites:
                                dims = [dim1, dim2, dim3]
                                d = {"metabolite" : metabolite, 
                                     "dims" : dims}
                                simulation = mrs_simulation.Simulation(d)
                                e_e.append(simulation.deflate(flavor))
                    
            return e_e
        elif flavor == Deflate.DICTIONARY:
            d = { }
            attributes = ("id", "name", "created", "investigator", "comment",
                          "isotope", "b0", 
                          "peak_search_ppm_low", "peak_search_ppm_high", 
                          "blend_tolerance_ppm", "blend_tolerance_phase",
                         )
                         
            for attribute in attributes:
                d[attribute] = getattr(self, attribute)
            
            d["pulse_sequence"] = None
            if self.pulse_sequence:
                d["pulse_sequence"] = self.pulse_sequence.deflate(flavor)
            
            d["user_static_parameters"] = copy.copy(self.user_static_parameters)

            d["metabolites"] = [metabolite.deflate(flavor) for metabolite 
                                                           in self.metabolites]

            d["simulations"] = [simulation.deflate(flavor) for simulation
                                                           in self.simulations]
            
            return d            
            

    def inflate(self, source):
        if hasattr(source, "makeelement"):
            # Quacks like an ElementTree.Element
            self.id = source.get("id")
            # Vespa < 0.1.6 didn't write the version attribute on the 
            # experiment element, so some XML files in the wild might not 
            # have a version attr. In that case, version 1.0.0 is implied.
            xml_version = source.get("version", "1.0.0")
            self.name = source.findtext("name")
            self.investigator = source.findtext("investigator")
            self.created = source.findtext("created")
            self.created = util_time.datetime_from_iso(self.created)
            self.comment = source.findtext("comment")
            self.b0 = float(source.findtext("main_field_frequency"))
            self.isotope = source.findtext("isotope")
            self.peak_search_ppm_low = \
                                float(source.findtext("peak_search_ppm_low"))
            self.peak_search_ppm_high = \
                                float(source.findtext("peak_search_ppm_high"))
            self.blend_tolerance_ppm = \
                                float(source.findtext("blend_tolerance_ppm"))
            self.blend_tolerance_phase = \
                                float(source.findtext("blend_tolerance_phase"))
            
            pulse_sequence_element = source.find("pulse_sequence")
            # Explicit test for None necessary here. See:
            # http://scion.duhs.duke.edu/vespa/project/ticket/35
            if pulse_sequence_element is not None:
                self.pulse_sequence = \
                    mrs_pulse_sequence.PulseSequence(pulse_sequence_element)

            parameter_elements = source.findall("user_static_parameter_value")
            self.user_static_parameters = [parameter_element.text 
                                                        for parameter_element
                                                        in parameter_elements]
            # Ensure parameters are all strings
            f = lambda the_param: "" if the_param is None else the_param
            self.user_static_parameters = [f(parameter) 
                                            for parameter 
                                            in self.user_static_parameters]
                    
            # I need to build a dict of metabs for the benefit of the
            # code that builds the simulations.
            metabolite_elements = source.findall("metabolite")
            metabolites = { }
            for metabolite_element in metabolite_elements:
                metabolite = mrs_metabolite.Metabolite(metabolite_element)
                metabolites[metabolite.id] = metabolite
                
            self.metabolites = list(metabolites.values())

            dims = [ ]
            if xml_version == "1.0.0":
                # In experiment XML version 1.0.0, simulations were direct 
                # children of the experiment element. The simulation XML
                # version was also 1.0.0.
                simulation_elements = source.findall("simulation")
                simulation_xml_version = "1.0.0"
            else:
                # In experiment XML versions > 1.0.0, there's a <simulations> 
                # element that's the child of the <experiment> element, and 
                # all individual <simulation> elements live in that. The 
                # version number for all simulations hangs off of the 
                # simulations element.
                e = source.find("simulations")
                simulation_xml_version = e.get("version")
                simulation_elements = e.findall("simulation")
                
            for simulation_element in simulation_elements:
                simulation = mrs_simulation.Simulation()
                simulation.inflate(simulation_element, simulation_xml_version, 
                                   metabolites)
                self.simulations.append(simulation)
                if simulation.dims not in dims:
                    dims.append(simulation.dims)
                
            # Here we use the nifty Python idiom zip(*dims) for matrix 
            # transposition as demonstrated in the example below:
            # >>> dims = [ (1, 2, 3), (1, 2, 3), (1, 2, 3), (1, 2, 3), ]
            # >>> print zip(*dims)
            # [(1, 1, 1, 1), (2, 2, 2, 2), (3, 3, 3, 3)]
            dims = list(zip(*dims))
            # Filter out duplicates & sort what's left
            self.dims = [sorted(list(set(dim))) for dim in dims]            

        elif hasattr(source, "keys"):
            # Quacks like a dict
            for key in list(source.keys()):
                if hasattr(self, key):
                    setattr(self, key, source[key])
                    
            if "metabolites" in source:
                self.metabolites = [mrs_metabolite.Metabolite(metabolite) \
                                                    for metabolite 
                                                    in source["metabolites"]]

            if "user_static_parameters" in source:
                self.user_static_parameters = source["user_static_parameters"]

            if "pulse_sequence" in source:
                self.pulse_sequence = \
                    mrs_pulse_sequence.PulseSequence(source["pulse_sequence"])
            

class Loop(object):
    """Metadata about a loop"""
    def __init__(self, dim=None):
        self.start = 0
        self.step = 0
        self.length = 0
        
        if dim:
            self.start, self.step, self.length = deconstruct_loop(dim)
            

    def __str__(self):
        return self.__unicode__()


    def __unicode__(self):
        # __unicode__() must return a Unicode object. 
        return "Start = %f, step = %f, length = %d" % \
                                    (self.start, self.step, self.length)
        
        
        
