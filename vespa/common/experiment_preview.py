# Python modules


# 3rd party modules


# Our modules
import vespa.common.mrs_metabolite as mrs_metabolite


class ExperimentPreview(object):
    """A lightweight version of an experiment that's good for populating
       lists of experiments (as in the experiment browser dialog).
    """    
    def __init__(self, attributes=None):
        self.id = ""
        self.name = ""
        self.is_public = False
        self.comment = ""
        self.b0 = None
        self.isotope = ""
        self.pulse_sequence_id = ""
        self.pulse_sequence_name = ""
        self.metabolites = [ ]

        if attributes is not None:
            self.inflate(attributes)
            
        if self.comment is None:
            self.comment = ""


    def __str__(self):
        return self.__unicode__()


    def __unicode__(self):
        lines = [ ]
        lines.append("--- Preview of Experiment %s ---" % self.id) 
        lines.append("Name: %s" % self.name) 
        lines.append("Public: %s" % ("True" if self.is_public else "False"))
        lines.append("Comment: %s" % self.comment[:100]) 
        b0 = ("%f" % self.b0) if (self.b0 is not None) else ""
        lines.append("B0: %s" % b0) 
        lines.append("Isotope: %s" % self.isotope)
        lines.append("Pulse Seq: %s (%s)" % (self.pulse_sequence_name, 
                                             self.pulse_sequence_id))
        metabolites = ["%s (%s)" % (m.name, m.id) for m in self.metabolites]
        lines.append("Metabolites: %s" % (", ".join(metabolites)))

        # __unicode__() must return a Unicode object. In practice the code
        # above always generates Unicode, but we ensure it here.
        return '\n'.join(lines)
        
        
    def inflate(self, source):
        if hasattr(source, "makeelement"):
            # Quacks like an ElementTree.Element
            
            # Experiment previews are never deflated to XML, so there's no
            # support for inflating them from XML
            raise NotImplementedError
            
        elif hasattr(source, "keys"):
            # Quacks like a dict
            for key in list(source.keys()):
                if hasattr(self, key):
                    setattr(self, key, source[key])

            self.metabolites = [mrs_metabolite.Metabolite(metabolite) \
                            for metabolite in source.get("metabolites", [ ])]
