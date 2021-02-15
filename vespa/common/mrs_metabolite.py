# Python modules

import xml.etree.cElementTree as ElementTree
import copy

# 3rd party modules


# Our modules
from vespa.common.constants import Deflate
import vespa.common.mrs_spin as mrs_spin
import vespa.common.mrs_j_coupling as mrs_j_coupling
import vespa.common.util.time_ as util_time
import vespa.common.util.xml_ as util_xml


class Metabolite(object):
    """A class describing metabolite, one of the fundamental objects in 
    Simulation. Metabolites include spins and J couplings.
    """
    # The XML_VERSION enables us to change the XML output format in the future 
    XML_VERSION = "1.0.0"

    def __init__(self, attributes=None):
        self.id = ""
        self.name = ""
        self.is_public = False
        self.creator = ""
        self.created = None
        self.comment = ""
        self.deactivated = None
        self.spins = [ ]
        self.j_couplings = [ ]
        
        if attributes is not None:
            self.inflate(attributes)

        if self.comment is None:
            self.comment = ""
        
        if self.created is None:
            self.created = util_time.now()
        
    __doc = """True if the metab is frozen, False otherwise. Note that this
    is only accurate in the context of the metab management dialog and may
    be inaccurate elsewhere or even raise an error. Use only within the 
    metabolite management dialog and its children!
    """
    def __get_is_frozen(self):
        # self.experiment_names exists only within the metab management
        # dialogs.
        return self.is_public or bool(self.experiment_names)
    is_frozen = property(__get_is_frozen, doc=__doc)


    def __lt__(self, other):
        """Python docs say that, "The sort routines are guaranteed to use  
        __lt__() when making comparisons between two objects. So, it is easy  
        to add a standard sort order to a class by defining an __lt__() 
        method.
        ref: http://docs.python.org/howto/sorting.html#odd-and-ends
        
        In the Metabolite class, we sort by the name, case-insensitively. 

        Note that metabs are also sorted in db.py. Keep that sorting in sync
        with this sorting.
        """
        return self.name.lower() < other.name.lower()


    def __str__(self):
        return self.__unicode__()


    def __unicode__(self):
        lines = [ ]
        lines.append("--- Metabolite %s ---" % self.id)
        lines.append("Name: %s" % self.name)
        lines.append("Public: %s" % ("True" if self.is_public else "False"))
        s = self.created.isoformat() if self.created else "[empty]"
        lines.append("Created: %s" % s)
        lines.append("Creator: %s" % self.creator)
        lines.append("Comment (abbr.): %s" % self.comment[:40])
        s = self.deactivated.isoformat() if self.deactivated else "[empty]"
        lines.append("Out of service: %s" % s)

        spins = [ "(%s, %f)" % (spin.isotope, spin.chemical_shift) \
                                                    for spin in self.spins]
        s = ", ".join(spins)
        lines.append("%d spins: %s" % (len(spins), s))

        # Here I format the J Couplings as they're formatted on the 
        # metabolite editing dialog.
        # FORMAT is the way I format each value. 
        FORMAT = "%.3f  "
        # COLUMN_WIDTH is the space the largest value occupies.
        COLUMN_WIDTH = len(FORMAT % -999)

        s = "%d J couplings: " % len(self.j_couplings)
        INDENT = len(s) * " "
        row = 1
        column = 1
        # elements_per_row e.g. with 6 spins: [5, 4, 3, 2, 1]
        elements_per_row = list(range(len(self.spins) - 1, 0, -1))
        for i, j_coupling in enumerate(self.j_couplings):
            # A j coupling value can't be saved to the database as None 
            # (NULL), but a newly-created JCoupling object has a value of 
            # None so I allow for that here.
            if j_coupling.value is None:
                value = "[empty]  "
            else:
                value = FORMAT % j_coupling.value
            s += value.rjust(COLUMN_WIDTH)
                
            if (i + 1) == sum(elements_per_row[:row]):
                # end of a row
                lines.append(s)
                row += 1
                column = row + 1
                # Each row is indented one more column than the previous one.
                s = INDENT + ((COLUMN_WIDTH * " ") * (row - 1))
                

        # __unicode__() must return a Unicode object. In practice the code
        # above always generates Unicode, but we ensure it here.
        return '\n'.join(lines)

        
    def clone(self):
        """Creates & returns a new metab that looks just like this one, but
           with a different id.
        """
        metabolite = copy.deepcopy(self)
        metabolite.id = None
        metabolite.created = util_time.now()
        metabolite.is_public = False
        
        return metabolite
        
        
    def deflate(self, flavor=Deflate.ETREE):
        if flavor == Deflate.ETREE:
            # m_e = metabolite_element
            m_e = ElementTree.Element("metabolite", {"id" : self.id,
                                                     "version" : self.XML_VERSION})
                                            
            for attribute in ("name", "created", "creator", "comment",
                              "deactivated", ):
                util_xml.TextSubElement(m_e, attribute, getattr(self, attribute))

            for spin in self.spins:
                m_e.append(spin.deflate(flavor))
            
            for j_coupling in self.j_couplings:
                m_e.append(j_coupling.deflate(flavor))

            return m_e
            
        elif flavor == Deflate.DICTIONARY:
            d = self.__dict__.copy()
            d["spins"] = [spin.deflate(flavor) for spin in d["spins"]]
            d["j_couplings"] = [j_coupling.deflate(flavor) for j_coupling \
                                                        in d["j_couplings"]]                                    
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
            deactivated = source.findtext("deactivated")
            if deactivated:
                self.deactivated = util_time.datetime_from_iso(deactivated)
                
            self.spins = [mrs_spin.Spin(spin) for spin in source.findall("spin")]

            self.j_couplings = [mrs_j_coupling.JCoupling(j_coupling) \
                                for j_coupling in source.findall("j_coupling")]

        elif hasattr(source, "keys"):
            # Quacks like a dict
            for key in list(source.keys()):
                if hasattr(self, key):
                    setattr(self, key, source[key])

            self.spins = [mrs_spin.Spin(spin) \
                                        for spin in source.get("spins", [ ])]
            self.j_couplings = [mrs_j_coupling.JCoupling(j_coupling) \
                            for j_coupling in source.get("j_couplings", [ ])]
            
        # Associate j couplings with spins
        spin1 = 0
        spin2 = 1
        for j_coupling in self.j_couplings:
            j_coupling.spin1 = self.spins[spin1]
            j_coupling.spin2 = self.spins[spin2]
            
            if spin2 == (len(self.spins) - 1):
                spin1 += 1
                spin2 = spin1 + 1
            else:
                spin2 += 1
                

