# Python modules

# 3rd party modules
import xml.etree.cElementTree as ElementTree

# Our modules
from vespa.analysis.mrs_metinfo import MetInfo
from vespa.common.constants import Deflate
from vespa.common.util.xml_ import TextSubElement
from vespa.common.mrs_generic_basis import GenericBasis



class UserPrior(object):
    """
    Contains settings for automated B0 and Phase correction methods, including
    a user defined basis object from which an ideal spectrum can be created.
    It also contains 'literature' information about various metabolites in the
    'metinfo' attribute.

    Versions:
      1.0.0 had self.spectrum = UserPriorSpectrum which was a thin skin over GenericBasis
      1.1.0 has self.basis = GenericSpectrum object for more clarity/simplicity

    """
    XML_VERSION = "1.1.0"
    
    def __init__(self, attributes=None):

        self.auto_b0_range_start      = 1.7
        self.auto_b0_range_end        = 3.4
        self.auto_phase0_range_start  = 1.85
        self.auto_phase0_range_end    = 2.25
        self.auto_phase1_range_start  = 2.75
        self.auto_phase1_range_end    = 3.60
        self.auto_phase1_pivot        = 2.01

        # Basic information input objects
        self.metinfo = MetInfo()
        self.basis = GenericBasis()

        if attributes is not None:
            self.inflate(attributes)


    def __str__(self):
        return self.__unicode__()


    def __unicode__(self):
        lines = [ ]

        lines.append("-----------  User Prior  -----------")
        lines.append("auto_b0_range_start     = %f" % self.auto_b0_range_start)
        lines.append("auto_b0_range_end       = %f" % self.auto_b0_range_end)
        lines.append("auto_phase0_range_start = %f" % self.auto_phase0_range_start)
        lines.append("auto_phase0_range_end   = %f" % self.auto_phase0_range_end)
        lines.append("auto_phase1_range_start = %f" % self.auto_phase1_range_start)
        lines.append("auto_phase1_range_end   = %f" % self.auto_phase1_range_end)
        lines.append("auto_phase1_pivot       = %f" % self.auto_phase1_pivot)
        
        metinfo = str(self.metinfo).split('\n')
        lines += ['   ' + line for line in metinfo]

        basis = str(self.basis).split('\n')
        lines += ['   ' + line for line in basis]

        return '\n'.join(lines)

    def deflate(self, flavor=Deflate.ETREE):
        if flavor == Deflate.ETREE:
            e = ElementTree.Element("user_prior",  {"version" : self.XML_VERSION})

            # These atttributes are all scalars and map directly to 
            # XML elements of the same name.
            for attribute in ("auto_b0_range_start", "auto_b0_range_end",
                    "auto_phase0_range_start", "auto_phase0_range_end",
                    "auto_phase1_range_start", "auto_phase1_range_end",
                    "auto_phase1_pivot", 
                    ):
                TextSubElement(e, attribute, getattr(self, attribute))

            e.append(self.metinfo.deflate())
            e.append(self.basis.deflate(flavor, tag="user_prior_basis", version=self.XML_VERSION))
            
            return e
            
        elif flavor == Deflate.DICTIONARY:
            return self.__dict__.copy()


    def inflate(self, source):
        if hasattr(source, "makeelement"):
            # Quacks like an ElementTree.Element
            
            # floats
            for attribute in ("auto_b0_range_start",
                    "auto_b0_range_end", "auto_phase0_range_start",
                    "auto_phase0_range_end", "auto_phase1_range_start",
                    "auto_phase1_range_end", "auto_phase1_pivot",
                             ):
                setattr(self, attribute, float(source.findtext(attribute)))

            self.metinfo = MetInfo(source.find("metinfo"))

            if self.XML_VERSION == '1.0.0':
                e = source.find("user_prior_spectrum")
            elif self.XML_VERSION == '1.1.0':
                e = source.find("user_prior_basis")
                if e is None:
                    e = source.find("user_prior_spectrum")      # just in case
            self.basis = GenericBasis(e)
            

        elif hasattr(source, "keys"):
            # Quacks like a dict
            for key in list(source.keys()):
                if hasattr(self, key):
                    setattr(self, key, source[key])


