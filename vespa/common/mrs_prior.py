# Python modules

import collections

# 3rd party modules
import numpy as np
import xml.etree.cElementTree as ElementTree

# Our modules
import vespa.common.util.ppm as util_ppm
import vespa.common.util.generic_spectral as util_spectral
import vespa.common.mrs_prior_metabolite as mrs_prior_metabolite
import vespa.common.util.xml_ as util_xml

from vespa.common.constants import Deflate

class BasisSetItem(object):
    """Represents the basis set for one metab as calculated by 
    calculate_full_basis_set().
    """
    def __init__(self):
        self.fids = [ ]
        self.peak_ppm = [ ]
        self.all_ppms = [ ]
    

class Prior(object):
    """ 
    This is the fundamental object that represents the data being
    manipulated in the Analysis program. 
    
    """
    # The XML_VERSION enables us to change the XML output format in the future 
    XML_VERSION = "1.0.0"

    def __init__(self, attributes=None):
        """     
        Prior Specific Attributes
        -------------------------------------
        
        source should be 'file', 'experiment' or 'default', but there's
          no enforcement of this rule. Similarly, source_id is usually
          a filename or UUID but can also be "default".
        
        basis set is keyed by metab names. It holds the output from 
          the most recent call to calculate_full_basis_set(). The output
          is one BasisSetItem instance for each metab in the PPM range.
        
        """
        
        self.source      = ''       
        self.source_id   = ''
        self.comment     = ''
                       
        self.nucleus     = ''       # FIXME bjs, this may be obsolete
        self.seqte       = 0.0      # FIXME bjs, this may be obsolete
                       
        self.metabolites = collections.OrderedDict()
        
        self.basis_set = collections.OrderedDict()
        
        if attributes is not None:
            self.inflate(attributes)


    ##### Standard Methods and Properties #####################################
    @property
    def names(self):
        """Returns a list of metabolite names. This property is read only."""
        return list(self.metabolites.keys())

    @property
    def all_ppms(self):
        """Returns a list of all ppm values sorted by metabolite name.
        This property is read only."""
        return [metabolite.ppms for metabolite in list(self.metabolites.values())]

    @property
    def basis_set_names(self):
        """Returns a list of the names of the metabolites in the basis set
        (which might be an empty list). This property is read only."""
        return list(self.basis_set.keys())


    def __str__(self):

        lines = [ ]
        lines.append("=============================================")
        lines.append("Prior Object")
        lines.append("=============================================")
        lines.append("Source: "+str(self.source))
        lines.append("Source ID: "+str(self.source_id))
        lines.append("Comment: "+str(self.comment))
        lines.append("Nucleus: "+str(self.nucleus))
        lines.append("Sequence TE: %f" % self.seqte)
        lines.append("Metabolite Names: "+str(self.names))
        lines.append(" ")
        for metabolite in list(self.metabolites.values()):
            lines.append(str(metabolite))

        # __unicode__() must return a Unicode object. In practice the code
        # above always generates Unicode, but we ensure it here.
        return '\n'.join(lines)


    def deflate(self, flavor=Deflate.ETREE):
        if flavor == Deflate.ETREE:
            e = ElementTree.Element("prior", {"version" : self.XML_VERSION})
                                            
            util_xml.TextSubElement(e, "source",    self.source)
            util_xml.TextSubElement(e, "source_id", self.source_id)
            util_xml.TextSubElement(e, "comment",   self.comment)
            util_xml.TextSubElement(e, "nucleus",   self.nucleus)
            util_xml.TextSubElement(e, "seqte",     self.seqte)

            for key in list(self.metabolites.keys()):
                met = self.metabolites[key]
                e.append(met.deflate(flavor))

            return e
            
        elif flavor == Deflate.DICTIONARY:
            raise NotImplementedError


    def inflate(self, source):
        if hasattr(source, "makeelement"):
            # Quacks like an ElementTree.Element
            
            self.source    = source.findtext("source")
            self.source_id = source.findtext("source_id")
            self.comment   = source.findtext("comment")
            self.nucleus   = source.findtext("nucleus")
            self.seqte     = float(source.findtext("seqte"))
            
            self.metabolites = collections.OrderedDict()
            metabolite_elements = source.findall("prior_metabolite")
            for metabolite_element in metabolite_elements:
                met = mrs_prior_metabolite.PriorMetabolite(metabolite_element)
                self.metabolites[met.name] = met

        elif hasattr(source, "keys"):
            # Quacks like a dict
            for key in list(source.keys()):
                if hasattr(self, key):
                    setattr(self, key, source[key])
                    
            for metabolite in source['prior_metabolites']:
                met = mrs_prior_metabolite.PriorMetabolite(metabolite)
                self.metabolites[met.name] = met
                


    ##### Object Specific Methods and Properties #####################################

    def calculate_full_basis_set(self, ppm_start, ppm_end, dataset ):
        '''
        Prior contains two representations of the same basis set data: 
        
        (1) A complete set of all metabolites and all spectral lines. This is
            stored in the self.metabolite dictionary as PriorMetabolite objects
             
        (2) A set of metabolite FID arrays at the resolution of the Dataset in
            which the Prior object is contained. These basis FIDs are created 
            from a truncated range of of spectral lines that are between a 
            start and end PPM range specified by the user/program. 
            
        Finally, the voigt chain stores a subset of the above FIDs array that 
        are actually selected to be part of the optimization model.
        
        This routine accepts the start/end of the ppm range and a Dataset 
        object which provides the spectral resolution information and populates
        the basis_set attribute.
        
        This method is typically only called on init and when the ppm start/end 
        value is changed.        
        
        '''
        # This function could possibly benefit from some results caching. 
        
        
        # calculate ppm start/end if not provided
        if ppm_start is None:
            ppm_start = dataset.pts2ppm(dataset.spectral_dims[0]-1)
        if ppm_end is None:
            ppm_end = dataset.pts2ppm(0)
        
        self.basis_set = collections.OrderedDict()

        for metabolite in list(self.metabolites.values()):
            ppms   = np.array(metabolite.ppms)
            areas  = np.array(metabolite.areas)
            phases = np.array(metabolite.phases)
            # filter for ppm range constraints
            index = ((ppms > ppm_start) & (ppms < ppm_end)).nonzero()  
            if len(ppms[index]) > 0:
                # Create a basis set item for this metab
                fid, peak_ppm = util_spectral.create_fid(areas[index], 
                                                         ppms[index], 
                                                         phases[index], 
                                                         dataset, 
                                                         calc_peakppm=True)
                basis_set_item = BasisSetItem()
                basis_set_item.fid = fid
                basis_set_item.peak_ppm = peak_ppm
                basis_set_item.all_ppms = list(ppms[index])

                self.basis_set[metabolite.name] = basis_set_item




#--------------------------------------------------------------------
# test code

def _test():

    import util_import
    import vespa.common.util.export as util_export
    import vespa.common.wx_gravy.common_dialogs as common_dialogs
    
    the_object = Prior()
    class_name = the_object.__class__.__name__
    filename = "_test_output_"+class_name+".xml"
    
    util_export.export(filename, [the_object])
    
    msg = ""
    try:
        importer = util_import.PriorImporter(filename)
    except IOError:
        msg = """I can't read the file "%s".""" % filename
    except SyntaxError:
        msg = """The file "%s" isn't a valid Simulation export file.""" % filename
        
    if msg:
        common_dialogs.message(msg, style=common_dialogs.E_OK)
    else:
        # Time to rock and roll!
        priors = importer.go()  
    
    print(priors[0])
    
    tom = 10


if __name__ == '__main__':
    _test()    
    
