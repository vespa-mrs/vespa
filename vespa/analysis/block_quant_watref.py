# Python modules
import os


# 3rd party modules
import numpy as np
from xml.etree.cElementTree import Element

# Our modules
import vespa.analysis.block_quant_identity as block_quant_identity
import vespa.analysis.chain_quant_watref as chain_quant_watref

import vespa.common.util.xml_ as util_xml
from vespa.common.constants import Deflate




# _CSS is the style sheet for the HTML that this code builds. Note that 
# wxPython's HTML control doesn't understand CSS. This is only for printing
# from the browser.
_CSS = """
    @page { margin: 20mm;          
            size: landscape;  
          }

    @media print {
        /* a4 = 210 x 297mm, US letter = 216 x 280mm. Two 20mm margins plus
        two 110mm block elements = 40 + 220 = 260mm which fits with room to
        spare on both paper sizes as long as the 'landscape' directive is
        respected.
        */
        div#table { width: 110mm; }
        div#img   { width: 110mm; }
    }
"""


# The 3 functions below are used for building HTML
def _format_column(column, places=4):
    """
    Given a column from the table, converts it to a nicely formatted string
    and returns it. The column can be float, int, bool or a string. Strings
    are returned untouched.

    """
    if isinstance(column, float):
        column = ("%%.%dg" % places) % column
    elif isinstance(column, int) or isinstance(column, bool):
        column = str(column)
    #else:
        # It's a string, leave it alone.

    return column


def _get_max_width(table, index):
    """Get the maximum width of the given column index"""
    return max([len(_format_column(row[index])) for row in table])    


def _pretty_space_table(table, places):
    """
    Returns a table of data, padded for alignment

    table (list of lists): The table to print. Each row must have
        the same number of columns.
    """
    col_paddings = []

    for i in range(len(table[0])):
        col_paddings.append(_get_max_width(table, i))

    lines = []
    for row in table:
        # left col
        line = row[0].center(col_paddings[0] + 2)
        # rest of the cols
        for i in range(1, len(row)):
            col = _format_column(row[i], places).center(col_paddings[i] + 2)
            line += col
        lines.append(line)
    
    return lines





class _Settings(object):
    """
    Settings object contains the parameter inputs used for processing in the 
    Chain object in this Block. Having a separate object helps to delineate 
    inputs/outputs and to simplify load/save of preset values.

    This object can also save/recall these values to/from an XML node.

    """
    # The XML_VERSION enables us to change the XML output format in the future
    XML_VERSION = "1.0.0"

    def __init__(self, attributes=None):
        """
        Most of these values appear in a GUI when the app first starts.

        Watref Processing Parameters
        -----------------------------------------

        xxx - xxxxxxxxxxxxxxxx
        xxx - xxxxxxxxxxxxxxxx
        xxx - xxxxxxxxxxxxxxxx
        xxx - xxxxxxxxxxxxxxxx
        xxx - xxxxxxxxxxxxxxxx
        
        GM, WM, CSF percentages Ernst et.al. 1993 Absolute quantitation of
        water and metabolites in the human brain. I. Compartments and water.
        J Magn Reson B 102:1-8

        """
        #------------------------------------------------------------
        # Watref Processing Variables
        #------------------------------------------------------------

        self.watref_dataset_id              = ''
        self.watref_filename                = ''
        self.sequence_te                    = 30.0
        self.water_averages                 = 1
        self.metabolite_averages            = 1
        self.apply_water_correction         = True
        self.tissue_content_gm              = 60.0
        self.tissue_content_wm              = 40.0
        self.tissue_content_csf             = 0.0
        self.water_content_gm               = 0.78
        self.water_content_wm               = 0.65
        self.water_content_csf              = 0.97
        self.water_t2_gm                    = 110.0
        self.water_t2_wm                    = 80.0
        self.water_t2_csf                   = 350.0
        self.apply_metabolite_correction    = True
        self.metabolite_t2                  = 160.0

        self.watref_dataset                 = None

        if attributes is not None:
            self.inflate(attributes)


    ##### Standard Methods and Properties #####################################

    def __str__(self):
        lines = []
        lines.append("------- {0} Object -------".format(self.__class__.__name__))
        lines.append("watref_dataset_id           : " + str(self.watref_dataset_id))
        lines.append("watref_filename             : " + str(self.watref_filename))
        lines.append("sequence_te                 : " + str(self.sequence_te))
        lines.append("water_averages              : " + str(self.water_averages))
        lines.append("metabolite_averages         : " + str(self.metabolite_averages))
        lines.append("apply_water_correction      : " + str(self.water_t2_csf))
        lines.append("tissue_content_gm           : " + str(self.tissue_content_gm))
        lines.append("tissue_content_wm           : " + str(self.tissue_content_wm))
        lines.append("tissue_content_csf          : " + str(self.tissue_content_csf))
        lines.append("water_content_gm            : " + str(self.water_content_gm))
        lines.append("water_content_wm            : " + str(self.water_content_wm))
        lines.append("water_content_csf           : " + str(self.water_content_csf))
        lines.append("water_t2_gm                 : " + str(self.water_t2_gm))
        lines.append("water_t2_wm                 : " + str(self.water_t2_wm))
        lines.append("water_t2_csf                : " + str(self.water_t2_csf))
        lines.append("apply_metabolite_correction : " + str(self.water_t2_csf))
        lines.append("metabolite_t2               : " + str(self.metabolite_t2))
        return '\n'.join(lines)


    def deflate(self, flavor=Deflate.ETREE):
        if flavor == Deflate.ETREE:
            e = Element("settings", {"version" : self.XML_VERSION})

            util_xml.TextSubElement(e, "watref_filename",             self.watref_filename)
            util_xml.TextSubElement(e, "water_averages",              self.water_averages)
            util_xml.TextSubElement(e, "sequence_te",                 self.sequence_te)
            util_xml.TextSubElement(e, "apply_water_correction",      self.apply_water_correction)
            util_xml.TextSubElement(e, "metabolite_averages",         self.metabolite_averages)
            util_xml.TextSubElement(e, "tissue_content_gm",           self.tissue_content_gm)
            util_xml.TextSubElement(e, "tissue_content_wm",           self.tissue_content_wm)
            util_xml.TextSubElement(e, "tissue_content_csf",          self.tissue_content_csf)
            util_xml.TextSubElement(e, "water_content_gm",            self.water_content_gm)
            util_xml.TextSubElement(e, "water_content_wm",            self.water_content_wm)
            util_xml.TextSubElement(e, "water_content_csf",           self.water_content_csf)
            util_xml.TextSubElement(e, "water_t2_gm",                 self.water_t2_gm)
            util_xml.TextSubElement(e, "water_t2_wm",                 self.water_t2_wm)
            util_xml.TextSubElement(e, "water_t2_csf",                self.water_t2_csf)
            util_xml.TextSubElement(e, "apply_metabolite_correction", self.apply_metabolite_correction)
            util_xml.TextSubElement(e, "metabolite_t2",               self.metabolite_t2)
            
            # In the next line, we *have* to save the uuid values from the
            # actual object rather than from the attribute above, in
            # order for the associated dataset uuid to reflect the new id
            # that is given in the top level dataset. Associated datasets are
            # given new temporary uuid values so that if the main dataset is
            # saved and immediately loaded back in, we do not get collisions
            # between the newly opened datasets and already existing ones.
            if self.watref_dataset is not None:
                util_xml.TextSubElement(e, "watref_dataset_id",        self.watref_dataset.id)

            return e

        elif flavor == Deflate.DICTIONARY:
            return self.__dict__.copy()


    def inflate(self, source):
        if hasattr(source, "makeelement"):
            # Quacks like an ElementTree.Element

            for name in ("apply_water_correction",
                         "apply_metabolite_correction" ):
                val = source.findtext(name)
                if val:
                    setattr(self, name, util_xml.BOOLEANS[val])

            for name in ("watref_filename",
                         "watref_dataset_id", ):
                val = source.findtext(name)
                if val is not None: setattr(self, name, val)

            for name in ("sequence_te",
                         "tissue_content_gm",
                         "tissue_content_wm",
                         "tissue_content_csf",
                         "water_content_gm",
                         "water_content_wm",
                         "water_content_csf",
                         "water_t2_gm",
                         "water_t2_wm",
                         "water_t2_csf",
                         "metabolite_t2",
                        ):
                val = source.findtext(name)
                if val:  setattr(self, name, float(val))

            for name in ("water_averages",
                         "metabolite_averages", 
                        ):
                val = source.findtext(name)
                if val:  setattr(self, name, int(val))

        elif hasattr(source, "keys"):
            # Quacks like a dict
            for key in list(source.keys()):
                if hasattr(self, key):
                    setattr(self, key, source[key])




class BlockQuantWatref(block_quant_identity.BlockQuantIdentity):
    """
    Building block to hold the state of a step in an MRS processing chain.
    Includes the functionality to save/recall this object to/from an XML node.

    Contains inputs/results for converting the fitted metabolite areas to 
    concentration values.

    """
    XML_VERSION = "1.0.0"

    def __init__(self, attributes=None):
        """
        Parameters:
            id (str): A permanent, unique identifying string for this object.
                Typically serves as a "source_id" for some other object. It
                is part of the provenance for this processing functor chain

            source_id (list): The unique identifier used to find the input data
                for this object. It may refer to one whole object that has only
                one result, OR it could refer to results inside an object that
                has multiple results.

            set (object): a Settings object

        """
        super().__init__(attributes)

        # processing parameters
        self.set = _Settings()

        # results storage
        self.watref_results = None

        if attributes is not None:
            self.inflate(attributes)

        self.chain = None


    ##### Standard Methods and Properties #####################################

    def __str__(self):
        lines = []
        lines.append("------- {0} Object -------".format(self.__class__.__name__))
        lines.append("\n")
        lines += _Settings.__str__(self).split('\n')
        lines.append("\n")
        lines.append("------- Main Object -------")
        lines.append("No additional data")

        return '\n'.join(lines)


    def create_chain(self, dataset):
        self.chain = chain_quant_watref.ChainQuantWatref(dataset, self)


    def check_parameter_dimensions(self, dataset):
        """
        Checks the "nparam" dimension in the results to see if the number of 
        parameters in the model has changed. Only resets results if this
        dimension has changed.
        
        """
        fit = dataset.blocks['fit']
        nparam = fit.nparam 
            
        if self.watref_results.shape[0] != nparam:
            self._reset_dimensional_data(dataset)



    def _reset_dimensional_data(self, dataset):
        """
        Resets (to zero) and resizes dimensionally-dependent data
        
        watref_results  - fit parameters from optimization
        
        """
        dims = dataset.spectral_dims
        fit  = dataset.blocks['fit']
        
        nparam = fit.nparam 
        
        if self.watref_results is None:
        
            self.watref_results = np.zeros((nparam, dims[1], dims[2], dims[3]))      
        
        else:
            param_dims = list(dims)
            param_dims[0] = nparam

            # maintain results if no dimension has changed
            if self.watref_results.shape[::-1] != param_dims:
                self.watref_results = np.zeros((nparam, dims[1], dims[2], dims[3]))      


    def get_associated_datasets(self, is_main_dataset=True):
        """
        Returns a list of datasets associated with this object

        'is_main_dataset' signals that this is the top level dataset gathering 
        associated datasets, and is used to stop circular references

        """
        datasets = block_quant_identity.BlockQuantIdentity.get_associated_datasets(self, is_main_dataset)

        if self.set.watref_dataset:
            # watref may have some ECC or CoilCombine processing that uses 
            # other datasets, so we need to return these as well as itself
            datasets += self.set.watref_dataset.get_associated_datasets(is_main_dataset=False)
            datasets += [self.set.watref_dataset,]
        else:
            return []
        
        return datasets


    def set_associated_datasets(self, datasets):
        """
        When we open a VIFF format file, main._import_file() calls this method
        to parse/store any datasets associated with this one as described below.
        
        """
        for dataset in datasets:
            if dataset.id == self.set.watref_dataset_id:
                self.set.watref_dataset = dataset


    def deflate(self, flavor=Deflate.ETREE):
        if flavor == Deflate.ETREE:
            e = Element("block_quant_watref", { "id" : self.id,
                                                "version" : self.XML_VERSION})

            util_xml.TextSubElement(e, "behave_as_preset", self.behave_as_preset)

            e.append(self.set.deflate())

            if not self.behave_as_preset:
                if self.watref_results is not None:
                    e.append(util_xml.numpy_array_to_element(self.watref_results,'watref_results'))

            return e

        elif flavor == Deflate.DICTIONARY:
            return self.__dict__.copy()


    def inflate(self, source):
        if hasattr(source, "makeelement"):
            # Quacks like an ElementTree.Element

            self.id = source.get("id")

            val = source.findtext("behave_as_preset")   # default is False
            if val is not None:
                self.behave_as_preset = util_xml.BOOLEANS[val]

            self.set = util_xml.find_settings(source, "")
            self.set = _Settings(self.set)

            if not self.behave_as_preset:

                # Explicit tests for None necessary in the code below.
                temp = source.find("watref_results")
                if temp is not None:
                    self.watref_results = util_xml.element_to_numpy_array(temp)

        elif hasattr(source, "keys"):
            # Quacks like a dict
            for key in list(source.keys()):
                if hasattr(self, key):
                    setattr(self, key, source[key])



######    Private methods    ################







#--------------------------------------------------------------------
# test code

def _test():

    pass


if __name__ == '__main__':
    _test()


