# Python modules


# 3rd party modules


# Our modules
import vespa.analysis.block_raw as block_raw
import vespa.common.mrs_data_raw as mrs_data_raw
import vespa.common.util.xml_ as util_xml
from vespa.common.constants import Deflate



class BlockRawProbep(block_raw.BlockRaw):
    """ 
    Building block to hold the state of a step in an MRS processing chain.
    Includes the functionality to save/recall this object to/from an XML node.

    Raw Blocks hold data loaded from file. They don't have 'inputs' for a
    Chain object. They do have the attributes inheirited from DataRaw. 

    This object is a sub-class of BlockRaw->DataRaw base classes

    We sub-class a 'Probep' version of this class so we can add an methods
    for get_/set_associated_datasets because PROBE-P Pfiles have both water
    and water suppressed data inside them.

    """
    XML_VERSION = "1.0.0"
    
    def __init__(self, attributes=None):
        super().__init__(attributes)

        self.metab_dataset = None
        self.water_dataset = None
        self.metab_dataset_id = ''
        self.water_dataset_id = ''


    ##### Standard Methods and Properties #####################################

    def __str__(self):
        lines = mrs_data_raw.DataRaw.__str__(self).split('\n')
        # Replace the heading line
        lines[0] = "------- {0} Object -------".format(self.__class__.__name__)
        lines.append("No printable data ")
        return '\n'.join(lines)


    def clear_associated_datasets(self):
        self.metab_dataset = None
        self.water_dataset = None
        self.metab_dataset_id = ''
        self.water_dataset_id = ''


    def get_associated_datasets(self, is_main_dataset=True):
        """
        Return list of datasets associated with this object
        - is_main_dataset: flag for top dataset, used to stop circular references
        
        """
        # order matters here - water then metab
        return [] if not is_main_dataset else [self.water_dataset, self.metab_dataset]


    def set_associated_datasets(self, datasets): 
        """ 'datasets' order matters here - water[0] then metab[1] """
        if self.metab_dataset_id == '':
            self.metab_dataset = datasets[1]
        else:
            for dataset in datasets:
                if self.metab_dataset_id == dataset.id:
                    self.metab_dataset = dataset

        if self.water_dataset_id == '':
            self.water_dataset = datasets[0]
        else:
            for dataset in datasets:
                if self.water_dataset_id == dataset.id:
                    self.water_dataset = dataset


    def deflate(self, flavor=Deflate.ETREE):

        e = super().deflate(flavor)
        if flavor == Deflate.ETREE:
            # Alter the tag name & XML version info   
            e.tag = "block_raw_probep"
            e.set("version", self.XML_VERSION)

            # Deflate the attribs specific to this sub-class
            if not self.behave_as_preset:
            
                if self.water_dataset:
                    util_xml.TextSubElement(e, "water_dataset_id", self.water_dataset.id)
            
            return e

        elif flavor == Deflate.DICTIONARY:
            return e


    def inflate(self, source):

        super().inflate(source)

        # Now I inflate the attribs that are specific to this class
        if hasattr(source, "makeelement"):
            # Quacks like an ElementTree.Element

            if not self.behave_as_preset:
                self.water_dataset_id = source.findtext("water_dataset_id")

            self.set = util_xml.find_settings(source, "block_raw_settings")
            self.set = block_raw._Settings(self.set)

        elif hasattr(source, "keys"):
            # Quacks like a dict
            for key in list(source.keys()):
                if not self.behave_as_preset:
                    if key == "water_dataset_id":
                        self.water_dataset_id = source["water_dataset_id"]
                    
                if key == "set":
                    setattr(self, key, source[key])


        


