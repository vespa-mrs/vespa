# Python modules


# 3rd party modules

# Our modules
import vespa.analysis.block_raw as block_raw
import vespa.common.mrs_data_raw as mrs_data_raw
import vespa.common.util.xml_ as util_xml
from vespa.common.constants import Deflate



class BlockRawCmrrSlaser(block_raw.BlockRaw):
    """ 
    Building block to hold the state of a step in an MRS processing chain.
    Includes the functionality to save/recall this object to/from an XML node.

    Raw Blocks hold data loaded from file. They don't have 'inputs' for a
    Chain object. They do have the attributes inheirited from DataRaw. 

    This object is a sub-class of BlockRaw->DataRaw base classes
    
    CMRR sLASER sequence has multiple datasets in one file. This object holds
    references to all of them. We have also added methods for get_/set_
    associated datasets based on the multiple file relationships

    """
    XML_VERSION = "1.0.0"
    
    def __init__(self, attributes=None):
        super().__init__(attributes)

        self.data_coil_combine      = None
        self.data_ecc1              = None
        self.data_water1            = None
        self.data_metab64           = None
        self.data_ecc2              = None
        self.data_water2            = None
        self.data_coil_combine_id   = ''
        self.data_ecc1_id           = ''
        self.data_water1_id         = ''
        self.data_metab64_id        = ''
        self.data_ecc2_id           = ''
        self.data_water2_id         = ''


    ##### Standard Methods and Properties #####################################

    def __str__(self):
        lines = mrs_data_raw.DataRaw.__str__(self).split('\n')
        # Replace the heading line
        lines[0] = "------- {0} Object -------".format(self.__class__.__name__)
        lines.append("No printable data ")
        return '\n'.join(lines)


    def clear_associated_datasets(self):
        self.data_coil_combine      = None
        self.data_ecc1              = None
        self.data_water1            = None
        self.data_metab64           = None
        self.data_ecc2              = None
        self.data_water2            = None
        self.data_coil_combine_id   = ''
        self.data_ecc1_id           = ''
        self.data_water1_id         = ''
        self.data_metab64_id        = ''
        self.data_ecc2_id           = ''
        self.data_water2_id         = ''


    def get_associated_datasets(self, is_main_dataset=True):
        """
        Returns a list of datasets associated with this object

        'is_main_dataset' signals that this is the top level dataset gathering 
        associated datasets, and is used to stop circular references
        
        """
        datasets = block_raw.BlockRaw.get_associated_datasets(self, is_main_dataset)

        # In dataset.deflate(), dataset.id is changed, deflated and then restored 
        # to original value to avoid id collisions if the VIFF file is loaded 
        # right back in. Here we return an unaltered list of associated datasets
        
        if is_main_dataset:
            assoc = []
            for item in [self.data_coil_combine, self.data_ecc1, self.data_water1, self.data_metab64, self.data_ecc2, self.data_water2]:
                if item is not None:
                    assoc.append(item)
            return assoc
        else:
            return []


    def set_associated_datasets(self, datasets): 
        """
        When we open a VIFF format file, main._import_file() calls this method
        to parse/store any datasets associated with this one as described below.
        
        """
        if self.data_coil_combine_id == '':
            if len(datasets) > 4:                       # 3 from twix, and maybe a MMol externally added later
                self.data_coil_combine = datasets[0]
        else:
            for dataset in datasets:
                if self.data_coil_combine_id == dataset.id:
                    self.data_coil_combine = dataset

        if self.data_ecc1_id == '':
            if len(datasets) > 4:
                self.data_ecc1 = datasets[1]
            else:    
                self.data_ecc1 = datasets[0]
        else:
            for dataset in datasets:
                if self.data_ecc1_id == dataset.id:
                    self.data_ecc1 = dataset

        if self.data_water1_id == '':
            if len(datasets) > 4:
                self.data_water1 = datasets[2]
            else:
                self.data_water1 = datasets[1]
        else:
            for dataset in datasets:
                if self.data_water1_id == dataset.id:
                    self.data_water1 = dataset

        if self.data_metab64_id == '':
            if len(datasets) > 4:
                self.data_metab64 = datasets[3]
            else:
                self.data_metab64 = datasets[2]
        else:
            for dataset in datasets:
                if self.data_metab64_id == dataset.id:
                    self.data_metab64 = dataset
                    
        if self.data_ecc2_id == '':
            if len(datasets) > 4:               # some older CMRR data did not take 2 ecc and water sets
                self.data_ecc2 = datasets[4]
        else:
            for dataset in datasets:
                if self.data_ecc2_id == dataset.id:
                    self.data_ecc2 = dataset

        if self.data_water2_id == '':
            if len(datasets) > 5:              # some older CMRR data did not take 2 ecc and water sets
                self.data_water2 = datasets[5]
        else:
            for dataset in datasets:
                if self.data_water2_id == dataset.id:
                    self.data_water2 = dataset


    def deflate(self, flavor=Deflate.ETREE):
        
        # Make my base class do its deflate work
        e = block_raw.BlockRaw.deflate(self, flavor)

        if flavor == Deflate.ETREE:
            # Alter the tag name & XML version info   
            e.tag = "block_raw_cmrr_slaser"
            e.set("version", self.XML_VERSION)
            
            # Deflate the attribs specific to this sub-class
            
            if not self.behave_as_preset:
            
                # These *have* to save object uuids rather than the 'id' attributes
                # herein so that the associated dataset uuids reflect the new ids
                # given in the top level dataset. This avoids id collisions between
                # (possibly) newly opened datasets and existing ones.

                if self.data_coil_combine is not None:
                    util_xml.TextSubElement(e, "data_coil_combine_id",  self.data_coil_combine.id)
                if self.data_ecc1 is not None:
                    util_xml.TextSubElement(e, "data_ecc1_id", self.data_ecc1.id)
                if self.data_water1 is not None:
                    util_xml.TextSubElement(e, "data_water1_id", self.data_water1.id)
                if self.data_metab64 is not None:
                    util_xml.TextSubElement(e, "data_metab64_id", self.data_metab64.id)
                if self.data_ecc2 is not None:
                    util_xml.TextSubElement(e, "data_ecc2_id", self.data_ecc2.id)
                if self.data_water2 is not None:
                    util_xml.TextSubElement(e, "data_water2_id", self.data_water2.id)

            return e

        elif flavor == Deflate.DICTIONARY:
            return e


    def inflate(self, source):

        # Make my base class do its inflate work
        block_raw.BlockRaw.inflate(self, source)

        # Now I inflate the attribs that are specific to this class
        if hasattr(source, "makeelement"):
            # Quacks like an ElementTree.Element

            if not self.behave_as_preset:
            
                self.data_coil_combine_id   = source.findtext("data_coil_combine_id")
                self.data_ecc1_id           = source.findtext("data_ecc1_id")
                self.data_water1_id         = source.findtext("data_water1_id")
                self.data_metab64_id        = source.findtext("data_metab64_id")
                self.data_ecc2_id           = source.findtext("data_ecc2_id")
                self.data_water2_id         = source.findtext("data_water2_id")
            
            self.set = block_raw._Settings(self.set)

        elif hasattr(source, "keys"):
            # Quacks like a dict
            for key in list(source.keys()):
                if not self.behave_as_preset:
                    if key == "data_coil_combine_id":
                        self.data_coil_combine_id = source["data_coil_combine_id"]
                    if key == "data_ecc1_id":
                        self.data_ecc1_id = source["data_ecc1_id"]
                    if key == "data_water1_id":
                        self.data_water1_id = source["data_water1_id"]
                    if key == "data_metab64_id":
                        self.data_metab64_id = source["data_metab64_id"]
                    if key == "data_ecc2_id":
                        self.data_ecc2_id = source["data_ecc2_id"]
                    if key == "data_water2_id":
                        self.data_water2_id = source["data_water2_id"]

                if key == "set":
                    setattr(self, key, source[key])


        


