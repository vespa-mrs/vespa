# Python modules


# 3rd party modules


# Our modules
import vespa.analysis.block_raw as block_raw
import vespa.common.mrs_data_raw as mrs_data_raw
import vespa.common.util.xml_ as util_xml
from vespa.common.constants import Deflate



class BlockRawEdit(block_raw.BlockRaw):
    """ 
    Building block to hold the state of a step in an MRS processing chain.
    Includes the functionality to save/recall this object to/from an XML node.

    Raw Blocks hold data loaded from file. They don't have 'inputs' for a
    Chain object. They do have the attributes inheirited from DataRaw. 

    This object is a sub-class of BlockRaw->DataRaw base classes

    'Edited' pulse sequence data sets usually have a pair of 'on'/'off' data,
    and likely a pair of processed 'sum'/'diff' data. This object holds
    references to all of them. We have also added methods for get_/set_
    associated datasets based on the multiple file relationships

    """
    XML_VERSION = "1.0.0"
    
    def __init__(self, attributes=None):

        self.data_on  = None
        self.data_off = None
        self.data_sum = None
        self.data_dif = None
        self.data_sum_indiv = None
        self.data_dif_indiv = None
        self.data_on_id  = ''
        self.data_off_id = ''
        self.data_sum_id = ''
        self.data_dif_id = ''
        self.data_sum_indiv_id = ''
        self.data_dif_indiv_id = ''

        # called last because inflate called inside here and values above would reset

        super().__init__(attributes)


    ##### Standard Methods and Properties #####################################

    def __str__(self):
        lines = mrs_data_raw.DataRaw.__str__(self).split('\n')
        # Replace the heading line
        lines[0] = "------- {0} Object -------".format(self.__class__.__name__)
        lines.append("No printable data ")
        return '\n'.join(lines)


    def clear_associated_datasets(self):
        self.data_on  = None
        self.data_off = None
        self.data_sum = None
        self.data_dif = None
        self.data_sum_indiv = None      # for derived dataset completeness
        self.data_dif_indiv = None      # for derived dataset completeness
        self.data_on_id  = ''
        self.data_off_id = ''
        self.data_sum_id = ''
        self.data_dif_id = ''
        self.data_sum_indiv_id = ''      # for derived dataset completeness
        self.data_dif_indiv_id = ''      # for derived dataset completeness


    def get_associated_datasets(self, is_main_dataset=True):
        """
        Return list of datasets associated with this object
        - is_main_dataset: flag for top dataset, used to stop circular references

        Order matters here
        - ON, OFF, then SUM, DIF, then SUM_INDIV, DIF_INDIV
        - because NOT a Fidsum object
        - Note. this is different in Fidsum object

        """
        r = []
        if is_main_dataset:
            # None values dealt with in notebook_datasets.global_poll_associated_tabs
            r = [self.data_on, self.data_off, self.data_sum, self.data_dif, self.data_sum_indiv, self.data_dif_indiv]

        return r


    def set_associated_datasets(self, datasets): 
        """
        When we open a VIFF format file, main._import_file() calls this method
        to parse/store any datasets associated with this one as described below.

        Order matters here
        - ON, OFF, then SUM, DIF, then SUM_INDIV, DIF_INDIV
        - because NOT a Fidsum object
        - Note. this is different in Fidsum object

        """
        ndatasets = len(datasets)

        # pass1 - use defined IDs to set
        for dataset in datasets:
            if self.data_on_id == dataset.id:
                self.data_on = dataset
            if self.data_off_id == dataset.id:
                self.data_off = dataset
            if self.data_sum_indiv_id == dataset.id:
                self.data_sum_indiv = dataset
            if self.data_dif_indiv_id == dataset.id:
                self.data_dif_indiv = dataset
            if self.data_sum_id == dataset.id:
                self.data_sum = dataset
            if self.data_dif_id == dataset.id:
                self.data_dif = dataset

        # pass2 - define by position if necessary
        if self.data_on_id == '':
            self.data_on = datasets[0]
            self.data_on_id = datasets[0].id

        if self.data_off_id == '':
            self.data_off = datasets[1]
            self.data_off_id = datasets[1].id

        if self.data_sum_indiv_id == '' and ndatasets > 2:
            self.data_sum_indiv = datasets[2]
            self.data_sum_indiv_id = datasets[2].id

        if self.data_dif_indiv_id == '' and ndatasets > 3:
            self.data_dif_indiv = datasets[3]
            self.data_dif_indiv_id = datasets[3].id

        if self.data_sum_id == '' and ndatasets > 4:
            self.data_sum = datasets[4]
            self.data_sum_id = datasets[4].id

        if self.data_dif_id == '' and ndatasets > 5:
            self.data_dif = datasets[5]
            self.data_dif_id = datasets[5].id



    def deflate(self, flavor=Deflate.ETREE):
        
        # Make my base class do its deflate work
        e = block_raw.BlockRaw.deflate(self, flavor)

        if flavor == Deflate.ETREE:
            # Alter the tag name & XML version info   
            e.tag = "block_raw_edit"
            e.set("version", self.XML_VERSION)
            
            # Deflate the attribs specific to this sub-class
            
            if not self.behave_as_preset:
            
                # In the next few lines, we *have* to save the uuid values from 
                # the actual objects rather than from the attributes above, in 
                # order for the associated dataset uuids to reflect the new ids
                # that are given in the top level dataset. Associated datasets are
                # given new temporary uuid values so that if the main dataset is 
                # saved and immediately loaded back in, we do not get collisions
                # between the newly opened datasets and already existing ones.
                if self.data_on:
                    util_xml.TextSubElement(e, "data_on_id",  self.data_on.id)
                if self.data_off:
                    util_xml.TextSubElement(e, "data_off_id", self.data_off.id)
                if self.data_sum:
                    util_xml.TextSubElement(e, "data_sum_id", self.data_sum.id)
                if self.data_dif:
                    util_xml.TextSubElement(e, "data_dif_id", self.data_dif.id)
                if self.data_sum_indiv:
                    util_xml.TextSubElement(e, "data_sum_indiv_id", self.data_sum_indiv.id)
                if self.data_dif_indiv:
                    util_xml.TextSubElement(e, "data_dif_indiv_id", self.data_dif_indiv.id)

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
            
                self.data_on_id  = source.findtext("data_on_id")
                self.data_off_id = source.findtext("data_off_id")
                self.data_sum_id = source.findtext("data_sum_id")
                self.data_dif_id = source.findtext("data_dif_id")
                self.data_sum_indiv_id = source.findtext("data_sum_indiv_id")
                self.data_dif_indiv_id = source.findtext("data_dif_indiv_id")

            # Look for settings under the old name as well as the standard name.
            self.set = util_xml.find_settings(source, "block_raw_edit_settings")
            self.set = block_raw._Settings(self.set)

        elif hasattr(source, "keys"):
            # Quacks like a dict
            for key in list(source.keys()):
                if not self.behave_as_preset:
                    if key == "data_on_id":
                        self.data_on_id = source["data_on_id"]
                    if key == "data_off_id":
                        self.data_off_id = source["data_off_id"]
                    if key == "data_sum_id":
                        self.data_sum_id = source["data_sum_id"]
                    if key == "data_dif_id":
                        self.data_dif_id = source["data_dif_id"]
                    if key == "data_sum_indiv_id":
                        self.data_sum_indiv_id = source["data_sum_indiv_id"]
                    if key == "data_dif_indiv_id":
                        self.data_dif_indiv_id = source["data_dif_indiv_id"]

                if key == "set":
                    setattr(self, key, source[key])


        


