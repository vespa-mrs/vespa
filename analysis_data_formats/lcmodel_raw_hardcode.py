"""
Routines for reading an LCModel *.RAW file. Note that this can ONLY happen
if there is already a dataset loaded from which we can 'steal' the SW, and 
DIM0 and other parameters.
"""


# Python modules
from __future__ import division
import os.path

# 3rd party modules
import numpy as np

# Our modules
import vespa.common.constants as constants
import vespa.common.util.misc as util_misc
import vespa.common.mrs_data_raw as mrs_data_raw
import vespa.common.configobj as configobj
import vespa.analysis.fileio.raw_reader as raw_reader 


# data is complex64 per Philips documentation for SDAT
NUMPY_DATA_TYPE = np.complex64




class LCModelRawHardcode(raw_reader.RawReader):
    # This inherits from raw_reader.RawReader (q.v.). The only methods you
    # need to implement are __init__() and read_raw(). You *may* want to 
    # override or supplement some of raw_reader.RawReader's other methods.

    def __init__(self):
        raw_reader.RawReader.__init__(self)
        
        self.filetype_filter = "LCMode RAW (*.RAW)|*.raw;*.RAW;"
        self.multiple = False
        
        
    def read_raw(self, filename, ignore_data=False, open_dataset=None):
        """
        Given the name of a .spar or .sdat file, returns a DataRawFidSum object
        populated with the parameters and data represented by the file pair.

        When ignore_data is True, this function only reads the parameters file
        which can be much faster than reading both params & data.
        """
        
        # must have an open dataset to borrow parameters from
        if open_dataset is None:
            raise raw_reader.OpenFileUserReadRawError("No open datasets, can not import/convert LCModel RAW file")
        
        # Read the SPAR file and extract the stuff I need.
        lines = open(filename, "rb").read()
        lines = lines.decode()

        # Now we find the start/end of the header section based on the standard
        # ">>> Begin of header <<<" and ">>> End of header <<<" marker lines

        istr = lines.find("$NMID")
        iend = lines.find("$END")

        head_only = lines[istr+5:iend]    # save this to return to user
        data_only = lines[iend+4:]

        data = util_misc.normalize_newlines(data_only)   
        data = data.replace("\n"," ")
        data = data.replace("  "," ")
        data = data.split(" ")

        data = [float(val)for val in data if util_misc.is_floatable(val) ]
        data = np.array(data)
        data = data.astype(dtype=np.float32).view(np.complex64)

        if data.shape[0] != open_dataset.raw_dims[0]:
            data = data[0:open_dataset.raw_dims[0]]

        d = { }
        
        d['headers']         = head_only
        d["data_source"]     = filename
        d["sequence_type"]   = "semi-laser"
        d["frequency"]       = open_dataset.frequency  # 123.244349
        d["sw"]              = open_dataset.sw         # 6002.44

        d["dims"]            = [1024,1,1,1]
        d["dims"][0]         = data.shape[-1]
        d["nucleus"]         = open_dataset.nucleus    # '1H'
        d["seqte"]           = open_dataset.seqte      # 35.0


        # Create a DataRaw out of the first set of data in the list.
        d["data"] = data #/ 10000.0      # empirical to allow area param to be larger
        raw = mrs_data_raw.DataRaw(d)

        return [raw,]





####################    Internal functions start here     ###############

