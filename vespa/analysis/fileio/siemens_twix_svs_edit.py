"""
Routine for reading a Siemens Twix meas.dat format file created by Siemens
SVS sequence svs_edit. This is WIP#529. This sequence take a number of 
FIDs for one SVS data acquisition inside. Every other FID is phase inverted
and needs to be combined in an (a-b) + (a-b) + ... fashion.

Data is returned in a DataRawFidsumEdit object populated with the file's data.

NB. The previous Twix svs_edit parser combined coil data inside the parser
    code where all four (On, Off, Sum, Diff) datasets had access to large
    signals by which to set coil weightings. In this parser, coils are combined
    after we create the Diff dataset. And often, we subtract out most of the
    signal for the DIFF dataset, so the coil combine algorithm in the Preprocess
    tab does not always work well (or at least give the same results as for the
    On, Off and Sum datasets). The WORKAROUND for this is to calculate coil
    combine weights in either On, Off or Sum, and then import them into the
    Diff dataset Preprocess tab using the 'External Dataset' method

"""

# Python modules
from __future__ import division

# 3rd party modules
import numpy as np

# Our modules
from vespa.analysis.fileio.siemens_twix import RawReaderSiemensTwix
from vespa.common.mrs_data_raw import DataRawEditFidsum



class RawReaderSiemensTwixSvsEdit(RawReaderSiemensTwix):
    """ Read a single Siemens Twix file into an DataRawEditFidsum object. """

    def read_raw(self, filename, ignore_data=False, open_dataset=None):
        """
        Given the name of a twix file, returns a DataRawEditFidsum object
        populated with the parameters and data represented by the file. 
        - ignore_data option is not implemented for this reader.
        - open_dataset attribute is not used in this reader.
        
        """
        twix, version_flag = self.get_twix(filename)
        d = self._get_parameters(twix.current, index='echo')

        data = d['data']

        if d["remove_os"]:
            data = self._remove_oversampling_basic(data)
            d["sw"] = d["sw"] / 2.0

        data = np.conjugate(data)

        dat_on  = data[0, :, :, :].copy()
        dat_off = data[1, :, :, :].copy()
        dat_sum = (dat_on + dat_off).copy()
        dat_dif = (dat_on - dat_off).copy()

        # Create a DataRawEditFidsum objects for the four different states
        # of the data ON, OFF, SUM and DIFFERENCE.
        
        d["data"] = dat_on
        d["data_source"] = filename+'.EditON'
        raw1 = DataRawEditFidsum(d)

        d["data"] = dat_off
        d["data_source"] = filename+'.EditOFF'
        raw2 = DataRawEditFidsum(d)

        d["data"] = dat_sum
        d["data_source"] = filename+'.Sum'
        raw3 = DataRawEditFidsum(d)

        d["data"] = dat_dif
        d["data_source"] = filename+'.Diff'
        raw4 = DataRawEditFidsum(d)

        return [raw1, raw2, raw3, raw4]






    
    
