"""
Routine for reading a Siemens Twix meas.dat format file created by Siemens
SVS sequence svs_edit. This is WIP#529. This sequence take a number of 
FIDs for one SVS data acquisition inside. Every other FID is phase inverted
and needs to be combined in an (a-b) + (a-b) + ... fashion.

Data is returned in a DataRaw object populated with the file's data.

"""


# Python modules
from __future__ import division
import os.path


# 3rd party modules


# Our modules
import vespa.common.constants as constants
from vespa.common.mrs_data_raw import DataRawEditFidsum
from vespa.analysis.fileio.dicom_siemens import RawReaderDicomSiemens



class RawReaderSiemensDicomFidsumMegaPress(RawReaderDicomSiemens):

    def __init__(self):
        """
        Reads multiple Siemens DICOMs file into an DataRawEditFidsum object.
        It implements the interface defined by raw_reader.RawReader (q.v.).

        """
        super().__init__()
        
        
    def read_raw(self, filename, ignore_data=False, *args, **kwargs):
        """
        Given twix filename, return a DataRawEditFidsum object with parameters
        and data. The 'Fidsum' flavor of raw data object has individual FIDs
        returned in it rather than one summed FID.
        - ignore_data option is not implemented for this reader.
        - open_dataset attribute is not used in this reader.
        
        """
        raw = super().read_raw(filename, ignore_data, kwargs['open_dataset'])

        # Change that DataRaw object into a DataRawEdit
        d = raw.deflate(constants.Deflate.DICTIONARY)

        return [DataRawEditFidsum(d),]
        

    def read_raws(self, ignore_data=False, *args, **kwargs):
        """
        Calls read_raw() once for each filename in self.filenames. It returns
        either a single DataRaw object or a list of them, depending on what
        read_raw() returns. (The implementation of read_raw() is 
        format-specific.)

        There are five possible cases --

        Cases 1,2,3,4 are listed in raw_reader.py
        
        This is an example of case 5

        5) There are multiple files that contain one/multiple FIDs which represent
        multiple datasets. This is the same as 4, just more complicated. This
        method doesn't try to handle this case; it's too format-specific. If
        you want to implement this, you'll have to write your own version of 
        this method.

        If this is being called when one or more datasets is already open,
        then one of those open datasets must be passed in the open_dataset 
        param. The attributes of the raw files that this method opens must be 
        consistent with the ones that are already open. If they're not,
        this code raises MultifileAttributeMismatchError.

        Callers should be ready to handle MultifileAttributeMismatchError,
        UnsupportedDimensionalityError, IOError, SIDataError, and 
        FileNotFoundError (for VASF and possibly other formats).
        
        In addition, DataRaw.concatenate() (q.v.) can raise ValueError
        which is not trapped by this method.
        """
        self.raws = [ ]

        for filename in self.filenames:
            try:
                # This try/except is here because we changed the API from
                # fixed attribute/keywords to variable entries and we don't
                # want to force users to rewrite their code. So, here we
                # first try the new way of calling the user method.
                raw = self.read_raw(filename, ignore_data, *args, **kwargs)
            except:
                # and this is the old way of calling the user method.
                raw = self.read_raw(filename, ignore_data)
            self.raws.append(raw)

        # Parse scan data into four arrays 
        dat_on, dat_off, dat_sum, dat_dif = _sort_data(self.raws)

        # Create a DataRawEditFidsum objects for the four different states
        # of the data ON, OFF, SUM and DIFFERENCE.

        # Change that DataRawFidsum object into a DataRawEditFidsum
        #
        # Concatenate the remaining FID data for this state into the first
        # raw object. This will result in one object with a 2D data set
        # inside it where the second dimension is the FID index. This is 
        # what Vespa expects a DataRawEditFidsum object to contain.

        raw1 = dat_on[0]
        raw1.data_sources[0] = raw1.data_sources[0]+'.EditON'
        for item in dat_on[1:]:
            raw1.concatenate(item)

        raw2 = dat_off[0]
        raw2.data_sources[0] = raw2.data_sources[0]+'.EditOFF'
        for data in dat_off[1:]:
            raw2.concatenate(data)

        raw3 = dat_sum[0]
        raw3.data_sources[0] = raw3.data_sources[0]+'.Sum'
        for data in dat_sum[1:]:
            raw3.concatenate(data)
            
        raw4 = dat_dif[0]
        raw4.data_sources[0] = raw4.data_sources[0]+'.Diff'
        for data in dat_dif[1:]:
            raw4.concatenate(data)

        self.raws = [raw1, raw2, raw3, raw4]
        
        self._check_consistency(fidsum=True)
        self._check_for_si()
        if kwargs['open_dataset']:
            self._check_compatibility(self.raws, kwargs['open_dataset'], fidsum=True)

        master = self.raws

        # TODO - bjs
        # - check scaling to get these numbers up
        # - check if ECC associate works in concert with other associations 

        return master


####################    Internal functions start here     ###############

def _sort_data(raws):
    """
    Given a list of raw data objects we process we sort the data into ON and 
    OFF lists based on the assumption that they alternate. If there are an
    odd number of raw data sets, the last one is ignored.

    When we are done parsing the ON and OFF state FIDs into arrays, we create
    SUM and DIF arrays to also send back.
    
    """
    nraw  = 2 * int(len(raws) / 2)       # get rid of odd number of items
    nraw2 = int(nraw/2)

    dat_on  = []
    dat_off = []

    for i in range(nraw2):             # this if data is arranged on/on/on/off/off/off
        dat_on.append(raws[i])
        dat_off.append(raws[i+nraw2])

    # This is an editing sequence, we process data into SUM and DIFF states

    dat_sum = []
    dat_dif = []
         
    for on,off in zip(dat_on, dat_off):
        d = on.deflate(constants.Deflate.DICTIONARY)
        add = DataRawEditFidsum(d)
        dif = DataRawEditFidsum(d)
        
        add.data = add.data + off.data
        dat_sum.append(add)
        dif.data = dif.data - off.data
        dat_dif.append(dif)

    return dat_on, dat_off, dat_sum, dat_dif


