"""
Reads a Siemens Twix meas.dat format file created by Siemens SVS sequence
svs_se. This is standared PRESS. A DataRawFidsum object is returned and may
have mult-averages, and/or multi-coil dimensions depeding on acquisition.

This should work for VB - VE software versions.

"""
# Python modules

# 3rd party modules
import numpy as np

# Our modules
from vespa.common.mrs_data_raw import DataRawFidsum
from vespa.analysis.fileio.siemens_twix import RawReaderSiemensTwix

# RAWDATA_SCALE is based on ICE_RAWDATA_SCALE in SpecRoFtFunctor.cpp
RAWDATA_SCALE = 131072.0 * 256.0


class RawReaderSiemensTwixSvsSe(RawReaderSiemensTwix):
    """
    Read a single Siemens Twix file from svs_se sequence into an DataRawFidsum
    object. This should work for VB - VE software versions.

    """
    def read_raw(self, filename, ignore_data=False, open_dataset=None):
        """
        Given the name of a twix file, returns a DataRawFidsum object
        populated with the parameters and data represented by the file.
        - ignore_data option is not implemented for this reader.
        - open_dataset attribute is not used in this reader.

        """
        twix, version_flag = self.get_twix(filename)
        d = self._get_parameters(twix, version_flag)

        navg = d["lAverages"]
        
        data = d['data']
        if d["remove_os"]:
            data = self._remove_oversampling_basic(data)
            d["sw"] = d["sw"] / 2.0

        data = np.conjugate(data)

        # there may be multiple repetitions of svs_se in this data
        nrep = data.shape[0]
        raws = []
        for i in range(nrep):
            d["data"] = data[i,:,:,:] * RAWDATA_SCALE / float(navg)
            d["data_source"] = filename+'_rep'+str(i).zfill(2) if nrep > 1 else filename
            raws.append(DataRawFidsum(d))

        return raws


