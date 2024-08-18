"""
Routine for reading a Siemens Twix meas.dat format file created by CMRR sLASER.
Data is returned in a DataRawCmrrSlaser object

These are for older data in Gulin Oz's 'Long' dataset from the BRP 2015

When we moved to using the PyMapVBVD library to read Siemens Twix files, We no
longer had access to the 'prep' scans data. Prep data is used in this data in
lieu of actual ECC and Water quant FIDs. So we created a separate module to
read the 'Long' data that still uses the original BJS Twix parser that has now
been renamed to RawReaderSiemensTwixBjs in siemens_twix_bjs.py module

"""


# Python modules

# 3rd party modules
import numpy as np

# Our modules
from vespa.common.mrs_data_raw import DataRawCmrrSlaser
from vespa.analysis.fileio.siemens_twix_bjs import RawReaderSiemensTwixBjs


# RAWDATA_SCALE is based on ICE_RAWDATA_SCALE in SpecRoFtFunctor.cpp
RAWDATA_SCALE = 131072.0 * 256.0


class RawReaderSiemensTwixSlaserCmrrVbGulinLong(RawReaderSiemensTwixBjs):

    def read_raw(self, filename, ignore_data=False, open_dataset=None):
        """
        Given filename, returns a DataRawCmrrSlaser object.
        - The ignore_data option is not implemented for this reader.
        - The open_dataset attribute is not used in this reader.

        """
        twix, version_flag = self.get_twix(filename)

        d = self._get_parameters(twix.current)
        d["data_source"] = filename

        data = d['data']
        prep = d['prep']

        _, _, navg, npts = data.shape

        if d["remove_os"]:
            data = self._remove_oversampling_basic(data)
            if prep is not None:
                prep = self._remove_oversampling_basic(prep)
            d["sw"] = d["sw"] / 2.0

        data = np.conjugate(data)
        if prep is not None:
            prep = np.conjugate(prep)

        raws = []

        # dataset1 - scan 0, water unsuppressed for coil combine and (maybe) ECC

        d["data"] = prep[:,:,0:1,:].copy() * RAWDATA_SCALE / float(1.0)
        d["data_source"] = filename + '.combine'
        raws.append(DataRawCmrrSlaser(d))

        # dataset2 - scan 1-2, water unsuppressed for water scale

        d["data"] = prep[:,:,1:2,:].copy() * RAWDATA_SCALE / float(1.0)
        d["data_source"] = filename + '.water1'
        raws.append(DataRawCmrrSlaser(d))

        # dataset3 - scans 5-68 (64 total), metabolite data with WS

        d["data"] = data.copy() * RAWDATA_SCALE / float(navg)
        d["data_source"] = filename + '.metab64'
        raws.append(DataRawCmrrSlaser(d))

        return raws


class RawReaderSiemensTwixSlaserCmrrVbGulinLongWater(RawReaderSiemensTwixBjs):

    def read_raw(self, filename, ignore_data=False, open_dataset=None):
        """
        Given filename, returns a DataRawCmrrSlaser object.
        - The ignore_data option is not implemented for this reader.
        - The open_dataset attribute is not used in this reader.

        """
        twix, version_flag = self.get_twix(filename)

        # 'scan' returns data in [ncha, navg, npts] order - easy to convert to mrs_data_raw dims
        d = self._get_parameters(twix.current, index='scan')
        d["data_source"] = filename

        data = d['data']
        prep = d['prep']
        nrep, _, navg, npts = data.shape

        if d["remove_os"]:
            data = self._remove_oversampling_basic(data)
            if prep is not None:
                prep = self._remove_oversampling_basic(prep)
            d["sw"] = d["sw"] / 2.0

        data = np.conjugate(data)
        if prep is not None:
            prep = np.conjugate(prep)

        raws = []

        # dataset1 - scan 0, water unsuppressed for coil combine and (maybe) ECC

        # tmp = prep[:,:,0,:].copy()
        # tmp.shape = (tmp.shape[0], tmp.shape[1], 1, tmp.shape[2])
        d["data"] = prep[:,:,0:1,:].copy() * RAWDATA_SCALE / float(1.0)
        d["data_source"] = filename + '.combine'
        raws.append(DataRawCmrrSlaser(d))

        # dataset2 - scan 1-2, water unsuppressed for water scale

        if prep.shape[2] > 1:
            d["data"] = prep[:,:,1:2,:].copy() * RAWDATA_SCALE / float(1.0)
        else:
            d["data"] = prep[:,:,0:1,:].copy() * RAWDATA_SCALE / float(1.0)
        d["data_source"] = filename + '.water1'
        raws.append(DataRawCmrrSlaser(d))

        # dataset3 - scans 5-68 (64 total), metabolite data with WS

        d["data"] = data.copy() * RAWDATA_SCALE / float(navg)
        d["data_source"] = filename + '.metab64'
        raws.append(DataRawCmrrSlaser(d))

        return raws

