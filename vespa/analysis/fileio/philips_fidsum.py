"""
Routines for reading a Philips .spar/.sdat formats and returning an 
DataRawFidsum object populated with the file's data.
"""

# Python modules

import os.path

# 3rd party modules
import numpy as np

# Our modules
import vespa.analysis.fileio.raw_reader as raw_reader
import vespa.analysis.fileio.util_exceptions as util_exceptions

from vespa.common.mrs_data_raw import DataRawFidsum
from vespa.analysis.fileio.philips_spar import RawReaderPhilipsSpar, \
                                               get_filename_pair, \
                                               _get_parameters




class RawReaderPhilipsFidsum(RawReaderPhilipsSpar):

    def __init__(self):

        raw_reader.RawReader.__init__(self)
        self.filetype_filter = "Spectra (*.spar)|*.spar;*.SPAR"
        self.multiple = True
        

    def read_raw(self, filename, ignore_data=False, *args, **kwargs):
        """
        Given Philips SPAR or SDAT filename, return populated DataRaw object.
        - 'ignore_data' has no effect on this parser

        """
        spar_filename, sdat_filename = get_filename_pair(filename)

        msg = ''
        if not os.path.isfile(spar_filename): msg += "No SPAR parameters file '%s'" % spar_filename
        if not os.path.isfile(sdat_filename): msg += "No SDAT data file '%s'" % sdat_filename
        if msg:
            raise util_exceptions.FileNotFoundError(msg)

        # Read the SPAR file and extract the stuff I need.
        header = open(spar_filename, "r").read()

        d = _get_parameters(header, sdat_filename)
        d["data_source"] = filename

        if not (d["data"].shape[0:2] == (1, 1)):
            msg = "SDAT data is unknown size, shape = %i" % d["data"].shape
            raise util_exceptions.IncorrectDimensionalityError(msg)

        return [DataRawFidsum(d),]


    
