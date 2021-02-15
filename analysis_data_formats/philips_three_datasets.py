"""
Routines for reading a Philips .spar/.sdat format file that has multiple
FIDs for one SVS data acquisition inside and returning a DataRaw object 
populated with the file's data.
"""


# Python modules
from __future__ import division
import os.path

# 3rd party modules
import numpy as np

# Our modules
import vespa.common.constants as constants
import vespa.common.util.misc as util_misc
import vespa.common.configobj as configobj
import vespa.analysis.fileio.util_philips as fileio_util
import vespa.analysis.fileio.raw_reader as raw_reader 

from vespa.common.mrs_data_raw import DataRawFidsum

# data is complex64 per Philips documentation for SDAT
NUMPY_DATA_TYPE = np.complex64
# BYTES_PER_ELEMENT expresses how many bytes each element occupies. You
# shouldn't need to change this definition.
BYTES_PER_ELEMENT = np.zeros(1, dtype=NUMPY_DATA_TYPE).nbytes


def get_filename_pair(filename):
    """
    Given the name of a SPAR data file (e.g. /home/me/foo.sdat) or
    parameters file (e.g. c:/stuff/xyz.spar), returns a tuple of
    (parameters_filename, data_filename). It doesn't matter if the
    filename is a fully qualified path or not.

    This is a little shaky on case sensitive file systems since I assume
    that the file extensions are either all upper or all lower case. If, for
    instance, the data file is foo.sdat and the param file is FOO.SPAR, this
    code won't generate the correct name.
    """
    # filenames are the same except for the last three letters.
    parameters_filename = data_filename = filename[:-3]

    if filename[-1:].isupper():
        data_filename += 'DAT'
        parameters_filename += 'PAR'
    else:
        data_filename += 'dat'
        parameters_filename += 'par'

    return (parameters_filename, data_filename)



class RawReaderPhilipsThreeDatasets(raw_reader.RawReader):
    # This inherits from raw_reader.RawReader (q.v.). The only methods you
    # need to implement are __init__() and read_raw(). You *may* want to 
    # override or supplement some of raw_reader.RawReader's other methods.

    def __init__(self):
        raw_reader.RawReader.__init__(self)
        
        # The sample files given to us all had uppercase extensions, so my guess
        # is that all Philips files are created this way. The GTK file open 
        # dialog takes case sensitivity seriously and won't show foo.SPAR if we
        # use a filetype filter of '*.spar'. Hence we specify both upper and
        # lower case extensions.
        self.filetype_filter = "Spectra (*.spar)|*.spar;*.SPAR;"
        self.multiple = False
        
        
    def read_raw(self, filename, ignore_data=False, open_dataset=None):
        """
        Given the name of a .spar or .sdat file, returns a DataRawFidSum object
        populated with the parameters and data represented by the file pair.

        When ignore_data is True, this function only reads the parameters file
        which can be much faster than reading both params & data.
        
        The open_dataset attribute is not used in this reader. 
        
        """
        parameters_filename, data_filename = get_filename_pair(filename)
        
        if not os.path.isfile(parameters_filename):
            raise raw_reader.FileNotFoundError("I can't find the parameters file '%s'" % parameters_filename)

        if not ignore_data and not os.path.isfile(data_filename):
            raise raw_reader.FileNotFoundError("I can't find the data file '%s'" % data_filename)
        
        # Read the SPAR file and extract the stuff I need.
        header = open(parameters_filename, "rb").read()

        d = _extract_parameters(header)

        d["data_source"] = filename

        # Read data, too, if the caller wants me to do so.
        dims = d["dims"]
        shape = dims[::-1]
        del d["dims"]

        if ignore_data:
            # Create zero data
            datas = [np.zeros(shape, NUMPY_DATA_TYPE)]
        else:
            datas = _read_data(data_filename, dims[0], d["rows"])

        # Create a DataRawFidsum out of the first set of data in the list.
        d["data"] = datas[0]
        raw1 = DataRawFidsum(d)
        # Concatenate the remainder onto the first.
        for data in datas[1:8]:
            data.shape = shape

            d["data"] = data

            raw1.concatenate(DataRawFidsum(d))

        # Create a DataRawFidsum out of the second set of data in the list.
        d["data"] = datas[8]
        raw2 = DataRawFidsum(d)
        # Concatenate the remainder onto the first.
        for data in datas[9:136]:
            data.shape = shape

            d["data"] = data

            raw2.concatenate(DataRawFidsum(d))

        # Create a DataRawFidsum out of the third set of data in the list.
        d["data"] = datas[136]
        raw3 = DataRawFidsum(d)
        # Concatenate the remainder onto the first.
        for data in datas[136:]:
            data.shape = shape

            d["data"] = data

            raw3.concatenate(DataRawFidsum(d))

        return [raw1,raw2,raw3]





####################    Internal functions start here     ###############

def _read_data(filename, npoints, nfids):
    """Given a filename, the number of points in each FID and the number of
    FIDs, reads the data from the file and returns a list of numpy arrays.
    
    The list will contain nfids arrays; the arrays will be 1D arrays with 
    npoints elements.
    """
    f = open(filename, "rb")

    datas = [ ]

    while nfids:
        data = f.read(npoints * BYTES_PER_ELEMENT)

        data = _decode_raw(data)

        data = np.fromiter(data, NUMPY_DATA_TYPE)
        data = np.conjugate(data)

        datas.append(data)

        nfids -= 1

    return datas


def _decode_raw(data):
    """
    Given a string of data in raw format, returns an iterable (tuple or list) 
    of Python complexes.
    """
    data = _vax_to_ieee_single_float(data)

    # Complex numbers are written as pairs of floats (or doubles). Here
    # I knit the (real, imaginary) pairs back into complex numbers.
    data = fileio_util.collapse_complexes(data)

    return data
    
    
def _extract_parameters(header):
    """
    Given the contents of an SPAR file as a string, extracts a few specific
    parameters and returns a flat dict containing those parameters and their 
    value. 
    The returned dict is appropriate for passing to DataRaw.inflate().
    """
    d = { }

    header = util_misc.normalize_newlines(header)

    # A copy of this goes into the dict.
    d["header"] = header

    # The header is in a proprietary format that we can massage into something
    # that looks enough like an INI file to make ConfigObj accept it. I'm 
    # not sure it's a great idea to abuse ConfigObj this way, but it seems
    # to work.
    header = header.split("\n")

    # The proprietary format uses exclamation points to indicate comments
    header = [line.replace('!', '#',   1) for line in header]

    # The proprietary format uses colon as the key/value delimiter.
    header = [line.replace(':', ' = ', 1) for line in header]

    header = configobj.ConfigObj(header)

    d["sequence_type"]  = header["scan_id"]
    d["frequency"]      = float(header["synthesizer_frequency"])/1000000.0
    d["sw"]             = float(header["sample_frequency"])
    
    d["dims"]           = [1,1,1,1024]
    d["dims"][0]        = int(header["samples"])
    d["dims"][1]        = int(header["dim2_pnts"])
    d["dims"][2]        = int(header["dim3_pnts"])
    d["dims"][3]        = 1 # FIXME bjs is there a third spatial dim?

    d["rows"]           = int(header["rows"])

    d["nucleus"]                = header["nucleus"]
    d["seqte"]                  = float(header["echo_time"])
    d["slice_thickness"]        = float(header["cc_size"])

    return d    
    

def _vax_to_ieee_single_float(data):
    """Converts a float in Vax format to IEEE format.

    data should be a single string of chars that have been read in from 
    a binary file. These will be processed 4 at a time into float values.
    Thus the total number of byte/chars in the string should be divisible
    by 4.
    
    Based on VAX data organization in a byte file, we need to do a bunch of 
    bitwise operations to separate out the numbers that correspond to the
    sign, the exponent and the fraction portions of this floating point
    number
    
    role :      S        EEEEEEEE      FFFFFFF      FFFFFFFF      FFFFFFFF
    bits :      1        2      9      10                               32
    bytes :     byte2           byte1               byte4         byte3    
    
    """
    f = []
    nfloat = int(len(data) / 4)
    for i in range(nfloat):

        byte2 = data[0 + i*4]
        byte1 = data[1 + i*4]
        byte4 = data[2 + i*4]
        byte3 = data[3 + i*4]
        
        # hex 0x80 = binary mask 10000000
        # hex 0x7f = binary mask 01111111
        
        sign  =  (ord(byte1) & 0x80) >> 7
        expon = ((ord(byte1) & 0x7f) << 1 )  + ((ord(byte2) & 0x80 ) >> 7 )
        fract = ((ord(byte2) & 0x7f) << 16 ) +  (ord(byte3) << 8 ) + ord(byte4)
        
        if sign == 0:
            sign_mult = 1.0
        else:
            sign_mult = -1.0;
        
        if 0 < expon:
            # note 16777216.0 == 2^24  
            val = sign_mult * (0.5 + (fract/16777216.0)) * pow(2.0, expon - 128.0)   
            f.append(val)
        elif expon == 0 and sign == 0:
            f.append(0)
        else: 
            f.append(0)
            # may want to raise an exception here ...
        
    return f