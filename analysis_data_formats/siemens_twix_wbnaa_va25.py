"""
Routines for reading a Siemens Twix meas.dat format file that has multiple
FIDs for one SVS data acquisition inside. Every other FID is phase inverted
and needs to be combined in an (a-b) + (a-b) + ... fashion.

Data is returned in a DataRaw object populated with the file's data.
"""


# Python modules
from __future__ import division
import os.path
import collections

# 3rd party modules
import numpy as np

# Our modules
import vespa.common.twix_parser as twix_parser
import vespa.common.constants as constants
import vespa.common.util.misc as util_misc
import vespa.common.mrs_data_raw_wbnaa as mrs_data_raw_wbnaa
import vespa.analysis.fileio.raw_reader as raw_reader 


# Change to True to enable the assert() statements sprinkled through the code
ASSERTIONS_ENABLED = False

# data is complex64 per Siemens documentation for Twix
NUMPY_DATA_TYPE = np.complex64
# BYTES_PER_ELEMENT expresses how many bytes each element occupies. You
# shouldn't need to change this definition.
BYTES_PER_ELEMENT = np.zeros(1, dtype=NUMPY_DATA_TYPE).nbytes


def get_filename_pair(filename):
    """
    Given the name of a twix VA25 data file (e.g. /home/me/foo.out) or
    parameters file (e.g. c:/stuff/xyz.asc), returns a tuple of
    (parameters_filename, data_filename). It doesn't matter if the
    filename is a fully qualified path or not.

    This is a little shaky on case sensitive file systems since I assume
    that the file extensions are either all upper or all lower case. If, for
    instance, the data file is foo.rsd and the param file is FOO.RSP, this
    code won't generate the correct name.
    """
    # filenames are the same except for the last letter.
    parameters_filename = data_filename = filename[:-3]

    if filename[-1:].isupper():
        data_filename += 'OUT'
        parameters_filename += 'ASC'
    else:
        data_filename += 'out'
        parameters_filename += 'asc'

    return (parameters_filename, data_filename)


class RawReaderSiemensTwixWbnaaVa25(raw_reader.RawReader):
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
        self.filetype_filter = "Twix File (*.out/*.asc)|*.dat;*.DAT;*.asc;*.ASC;"
        self.multiple = False
        
        
    def read_raw(self, filename, ignore_data=False, open_dataset=None):
        """
        Given the name of a .spar or .sdat file, returns a DataRawWbnaa object
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

        scans, _ = twix_parser.read(data_filename)

        # Read the RSP file and extract the stuff I need.
        #
        # The with statement is handy when you have two related operations
        # which you'd like to execute as a pair with a block of code in 
        # between. The classic example is opening a file, manipulating the
        # file, then closing it.
        with open(parameters_filename, "rb") as f:
            evps = f.read()
        
        d = _extract_parameters(evps)

        d["data_source"] = filename

        # Read data, too, if the caller wants me to do so. 
        dims = d["dims"]
        shape = dims[::-1]
        del d["dims"]

        datas = _read_data(scans)

        # Create a DataRawWbnaa out of the first set of data in the list.
        d["data"] = datas[0]
        raw = mrs_data_raw_wbnaa.DataRawWbnaa(d)
        
        # Concatenate the remainder onto the first.
        for data in datas[1:]:
            data.shape = shape
            d["data"] = data
            raw.concatenate(mrs_data_raw_wbnaa.DataRawWbnaa(d))

        return [raw,]





####################    Internal functions start here     ###############

def _read_data(scans):
    """
    Given a filename, the number of points in each FID and the number of
    FIDs, reads the data from the file and returns a list of numpy arrays.
    
    The list will contain nfids arrays; the arrays will be 1D arrays with 
    npoints elements.
    
    """
    nscans = int(len(scans)/2)
    datas = []

    for i in range(nscans):
        dat1 = np.array(scans[i*2].data)
        dat2 = np.array(scans[i*2+1].data)
        data = np.conjugate(dat1-dat2)
        datas.append(data)

## comment out the code above and uncomment below to see all scans
## rather than doing difference paired data

#    for scan in scans:
#        datas.append(np.conjugate(np.array(scan.data)))

    return datas


    
def _extract_parameters(evps):
    """
    Given the contents of an SPAR file as a string, extracts a few specific
    parameters and returns a flat dict containing those parameters and their 
    value. 
    The returned dict is appropriate for passing to DataRaw.inflate().
    
    """
    d = { }

    header, clean_header = _parse_protocol_data(evps)

    # A copy of this goes into the dict.
    d["header"] = clean_header

    remove_oversample_flag = header.get("sSpecPara.ucRemoveOversampling", "0x0")
    remove_oversample_flag = (remove_oversample_flag.strip() == "0x1")
    d["sw"]             = 1.0 / (float(header.get("sRXSPEC.alDwellTime[0]", 1.0)) * 1e-9)

    d["readout_os"]     = float(header.get("m_flReadoutOSFactor", "1"))
    d["sequence_type"]  = header.get("tSequenceFileName", "wbnaa")
    d["frequency"]      = float(header["sTXSPEC.asNucleusInfo[0].lFrequency"])/1000000.0
    d["dims"]           = mrs_data_raw_wbnaa.DataRawWbnaa.DEFAULT_DIMS
    d["dims"][0]        = int(header["sSpecPara.lVectorSize"]) * int(d["readout_os"])
    d["dims"][1]        = 1 # concat will take care of header["lAverages"]
    d["dims"][2]        = 1
    d["dims"][3]        = 1 
    d["seqte"]          = float(header["alTE[0]"])*1e-6
    
    d["nucleus"] = header["sTXSPEC.asNucleusInfo[0].tNucleus"].replace('"',' ').strip()

    return d    
    
    
    
def _parse_protocol_data(protocol_data):
    """
    Returns a dictionary containing the name/value pairs inside the
    "ASCCONV" section of the MrProtocol or MrPhoenixProtocol elements
    of a Siemens CSA Header tag.
    """
    # Protocol_data is a large string (e.g. 32k) that lists a lot of
    # variables in a JSONish format with which I'm not familiar. Following
    # that there's another chunk of data delimited by the strings you see
    # below.
    # That chunk is a list of name=value pairs, INI file style. We
    # ignore everything outside of the ASCCONV delimiters. Everything inside
    # we parse and return as a dictionary. 
    #
    # As of the Siemens VD scanner software version the starting string is
    # no longer ### ASCCONV BEGIN ### rather it seems to have some other
    # info about what was converted inserted after the BEGIN and before 
    # the ### delimiter. To get around this for now, we search just for the  clean_header[-1] len(protocol_data)
    # beginning of the string ### ASCONV BEGIN, and then throw away the
    # first line after we split the string into lines.
    
    # Parse First Section ...
    
    start = protocol_data.find("### ASCCONV BEGIN")
    end = protocol_data.find("### ASCCONV END ###")

    _my_assert(start != -1)
    _my_assert(end != -1)

    clean_start = start
    clean_end   = end + len("### ASCCONV END ###")
    clean_header1 = protocol_data[clean_start:clean_end]

    start += len("### ASCCONV BEGIN ###")
    protocol_data1 = protocol_data[start:end]

    lines1 = protocol_data1.split('\n')
    lines1 = lines1[1:]

    # The two lines of code below turn the 'lines' list into a list of
    # (name, value) tuples in which name & value have been stripped and
    # all blank lines have been discarded.
    f = lambda pair: (pair[0].strip(), pair[1].strip())
    lines1 = [f(line.split('=')) for line in lines1 if line]

    # Parse Second Section ... YAPS buffer

    start2 = protocol_data.find("The YAPS buffer:")
    _my_assert(start2 != -1)

    clean_start2 = start2
    clean_header2 = protocol_data[clean_start2:]

    start2 += len("The YAPS buffer:")
    protocol_data2 = protocol_data[start2:]

    _lines2 = protocol_data2.split('\n')
    _lines2 = _lines2[1:]
    lines2 = []

    for line in _lines2:
        if line:
            akey, aval = line.split('=')
            aval = aval.translate(None,"[]")
            lines2.append((akey.strip(),aval.strip()))


    # Combine both sections and return dictionary

    lines = lines1 + lines2
    clean_header = clean_header1 + "\n" + clean_header2

    return dict(lines), clean_header


def _my_assert(expression):
    if ASSERTIONS_ENABLED:
        assert(expression)
        
 

    
        

            