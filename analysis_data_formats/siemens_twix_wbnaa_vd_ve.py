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
import vespa.common.twix_parser_multi_raid as twix_parser_multi_raid
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



class RawReaderSiemensTwixWbnaaVDVE(raw_reader.RawReader):
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
        self.filetype_filter = "Twix File (*.dat)|*.dat;*.DAT;"
        self.multiple = False
        
        
    def read_raw(self, filename, ignore_data=False, open_dataset=None):
        """
        Given the name of a .spar or .sdat file, returns a DataRawWbnaa object
        populated with the parameters and data represented by the file pair.

        When ignore_data is True, this function only reads the parameters file
        which can be much faster than reading both params & data.
        
        The open_dataset attribute is not used in this reader. 
        
        """
        #scans, evps = twix_parser.read(filename)

        twix = twix_parser_multi_raid.TwixMultiRaid()
        twix.populate_from_file(filename)

        evps = twix.measurements[-1].evps
        
        d = _extract_parameters(evps)

        d["data_source"] = filename

        # Read data, too, if the caller wants me to do so. 
        dims = d["dims"]
        shape = dims[::-1]
        del d["dims"]

        #datas = _read_data(scans)
        raw = twix.get_data_numpy_channel_scan_order()
        nscans = int(raw.shape[-2]/2)
        datas = []
    
        for i in range(nscans):
            dat = np.conjugate(raw[0,i*2,:] - raw[0,i*2+1,:])
            datas.append(dat)
        
        # Create a DataRawWbnaa out of the first set of data in the list.
        d["data"] = datas[0]
        raw = mrs_data_raw_wbnaa.DataRawWbnaa(d)
        
        # Concatenate the remainder onto the first.
        for data in datas[1:]:
            data.shape = shape
            d["data"] = data
            raw.concatenate(mrs_data_raw_wbnaa.DataRawWbnaa(d))

        return [raw,]





################    Internal helper functions start here     #############

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

    header, clean_header = _parse_protocol_data(evps[3][1])

    # A copy of this goes into the dict.
    d["header"] = clean_header

    remove_oversample_flag = header.get("sSpecPara.ucRemoveOversampling", "0x0")
    remove_oversample_flag = (remove_oversample_flag.strip() == "0x1")
    d["sw"]             = 1.0 / (float(header.get("sRXSPEC.alDwellTime[0]", 1.0)) * 1e-9)

    d["readout_os"]     = float(_get_siemens_xprotocol(evps[0][1], "ReadoutOS", 1.0))
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
    #
    start = protocol_data.find("### ASCCONV BEGIN")
    end = protocol_data.find("### ASCCONV END ###")

    _my_assert(start != -1)
    _my_assert(end != -1)

    clean_start = start
    clean_end   = end + len("### ASCCONV END ###")
    clean_header = protocol_data[clean_start:clean_end]

    start += len("### ASCCONV BEGIN ###")
    protocol_data = protocol_data[start:end]

    lines = protocol_data.split('\n')
    lines = lines[1:]

    # The two lines of code below turn the 'lines' list into a list of
    # (name, value) tuples in which name & value have been stripped and
    # all blank lines have been discarded.
    f = lambda pair: (pair[0].strip(), pair[1].strip())
    lines = [f(line.split('=')) for line in lines if line]

    return dict(lines), clean_header


def _my_assert(expression):
    if ASSERTIONS_ENABLED:
        assert(expression)
        
 
def _get_siemens_xprotocol(head_only, key, default):
    
    head = util_misc.normalize_newlines(head_only)   
    head = head.split("\n")
    
    # find substring 'key' in list even if item in list is not iterable
    items = [el for el in head if isinstance(el, collections.Iterable) and (key in el)]
    
    for item in items:
        start = item.find("{")
        end = item.find("}")
        if start != -1 and end != -1:
            temp = item[start+1:end]
            # remove duplicate white space
            temp = " ".join(temp.split())
            temp = temp.split()
            
            if len(temp) == 1:
                return temp[0]
            elif len(temp) == 3:
                return temp[2]
            else:
                return temp
    
    return default
    
        

            