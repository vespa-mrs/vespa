"""
Read a Siemens Twix meas.dat file created by Siemens svs_se SVS MRS sequence. 
This is a standared PRESS sequence. It acquires a number of FIDs for one data 
acquisition. If multi-Rx coil is used each scan has multiple FIDs associated
with it, one for each receive channel. 

Data is returned in a DataRawUncomb object populated with each of the FIDs
acquired. The data array has [nfids, nchans, npts] dimensions. DataRawUncomb
is a 'Fidsum' child class. It has a Pre-processing tab inserted into the 
processing chain. This is where the multi-coil data is combined, now, rather
than here in the data parsing "Raw" step.

"""


# Python modules
from __future__ import division
import collections

# 3rd party modules
import numpy as np

# Our modules
import vespa.common.twix_parser as twix_parser
import vespa.common.constants as constants
import vespa.common.util.misc as util_misc
import vespa.common.mrs_data_raw_uncomb as mrs_data_raw_uncomb
import vespa.analysis.src.fileio.raw_reader as raw_reader 
import vespa.common.twix_sort as twix_sort


# Change to True to enable the assert() statements sprinkled through the code
ASSERTIONS_ENABLED = False

# data is complex64 per Siemens documentation for Twix
NUMPY_DATA_TYPE = np.complex64
# BYTES_PER_ELEMENT expresses how many bytes each element occupies. You
# shouldn't need to change this definition.
BYTES_PER_ELEMENT = np.zeros(1, dtype=NUMPY_DATA_TYPE).nbytes

# RAWDATA_SCALE is based on ICE_RAWDATA_SCALE in SpecRoFtFunctor.cpp
RAWDATA_SCALE = 131072.0 * 256.0


class RawReaderSiemensTwixSvsSeUncomb(raw_reader.RawReader):
    # This inherits from raw_reader.RawReader (q.v.). The only methods you
    # need to implement are __init__() and read_raw(). You *may* want to 
    # override or supplement some of raw_reader.RawReader's other methods.

    def __init__(self):
        raw_reader.RawReader.__init__(self)
        
        # The twix files all have an extension of *.dat, thus the filter below.
        self.filetype_filter = "Twix File (*.dat)|*.dat;*.DAT;"
        self.multiple = False
        
        
    def read_raw(self, filename, ignore_data=False, open_dataset=None):
        """
        Given the name of a twix file, returns a DataRawUncomb object
        populated with the parameters and data represented by the file. 
        The 'Fidsum' flavor of raw data object has individual FIDs returned in
        it rather than one summed FID.

        The ignore_data option is not implemented for this reader. 
        
        The open_dataset attribute is not used in this reader. 
        
        """
        measurement = twix_parser.TwixRaid()
        measurement.populate_from_file(filename)

        dat = twix_sort.TwixSort(measurement).get_data_numpy_scan_channel_order()
        
        # Fill a dictionary with the minimum necessary parameters needed from
        # the twix file in order to create the DataRawUncomb object.
        d = _extract_parameters(measurement.evps)

        coil_ids = sorted(set([item.channel_id for item in measurement.scans]))

        d["filename"]    = filename
        d["nscans"]      = len(measurement.scans)
        d["ncoils"]      = len(coil_ids)
        d["nfids"]       = int(d["nscans"]/d["ncoils"])
        d["dims"][2]     = d["ncoils"]  
        d["start_point"] = measurement.scans[0].free_parameters[0]
        
        # scale, remove os, and apply complex conj to data
        dat = _siemens_process(dat, d)

        if d["remove_os"]:
            # update parameter dictionary if oversampling removed
            d["sw"]   = d["sw"] / 2.0
        
        del d["dims"]

        # Create a DataRawUncomb object.
        
        d["data"] = dat
        d["data_source"] = filename
        raw = mrs_data_raw_uncomb.DataRawUncomb(d)
        
        return [raw, ]





####################    Internal functions start here     ###############

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
    d["remove_os"]      = remove_oversample_flag
    d["readout_os"]     = float(_get_siemens_xprotocol(evps[0][1], "ReadoutOS", 1.0))
    d["sequence_type"]  = header.get("tSequenceFileName", "wbnaa")
    d["frequency"]      = float(header["sTXSPEC.asNucleusInfo[0].lFrequency"])/1000000.0
    d["dims"]           = mrs_data_raw_uncomb.DataRawUncomb.DEFAULT_DIMS
    d["dims"][0]        = int(header["sSpecPara.lVectorSize"]) 
    d["dims"][1]        = int(header["lAverages"])*(int(header["lRepetitions"])+1)
    d["dims"][2]        = 1   # will be set at read_raw level using scans.channel_ids
    d["dims"][3]        = 1 
    d["seqte"]          = float(header["alTE[0]"])*1e-6
    
    nuc = header["sTXSPEC.asNucleusInfo[0].tNucleus"].replace('"',' ').strip()
    if nuc == '1H':
        d["midppm"] = constants.DEFAULT_PROTON_CENTER_PPM
    else:
        d["midppm"] = constants.DEFAULT_XNUCLEI_CENTER_PPM

    d["nucleus"] = nuc

    return d  
    

def _siemens_process(scans, d):
    """
    This Import method was originally developed to read/process from Siemens
    IceSpectroEdit, so I leave in the documentation below describing where the
    processing code came from.
    
    Given a numpy array of [fids,scans,pts] and the extracted parameter dictionary 
    we apply the standard Siemens processing as for ICE program IceSpectroEdit. 
    My way of removing oversampling seems to not match 100% but the
    other steps are in the order and perform very similarly to the ICE program.

    """
    dims      = d["dims"]
    dim0      = dims[0]
    ncoils    = d["ncoils"]
    nfids     = d["nfids"]

    # lead to typical integral values of 1 or larger which are nicely displayed
    scale = RAWDATA_SCALE / float(nfids) 

    dat = np.empty([ncoils,nfids,dim0], dtype=np.complex128)

    # for each FID/Channel:
    #   scale for range that is easily displayed
    #   remove pre-/post- points and oversampling, if any
    #   apply complex conjugate to swap x-axis for display
    
    for i in range(nfids):
        
        for j in range(ncoils):
           
            chan = scans[i,j,:] * scale
            chan = _remove_oversampling_basic(chan, d)
            chan = np.conjugate(chan)

            dat[j,i,:] = chan       # index coils on outside so eventually collapse to 1,1,nfid,dim0
            
    return dat
    
    
def _parse_protocol_data(protocol_data):
    """
    Returns a dictionary containing the name/value pairs inside the
    "ASCCONV" section of the MrProtocol or MrPhoenixProtocol elements
    of a Siemens CSA Header tag.
    
    Protocol_data is a large string (e.g. 32k) that lists a lot of
    variables in a JSONish format with which I'm not familiar. Following
    that there's another chunk of data delimited by the strings you see
    below.
    
    That chunk is a list of name=value pairs, INI file style. We
    ignore everything outside of the ASCCONV delimiters. Everything inside
    we parse and return as a dictionary. 
    
    As of the Siemens VD scanner software version the starting string is
    no longer ### ASCCONV BEGIN ### rather it seems to have some other
    info about what was converted inserted after the BEGIN and before 
    the ### delimiter. To get around this for now, we search just for the  
       clean_header[-1] len(protocol_data)
    beginning of the string ### ASCONV BEGIN, and then throw away the
    first line after we split the string into lines.

    """
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



def _remove_oversampling_basic(chan, d):   
    """ 
    Basic FFT algorithm to remove data acquisition oversampling. This is a 
    simple OS removal, not the full Siemens algorithm.
    
    There are five steps to remove oversampling:
    1. forward transform data at acqdim0 resolution
    2. shift data by dim0/2 so we have a FID centered at dim0/2 before snip
    3. snip out first dim0 points and roll back by -dim0/2
    4. inverse transform data at dim0 resolution
    5. (omitted) multiply by 0.5 to maintain same scale as original data
    
    (based on the results in the actual IceSpecEdit program, but still these 
    results are a few percent off of Siemens results, ah well) 
    
    """
    dims        = d["dims"]
    dim0        = dims[0]
    acqdim0     = dims[0] * int(d["readout_os"])
    remove_os   = d["remove_os"]
    start_point = d['start_point']
    end_point   = start_point + acqdim0
    

    # remove extra points at the beginning/end of FID
    chan = np.array(chan[start_point:end_point])

    if remove_os:
        chan = np.fft.fft(chan)
        chan = np.roll(chan, int(dim0/2))
        chan = np.fft.ifft(np.roll(chan[:dim0], int(-dim0/2))) 

    return chan


def _remove_oversampling(scan, d):   
    """ not currently in use as incomplete re. Siemens ICE code """
    
    vector_size = d["vector_size"]
    remove_os   = d["remove_os"]
    left_points = scan.pre
    right_points = scan.post
    reduced_points = scan.samples_in_scan - scan.pre - scan.post
    half_vector_size = int(vector_size / 2)
    
    if (reduced_points % vector_size) != 0:
        raise ValueError('remove_oversampling: final data size not multiple of vector size.')


    if not remove_os:
        # keep oversampled points but remove extra points
        start_point = scan.pre
        scan.data = scan.data[start_point:start_point+reduced_points]
    else:
        # remove oversampled data
        shift_points = scan.post if scan.post < scan.pre else scan.pre
        
        if shift_points == 0:
            # no extra pts available, final data will show truncation artifact
            start_point = scan.pre
            data = np.array(scan.data[start_point:start_point+reduced_points])
            data = np.fft.fft(data, n=vector_size)
            data = np.fft.ifft(data) * 0.5
            scan.data = data.tolist()
            
        else:
            # Extra pts available to use for removing truncation artifact.
            # Process data twice, centering signal to the left and right of kSpaceCentreColumn (TE)
            # Retrieve half of final signal from each set.
            pass