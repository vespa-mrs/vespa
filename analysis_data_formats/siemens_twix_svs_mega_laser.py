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
import collections

# 3rd party modules
import numpy as np

# Our modules
import vespa.common.twix_parser as twix_parser
import vespa.common.constants as constants
import vespa.common.util.misc as util_misc
from vespa.common.mrs_data_raw import DataRawEditFidsum
import vespa.analysis.fileio.raw_reader as raw_reader 

from vespa.common.constants import DEGREES_TO_RADIANS, RADIANS_TO_DEGREES

# Change to True to enable the assert() statements sprinkled through the code
ASSERTIONS_ENABLED = False

# data is complex64 per Siemens documentation for Twix
NUMPY_DATA_TYPE = np.complex64
# BYTES_PER_ELEMENT expresses how many bytes each element occupies. You
# shouldn't need to change this definition.
BYTES_PER_ELEMENT = np.zeros(1, dtype=NUMPY_DATA_TYPE).nbytes

# RAWDATA_SCALE is based on ICE_RAWDATA_SCALE in SpecRoFtFunctor.cpp
RAWDATA_SCALE = 131072.0 * 256.0


class RawReaderSiemensTwixSvsMegaLaser(raw_reader.RawReader):
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
        Given the name of a twix file, returns a DataRawEditFidsum object
        populated with the parameters and data represented by the file. 
        The 'Fidsum' flavor of raw data object has individual FIDs returned in
        it rather than one summed FID.

        The ignore_data option is not implemented for this reader. 
        
        The open_dataset attribute is not used in this reader. 
        
        """
        scans, evps = twix_parser.read(filename)
        
        # Fill a dictionary with the minimum necessary parameters needed from
        # the twix file in order to create the DataRawEditFidsum objects.
        d = _extract_parameters(evps)

        # Parse scan data into four arrays 
        dat_on, dat_off, dat_sum, dat_dif = _read_data(scans, d)

        if d["remove_os"]:
            #pass
            # update parameter dictionary if oversampling removed
            #
            # Note. we are in Twix file here so OS has not yet been
            #  removed. So no divide by 2 here.
            #
            d["sw"]   = d["sw"] / 2.0
        
        dims = d["dims"]
        shape = dims[::-1]
        del d["dims"]

        # Create a DataRawEditFidsum objects for the four different states
        # of the data ON, OFF, SUM and DIFFERENCE.
        
        d["data"] = dat_on[0]
        d["data_source"] = filename+'.EditON'
        raw1 = DataRawEditFidsum(d)
        # Concatenate the remaining FID data for this state into the first
        # raw object. This will result in one object with a 2D data set
        # inside it where the second dimension is the FID index. This is 
        # what Vespa expects a DataRawEditFidsum object to contain.
        for data in dat_on[1:]:
            data.shape = shape
            d["data"] = data
            raw1.concatenate(DataRawEditFidsum(d))

        d["data"] = dat_off[0]
        d["data_source"] = filename+'.EditOFF'
        raw2 = DataRawEditFidsum(d)
        for data in dat_off[1:]:
            data.shape = shape
            d["data"] = data
            raw2.concatenate(DataRawEditFidsum(d))

        d["data"] = dat_sum[0]
        d["data_source"] = filename+'.Sum'
        raw3 = DataRawEditFidsum(d)
        for data in dat_sum[1:]:
            data.shape = shape
            d["data"] = data
            raw3.concatenate(DataRawEditFidsum(d))

        d["data"] = dat_dif[0]
        d["data_source"] = filename+'.Diff'
        raw4 = DataRawEditFidsum(d)
        for data in dat_dif[1:]:
            data.shape = shape
            d["data"] = data
            raw4.concatenate(DataRawEditFidsum(d))


        # TODO - bjs
        # - check scaling to get these numbers up
        # - check if ECC associate works in concert with other associations 

        return [raw1, raw2, raw3, raw4]





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

    d["remove_os"]   = remove_oversample_flag

    d["readout_os"]     = float(_get_siemens_xprotocol(evps[0][1], "ReadoutOS", 1.0))
    d["sequence_type"]  = header.get("tSequenceFileName", "wbnaa")
    d["frequency"]      = float(header["sTXSPEC.asNucleusInfo[0].lFrequency"])/1000000.0
    d["dims"]           = [1,1,1,2048]
    d["dims"][0]        = int(header["sSpecPara.lVectorSize"]) 
    d["dims"][1]        = 1 # concat will take care of header["lAverages"]
    d["dims"][2]        = 1
    d["dims"][3]        = 1 
    d["seqte"]          = float(header["alTE[0]"])*1e-6
    
    d["nucleus"] = header["sTXSPEC.asNucleusInfo[0].tNucleus"].replace('"',' ').strip()

    return d  
    
    
    
def _read_data(scans, d):
    """
    Given a list of scans and the extracted parameter dictionary we process
    the data in as similar a way to Siemens ICE program IceSpectroEdit as
    we can. My way of removing oversampling seems to not match 100% but the
    other steps are in the order and perform very similarly to the ICE program.

    When we are done parsing the ON and OFF state FIDs into arrays, we create
    SUM and DIF arrays to also send back.
    
    """
    dims    = d["dims"]
    dim0    = dims[0]
    acqdim0 = dims[0] * int(d["readout_os"])
    remove_os = d["remove_os"]
    
    
    coil_ids = sorted(set([item.channel_id for item in scans]))
    nscans   = len(scans)
    ncoils   = len(coil_ids)

    # determine if there are any extra points at the beginning or end of FID
    start_point = scans[0].free_parameters[0]
    end_point   = start_point + acqdim0

    # lead to typical integral values of 1 or larger which are nicely displayed
    scale = RAWDATA_SCALE / float(nscans/ncoils) 

    
    # Parse each group of FIDs for N channels for each average as separate
    # from all other scans. Perform the following steps:
    # - scale, all channels
    # - remove OS, all channels
    # - if needed, we remove oversampling, all channels
    # - convert to complex conjugate for proper plotting
    # - accumulate weighting factors and phases, all channels
    # - collate all channels as numpy arrays
    # - calc final weights for each channel 
    # - apply weight and phase corrections to all channels
    # - sum channels to one array, sort into ON or OFF states
    dat_on  = []
    dat_off = []

    for i in range(int(nscans/ncoils)):
        
        # for each average, remove pre-/post- points and apply phase to  
        # correct for coil geometry
        chans = []
        weight = []
        phases = []
        lamda  = 0.0
        lamda2 = 0.0
        
        for j in range(ncoils):
            chan = scans[i*ncoils+j].data
            chan = np.array(chan[start_point:end_point])*scale

            if remove_os:
                #
                # Note that this is a simple OS removal, not the full Siemens
                # algorithm. I have not had time to implement and test that.
                #
                # There are five steps to remove oversampling:
                # 1. forward transform data at acqdim0 resolution
                # 2. shift data by dim0/2 so we have a FID centered at dim0/2 before snip
                # 3. snip out first dim0 points and roll back by -dim0/2
                # 4. inverse transform data at dim0 resolution
                # 5. (omitted) multiply by 0.5 to maintain same scale as original data
                #    (based on the results in the actual IceSpecEdit program,
                #     but still these results are a few percent off of Siemens 
                #     results, ah well) 
                chan = np.fft.fft(chan)
                chan = np.roll(chan, int(dim0/2))
                chan = np.fft.ifft(np.roll(chan[:dim0], int(-dim0/2))) 

            # apply complex conjugate to swap x-axis for proper display of data
            chan = np.conjugate(chan)

            magn = np.abs(chan[0])
            phas = np.conjugate(chan[0])/magn        # normalized complex conj to cancel phase
            phases.append(phas)
            weight.append(magn)
            lamda2 += magn
            lamda  += magn * magn
            chans.append(chan) 
            
        if lamda == 0: lamda = 1.0
        
        flag_norm_to_sum = False    # default for now
        if flag_norm_to_sum:
            # spectro data based normalization to sum of sensitivities
            lamda = lamda2/lamda
        else:
            # spectro data based normalization to square root of sum of squared sensitivities
            lamda = 1.0/np.sqrt(lamda)
        
        weight = [wt*lamda for wt in weight]
        
        # apply weighting and phase corrections
        for j,chan in enumerate(chans):
            chans[j] = chan * weight[j] * phases[j]
        
        # sum corrected FIDs from each coil into one combined FID
        data = np.sum(chans, axis=0)    
              
        # sort into ON or OFF states
        if i % 2 == 0:
            dat_on.append(data)
        else:
            dat_off.append(data)

    dat_sum = []
    dat_dif = []
         
    # Final step, since this is an editing sequence we process data 
    # into ON and OFF states as well. 
    for on,off in zip(dat_on, dat_off):
        dat_sum.append(on + off)
        dat_dif.append(on - off)

    return dat_on, dat_off, dat_sum, dat_dif

    
    
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



def _remove_oversampling(scan, d):   
    """ not currently in use as incomplete """
    
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
