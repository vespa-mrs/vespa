"""
Routines for reading a Siemens *.rda format and returning an 
DataRaw object populated with the file's data.
"""

# Python modules
from __future__ import division

import os
import struct

from ast import literal_eval
from collections import OrderedDict

# 3rd party modules
import numpy as np

# Our modules
import vespa.common.constants               as constants
import vespa.common.util.misc               as util_misc
import vespa.common.configobj               as configobj
import vespa.common.util.fileio             as util_fileio
import vespa.analysis.fileio.raw_reader     as raw_reader 
import vespa.analysis.fileio.util_philips   as fileio_util_philips
import vespa.analysis.fileio.util_exceptions as util_exceptions

from vespa.common.mrs_data_raw import DataRawFidsum


# complex data are binary data in little endian single precision IEEE float format as per Philips
NUMPY_DATA_TYPE = np.complex64

def _get_filename_pair(filename):
    """
    Given the name of either a Philips .data/.list file returns a tuple of
    (list_filename, data_filename). It doesn't matter if the
    filename is a fully qualified path or not.

    """
    base = os.path.splitext(os.path.basename(filename))[0]
    path = os.path.dirname(os.path.abspath(filename))
    
    fname_data = os.path.join(path, base+'.data')
    fname_list = os.path.join(path, base+'.list')

    return (fname_list, fname_data)


def _convert_string(s):
    '''
    Try to convert string literal to a number or a bool such that '1.0' and '1' 
    will both return 1 and that 'False' and 'True' return as bools not numbers.

    '''
    if isinstance(s, str):      # It's a string.  Does it represnt a literal?
        try:
            val = literal_eval(s)
        except:                 # s doesn't represnt any sort of literal so no conversion done.
            val = s
    else:                       # It's already something other than a string
        val = s

    if isinstance(val, float):  # Is the float actually an int? (i.e. is the float 1.0 ?)
        if val.is_integer(): 
            return int(val)
        return val              # It really is a float

    return val



class RawReaderPhilipsDataList(raw_reader.RawReader):
    
    def __init__(self):
        raw_reader.RawReader.__init__(self)
        
        self.filetype_filter = "Philips RAW data file (*.data)|*.data"
        self.multiple = False
        

    def read_raw(self, filename, ignore_data=False, default_params=None): 
        """
        Given the name of a .data file, returns a DataRaw object
        populated with the parameters and data represented by the file pair.

        When ignore_data is True, this function only reads the parameters file
        which might be faster than reading both params & data.
        """
        list_filename, data_filename = _get_filename_pair(filename)

        if not os.path.isfile(list_filename):
            msg = "I can't find the .list file '%s'" % list_filename
            raise util_exceptions.FileNotFoundError(msg)

        if not ignore_data and not os.path.isfile(data_filename):
            msg = "I can't find the .data file '%s'" % data_filename
            raise util_exceptions.FileNotFoundError(msg)
        
        # Read the .list file and extract the stuff I need.
        header = open(list_filename, "rb").read()
    
        d, data_vector_index = _extract_parameters(header)
        
        data, noise, dynamic = _parse_data(d, data_vector_index, data_filename)

        base_params = {"frequency" : 127.748488, "sw" : 2500.0 }
        if default_params is None:
            default_params = base_params
        
        for key in ["frequency", "sw"]:
            d[key] = default_params[key] if key in default_params.keys() else base_params[key]
        
        else:
            d["frequency"] = 127.748488
            d["sw"]        = 2500.0
            
        d["dims"] = list(d["dims_0"])        
        d['data_source'] = data_filename

        d["data"] = data

        return [DataRawFidsum(d),]
    
    
####################    Internal functions start here     ###############

    
def _extract_parameters(header, frequency=127.9):
    """
    Given the contents of a .list file as a string, extracts a few specific
    parameters and returns a flat dict containing those parameters and their 
    value. 

    """
    d = OrderedDict()
    
    header = util_misc.normalize_newlines(header)
    header = header.split("\n")     # now have list of strings
    
    # header = [item.split() for item in header]
    
    indx_general   = [idx for idx, s in enumerate(header) if 'GENERAL INFORMATION' in s][0]
    indx_data_str  = [idx for idx, s in enumerate(header) if 'START OF DATA VECTOR INDEX' in s][0]
    indx_data_end  = [idx for idx, s in enumerate(header) if 'END OF DATA VECTOR INDEX' in s][0]
    indx_head_only = [idx for idx, s in enumerate(header) if '# Complex data vector types:' in s][0]
    
    hdr_top     = header[0:indx_general]
    hdr_general = header[indx_general:indx_data_str]
    hdr_data    = header[indx_data_str:indx_data_end]
    hdr_last    = header[indx_data_end:]
    
    hdr_head_only = "\n".join(header[0:indx_head_only])
    
    # Parse software/system identification ------------------------------------
    
    items = [['Scan name   :', 'scan_name'],
             ['Dataset name:', 'dataset_name'],
             ['Gyroscan SW release            :', 'gyroscan_sw_release'],
             ['Reconstruction Host SW Version :', 'recon_host_sw_version'],
             ['Reconstruction   AP SW Version :', 'recon_ap_sw_version'],
            ]
    
    for substr, label in items:
        for line in hdr_top:
            if substr in line : 
                d[label] = (line.split(substr,1)[1]).strip()
                break
    
    # Parse general information section ---------------------------------------
    
    # Have to know number_of_mixes to process other parameter
    
    for line in hdr_general:
        if 'number_of_mixes' in line: 
            d['number_of_mixes'] = _convert_string(line.split()[-1])
            nmix = d['number_of_mixes']
            break
    
    if 'number_of_mixes' not in d.keys():
        msg = "No 'number_of_mixes' parameter in header."
        raise util_exceptions.IncompleteHeaderParametersError(msg)
    
    for line in hdr_general:
        _parse_parameter_line(d, line, nmix, exclude_labels=['number_of_mixes',])
        
    # Parse data vector index section -----------------------------------------
    #
    # (typ, mix, dyn, card, echo, loca, chan, extr1,extr2, kx, ky, kz, aver, sign, rf, grad, enc, rtop, rr, size, offset)

    data_vector_index = []
    ncoils = set()
    
    for line in hdr_data:
        line = [_convert_string(item) for item in line.split()]
        if line[0] != '#':
            data_vector_index.append(line)
            ncoils.add(line[6])
            
    d['ncoil'] = len(ncoils)

    # Parse last information section ------------------------------------------

    for line in hdr_last:
        _parse_parameter_line(d, line, nmix, exclude_labels=[])


    # calculate a few more parameters -----------------------------------------

    nextra1 = d['number_of_extra_attribute_1_values']
    nave    = d['number_of_signal_averages']

    d["frequency"]   = frequency             # default value
    d["sw_0"]        = float(d['F_range'][0][1] - d['F_range'][0][0] + 1)
    d["dims_0"]      = [1,1,1,1024]
    d["dims_0"][0]   = int(d['t_range'][0][1] - d['t_range'][0][0] + 1)
    d["dims_0"][1]   = int(nextra1[0] * nave[0])
    d["dims_0"][2]   = int(d['ncoil'])
    d["dims_0"][3]   = int(1)

    if nmix == 2: 
        d["sw_1"]      = float(d['F_range'][1][1] - d['F_range'][1][0] + 1)
        d["dims_1"]    = [1,1,1,1024]
        d["dims_1"][0] = int(d['t_range'][1][1] - d['t_range'][1][0] + 1)
        d["dims_1"][1] = int(nextra1[1] * nave[1])
        d["dims_1"][2] = int(d['ncoil'])
        d["dims_1"][3] = int(1)

    d["header"] = hdr_head_only

    return d, data_vector_index
    


def _parse_parameter_line(d, line, nmix, exclude_labels=''):    
    """
    This parses a *parameter* line, not a data vector index line
    
    """
    line = line.split()
    if len(line) >= 1 and line[0] == '.':
        label = line[4]
        imix  = _convert_string(line[1])
        
        if label in exclude_labels:
            return
        
        if label not in d.keys():
            d[label] = ['',] * nmix
            
        if len(line) == 7:
            d[label][imix] = _convert_string(line[6])
        elif len(line) == 8: 
            d[label][imix] = [_convert_string(line[6]),_convert_string(line[7])]
        else:
            msg = "More than 2 values in header parameter = "+label
            raise util_exceptions.IncompleteHeaderParametersError(msg)

    return


def _parse_data(d, data_vector_index, data_filename):
    """
    Reads data from a .data file and returns it in a numpy array. Note that
    the array is unshaped (one dimensional).

    From .data/.list header -  data vector indices
       0    1    2    3     4     5     6      7      8     9  10  11   12    13   14   15    16   17   18   19     20
     (typ, mix, dyn, card, echo, loca, chan, extr1, extr2, kx, ky, kz, aver, sign, rf, grad, enc, rtop, rr, size, offset)


    """
    format = "f"
    fsize  = struct.calcsize(format)
    ncoil  = d['ncoil']
    
    indx_noise    = []
    indx_dynamic  = []
    indx_spectrum = []
    
    for line in data_vector_index:
        if line[0] == 'NOI':
            indx_noise.append(line)
        if line[0] == 'STD':
            if line[1] == 1:
                indx_dynamic.append(line)
            elif line[1] == 0:
                indx_spectrum.append(line)
   
    raw = open(data_filename, "rb").read()
   
    # Read in STD SPECTRUM FIDS

    data_spectrum = np.ndarray(d["dims_0"][::-1], dtype=np.complex64)
    for i, line in enumerate(indx_spectrum):
        offset = line[20]
        dsize  = line[19]
        count = dsize // fsize                  # element count
        data = raw[offset:offset+dsize]
        data = struct.unpack(format * count, data)
        data = util_fileio.collapse_complexes(data)
        data_spectrum[0, line[6], int(i//ncoil), :] = data

    # Read in NOISE FIDS
    
    data_noise = None
    if indx_noise != []:
        
        shape0_noise = list(d["dims_0"][::-1])
        shape0_noise[2] = len(indx_noise) // d['ncoil']        
        data_noise = np.ndarray(shape0_noise, dtype=np.complex64)
        for i, line in enumerate(indx_noise):
            offset = line[20]
            dsize  = line[19]
            count = dsize // fsize                  # element count
            data = raw[offset:offset+dsize]
            data = struct.unpack(format * count, data)
            data = util_fileio.collapse_complexes(data)
            data_noise[0, line[6], int(i//ncoil), :] = data

    # Read in DYNAMIC FIDS
    
    data_dynamic = None    
    if indx_dynamic != []:

        data_dynamic = np.ndarray(d["dims_1"][::-1], dtype=np.complex64)
        for i,line in enumerate(indx_dynamic):
            offset = line[20]
            dsize  = line[19]
            count = dsize // fsize                  # element count
            data = raw[offset:offset+dsize]
            data = struct.unpack(format * count, data)
            data = util_fileio.collapse_complexes(data)
            data_dynamic[0, line[6], int(i//ncoil), :] = data
        
    data_spectrum = np.conjugate(data_spectrum)

    return data_spectrum, data_noise, data_dynamic


    
    
#------------------------------------------------------------

def _test():

    fpath = "C:\\Users\\bsoher\\projects\\2015_gulin_BRP\\cross_platform_development\\philips_platform\\testdata\\2019_04_18_sandeep_press_braino_various_formats\\RAWDATA_DATALIST\\"
    #fpath = "D:\\Users\\bsoher\\projects\\2015_gulin_BRP\\cross_platform_development\\philips_platform\\testdata\\2019_04_18_sandeep_press_braino_various_formats\\RAWDATA_DATALIST\\"
    #fpath = "D:\Users\\bsoher\\code\\repository_svn\\sample_data\\philips_multiple_formats_braino\\RAWDATA_DATALIST\\"
    fbase = "raw_008.data"
    
    fname = fpath+fbase

    fobj = RawReaderPhilipsDataList()

    r = fobj.read_raw(fname)
    
    bob = 10


if __name__ == '__main__':
    _test()    