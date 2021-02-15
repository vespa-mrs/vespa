"""
Routines for reading a Philips .spar/.sdat formats and returning an 
DataRaw object populated with the file's data.

"""

# Python modules
import re
import os.path

# 3rd party modules
import numpy as np
from scipy.spatial.transform import Rotation

# Our modules
import vespa.common.mrs_data_raw as mrs_data_raw
import vespa.common.configobj as configobj
import vespa.analysis.fileio.raw_reader as raw_reader

import vespa.analysis.fileio.util_exceptions as util_exceptions


_EOL_REGEX = re.compile(r"(?:\r\n)|\r|\n")


def get_filename_pair(filename):
    """
    Given one *.spar/*.sdat filename, returns tuple with both filenames
    It doesn't matter if the filename is a fully qualified path or not.
    - one assumption, the extension are either all caps or all lower

    """
    spar_filename = sdat_filename = filename[:-3]
    if filename[-1:].isupper():
        sdat_filename += 'DAT'
        spar_filename += 'PAR'
    else:
        sdat_filename += 'dat'
        spar_filename += 'par'

    return (spar_filename, sdat_filename)



class RawReaderPhilipsSpar(raw_reader.RawReader):

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
        if not os.path.isfile(spar_filename):
            msg += "SPAR parameter file not found - '%s'\n" % spar_filename
        if  not os.path.isfile(sdat_filename):
            msg += "SDAT data file not found - '%s'\n" % sdat_filename
        if msg:
            raise util_exceptions.FileNotFoundError(msg)
        
        spar_str = open(spar_filename, "r").read()
        d = _get_parameters(spar_str, sdat_filename)

        d["data_source"] = filename

        if not (d["data"].shape[0:2] == (1,1)):
            msg = "SDAT data is unknown size, shape = %i" % d["data"].shape
            raise util_exceptions.IncorrectDimensionalityError(msg)

        if d["data"].shape[2] != 1:         # not Fidsum, but multiple rows,
            data = d["data"]                # so 'load into screen'
            raws = []
            for i in range(d["data"].shape[2]):
                d["data"] = data[0,0,i,:]
                d["data"].shape = 1,1,1,d["data"].shape[0]
                raws.append(mrs_data_raw.DataRaw(d))
        else:
            raws = [mrs_data_raw.DataRaw(d),]

        return raws
    
    
####################    Internal functions start here     ###############


def _get_parameters(header, sdat_filename):
    """ Given a SPAR string and SDAT filename, extract parameters needed into a dict. """

    # Massage SPAR format into pseudo-INI file for ConfigObj. SPAR uses
    # exclamation points for comments, colon for key/value delimiter
    header = _EOL_REGEX.sub("\n", header)  # normalize 'newlines'
    header = header.split("\n")
    header = [line.replace('!', '#',   1) for line in header]
    header = [line.replace(':', ' = ', 1) for line in header]
    header = configobj.ConfigObj(header)

    # parse data -----------------------------------------------

    data = open(sdat_filename, "rb").read()
    data = _vax_to_ieee_single_float(data)
    data_iter = iter(data)
    data = [complex(r, i) for r, i in zip(data_iter, data_iter)]
    complex_data = np.fromiter(data, np.complex64)
    complex_data = np.conjugate(complex_data)
    complex_data.shape = [1, 1, int(header["rows"]), int(header["samples"])]

    # parse transform and other MRS parameters -----------------

    try:
        angle_lr = float(header['lr_angulation'])
        angle_ap = float(header['ap_angulation'])
        angle_hf = float(header['cc_angulation'])
        dim_lr   = float(header['lr_size'])
        dim_ap   = float(header['ap_size'])
        dim_hf   = float(header['cc_size'])
        shift_lr = float(header['lr_off_center'])
        shift_ap = float(header['ap_off_center'])
        shift_hf = float(header['cc_off_center'])

        voxel_size = np.array([dim_lr, dim_ap, dim_hf])

        tform = np.zeros((4,4))
        scaling_mat = np.diag([dim_lr,dim_ap,dim_hf])
        rot = Rotation.from_euler('xyz', [-angle_lr,-angle_ap,angle_hf], degrees=True)
        tform[0:3,0:3] = rot.as_matrix() @ scaling_mat
        tform[3,3] = 1.0
        tform[0:3,3] = [-shift_lr, -shift_ap, shift_hf]
    except:
        # this will trigger default
        voxel_size = [20.0, 20.0, 20.0]
        tform = None

    hdr = str(header)
    hdr = hdr.split(',')
    hdr = '\n'.join(hdr)

    params = {'sw'              : float(header["sample_frequency"]),
              'frequency'       : float(header["synthesizer_frequency"])/1000000.0,
              'resppm'          : 4.7,
              'echopeak'        : 0.0,
              'nucleus'         : header["nucleus"],
              'seqte'           : float(header["echo_time"]),
              'seqtr'           : float(header["repetition_time"]),
              'voxel_dimensions' : voxel_size,
              'header'          : hdr,
              'transform'       : tform,
              'data'            : complex_data}

    return params


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

        byte2 = data[0 + i * 4]
        byte1 = data[1 + i * 4]
        byte4 = data[2 + i * 4]
        byte3 = data[3 + i * 4]

        # hex 0x80 = binary mask 10000000
        # hex 0x7f = binary mask 01111111

        sign = (byte1 & 0x80) >> 7
        expon = ((byte1 & 0x7f) << 1) + ((byte2 & 0x80) >> 7)
        fract = ((byte2 & 0x7f) << 16) + (byte3 << 8) + byte4

        if sign == 0:
            sign_mult = 1.0
        else:
            sign_mult = -1.0

        if 0 < expon:
            # note 16777216.0 == 2^24
            val = sign_mult * (0.5 + (fract / 16777216.0)) * pow(2.0, expon - 128.0)
            f.append(val)
        elif expon == 0 and sign == 0:
            f.append(0)
        else:
            f.append(0)
            # may want to raise an exception here ...

    return f