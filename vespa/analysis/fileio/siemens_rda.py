"""
Routines for reading a Siemens *.rda format and returning an 
DataRaw object populated with the file's data.

"""
# Python modules
import re
import struct

# 3rd party modules
import numpy as np

# Our modules
import vespa.common.mrs_data_raw as mrs_data_raw
import vespa.common.configobj as configobj
import vespa.common.util.fileio as util_fileio
import vespa.analysis.fileio.raw_reader as raw_reader 

from vespa.common.base_transform import transformation_matrix


_EOL_REGEX = re.compile(r"(?:\r\n)|\r|\n")




class RawReaderSiemensRda(raw_reader.RawReader):

    def __init__(self):

        raw_reader.RawReader.__init__(self)
        self.filetype_filter = "Spectra (*.rda)|*.rda"
        self.multiple = True
        

    def read_raw(self, filename, ignore_data=False, *args, **kwargs):
        """
        Given Siemens .RDA filename, return populated DataRaw object.
        - 'ignore_data' has no effect on this parser

        """
        rda_str = open(filename, "rb").read()
        d = _get_parameters(rda_str)
        d["data_source"] = filename

        return [mrs_data_raw.DataRaw(d),]
    
    
####################    Internal functions start here     ###############

def _get_parameters(rda_string):
    """
    Given a RDA string, extract parameters needed into a dict.
    Header starts/ends with known text markers. Data is right after header.

    """

    # parse header, decode to a string ----------------------------------------

    rda_str = rda_string.decode('latin-1')
    istr = rda_str.find(">>> Begin of header <<<")
    iend = rda_str.find(">>> End of header <<<")
    head_only = rda_str[istr+25:iend]    # save this to return to user

    # Massage RDA format into a pseudo-INI file for ConfigObj()
    head = _EOL_REGEX.sub("\n", head_only)  # normalize 'newlines'
    head = head.split("\n")
    head = [line.replace(':', ' = ', 1) for line in head]
    hdr  = configobj.ConfigObj(head)

    # parse data - has to be the original byte array --------------------------

    data = rda_string[iend+23:]
    element_count = len(data) // struct.calcsize("dd")
    try:
        data = struct.unpack("dd" * element_count, data)
    except struct.error:
        msg = "Unexpected input encountered while reading raw data"
        raise util_fileio.UnreadableDataError(msg)

    data_iter = iter(data)
    data = [complex(r, i) for r, i in zip(data_iter, data_iter)]
    complex_data = np.fromiter(data, np.complex64)
    complex_data.shape = [int(hdr["CSIMatrixSize[2]"]),
                          int(hdr["CSIMatrixSize[1]"]),
                          int(hdr["CSIMatrixSize[0]"]),
                          int(hdr["VectorSize"])]

    # parse transform and other MRS parameters -----------------

    try:
        voxel_size = [float(hdr["PixelSpacingCol"]),
                      float(hdr["PixelSpacingRow"]),
                      float(hdr["PixelSpacing3D"])]

        col_vector = np.array([float(hdr["image_normal_sagittal"]),
                               float(hdr["image_normal_coronal"]),
                               float(hdr["image_normal_transverse"])])

        row_vector = np.array([float(hdr['image_column_sagittal']),
                               float(hdr['image_column_coronal']),
                               float(hdr['image_column_transverse'])])

        voi_position = np.array([float(hdr['ColumnVector[0]']),
                                 float(hdr['ColumnVector[1]']),
                                 float(hdr['ColumnVector[2]'])])

        tform = transformation_matrix(row_vector, col_vector, voi_position, voxel_size)

    except:
        # this will trigger default
        voxel_size = [20.0, 20.0, 20.0]
        tform = None


    params = {'sw'          : 1.0/(float(hdr["DwellTime"])*0.000001),
              'frequency'   : float(hdr["MRFrequency"]),
              'resppm'      : 4.7,
              'echopeak'    : 0.0,
              'nucleus'     : hdr["Nucleus"],
              'seqte'       : float(hdr["TE"]),
              'seqtr'       : float(hdr["TR"]),
              'voxel_dimensions' : voxel_size,
              'header'      : head_only,
              'transform'   : tform,
              'data'        : complex_data}


    return params

