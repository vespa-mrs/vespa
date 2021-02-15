"""
Routines for reading a VASF format .rsp/.rsd file pair and returning an 
DataRaw object populated with the files' data.

"""
# Python modules
import re
import struct
import os.path

# 3rd party modules
import numpy as np

# Our modules
import vespa.analysis.fileio.raw_reader as raw_reader
import vespa.analysis.fileio.util_exceptions as util_exceptions
import vespa.common.constants as constants
import vespa.common.mrs_data_raw as mrs_data_raw
import vespa.common.configobj as configobj
import vespa.common.util.fileio as util_fileio
from vespa.common.base_transform import transformation_matrix
from vespa.common.mrs_data_raw import DataRawFidsum
from vespa.common.constants import Deflate

_EOL_REGEX = re.compile(r"(?:\r\n)|\r|\n")




def get_filename_pair(filename):
    """
    Given the name of a VASF data file (*.rsd) or parameter file (*.rsp) return
    a tuple of (parameters_filename, data_filename). It doesn't matter if the
    filename is a fully qualified path or not.
    - assumes extensions are all caps or all lower

    """
    param_filename = data_filename = filename[:-1]
    if filename[-1:].isupper():
        data_filename += 'D'
        param_filename += 'P'
    else:
        data_filename += 'd'
        param_filename += 'p'

    return (param_filename, data_filename)


class RawReaderVasf(raw_reader.RawReader):

    def __init__(self):

        raw_reader.RawReader.__init__(self)
        self.filetype_filter = "Spectra (*.rsd)|*.rsd;"
        self.multiple = True
        
        
    def read_raw(self, filename, ignore_data=False, *args, **kwargs):
        """
        Given the name of a .rsp or .rsd file, returns an DataRaw object
        populated with the parameters and data represented by the file pair.
        'ignore_data' has no effect on this parser
        'open_dataset' does not affect this parser

        """
        fparam, fdata = get_filename_pair(filename)

        if not os.path.isfile(fparam):
            msg = "Parameters file not found - '%s'" % fparam
            raise util_exceptions.FileNotFoundError(msg)

        if not ignore_data and not os.path.isfile(fdata):
            msg = "Data file not found - '%s'" % fdata
            raise util_exceptions.FileNotFoundError(msg)

        # Read the RSP file and extract the stuff I need.
        header = open(fparam, "rb").read()
        header = header.decode('utf-8')
        header = _EOL_REGEX.sub("\n", header)  # normalize 'newlines'
    
        d = _get_parameters(header, fdata)
    
        d["header"]      = header
        d["data_source"] = fdata

        return [mrs_data_raw.DataRaw(d),]



####################    Internal functions start here     ###############

def _get_parameters(header, fdata):
    """ Given the text of a .rsp file, extracts a subset of the parameters """

    # Massage header into pseudo-INI file for ConfigObj()
    header = header.split("\n")
    header = [line.replace(';', '#', 1) for line in header]
    header = configobj.ConfigObj(header)

    file_ = header['FILE INFORMATION']
    meas  = header['MEASUREMENT INFORMATION']
    data  = header['DATA INFORMATION']

    # parse data - has to be the original byte array --------------------------

    data_shape = [int(data["data_size_partition"]),
                  int(data["data_size_column"]),
                  int(data["data_size_line"]),
                  int(data["data_size_spectral"])]

    data_type = constants.DataTypes.any_type_to_numpy(file_["data_type"])
    data_type = constants.DataTypes.any_type_to_internal(data_type)
    is_xdr    = (file_.get("data_format", "xdr").lower() == "xdr")
    complex_data = open(fdata, "rb").read()
    element_count = len(complex_data) // constants.DataTypes.XDR_TYPE_SIZES[data_type]

    if is_xdr:
        complex_data = util_fileio.decode_xdr(complex_data, data_type, element_count)
    else:
        type_hash = {constants.DataTypes.COMPLEX64: "ff",
                     constants.DataTypes.COMPLEX128: "dd",
                     constants.DataTypes.FLOAT64: "d",
                     constants.DataTypes.FLOAT32: "f",
                     constants.DataTypes.INT32: "l",
                     constants.DataTypes.BYTE: "B", }
        try:
            complex_data = struct.unpack(type_hash[data_type] * element_count, complex_data)
        except struct.error:
            msg = "Unexpected input encountered while reading raw data"
            raise util_fileio.UnreadableDataError(msg)

        if constants.DataTypes.is_complex(data_type):
            data_iter = iter(complex_data)
            complex_data = [complex(r, i) for r, i in zip(data_iter, data_iter)]

    data_type = constants.DataTypes.any_type_to_numpy(data_type)
    complex_data = np.fromiter(complex_data, data_type)
    complex_data.shape = data_shape

    # parse transform and other MRS parameters -----------------

    freq = 127.9  if "frequency"   not in meas else float(meas["frequency"]) * 1e-6
    sw   = 2500.0 if "sweep_width" not in meas else float(meas["sweep_width"])
    tr   = 2000.0 if "repetition_time_1" not in meas else float(meas["repetition_time_1"])

    # seqte
    te = 0.0
    if "echo_time" in meas:
        te = float(meas["echo_time"])
    elif "echo_time_1" in meas:
        te = float(meas["echo_time_1"])

    try:
        voxel_size = [float(data["image_dimension_line"]),
                      float(data["image_dimension_column"]),
                      float(data["image_dimension_partition"])]

        voi_position = np.array([float(data["image_position_coronal"]),
                                 float(data["image_position_sagittal"]),
                                 float(data["image_position_transverse"])])

        row_vector = np.array([float(meas['image_normal_coronal']),
                               float(meas['image_normal_sagittal']),
                               float(meas['image_normal_transverse'])])

        col_vector = np.array([float(meas['image_column_coronal']),
                               float(meas['image_column_sagittal']),
                               float(meas['image_column_transverse'])])

        tform = transformation_matrix(row_vector, col_vector, voi_position, voxel_size)

    except:
        # this will trigger default
        voxel_size = [20.0, 20.0, 20.0]
        tform = None


    params = {'sw'          : sw,
              'frequency'   : freq,
              'resppm'      : 4.7,
              'echopeak'    : 0.0,
              'nucleus'     : meas["nucleus"],
              'seqte'       : te,
              'seqtr'       : tr,
              'voxel_dimensions' : voxel_size,
              'header'      : '',                       # set in read_raw
              'transform'   : tform,
              'data'        : complex_data}

    return params


class RawReaderVasfFidsum(RawReaderVasf):
    """ Read multiple VASF files into a DataRawFidsum object """

    def __init__(self):
        RawReaderVasf.__init__(self)

    def read_raw(self, filename, ignore_data=False, *args, **kwargs):
        raw = super().read_raw(filename, ignore_data)[0]
        raw = DataRawFidsum(raw.deflate(Deflate.DICTIONARY))
        return [raw, ]
