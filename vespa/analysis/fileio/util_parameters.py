# Python modules


# 3rd party modules


# Our modules
import vespa.common.constants as constants
import vespa.common.mrs_data_raw as mrs_data_raw


def map_parameters(file_, identification=None, measurement=None, data=None):
    """
    Given one or four parameter dicts (see below) that describe data (e.g. 
    from a VASF file), this function remaps a subset of the dict contents 
    into a new, single, flat dict suitable for passing to 
    mrs_data_raw.DataRaw().

    Callers from VASF code will probably want to pass four dicts; they represent
    the four sections of a VASF .rsp file ([FILE INFORMATION], 
    [IDENTIFICATION INFORMATION], [MEASUREMENT INFORMATION], 
    [DATA INFORMATION]).

    Callers from other code (e.g. dicom_siemens.py) will probably only have 
    one dict to pass that contains the same keys as the four VASF-style dicts.

    A better name for this function might be map_vasf-style_parameters() but
    I didn't use that for obvious reasons.
    """
    d = { }

    # If I was given only one dict, I use the same dict for all four "sections".
    if not identification:
        identification = file_
        measurement = file_
        data = file_

    # Data type (float32, complex64, etc.)
    if "data_type" in file_:
        d["data_type"] = file_["data_type"]
        d["data_type"] = constants.DataTypes.any_type_to_numpy(d["data_type"])

    # Whether or not data is stored in XDR format
    d["is_xdr"] = (file_.get("data_format", "xdr").lower() == "xdr")

    # Seq type (e.g. multislice_si)
    if "sequence_type" in identification:
        d["sequence_type"] = identification["sequence_type"]

    # Freq is assumed to be in Hz, we store it in Mhz
    if "frequency" in measurement:
        d["frequency"] = float(measurement["frequency"]) * 1e-6

    # sweep width
    if "sweep_width" in measurement:
        d["sw"] = float(measurement["sweep_width"])

    # Data dimensions
    d["dims"] = mrs_data_raw.DataRaw.DEFAULT_DIMS
    if "measure_size_spectral" in measurement:
        d["dims"][0] = measurement["measure_size_spectral"]
    if "measure_size_line" in measurement:
        d["dims"][1] = measurement["measure_size_line"]
    if "measure_size_column" in measurement:
        d["dims"][2] = measurement["measure_size_column"]
    if "measure_size_partition" in measurement:
        d["dims"][3] = measurement["measure_size_partition"]

    # If data_size_xxxx info is present, that trumps the dimensionality 
    # described by measure_size_xxxx. 
    # Our sample data file 4007_smh_spat.rsp is an example of this.
    if "data_size_spectral" in data:
        d["dims"][0] = data["data_size_spectral"]
    if "data_size_line" in data:
        d["dims"][1] = data["data_size_line"]
    if "data_size_column" in data:
        d["dims"][2] = data["data_size_column"]
    if "data_size_partition" in data:
        d["dims"][3] = data["data_size_partition"]

    d["dims"] = [int(dim) for dim in d["dims"]]

    # nucleus/isotope
    if "nucleus" in measurement:
        d["nucleus"] = measurement["nucleus"]

    # seqte
    if "echo_time" in measurement:
        d["seqte"] = float(measurement["echo_time"])
    elif "echo_time_1" in measurement:
        d["seqte"] = float(measurement["echo_time_1"])

    # Miscellany
    if "echo_position" in measurement:
        d["echopeak"] = float(measurement["echo_position"])

    if "flip_angle" in measurement:
        d["flip_angle"] = float(measurement["flip_angle"])

    if "slice_thickness" in data:
        d["slice_thickness"] = float(data["slice_thickness"])

    if "slice_orientation_pitch" in data:
        d["slice_orientation_pitch"] = data["slice_orientation_pitch"]

    if "slice_orientation_roll" in data:
        d["slice_orientation_roll"] = data["slice_orientation_roll"]


    # image position
    d["image_position"] = mrs_data_raw.DataRaw.DEFAULT_IMAGE_POSITION
    if "image_position_sagittal" in measurement:
        d["image_position"][0] = float(measurement["image_position_sagittal"])
    if "image_position_coronal" in measurement:
        d["image_position"][1] = float(measurement["image_position_coronal"])
    if "image_position_transverse" in measurement:
        d["image_position"][2] = float(measurement["image_position_transverse"])

    # image dimension
    d["image_dimension"] = mrs_data_raw.DataRaw.DEFAULT_IMAGE_DIMENSION
    if "image_dimension_sagittal" in data:
        d["image_dimension"][0] = float(data["image_dimension_sagittal"])
    if "image_dimension_coronal" in data:
        d["image_dimension"][1] = float(data["image_dimension_coronal"])
    if "image_dimension_transverse" in data:
        d["image_dimension"][2] = float(data["image_dimension_transverse"])

    # image orientation column
    d["image_orient_col"] = mrs_data_raw.DataRaw.DEFAULT_IMAGE_ORIENTATION_COLUMN
    if "image_column_sagittal" in measurement:
        d["image_orient_col"][0] = float(measurement["image_column_sagittal"])
    if "image_column_coronal" in measurement:
        d["image_orient_col"][1] = float(measurement["image_column_coronal"])
    if "image_column_transverse" in measurement:
        d["image_orient_col"][2] = float(measurement["image_column_transverse"])

    # image orientation row
    d["image_orient_row"] = mrs_data_raw.DataRaw.DEFAULT_IMAGE_ORIENTATION_ROW
    if "image_normal_sagittal" in measurement:
        if not measurement["image_normal_sagittal"] == '':
            d["image_orient_row"][0] = float(measurement["image_normal_sagittal"])
    if "image_normal_coronal" in measurement:
        if not measurement["image_normal_coronal"] == '':
            d["image_orient_row"][1] = float(measurement["image_normal_coronal"])
    if "image_normal_transverse" in measurement:
        if not measurement["image_normal_transverse"] == '':
            d["image_orient_row"][2] = float(measurement["image_normal_transverse"])

    measure_time = [0.0,]
    if "measure_time" in list(file_.keys()):
        if file_["measure_time"] != '':
            hms, frac = file_["measure_time"].split('.')
            measure_time = float(hms[0:2])*3600 + float(hms[2:4])*60 + float(hms[4:6]) + float('0.'+frac)
    elif "measure_time" in list(identification.keys()):
        if identification["measure_time"] != '':
            hms = identification["measure_time"]
            measure_time = float(hms[0:2])*3600 + float(hms[3:5])*60 + float(hms[6:8])
    d["measure_time"] = [measure_time,]

    return d    
