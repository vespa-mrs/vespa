# Python modules
import os
import struct

# 3rd party modules
import pydicom
import numpy as np
from scipy.spatial.transform import Rotation

# Our modules
import vespa.analysis.fileio.raw_reader as raw_reader
import vespa.common.util.config as util_config
import vespa.common.util.misc as util_misc
from vespa.analysis.fileio.util_exceptions import IncompleteHeaderParametersError
from vespa.common.mrs_data_raw import DataRaw, DataRawFidsum
from vespa.common.constants import Deflate

# need for inline processing - no wx
try:
    import wx       # just a test here
    from vespa.analysis.fileio.dicom_browser_dialog import PhilipsMrsBrowser
except:
    PhilipsMrsBrowser = None

"""
Here is a primer on how data is taken and what a 'dynamic' means ...

for kD = 1:nDynamics

    if start_cycles > 0 
        for kStart = 1:start_cycles
            do startup acq - no data - aka 'dummy scans' from GE
        end
    end
    
    if spectral_corr == YES
        for kSpectral = 1:spectral_corr_NSA
            do water acq - these are the so called 'dynamic water scans'
        end
    end
    water_data = average(do water acq)
    
    for kNSA = 1:NSA       # number of signal averages = NSA
        do water-suppressed acq
    end
    meta_data = average(do water-suppressed acq, maybe apply water_data ECC here too?)
end

So if NSA = 4 and nDynamics = 16, we would have 16 FIDs saved to file, but 
64 total averages for SNR, and for total acq time we add TR * (start_cycles + NSA*nDynamics)
"""




class RawReaderDicomPhilips(raw_reader.RawReader):
    """
    Reads a Philips DICOM file into an DataRaw object.
    It implements the interface defined by raw_reader.RawReader (q.v.).

    """
    def __init__(self):
        raw_reader.RawReader.__init__(self)


    def pickfile(self, default_path=""):
        """ 
        The default here is to allow multiple filenames. Each will be treated as
        a separate MRS file and 'loaded into the screen'.
        
        """
        if not os.path.exists(default_path):
            default_path = util_config.VespaConfig()["general"].get("last_dicom_browse_path", "")

        if not os.path.exists(default_path):
            default_path = util_misc.get_documents_dir()

        if PhilipsMrsBrowser is not None:

            dialog = PhilipsMrsBrowser(multi_select=self.multiple,
                                       default_path=default_path,
                                       show_tags=False,
                                       preview_size=None)
            dialog.ShowModal()
            self.filenames = dialog.filenames
        else:
            self.filenames = []

        if self.filenames:
            config = util_config.VespaConfig()
            config["general"]["last_dicom_browse_path"] = dialog.path
            config.write()

        return bool(self.filenames)


    def read_raw(self, filename, ignore_data=False, *args, **kwargs):
        """
        Given Philips DICOM filename, return populated DataRaw object
        Ignore data has no effect on this parser

        """
        dataset = pydicom.dcmread(filename)

        if not _is_mrs_dicom(dataset):
            raise ValueError("Dataset does not have MRS data.")

        d = _get_parameters_philips_proprietary(dataset)

        d["data_source"] = filename

        return [DataRaw(d),]



class RawReaderDicomPhilipsFidsum(RawReaderDicomPhilips):
    """ Read multiple Philips DICOMs file into a DataRawFidsum object """

    def __init__(self):
        RawReaderDicomPhilips.__init__(self)

    def read_raw(self, filename, ignore_data=False, *args, **kwargs):
        """ Call base class then convert DataRaw into DataRawFidsum object. """
        raw = super().read_raw(filename, ignore_data)[0]
        raw = DataRawFidsum(raw.deflate(Deflate.DICTIONARY))
        return [raw,]



####################    Internal functions start here     ###############

def _get_parameters_philips_proprietary(dataset):
    """ Returns a subset of the parameters from a Pydicom dataset """

    ds = dataset

    if (0x0028,0x9001) in ds:
        averages = ds[0x0028,0x9001].value
    elif (0x2005,0x140f) in ds:
        if "DataPointRows" in ds[0x2005,0x140f].value[0]:
            averages = ds[0x2005,0x140f].value[0].DataPointRows
    elif (0x2001,0x1081) in ds:
        averages = ds[0x2001,0x1081].value     # maybe?  number of dynamic scans
    else:
        averages = 1

    if (0x2005,0x1315) in ds:
        spectral_points = ds[0x2005,0x1315].value
    elif (0x2005,0x140f) in ds:
        if "DataPointColumns" in ds[0x2005,0x140f].value[0]:
            spectral_points = ds[0x2005,0x140f].value[0].DataPointColumns
    elif (0x0028,0x9002) in ds:
        spectral_points = ds[0x0028,0x9002].value
    else:
        raise IncompleteHeaderParametersError("spectral_points")

    data = ds[0x2005, 0x1270].value

    if isinstance(data, (bytes, bytearray)):
        # Big simplifying assumptions --
        # 0) Unpack byte stream into a series of float32 values
        # 1) Floats a series of complex numbers organized as ririri...
        # 2) Data is little endian.
        data = struct.unpack("<%df" % (len(data) / 4), data)

    data_iter = iter(data)
    data = [complex(r, i) for r, i in zip(data_iter, data_iter)]
    complex_data = np.fromiter(data, np.complex64)
    complex_data = complex_data.conjugate() # empirical
    complex_data.shape = (1, 1, int(averages), int(spectral_points))

    # parse transform and other MRS parameters -----------------

    tr = float(ds.RepetitionTime) if "RepetitionTime" in ds else float(ds[0x2005,0x1030][0])

    # Thoughts on how to find VoI 'orientation'
    # - (2005,1566) is string indicating svs 90 pulse, 'slice slab' orientation
    # - e.g. 'FH' means slice slab is Axial, 180s are Sag/Cor for PRESS
    # - (0018, 5100) is Patient Position, e.g. 'HFS' for how they lie on table

    try:
        section  = ds[0x2005,0x1085][0]
        angle_lr = float(section[0x02005, 0x1056].value)
        angle_ap = float(section[0x02005, 0x1054].value)
        angle_hf = float(section[0x02005, 0x1055].value)
        dim_lr   = float(section[0x02005, 0x1059].value)
        dim_ap   = float(section[0x02005, 0x1057].value)
        dim_hf   = float(section[0x02005, 0x1058].value)
        shift_lr = float(section[0x02005, 0x105c].value)
        shift_ap = float(section[0x02005, 0x105a].value)
        shift_hf = float(section[0x02005, 0x105b].value)

        voxel_size = np.array([dim_lr,dim_ap,dim_hf])

        tform = np.zeros((4,4))
        scaling_mat = np.diag([dim_lr,dim_ap,dim_hf])
        rot = Rotation.from_euler('xyz', [-angle_lr,-angle_ap,angle_hf], degrees=True)
        tform[0:3,0:3] = rot.as_matrix() @ scaling_mat
        tform[3,3] = 1.0
        tform[0:3,3] = [-shift_lr, -shift_ap, shift_hf]
    except:
        # this will trigger default
        voxel_size = [20.0,20.0,20.0]
        tform = None

    params = {'is_dicom_sop': False,
              'sw'          : float(ds[0x2005,0x1357].value),
              'frequency'   : float(ds[0x2001,0x1083].value),
              'resppm'      : 4.7,
              'echopeak'    : 0.0,
              'nucleus'     : ds[0x2001,0x1087].value,
              'seqte'       : float(ds[0x2001,0x1025].value),
              'seqtr'       : tr,
              'voxel_dimensions' : voxel_size,
              'header'      : str(dataset),
              'transform'   : tform,
              'data'        : complex_data}

    return params


def _is_mrs_dicom(dataset):
    """ returns True if all criteria (Sandeep) are met """

    if not isinstance(dataset, pydicom.dataset.Dataset):
        raise ValueError("Object passed in not a dicom Dataset.")

    if not "ProtocolName" in dataset:
        return False

    if dataset.ProtocolName == 'ExamCard':
        return False

    if not (0x2005, 0x10c0) in dataset:
        return False

    if dataset[0x2005, 0x10c0].value != 'SPECTRO':
        return False

    if not (0x2005, 0x1270) in dataset:
        return False

    return True

