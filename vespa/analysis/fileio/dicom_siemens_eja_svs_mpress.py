# Python modules
import os
import struct
import pprint
pp = pprint.PrettyPrinter(depth=2)

# 3rd party modules
import pydicom
import pydicom.dicomio
import numpy as np

from pydicom.values import convert_numbers

# Our modules
import vespa.analysis.fileio.raw_reader as raw_reader
import vespa.common.util.config as util_config
import vespa.common.util.misc as util_misc
from vespa.common.base_transform import transformation_matrix
from vespa.common.mrs_data_raw import DataRawEdit, DataRawEditFidsum

from vespa.common.constants import Deflate

# need for inline processing - no wx
try:
    import wx   # just a test here
    from vespa.analysis.fileio.dicom_browser_dialog import SiemensMrsBrowser
except:
    SiemensMrsBrowser = None



# DICOM standard tags
TAG_SOP_CLASS_UID = (0x0008, 0x0016)



class RawReaderDicomSiemensXaEjaSvsMpress(raw_reader.RawReader):
    """ Read a single Siemens DICOM file into an DataRaw object. """

    def __init__(self):
        raw_reader.RawReader.__init__(self)
        self.multiple = False


    def pickfile(self, default_path=""):
        """ Default is multiple filenames that will be 'loaded into the screen' """

        if not os.path.exists(default_path):
            default_path = util_config.VespaConfig()["general"].get("last_dicom_browse_path", "")
        if not os.path.exists(default_path):
            default_path = util_misc.get_documents_dir()

        if SiemensMrsBrowser is not None:
            dialog = SiemensMrsBrowser(multi_select=self.multiple,
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


    def read_raw(self, filename, ignore_data, *args, **kwargs):
        """
        Given Siemens DICOM filename, return a populated DataRaw object

        The sop_class_uid flags if this is the older proprietary Siemens hack
        format or the newer DICOM standard MR Spectroscopy Storage object.

        Ignore data has no effect on this parser

        """

        # the .IMA format is a DICOM standard, but Siemens stores a lot of
        # information inside a private and very complicated header with its own
        # data storage format, we have to get that information out along with
        # the data. we start by reading in the DICOM file completely
        dataset = pydicom.dicomio.read_file(filename)

        sop_class_uid = pydicom.uid.UID(str(dataset['SOPClassUID'].value.upper()))

        if sop_class_uid.name == 'MR Spectroscopy Storage':
            d = _get_parameters_dicom_sop(dataset)
        else:
            raise ValueError('Error (RawReaderDicomSiemensFidsumXaEjaSvsMpress): Old style Siemens DICOM should not be here.')

        d["data_source"] = filename

        data = d['data']

        dat_on  = data[1, :, :, :].copy()
        dat_off = data[0, :, :, :].copy()
        dat_sum = (dat_on+dat_off).copy()
        dat_dif = (dat_on-dat_off).copy()

        # Create a DataRawEditFidsum objects for the four different states
        # of the data ON, OFF, SUM and DIFFERENCE.

        d["data"] = dat_on
        d["data_source"] = filename+'.EditON'
        raw_on = DataRawEdit(attributes=d)

        d["data"] = dat_off
        d["data_source"] = filename+'.EditOFF'
        raw_off = DataRawEdit(attributes=d)

        d["data"] = dat_sum
        d["data_source"] = filename+'.Sum'
        raw_sum = DataRawEdit(attributes=d)

        d["data"] = dat_dif
        d["data_source"] = filename+'.Diff'
        raw_dif = DataRawEdit(attributes=d)

        raws = [raw_on, raw_off, raw_sum, raw_dif]

        return raws


class RawReaderDicomSiemensXaEjaSvsMpressOnOff(RawReaderDicomSiemensXaEjaSvsMpress):
    """ Read a single Siemens DICOM file into an DataRaw object. """

    def __init__(self):
        raw_reader.RawReader.__init__(self)
        self.multiple = False


    def read_raw(self, filename, ignore_data, *args, **kwargs):
        """
        Given Siemens DICOM filename, return a populated DataRaw object

        The sop_class_uid flags if this is the older proprietary Siemens hack
        format or the newer DICOM standard MR Spectroscopy Storage object.

        Ignore data has no effect on this parser

        """

        # the .IMA format is a DICOM standard, but Siemens stores a lot of
        # information inside a private and very complicated header with its own
        # data storage format, we have to get that information out along with
        # the data. we start by reading in the DICOM file completely
        dataset = pydicom.dicomio.read_file(filename)

        sop_class_uid = pydicom.uid.UID(str(dataset['SOPClassUID'].value.upper()))

        if sop_class_uid.name == 'MR Spectroscopy Storage':
            d = _get_parameters_dicom_sop(dataset)
        else:
            raise ValueError('Error (RawReaderDicomSiemensFidsumXaEjaSvsMpress): Old style Siemens DICOM should not be here.')

        d["data_source"] = filename

        data = d['data']

        dat_on  = data[1, :, :, :].copy()
        dat_off = data[0, :, :, :].copy()

        # Create a DataRawEditFidsum objects for the four different states
        # of the data ON, OFF, SUM and DIFFERENCE.

        d["data"] = dat_on
        d["data_source"] = filename+'.EditON'
        raw_on = DataRawEdit(attributes=d)

        d["data"] = dat_off
        d["data_source"] = filename+'.EditOFF'
        raw_off = DataRawEdit(attributes=d)

        raws = [raw_on, raw_off]

        return raws



class RawReaderDicomSiemensFidsumXaEjaSvsMpress(RawReaderDicomSiemensXaEjaSvsMpress):
    """ Read multiple Siemens DICOMs file into a DataRawFidsum object """

    def __init__(self):
        RawReaderDicomSiemensXaEjaSvsMpress.__init__(self)
        self.multiple = True

    def read_raw(self, filename, ignore_data=False, get_header=True, *args, **kwargs):
        """
        Given Siemens DICOM filename, return list with a single dict containing
        essential data to create an MrsDataRaw or MrsDataEdit object

        The sop_class_uid flags an error if this is the older proprietary Siemens
        format

        Ignore data has no effect on this parser

        SOPClassID for spectro from https://dicom.nema.org/dicom/2013/output/chtml/part04/sect_B.5.html
        - sop_class_uid == '1.2.840.10008.5.1.4.1.1.4.2' is 'MR Spectroscopy Storage'

        """

        dataset = pydicom.dicomio.dcmread(filename)

        sop_class_uid = pydicom.uid.UID(dataset['SOPClassUID'].value)

        if sop_class_uid.name == 'MR Spectroscopy Storage' or sop_class_uid == '1.2.840.10008.5.1.4.1.1.4.2':
            d = _get_parameters_dicom_sop(dataset, get_header=get_header)
        else:
            raise ValueError('Error (RawReaderDicomSiemensFidsumXaEjaSvsMpress): Old style Siemens DICOM should not be here.')

        d["data_source"] = filename

        return [d, ]


    def read_raws(self, ignore_data=False, *args, **kwargs):
        """
        Calls read_raw() once for each filename in self.filenames.
        - returns list with four DataRawEditFidsum objects for ON, OFF, SUM and DIFF
        - read_raw() method should return one FID data set per file.

        """
        self.raws_on = []
        self.raws_off = []

        nfiles = len(self.filenames)
        get_header = True

        for i, fname in enumerate(self.filenames):
            get_header = (i==0) or (i==int(nfiles/2))
            raw = self.read_raw(fname, ignore_data, get_header=get_header, *args, **kwargs)
            if i < int(nfiles/2):
                self.raws_off.append(raw[0])
            else:
                self.raws_on.append(raw[0])

        d = self.raws_on[0]
        d["data_source"] = d["data_source"]+'.EditON'
        raw_on = DataRawEditFidsum(d)
        for d in self.raws_on[1:]:
            raw_on.concatenate(DataRawEditFidsum(d))

        d = self.raws_off[0]
        d["data_source"] = d["data_source"]+'.EditOFF'
        raw_off = DataRawEditFidsum(d)
        for d in self.raws_off[1:]:
            raw_off.concatenate(DataRawEditFidsum(d))

        dat_on  = raw_on.data.copy()
        dat_off = raw_off.data.copy()

        dat_sum = dat_on + dat_off
        dat_dif = dat_on - dat_off

        d = self.raws_on[0]
        d["data"] = dat_sum
        d["data_source"] = d["data_source"]+'.SumIndiv'
        raw_sum = DataRawEditFidsum(d)

        d = self.raws_on[0]
        d["data"] = dat_dif
        d["data_source"] = d["data_source"]+'.DiffIndiv'
        raw_dif = DataRawEditFidsum(d)

        self.raws = [raw_on, raw_off, raw_sum, raw_dif]

        self._check_consistency(fidsum=True)
        self._check_for_si()
        if 'open_dataset' in list(kwargs.keys()):
            if kwargs['open_dataset'] is not None:
                self._check_compatibility(self.raws, kwargs['open_dataset'], fidsum=True)

        return self.raws


class RawReaderDicomSiemensFidsumXaEjaSvsMpressOnOff(RawReaderDicomSiemensFidsumXaEjaSvsMpress):
    """ Read multiple Siemens DICOMs file into a DataRawFidsum object """

    def __init__(self):
        RawReaderDicomSiemensXaEjaSvsMpress.__init__(self)
        self.multiple = True


    def read_raws(self, ignore_data=False, *args, **kwargs):
        """
        Calls read_raw() once for each filename in self.filenames.
        - returns list with four DataRawEditFidsum objects for ON, OFF, SUM and DIFF
        - read_raw() method should return one FID data set per file.

        """
        self.raws_on = []
        self.raws_off = []

        nfiles = len(self.filenames)
        get_header = True

        for i, fname in enumerate(self.filenames):
            get_header = (i==0) or (i==int(nfiles/2))
            raw = self.read_raw(fname, ignore_data, get_header=get_header, *args, **kwargs)
            if i < int(nfiles/2):
                self.raws_off.append(raw[0])
            else:
                self.raws_on.append(raw[0])

        d = self.raws_on[0]
        d["data_source"] = d["data_source"]+'.EditON'
        raw_on = DataRawEditFidsum(d)
        for d in self.raws_on[1:]:
            raw_on.concatenate(DataRawEditFidsum(d))

        d = self.raws_off[0]
        d["data_source"] = d["data_source"]+'.EditOFF'
        raw_off = DataRawEditFidsum(d)
        for d in self.raws_off[1:]:
            raw_off.concatenate(DataRawEditFidsum(d))

        self.raws = [raw_on, raw_off]

        self._check_consistency(fidsum=True)
        self._check_for_si()
        if 'open_dataset' in list(kwargs.keys()):
            if kwargs['open_dataset'] is not None:
                self._check_compatibility(self.raws, kwargs['open_dataset'], fidsum=True)

        return self.raws


class RawReaderDicomSiemensFidsumXaEjaSvsMpressOnOffAdv(RawReaderDicomSiemensFidsumXaEjaSvsMpress):
    """ Read multiple Siemens DICOMs file into a DataRawFidsum object """

    def __init__(self):
        RawReaderDicomSiemensXaEjaSvsMpress.__init__(self)
        self.multiple = True


    def read_raws(self, ignore_data=False, *args, **kwargs):
        """
        Calls read_raw() once for each filename in self.filenames.
        - returns list with four DataRawEditFidsum objects for ON, OFF, SUM and DIFF
        - read_raw() method should return one FID data set per file.

        """
        self.raws_on = []
        self.raws_off = []

        nfiles = len(self.filenames)
        get_header = True

        raws = []
        raw_avgs = 0

        for i, fname in enumerate(self.filenames):
            get_header = (i==0) or (i==int(nfiles/2))
            raw = self.read_raw(fname, ignore_data, get_header=get_header, *args, **kwargs)
            raw_avgs += raw[0]['data'].shape[0]
            raws.append(raw)

        tmp_raw = raws[0][0].copy()

        avg_cnt = 0
        for raw in raws:
            for j in range(raw[0]['data'].shape[0]):

                tmp_raw['data'] = raw[0]['data'][j,0,0,:].copy()
                tmp_raw['data'] = tmp_raw['data'].conjugate()
                tmp_raw['data'].shape = 1,1,1,len(tmp_raw['data'])

                if avg_cnt < int(raw_avgs/2):
                    self.raws_off.append(tmp_raw.copy())
                else:
                    self.raws_on.append(tmp_raw.copy())
                avg_cnt += 1

        d = self.raws_on[0]
        d["data_source"] = d["data_source"]+'.EditON'
        raw_on = DataRawEditFidsum(d)
        for d in self.raws_on[1:]:
            raw_on.concatenate(DataRawEditFidsum(d))

        d = self.raws_off[0]
        d["data_source"] = d["data_source"]+'.EditOFF'
        raw_off = DataRawEditFidsum(d)
        for d in self.raws_off[1:]:
            raw_off.concatenate(DataRawEditFidsum(d))

        self.raws = [raw_on, raw_off]

        self._check_consistency(fidsum=True)
        self._check_for_si()
        if 'open_dataset' in list(kwargs.keys()):
            if kwargs['open_dataset'] is not None:
                self._check_compatibility(self.raws, kwargs['open_dataset'], fidsum=True)

        return self.raws




####################    Internal functions start here     ###############

CSA1 = 0
CSA2 = 1

ima_types = {
    "floats": ["NumberOfAverages", "RSatPositionSag", "PercentPhaseFieldOfView",
               "RSatOrientationSag", "MixingTime", "PercentPhaseFieldOfView",
               "RSatPositionCor", "InversionTime", "RepetitionTime",
               "VoiThickness", "TransmitterReferenceAmplitude",
               "ImageOrientationPatient", "SliceThickness",
               "RSatOrientationTra", "PixelBandwidth", "SAR",
               "PixelSpacing", "ImagePositionPatient", "VoiPosition",
               "SliceLocation","FlipAngle", "VoiInPlaneRotation", "VoiPhaseFoV",
               "SliceMeasurementDuration", "HammingFilterWidth",
               "RSatPositionTra", "MagneticFieldStrength", "VoiOrientation",
               "PercentSampling", "EchoTime", "VoiReadoutFoV", "RSatThickness",
               "RSatOrientationCor", "ImagingFrequency", "TriggerTime",
               "dBdt", "TransmitterCalibration", "PhaseGradientAmplitude",
               "ReadoutGradientAmplitude", "SelectionGradientAmplitude",
               "GradientDelayTime", "dBdt_max", "t_puls_max", "dBdt_thresh",
               "dBdt_limit", "SW_korr_faktor", "Stim_lim", "Stim_faktor"],
    "integers": ["Rows", "Columns", "DataPointColumns",
                 "SpectroscopyAcquisitionOut-of-planePhaseSteps",
                 "EchoPartitionPosition", "AcquisitionMatrix",
                 "NumberOfFrames", "EchoNumbers", "RealDwellTime",
                 "EchoTrainLength", "EchoLinePosition",
                 "EchoColumnPosition", "SpectroscopyAcquisitionDataColumns",
                 "SpectroscopyAcquisitionPhaseColumns",
                 "SpectroscopyAcquisitionPhaseRows", "RfWatchdogMask",
                 "NumberOfPhaseEncodingSteps", "DataPointRows",
                 "UsedPatientWeight", "NumberOfPrescans",
                 "Stim_mon_mode", "Operation_mode_flag", "CoilId",
                 "MiscSequenceParam", "MrProtocolVersion",
                 "ProtocolSliceNumber"],
    "strings": ["ReferencedImageSequence", "ScanningSequence", "SequenceName",
                "ImagedNucleus", "TransmittingCoil", "PhaseEncodingDirection",
                "VariableFlipAngleFlag", "SequenceMask",
                "AcquisitionMatrixText", "MultistepIndex",
                "DataRepresentation", "SignalDomainColumns",
                "k-spaceFiltering", "ResonantNucleus",
                "ImaCoilString", "FrequencyCorrection",
                "WaterReferencedPhaseCorrection", "SequenceFileOwner",
                "CoilForGradient", "CoilForGradient2",
                "PositivePCSDirections", ],
}

def _read_csa_header(csa_header_bytes):
    # two possibilities exist here, either this is a CSA2 format beginning with
    # an SV10 string, or a CSA1 format which doesn't. in CSA2 after the "SV10"
    # are four junk bytes, then the number of tags in a uint32 and a delimiter
    # uint32 containing the value 77. in CSA1 there is just the number of tags
    # and the delimiter. after that the two formats contain the same structure
    # for each tag, but the definition of the size of the items in each tag is
    # different between the two versions

    if csa_header_bytes[:4] == "SV10".encode('latin-1'):
        num_tags, delimiter = struct.unpack("<II", csa_header_bytes[8:16])
        header_offset = 16
        header_format = CSA2
    else:
        num_tags, delimiter = struct.unpack("<II", csa_header_bytes[:8])
        header_offset = 8
        header_format = CSA1
    # now we can iteratively read the tags and the items inside them
    csa_header = {}
    for i in range(num_tags):
        the_bytes = csa_header_bytes[header_offset:(header_offset + 84)]
        r = struct.unpack("<64si4siii", the_bytes)
        name, vm, vr, syngo_dt, nitems, delimiter = r
        header_offset += 84

        # the name of the tag is 64 bytes long, but the string we want is
        # null-terminated inside, so extract the real name by taking only bytes
        # up until the first 0x00
        name = name.decode('latin-1')
        name = name.split("\x00", 1)[0]
        # read all the items inside this tag
        item_list = []
        for j in range(nitems):
            sizes = struct.unpack("<4L", csa_header_bytes[header_offset:(header_offset + 16)])
            header_offset += 16
            if header_format == CSA2:
                item_length = sizes[1]
                if (header_offset + item_length) > len(csa_header_bytes):
                    item_length = len(csa_header_bytes) - header_offset
            elif header_format == CSA1:
                item_length = sizes[0]
            item_bytes = csa_header_bytes[header_offset:(header_offset + item_length)]
            item, = struct.unpack("<%ds" % item_length, item_bytes)
            item = item.decode('latin-1')
            item = item.split("\x00", 1)[0]
            if item_length > 0:
                if name in ima_types["floats"]:
                    item = float(item)
                elif name in ima_types["integers"]:
                    item = int(item)
                elif name in ima_types["strings"]:
                    pass
                else:
                    pass
                    # warnings.warn("Unhandled name {0} with vr {1} and value {2}".format(name, vr, item))
                item_list.append(item)
            header_offset += item_length
            header_offset += (4 - (item_length % 4)) % 4  # move offset to next 4 byte boundary
        if len(item_list) == 1:
            item_list = item_list[0]
        csa_header[name] = item_list
    return csa_header


def _get_parameters_dicom_sop(dataset, get_header=True):
    """
    Returns a subset of the parameters from a Pydicom dataset

    """
    # get shape of the data (slices, rows, columns, fid_points)
    npts    = int(float(dataset[0x5200,0x9229][0][0x0018,0x9103][0]["SpectroscopyAcquisitionDataColumns"].value))
    navg    = int(float(dataset[0x5200,0x9229][0][0x0018,0x9119][0]['NumberOfAverages'].value))
    nframes = int(float(dataset[0x0028,0x0008].value))

    # decide if we have summed spectra or indiv fids
    avg_flag = True if navg > 1 else False

    navg = 1
    data_shape = nframes, 1, navg, npts

    dataf =  convert_numbers(dataset['SpectroscopyData'].value, True, 'f') # (0x5600, 0x0020)
    data_iter = iter(dataf)
    data = [complex(r, i) for r, i in zip(data_iter, data_iter)]
    complex_data = np.fromiter(data, dtype=np.complex64)
    complex_data.shape = data_shape
    complex_data = complex_data.conjugate()        # different in XA30 bjs

    try:
        if (0x0020,0x9116) in list(dataset[0x5200,0x9229][0].keys()):
            iorient = dataset[0x5200,0x9229][0][0x0020,0x9116][0]['ImageOrientationPatient'].value
        elif (0x0020,0x9116) in list(dataset[0x5200,0x9230][0].keys()):
            iorient = dataset[0x5200,0x9230][0][0x0020,0x9116][0]['ImageOrientationPatient'].value
        else:
            iorient = [-1, 0, 0, 0, 1, -0]
        row_vector    = np.array(iorient[0:3])
        col_vector    = np.array(iorient[3:6])
        voi_position  = dataset[0x5200,0x9230][0][0x0020,0x9113][0]['ImagePositionPatient'].value

        voxel_size = [dataset[0x0018,0x9126][0]['SlabThickness'].value,
                      dataset[0x0018,0x9126][1]['SlabThickness'].value,
                      dataset[0x0018,0x9126][2]['SlabThickness'].value]

        tform = transformation_matrix(row_vector, col_vector, voi_position, voxel_size)

    except:
        # this will trigger default
        voxel_size = np.array([20.0, 20.0, 20.0])
        tform = None

    hdr = str(dataset) if get_header else ''

    params = {'is_dicom_sop': True,
              'sw'          : dataset["SpectralWidth"].value,
              'frequency'   : dataset["TransmitterFrequency"].value,
              'resppm'      : 4.7,
              'echopeak'    : 0.0,
              'nucleus'     : dataset["ResonantNucleus"].value,
              'seqte'       : float(dataset[0x5200,0x9229][0][0x0018,0x9114][0]['EffectiveEchoTime'].value),
              'seqtr'       : float(dataset[0x5200,0x9229][0][0x0018,0x9112][0]['RepetitionTime'].value),
              'voxel_dimensions' : voxel_size,
              'header'      : hdr,
              'transform'   : tform,
              'data'        : complex_data}

    return params






