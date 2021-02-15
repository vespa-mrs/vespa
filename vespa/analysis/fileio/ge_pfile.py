# Python modules

# 3rd party modules
import numpy as np

# Our modules
import vespa.analysis.fileio.raw_reader as raw_reader 
import vespa.common.ge_read_pfile as read_pfile

from vespa.common.mrs_data_raw import DataRawProbep, DataRawCmrrSlaser
from vespa.common.base_transform import transformation_matrix
from vespa.analysis.fileio.util_exceptions import UnsupportedPulseSequenceError, \
                                                  OpenFileUserReadRawError, \
                                                  SIDataError




class RawReaderGePfile(raw_reader.RawReader):

    def __init__(self):
        raw_reader.RawReader.__init__(self)
        self.filetype_filter = "GE PROBE-P (*.7)|*.7;"
        self.multiple = False
        

class RawReaderGeProbep(RawReaderGePfile):
    """
    Given the name of a GE Pfile (*.7) file, the default is to return two
    DataRawProbep objects, one for water, one for metabolites
    - ignore_data has no effect in this method

    - the mapper() method returns data as (x, y, z, navg, ncoil, nspectralpts)
    - Analysis wants (1,ncoil,navg,nspectralpts) so we ditch a few initial
        dimensions and swap ncoil/navg axes in the numpy array.

    """

    def read_raw(self, filename, ignore_data=False, *args, **kwargs):
        pfile = read_pfile.Pfile(filename)
        psd = pfile.hdr.rhi_psdname.decode('utf-8').lower()

        if psd != 'probe-p':
            msg = 'Not PROBE-P data, contains data from (%s) sequence. Returning.' % (psd,)
            raise UnsupportedPulseSequenceError(msg)

        if not pfile.is_svs:
            raise SIDataError

        d = _get_parameters(pfile)

        # parse data for water and metab --------------------------------------

        if pfile.map.raw_unsuppressed is None:
            raise OpenFileUserReadRawError("No Water data in PROBE file.")

        d["data_source"] = filename+'.water'
        water = pfile.map.raw_unsuppressed                  # typically (1,1,1,navg,ncoil,npts)
        water.shape = water.shape[(len(water.shape) - 4):]
        d["data"] = np.transpose(water, [0,2,1,3])          # swap ncoil/navg axes
        wat = DataRawProbep(d)

        d["data_source"] = filename+'.metab'
        metab = pfile.map.raw_suppressed                    # typically (1,1,1,navg,ncoil,npts)
        metab.shape = metab.shape[(len(metab.shape)-4):]
        d["data"] = np.transpose(metab, [0,2,1,3])          # swap ncoil/navg axes
        met = DataRawProbep(d)

        return [wat, met]       # order matters here to auto-label in Notebook


class RawReaderGeOslaserCmrr(RawReaderGePfile):
    """
    Given the name of a GE Pfile (*.7) file, the default is to return six
    DataRawProbep objects for: coil combine, watquant1, ecc1, watquant2, ecc2
    and metabolites.
    - ignore_data has no effect in this method

    """

    def read_raw(self, filename, ignore_data=False, *args, **kwargs):
        pfile = read_pfile.Pfile(filename)
        psd = pfile.hdr.rhi_psdname.decode('utf-8').lower()

        if psd != 'oslaser':
            msg = 'Not osLASER data, contains data from (%s) sequence. Returning.' % (psd,)
            raise UnsupportedPulseSequenceError(msg)

        if not pfile.is_svs:
            raise SIDataError

        d = _get_parameters(pfile)

        # parse data for water and metab --------------------------------------

        raws = []

        water = pfile.map.raw_unsuppressed                  # typically (1,1,1,navg,ncoil,npts)
        metab = pfile.map.raw_suppressed
        water.shape = water.shape[(len(water.shape) - 4):]
        metab.shape = metab.shape[(len(metab.shape) - 4):]
        water = np.transpose(water, [0, 2, 1, 3])           # swap ncoil/navg axes
        metab = np.transpose(metab, [0, 2, 1, 3])

        if water.shape[2] != 8:
            msg = "There were not 8 prepare FIDs (%d) in oslaser data file, returning." % (water.shape[2],)
            raise OpenFileUserReadRawError(msg)

        # dataset1 - scan 2, water unsuppressed for coil combine - this is a fudge

        d["data"] = water[:, :, 2:3, :]
        d["data_source"] = filename + '.combine'
        raws.append(DataRawCmrrSlaser(d))

        # dataset2 - scan 2-3 (2 total), water unsuppressed for ECC

        d["data"] = water[:, :, 2:4, :].copy()
        d["data_source"] = filename + '.ecc1'
        raws.append(DataRawCmrrSlaser(d))

        # dataset3 - scans 0-1 (2 total), water unsuppressed for water scale

        d["data"] = water[:, :, 0:2, :].copy()
        d["data_source"] = filename + '.water1'
        raws.append(DataRawCmrrSlaser(d))

        # dataset5 - scans 6-7 (2 total), water unsuppressed for ecc

        d["data"] = water[:, :, 6:8, :].copy()
        d["data_source"] = filename + '.ecc2'
        raws.append(DataRawCmrrSlaser(d))

        # dataset6 - scans 4-5 (2 total), water unsuppressed for water scale

        d["data"] = water[:, :, 4:6, :].copy()
        d["data_source"] = filename + '.water2'
        raws.append(DataRawCmrrSlaser(d))

        # dataset4 - metab array

        d["data"] = metab
        d["data_source"] = filename + '.metab64'
        raws.append(DataRawCmrrSlaser(d))

        return raws

        # d["data_source"] = filename + '.water'
        # water = pfile.map.raw_unsuppressed  # typically (1,1,1,navg,ncoil,npts)
        # water.shape = water.shape[(len(water.shape) - 4):]
        # d["data"] = np.transpose(water, [0, 2, 1, 3])  # swap ncoil/navg axes
        # wat = DataRawProbep(d)
        #
        # d["data_source"] = filename + '.metab'
        # metab = pfile.map.raw_suppressed  # typically (1,1,1,navg,ncoil,npts)
        # metab.shape = metab.shape[(len(metab.shape) - 4):]
        # d["data"] = np.transpose(metab, [0, 2, 1, 3])  # swap ncoil/navg axes
        # met = DataRawProbep(d)
        #
        # return [wat, met]  # order matters here to auto-label in Notebook



####################    Internal functions start here     ###############

def _get_parameters(pfile):
    """ Given the pfile header, extracts a few parameters into a dict """

    hdr = pfile.hdr

    # A pfile contains many attributes. We're only interested in a few ------

    sw = float(hdr.rhr_rh_user0)  # in Hz
    freq = float(hdr.rhr_rh_ps_mps_freq) / 1e7  # in MHz
    te = float(hdr.rhi_te) / 1000  # in ms
    tr = float(hdr.rhi_tr) / 1000  # in ms

    #  Use the mps freq and field strength to determine gamma and thus isotope
    gamma = (hdr.rhr_rh_ps_mps_freq * 1e-7) / (hdr.rhe_magstrength / 10000.0)
    if abs(gamma - 42.57) < 0.3:
        nucleus = "1H"
    elif abs(gamma - 10.7) < 0.3:
        nucleus = "13C"
    elif abs(gamma - 17.2) < 0.3:
        nucleus = "31P"
    else:
        nucleus = "1H"
    resppm = 4.7 if nucleus == "1H" else 0.0

    pulseq = hdr.rhi_psdname.decode('utf-8').lower()

    # calculate directional cosines to get row/col vectors ---------------------
    dcos = pfile.map.get_dcos
    dcos[dcos == 0.0] = 0.0             # remove -0.0 values

    try:
        voxel_size = pfile.map.get_select_box_size
        col_vector = dcos[1][0:3]
        row_vector = dcos[0][0:3]
        voi_position = pfile.map.get_select_box_center
        tform = transformation_matrix(row_vector, col_vector, voi_position, voxel_size)
    except:
        # this will trigger default
        voxel_size = [20.0, 20.0, 20.0]
        tform = None

    hdr = "\n".join(pfile.dump_header_strarr())

    params = {'pulseq' : pulseq,
              'sw': sw,
              'frequency': freq,
              'resppm': resppm,
              'echopeak': 0.0,
              'nucleus': nucleus,
              'seqte': te,
              'seqtr': tr,
              'voxel_dimensions': voxel_size,
              'header': hdr,
              'transform': tform}

    return params




