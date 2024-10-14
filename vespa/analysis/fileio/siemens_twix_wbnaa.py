"""
Routines for reading a Siemens Twix meas.dat format file that has multiple
FIDs for one SVS data acquisition inside. Every other FID is phase inverted
and needs to be combined in an (a-b) + (a-b) + ... fashion.

Data is returned in a DataRaw object populated with the file's data.
"""


# Python modules
import struct
import collections

# 3rd party modules
import numpy as np

# Our modules
import vespa.analysis.fileio.raw_reader as raw_reader
import vespa.common.util.misc as util_misc
from vespa.common.mrs_data_raw_wbnaa import DataRawWbnaa
from vespa.common.base_transform import transformation_matrix, rotation_matrix
from vespa.common.twix_parser import TwixRaid
from vespa.common.twix_parser_multi_raid import TwixMultiRaid



# Change to True to enable the assert() statements sprinkled through the code
ASSERTIONS_ENABLED = False

# data is complex64 per Siemens documentation for Twix
NUMPY_DATA_TYPE = np.complex64
# BYTES_PER_ELEMENT expresses how many bytes each element occupies. You
# shouldn't need to change this definition.
BYTES_PER_ELEMENT = np.zeros(1, dtype=NUMPY_DATA_TYPE).nbytes



class RawReaderSiemensTwixWbnaa(raw_reader.RawReader):
    # This inherits from raw_reader.RawReader (q.v.). The only methods you
    # need to implement are __init__() and read_raw(). You *may* want to 
    # override or supplement some of raw_reader.RawReader's other methods.

    def __init__(self):
        raw_reader.RawReader.__init__(self)
        
        # The sample files given to us all had uppercase extensions, so my guess
        # is that all Philips files are created this way. The GTK file open 
        # dialog takes case sensitivity seriously and won't show foo.SPAR if we
        # use a filetype filter of '*.spar'. Hence we specify both upper and
        # lower case extensions.
        self.filetype_filter = "Twix File (*.dat)|*.dat;*.DAT;"
        self.multiple = False
        
        
    def read_raw(self, filename, ignore_data=False, open_dataset=None):
        """
        Given the name of a .spar or .sdat file, returns a DataRawWbnaa object
        populated with the parameters and data represented by the file pair.

        When ignore_data is True, this function only reads the parameters file
        which can be much faster than reading both params & data.
        
        The open_dataset attribute is not used in this reader. 
        
        """
        twix, version_flag = self.get_twix(filename)

        d = self._get_parameters(twix, version_flag)

        # NB. bjs 9/2024 - May be broken by update to twix_parser.py code, may
        #   need local _get_parameters() method to properly organize data

        # Create a DataRawWbnaa out of the first set of data in the list.
        datas = d["data"]       # this is a list of numpy arrays [1,ncoil,nfids,npts]

        d["data"] = datas[0]
        d["data_source"] = filename
        raw = DataRawWbnaa(d)
        
        # Concatenate the remainder onto the first.
        for data in datas[1:]:
            # data.shape = shape
            d["data"] = data
            raw.concatenate(DataRawWbnaa(d))

        return [raw,]


    def get_twix(self, filename):

        with open(filename, 'rb') as fin:
            # we can tell the type of file from the first two uints in the header
            uint1, uint2 = struct.unpack("II", fin.read(8))
            version_flag = 'vd' if uint1 == 0 and uint2 <= 64 else 'vb'

        if version_flag == 'vb':
            twix = TwixRaid()
            twix.populate_from_file(filename)
        else:
            twix = TwixMultiRaid()
            twix.populate_from_file(filename)

        return twix, version_flag


    def _remove_oversampling_basic(self, data):
        """ Basic data oversampling removal, not full Siemens algorithm """
        nrep, ncoil, nfid, npts = data.shape
        data = np.roll(np.fft.fft(data, axis=3), int(npts/4), axis=3)
        data = np.fft.ifft(np.roll(data[:,:,:,0:int(npts/2)], int(-npts/4), axis=3), axis=3)

        return data


    ################    Internal helper functions start here     #############

    def _get_parameters(self, twix, index='rep'):
        """ Return parameters and data from a Siemens Twix """

        msg = 'NB. bjs - broken by Aug 2024 twix update, fix please!'
        raise ValueError(msg)

        evps = twix.evps

        header, clean_header = self._parse_protocol_data(evps[3][1])

        software = header['sProtConsistencyInfo.tBaselineString'].lower()
        if   'n4_vb' in software: software_version = 'vb'
        elif 'n4_vd' in software: software_version = 'vd'
        elif 'n4_ve' in software: software_version = 've'

        wiplong, wipdouble = self._get_xprot_wipvars(evps[2][1])

        if software_version == 'vb':
            ref_nscans  = int(header['sSpecPara.lAutoRefScanNo']) if 'sSpecPara.lAutoRefScanNo' in header.keys() else 0
            ref_flag    = int(ref_nscans != 0)
            prep_nscans = 0
            if 'sSpecPara.lPreparingScans' in header.keys():
                prep_nscans = int(header['sSpecPara.lPreparingScans'])
        else:       # VD and VE
            ref_flag    = int(header['sSpecPara.lAutoRefScanMode'])
            ref_nscans  = int(header['sSpecPara.lAutoRefScanNo'])
            prep_nscans = 0
            if 'sSpecPara.lPreparingScans' in header.keys():
                prep_nscans = int(header['sSpecPara.lPreparingScans'])

        lAverages = int(header['lAverages'])     # number of metabolite scans

        # get data ------------------------------------------

        if index == 'rep':
            data, prep = twix.get_data_numpy_rep_coil_avg_npts_order(return_prep=True)
        elif index == 'echo':
            data, prep = twix.get_data_numpy_echo_coil_avg_npts_order(return_prep=True)
        elif index == 'scan':
            data = twix.get_data_numpy_scan_channel_order()     # differs from Siemens twix regular
            while len(data.shape) < 4:
                data = np.expand_dims(data, axis=0)
            prep = None
        else:
            data, prep = twix.get_data_numpy_rep_coil_avg_npts_order(return_prep=True)

        if software_version == 'vb':
            fid_str = twix.scans[0].free_parameters[0]            # extra points at begin of FID?
        else:
            fid_str = twix.scans[0].scan_header.free_parameters[0]

        acqdim0 = int(2 ** np.floor(np.log2(data.shape[3] - fid_str)))      # largest pow(2) given shape[3]
        data = data[:, :, :, fid_str:fid_str+acqdim0].copy()
        if prep is not None:
            prep = prep[:, :, :, fid_str:fid_str + acqdim0].copy()

        # these next bits are WBNAA specific ------------------------------

        nscans = int(data.shape[2] / 2)
        datas = []

        for i in range(nscans):
            dat1 = data[:,:,i*2,:]
            dat2 = data[:,:,i*2+1,:]
            tmp = np.conjugate(dat1 - dat2)
            while len(tmp.shape) < 4:
                tmp = np.expand_dims(tmp, axis=0)
            datas.append(tmp)

        ## comment out the code above and uncomment below to see all scans
        ## rather than doing difference paired data

        #    for scan in scans:
        #        datas.append(np.conjugate(np.array(scan.data)))

        # get header info -----------------------------------

        if software_version == 'vb':
            # get voxel size
            ro_fov = self._get_xprot(evps[0][1],"VoI_RoFOV", 20.0)
            pe_fov = self._get_xprot(evps[0][1],"VoI_PeFOV", 20.0)
            slice_thickness = self._get_xprot(evps[0][1],"VoI_SliceThickness", 20.0)

            # get position information
            pos_sag = self._get_xprot(evps[0][1],"VoI_Position_Sag", 0.0)
            pos_cor = self._get_xprot(evps[0][1],"VoI_Position_Cor", 0.0)
            pos_tra = self._get_xprot(evps[0][1],"VoI_Position_Tra", 0.0)

            # get orientation information
            in_plane_rot = self._get_xprot(evps[0][1],"VoI_InPlaneRotAngle", 0.0)
            normal_sag   = self._get_xprot(evps[0][1],"VoI_Normal_Sag", 1.0)
            normal_cor   = self._get_xprot(evps[0][1],"VoI_Normal_Cor", 0.0)
            normal_tra   = self._get_xprot(evps[0][1],"VoI_Normal_Tra", 0.0)

        else:
            # get voxel size
            ro_fov          = self._get_hdr_float_def(header, 'sSpecPara.sVoI.dReadoutFOV', 20.0)
            pe_fov          = self._get_hdr_float_def(header, 'sSpecPara.sVoI.dPhaseFOV',   20.0)
            slice_thickness = self._get_hdr_float_def(header, 'sSpecPara.sVoI.dThickness',  20.0)

            # get position information
            pos_sag = self._get_hdr_float_def(header, 'sSpecPara.sVoI.sPosition.dSag', 0.0)
            pos_cor = self._get_hdr_float_def(header, 'sSpecPara.sVoI.sPosition.dCor', 0.0)
            pos_tra = self._get_hdr_float_def(header, 'sSpecPara.sVoI.sPosition.dTra', 0.0)

            # get orientation information
            in_plane_rot = self._get_hdr_float_def(header, 'sSpecPara.sVoI.InPlaneRot',   0.0)
            normal_sag   = self._get_hdr_float_def(header, 'sSpecPara.sVoI.sNormal.dSag', 0.0)
            normal_cor   = self._get_hdr_float_def(header, 'sSpecPara.sVoI.sNormal.dCor', 0.0)
            normal_tra   = self._get_hdr_float_def(header, 'sSpecPara.sVoI.sNormal.dTra', 0.0)

        voxel_size = [ro_fov, pe_fov, slice_thickness]
        voxel_pos  = [pos_sag, pos_cor, pos_tra]

        # the orientation is stored in a somewhat strange way - a normal vector and
        # a rotation angle. To get the row vector, we first use Gram-Schmidt to
        # make [-1, 0, 0] (the default row vector) orthogonal to the normal, and
        # then rotate that vector by the rotation angle  - (from Suspect package)

        x_vector      = np.array([-1, 0, 0])
        normal_vector = np.array([normal_sag, normal_cor, normal_tra])
        orthogonal_x  = x_vector - np.dot(x_vector, normal_vector) * normal_vector
        orthonormal_x = orthogonal_x / np.linalg.norm(orthogonal_x)
        rot_matrix    = rotation_matrix(in_plane_rot, normal_vector)
        row_vector    = np.dot(rot_matrix, orthonormal_x)
        col_vector    = np.cross(row_vector, normal_vector)

        tform = transformation_matrix(row_vector, col_vector, voxel_pos, voxel_size)

        params = {'wiplong'     : wiplong,
                  'wipdouble'   : wipdouble,
                  'ref_flag'    : ref_flag,
                  'ref_nscans'  : ref_nscans,
                  'prep_nscans' : prep_nscans,
                  'lAverages'   : lAverages,
                  'remove_os'   : header.get("sSpecPara.ucRemoveOversampling", "0x0").strip() == "0x1",
                  'seqname'     : header.get("tSequenceFileName", "wbnaa"),
                  'sw'          : 1.0 / (float(header.get("sRXSPEC.alDwellTime[0]", 1.0)) * 1e-9),
                  'frequency'   : float(header["sTXSPEC.asNucleusInfo[0].lFrequency"])/1000000.0,
                  'resppm'      : 4.7,
                  'echopeak'    : 0.0,
                  'nucleus'     : header["sTXSPEC.asNucleusInfo[0].tNucleus"].replace('"',' ').strip(),
                  'seqte'       : float(header["alTE[0]"])*1e-6,
                  'seqtr'       : float(header["alTR[0]"])*1e-6,
                  'voxel_dimensions' : voxel_size,
                  'header'      : clean_header,
                  'transform'   : tform,
                  'data'        : datas,
                  'prep'        : prep  }

        # d["readout_os"] = float(self._get_xprot(evps[0][1], "ReadoutOS", 1.0))

        return params


    def _parse_protocol_data(self, prot):
        """
        Return dict of MrProtocol/MrPhoenixProtocol item from Siemens CSA Header.

        'Prot' (~ 32k string) is a JSONish format followed by a long string
        delimited by 'ASCCONV' markers that contains a list of name=value pairs. We
        ignore data outside of ASCCONV delimiters and return inside item as a dict.

        As of the Siemens VD software version the starting string is no longer just
        ### ASCCONV BEGIN ### but has other info about what was converted inserted
        after the BEGIN and before the ### delimiter. To get around this for now,
        we search just for ### ASCONV BEGIN, and then throw away the first line
        after we split the string into lines.

        """
        str, end = prot.find("### ASCCONV BEGIN"), prot.find("### ASCCONV END ###")

        clean = prot[str:end+len("### ASCCONV END ###")]
        prot  = prot[str+len("### ASCCONV BEGIN ###"):end]

        f = lambda pair: (pair[0].strip(), pair[1].strip())

        lines = prot.split('\n')
        lines = lines[1:]
        lines = [f(line.split('=')) for line in lines if line]

        return dict(lines), clean


    def _get_xprot(self, head_only, key, default):

        head = util_misc.normalize_newlines(head_only)
        head = head.split("\n")

        # find substring 'key' in list even if item in list is not iterable
        items = [el for el in head if isinstance(el, collections.Iterable) and (key in el)]

        for item in items:
            read_double = True if 'ParamDouble' in item else False
            read_int    = True if 'ParamLong'   in item else False
            read_string = True if 'ParamString' in item else False

            start = item.find("{")
            end = item.find("}")
            if start != -1 and end != -1:
                temp = item[start+1:end]
                temp = " ".join(temp.split())       # remove duplicate white space
                temp = temp.split()

                if read_double:
                    if   len(temp) == 1: return float(temp[0])
                    elif len(temp) == 2: return float(0.0)
                    elif len(temp) == 3: return float(temp[2])
                elif read_int:
                    if   len(temp) == 1: return int(temp[0])
                    elif len(temp) == 2: return int(0)
                    elif len(temp) == 3: return int(temp[2])
                else:
                    if   len(temp) == 1: return temp[0]
                    elif len(temp) == 2: return temp
                    elif len(temp) == 3: return temp[2]

        return default


    def _get_xprot_wipvars(self, hdr):

        wip = hdr[hdr.lower().find('<ParamMap."sWiPMemBlock">'.lower()):]

        alf = wip[wip.find('<ParamLong."alFree">'):]
        alf = alf[alf.find('{')+1:alf.find('}')]
        alf = [int(item) for item in alf.replace('\n', '').split()]

        adf = wip[wip.find('<ParamDouble."adFree">'):]
        adf = adf[adf.find('{')+1:adf.find('}')]
        adf = [float(item) for item in adf.replace('<Precision> 6', '').replace('\n', '').split()]

        return alf, adf


    def _get_hdr_float_def(self, header, label, default):
        return float(header[label]) if label in header.keys() else default


# def _read_data(scans):
#     """
#     Given a filename, the number of points in each FID and the number of
#     FIDs, reads the data from the file and returns a list of numpy arrays.
#
#     The list will contain nfids arrays; the arrays will be 1D arrays with
#     npoints elements.
#
#     """
#     nscans = int(len(scans) / 2)
#     datas = []
#
#     for i in range(nscans):
#         dat1 = np.array(scans[i * 2].data)
#         dat2 = np.array(scans[i * 2 + 1].data)
#         data = np.conjugate(dat1 - dat2)
#         datas.append(data)
#
#     ## comment out the code above and uncomment below to see all scans
#     ## rather than doing difference paired data
#
#     #    for scan in scans:
#     #        datas.append(np.conjugate(np.array(scan.data)))
#
#     return datas


            