"""
Routine for reading a Siemens Twix meas.dat format file for MRS data. 

Data is returned in a DataRawFidsum object populated with each of the FIDs
acquired. If multi-Rx coil used, the data is returned uncombined.

This was the original BJS parser for Twix. In 2024, we moved to using the
PyMapVBVD library.  However, in that module we do not have access to the
data from 'prep' scans.  So we keep this copy around in case it is needed.

Seems to have issues with VE data format, but does OK with VB.

"""

# Python modules
import struct
import collections

# 3rd party modules
import numpy as np

# Our modules
import vespa.analysis.fileio.raw_reader as raw_reader
import vespa.common.util.misc as util_misc
from vespa.common.mrs_data_raw import DataRawFidsum
from vespa.common.base_transform import transformation_matrix, rotation_matrix
from vespa.common.twix_parser import TwixRaid
from vespa.common.twix_parser_multi_raid import TwixMultiRaid


RAWDATA_SCALE = 131072.0 * 256.0    # re. ICE_RAWDATA_SCALE in SpecRoFtFunctor.cpp


class RawReaderSiemensTwixBjs(raw_reader.RawReader):
    """ Read a single Siemens Twix file into an DataRawFidsum object. """

    def __init__(self):
        raw_reader.RawReader.__init__(self)
        self.filetype_filter = "Twix File (*.dat)|*.dat;*.DAT;"
        self.multiple = False
        
        
    def read_raw(self, filename, ignore_data=False, open_dataset=None):
        """
        Given Siemens twix filename, return a populated DataRawFidsum object
        - Ignore data has no effect on this method
        - The open_dataset attribute is not used in this reader.
        
        """
        twix, version_flag = self.get_twix(filename)

        d = self._get_parameters(twix.current, version_flag)

        data = d['data']
        if d["remove_os"]:
            data = self._remove_oversampling_basic(data)
            d["sw"] = d["sw"] / 2.0

        data *= RAWDATA_SCALE / float(data.shape[2])    # for typical integral values ~1 for nice display
        d['data'] = np.conjugate(data)                  # conjugate swaps x-axis for proper display of data
        d["data_source"] = filename

        raw = DataRawFidsum(d)

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


    ####################    Internal functions start here     ###############

    def _get_parameters(self, twix, software_version, index='rep'):
        """ Return parameters and data from a Siemens Twix """

        evps = twix.evps

        header, clean_header = self._parse_protocol_data(evps[3][1])
        hdr_dicom = evps[1][1]
        # CMRR sLASER info - need to move later
        #
        # MEAS.sSpecPara.lAutoRefScanMode [aushFreePara2 for VB] >1, water refs are saved
        # MEAS.sSpecPara.lAutoRefScanNo   [aushFreePara3 for VB] = number of scans acquired for ecc and water scaling references at start and end of protocol.

        if 'sProtConsistencyInfo.tBaselineString' not in list(header.keys()):
            if "syngo MR XA" in hdr_dicom:
                software_version = 'nx_va'
            else:
                software_version = 'xx'
        else:
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

        #----------------------------------------------------------------------
        # The following boolean calculations are done to help us try to figure
        # out what specific sequence this is. This base parser does not act on
        # these flags, but they are available for use by drived parsers.

        seqstr = header['tSequenceFileName']

        # Ralf Mekle's SPECIAL, or CIBM SPECIAL sequence
        isSpecial = 'rm_special' in seqstr or 'vq_special' in seqstr

        # Jamie Near's SPECIAL, Masoumeh Dehghani's Adiabatic SPECIAL, SPECIAL or InvRecov SPECIAL
        isjnSpecial = 'jn_svs_special'   in seqstr or \
                      'md_Adiab_Special' in seqstr or \
                      'md_Special'       in seqstr or \
                      'md_Inv_special'   in seqstr

        isjnMP = 'jn_MEGA_GABA' in seqstr   # Jamie Near's MEGA-PRESS sequence

        # Is this another one of Jamie Near's or a sequence derived from Jamie Near's sequences (by Masoumeh Dehghani)
        isjnseq = 'jn_' in seqstr or 'md_' in seqstr

        isWIP529 = 'edit_529' in seqstr     # Is this WIP 529 (MEGA-PRESS)?
        isWIP859 = 'edit_859' in seqstr     # Is this WIP 859 (MEGA-PRESS)?

        # One of Eddie Auerbach or Dinesh Deelchand (CMRR, U Minnesota) sequences?
        isMinn   = 'eja_svs_' in seqstr or 'dkd_svs_' in seqstr

        # Siemens PRESS or STEAM and make sure it's not 'eja_svs_steam'
        isSiemens=('svs_se' in seqstr or 'svs_st' in seqstr) and not ('eja_svs' in seqstr)


        # get data ------------------------------------------

        if index == 'rep':
            data, prep = twix.get_data_numpy_rep_coil_avg_npts_order(return_prep=True)
        elif index == 'echo':
            data, prep = twix.get_data_numpy_echo_coil_avg_npts_order(return_prep=True)
        elif index == 'scan':
            data, prep = twix.get_data_numpy_channel_scan_order_with_prep()
            while len(data.shape) < 4:
                data = np.expand_dims(data, axis=0)
            if prep is not None:
                while len(prep.shape) < 4:
                    prep = np.expand_dims(prep, axis=0)
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
                  'isSpecial'   : isSpecial,
                  'isjnSpecial' : isjnSpecial,
                  'isjnMP'      : isjnMP,
                  'isjnseq'     : isjnseq,
                  'isWIP529'    : isWIP529,
                  'isWIP859'    : isWIP859,
                  'isMinn'      : isMinn,
                  'isSiemens'   : isSiemens,
                  'remove_os'   : header.get("sSpecPara.ucRemoveOversampling", "0x0").strip() == "0x1",
                  'seqname'     : header.get("tSequenceFileName", "svs_xxx"),
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
                  'data'        : data,
                  'prep'        : prep  }

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
        items = [el for el in head if isinstance(el, collections.abc.Iterable) and (key in el)]

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



#------------------------------------------------------------------------------
# test code

def _test():

    import numpy as np
    from matplotlib import pyplot as plt

    fname = r"D:\Users\bsoher\code\repository_svn\sample_data\siemens_twix_vb19_svs_se\meas_MID00165_FID50843_eja_svs_press.dat"

    twix = RawReaderSiemensTwix()
    twix.read_raw(fname)


    # plt.plot(savdat.real)
    # plt.plot(chain.data.real, 'b', linewidth=2)
    # plt.show()


    bob = 10
    bob += 1



if __name__ == '__main__':
    """ Just testing the read_raw method here """
    _test()


# """
#     Gleaned info from Jamie Near's FID-A io files - Siemens Twix
#
#     from io_loadspec_twix.m
#     -----------------------
#
#     (line 50)
#     % find out what sequence, the data were acquired with.  If this is a
#     % multi-raid file, then the header may contain multiple instances of
#     % 'tSequenceFileName' for different scans (including a pre-scan).
#     % Therefore, if multi-raid file, we will need to do a bit of extra digging
#     % to find the correct sequence name.
#
#     sequence=twix_obj.hdr.Config.SequenceFileName;
#
#     % Try to find out what sequnece this is:
#     isSpecial=~isempty(strfind(sequence,'rm_special')) ||...  %Is this Ralf Mekle's SPECIAL sequence?
#                 ~isempty(strfind(sequence,'vq_special'));  %or the CIBM SPECIAL sequence?
#     isjnSpecial=~isempty(strfind(sequence,'jn_svs_special')) ||...  %or Jamie Near's SPECIAL sequence?
#                 ~isempty(strfind(sequence,'md_Adiab_Special')) ||... %or Masoumeh Dehghani's Adiabatic SPECIAL sequence?
#                 ~isempty(strfind(sequence,'md_Special')) ||... %or another version of Masoumeh Dehghani's SPECIAL sequence?
#                 ~isempty(strfind(sequence,'md_Inv_special')); %or Masoumeh Dehghani's Inversion Recovery SPECIAL sequence?
#     isjnMP=~isempty(strfind(sequence,'jn_MEGA_GABA')); %Is this Jamie Near's MEGA-PRESS sequence?
#     isjnseq=~isempty(strfind(sequence,'jn_')) ||... %Is this another one of Jamie Near's sequences
#             ~isempty(strfind(sequence,'md_'));      %or a sequence derived from Jamie Near's sequences (by Masoumeh Dehghani)?
#     isWIP529=~isempty(strfind(sequence,'edit_529')); %Is this WIP 529 (MEGA-PRESS)?
#     isWIP859=~isempty(strfind(sequence,'edit_859')); %Is this WIP 859 (MEGA-PRESS)?
#     isMinn=~isempty(strfind(sequence,'eja_svs_')); %Is this one of Eddie Auerbach's (CMRR, U Minnesota) sequences?
#     isSiemens=(~isempty(strfind(sequence,'svs_se')) ||... %Is this the Siemens PRESS seqeunce?
#                 ~isempty(strfind(sequence,'svs_st'))) && ... % or the Siemens STEAM sequence?
#                 isempty(strfind(sequence,'eja_svs'));    %And make sure it's not 'eja_svs_steam'.
#
#
#     (line 359)
#     %Find the number of points acquired before the echo so that this
#     %information can be stored in the .pointsToLeftshfit field of the data
#     %structure.  Depending on the pulse sequence used to acquire the data, the
#     %header location of this parameter is different.  For product PRESS
#     %seqeunces, the value is located in twix_obj.image.freeParam(1).  For WIP
#     %sequences, the value is located in twix_obj.image.cutOff(1,1).  For CMRR
#     %sequences, the value is located in twix_obj.image.iceParam(5,1).  Special
#     %thanks to Georg Oeltzschner for decoding all of this and sharing the
#     %information with me:
#
#     if isWIP529 || isWIP859
#         leftshift = twix_obj.image.cutOff(1,1);       self.pre in Python twix_parser
#     elseif isSiemens
#         leftshift = twix_obj.image.freeParam(1);
#     elseif isMinn
#         leftshift = twix_obj.image.iceParam(5,1);
#     else
#         leftshift = twix_obj.image.freeParam(1);
#     end
# """

# def _remove_oversampling(scan, d):
#     """ not currently in use as incomplete """
#
#     vector_size = d["vector_size"]
#     remove_os   = d["remove_os"]
#     left_points = scan.pre
#     right_points = scan.post
#     reduced_points = scan.samples_in_scan - scan.pre - scan.post
#     half_vector_size = int(vector_size / 2)
#
#     if (reduced_points % vector_size) != 0:
#         raise ValueError('remove_oversampling: final data size not multiple of vector size.')
#
#
#     if not remove_os:
#         # keep oversampled points but remove extra points
#         start_point = scan.pre
#         scan.data = scan.data[start_point:start_point+reduced_points]
#     else:
#         # remove oversampled data
#         shift_points = scan.post if scan.post < scan.pre else scan.pre
#
#         if shift_points == 0:
#             # no extra pts available, final data will show truncation artifact
#             start_point = scan.pre
#             data = np.array(scan.data[start_point:start_point+reduced_points])
#             data = np.fft.fft(data, n=vector_size)
#             data = np.fft.ifft(data) * 0.5
#             scan.data = data.tolist()
#
#         else:
#             # Extra pts available to use for removing truncation artifact.
#             # Process data twice, centering signal to the left and right of kSpaceCentreColumn (TE)
#             # Retrieve half of final signal from each set.
#             pass


