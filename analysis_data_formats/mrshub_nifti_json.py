"""
Routine for reading a NIfTI/JSON format file returning an
DataRawFidsum object populated with the file's data.
Still in development.  So just a first pass here.

"""
# Python modules
import json

# 3rd party modules
import numpy as np
import nibabel as ni
import nibabel.orientations as nio

# Our modules
from vespa.common.mrs_data_raw import DataRawFidsum
from vespa.analysis.fileio.raw_reader import RawReader

from vespa.common.constants import DEGREES_TO_RADIANS, RADIANS_TO_DEGREES



class RawReaderNiftiJson(RawReader):

    def __init__(self):

        RawReader.__init__(self)
        self.filetype_filter = "NIfTI File (*.nii/*.nii.gz)|*.nii;*.NII;*.nii.gz;*.NII.GZ"
        self.multiple = False
        
        
    def read_raw(self, filename, ignore_data=False, open_dataset=None):
        """
        Given NIfTI *.nii filename, return populated DataRawFidsum object.
        - 'ignore_data' has no effect on this parser

        """
        nii_obj = ni.load(filename)
        d = _get_parameters(nii_obj)
        d["data_source"] = filename

        return [DataRawFidsum(d),]





####################    Internal functions start here     ###############

def _get_parameters(nii_obj):
    """
    Given the contents of a NIfTI file as a Nibabel object, extract a few
    specific parameters and returns a flat dict containing those parameters
    and their value.
    The returned dict is appropriate for passing to DataRaw.inflate().

    """
    d = { }

    hdr = nii_obj.get_header().copy()
    dat = nii_obj.dataobj[:,:,:]
    aff = nii_obj.get_affine().copy()
    ext = nii_obj.header.extensions[0].get_content().decode('utf8')
    jsn = json.loads(ext)

    # Massage JSON/NIfTI format into a readable list of strings
    head = json.dumps(jsn, indent=2)
    head = 'Nifti Header\n' + str(hdr) + '\n' + '\nJSON Header Extension\n' + head

    # parse data - has to be the original byte array --------------------------

    if "dim_5" not in jsn.keys():
        data = np.squeeze(dat.copy())
        data = np.pad(data, (0, 1), 'constant', constant_values=(0+0j, 0+0j))
        data = np.conjugate(data)
        data.shape = 1,1,1,data.shape[0]
    else:
        if jsn["dim_5"].upper() == "DIM_COIL":
            if jsn["dim_6"].upper() == "DIM_DYN":
                if jsn["dim_6_use"] == "Signal repetitions":
                    data = np.squeeze(dat.copy())
                    data = data.transpose(1,2,0)
                    data = np.conjugate(data)
                else:
                    raise ValueError('unknown dim_6_use field')
            else:
                raise ValueError('unknown dim_6 field')
        else:
            raise ValueError('unknown dim_5 field')

    # parse transform and other MRS parameters -----------------

    te = float(jsn["EchoTime"]) if "EchoTime" in jsn.keys() else 0.0
    tr = float(jsn["RepetitionTime"]) if "RepetitionTime" in jsn.keys() else 0.0

    try:
        nii_ornt   = nio.io_orientation(aff)
        dcm_ornt   = nio.axcodes2ornt(('L', 'P', 'S'))
        ornt_trans = nio.ornt_transform(nii_ornt, dcm_ornt)
        aff_trans  = nio.inv_ornt_aff(ornt_trans, [1,1,1])
        tform      = np.dot(aff, aff_trans)
        voxel_size = np.linalg.norm(tform, axis=0)[0:3]  # from common.base_transform.py

    except:
        # this will trigger default
        voxel_size = [20.0, 20.0, 20.0]
        tform = None


    params = {'sw'          : jsn["SpectralWidth"],
              'frequency'   : float(jsn["TransmitterFrequency"]/1e6),
              'resppm'      : 4.7,
              'echopeak'    : 0.0,
              'nucleus'     : jsn["ResonantNucleus"],
              'seqte'       : te,
              'seqtr'       : tr,
              'voxel_dimensions' : voxel_size,
              'header'      : head,
              'transform'   : tform,
              'data'        : data    }


    return params

    return d  
    
    
    
