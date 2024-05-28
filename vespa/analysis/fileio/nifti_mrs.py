"""
Routines for reading the NIfTI-MRS format and returning a
DataRaw object populated with the file's data.

"""

# Python modules
import os.path
import json

# 3rd party modules
import numpy as np
import nibabel as nib

# Our modules
import vespa.common.mrs_data_raw as mrs_data_raw
import vespa.analysis.fileio.raw_reader as raw_reader
import vespa.analysis.fileio.util_exceptions as util_exceptions





class RawReaderNiftiMrs(raw_reader.RawReader):

    def __init__(self):

        raw_reader.RawReader.__init__(self)
        self.filetype_filter = "Spectra (*.nii,*.nii.gz)|*.nii;*.nii.gz"
        self.multiple = False


    def read_raw(self, filename, ignore_data=False, *args, **kwargs):
        """
        Given NIfTI filename, return populated DataRaw object.
        - 'ignore_data' has no effect on this parser

        """
        msg = ''
        if not os.path.isfile(filename):
            msg += "NIfTI-MRS file not found - '%s'\n" % filename
        if msg:
            raise util_exceptions.FileNotFoundError(msg)

        img = nib.load(filename)
        d = _get_parameters(img, filename)

        d["data_source"] = filename

        if not (d["data"].shape[0:3] == (1,1,1)):
            # svs data should have x,y,z dims = 1
            msg = "NIfTI-MRS data is unknown size, shape = %s" % str(d["data"].shape)
            raise util_exceptions.IncorrectDimensionalityError(msg)

        if len(d["data"].shape) > 4:
            for item in d["data"].shape[4:]:
                if item != 1:
                    # not dealing with Navg or Ncoil yet
                    msg = "NIfTI-MRS data shape not a single FID = %s" % str(d["data"].shape)
                    raise util_exceptions.IncorrectDimensionalityError(msg)

        d["data"].shape = 1,1,1,d["data"].shape[3]

        raws = [mrs_data_raw.DataRaw(d),]

        return raws


####################    Internal functions start here     ###############


def _get_parameters(img, filename):
    """ Given a NIfTI-MRS object and filename, extract parameters needed into a dict. """

    # parse data -----------------------------------------------

    complex_data = img.get_fdata(dtype=np.complex64)
    complex_data = np.array(complex_data.tolist())

    # parse header -----------------------------------------------

    header = img.header
    dwelltime = header['pixdim'][4]
    intent_name = header.get_intent()[2]
    affine = header.get_best_affine()

    hdr_ext_codes = img.header.extensions.get_codes()
    extn = json.loads(img.header.extensions[hdr_ext_codes.index(44)].get_content())

    keys = extn.keys()
    nucleus = extn['ResonantNucleus'][0] if 'ResonantNucleus' in keys else '1H'
    frequency = extn['SpectrometerFrequency'][0] if 'SpectrometerFrequency' in keys else 128.0
    te = extn['EchoTime'] if 'EchoTime' in keys else 0.0
    tr = extn['RepetitionTime'] if 'RepetitionTime' in keys else 0.0

    #  Vespa comment is '\n' delineated string, extract for clarity in GUI
    if 'VespaComment' in keys:
        value = extn['VespaComment']['Value']
        extn['VespaComment']['Value'] = 'extracted for display'     # TODO bjs - why is this here overwriting?
        comment = '\nVespaComment - Value\n--------------------------------------\n'
        comment += value
    else:
        comment = '\n'

    resppm = 4.7
    echopeak = 0.0
    voxel_size = [10000.0, 10000.0, 10000.0]
    tform = affine

    # format header for display ------------------------------------------

    hdr1 = str(header)                  # the Nifti main header
    hdr2 = json.dumps(extn, indent=4)   # the JSON sidecar extension header
    hdr = hdr1+'\n\nNIfTI-MRS - Extension\n-----------------------------------------\n'+hdr2+'\n'+comment

    params = {'sw'               : float(1.0/dwelltime),
              'frequency'        : frequency,
              'resppm'           : resppm,
              'echopeak'         : echopeak,
              'nucleus'          : nucleus,
              'seqte'            : te,
              'seqtr'            : tr,
              'voxel_dimensions' : voxel_size,
              'header'           : hdr,
              'transform'        : tform,
              'data'             : complex_data}

    return params


