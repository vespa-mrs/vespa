"""
Routines for reading a Bruker file (ser or fid and acqus) and returning an
DataRaw object populated with the files' data.

Most of the work for reading Bruker files is done by modules we swiped from
the NMRGlue project by J. J. Helmus and C.P. Jaroniec, nmrglue,
http://code.google.com/p/nmrglue, The Ohio State University.

"""
# Python modules
import os.path

# 3rd party modules
import numpy as np
import vespa.analysis.fileio.nmrglue.bruker_for_vespa as nmrglue_bruker

# Our modules
import vespa.common.util.misc as util_misc
import vespa.common.mrs_data_raw as mrs_data_raw
import vespa.analysis.fileio.raw_reader as raw_reader 
import vespa.analysis.fileio.util_exceptions as util_exceptions




class RawReaderBruker(raw_reader.RawReader):
    def __init__(self):
        raw_reader.RawReader.__init__(self)
        self.filetype_filter = "Spectra (acqus,fid,ser)|acqus;fid;ser"
        self.multiple = False
        

    def read_raw(self, path, ignore_data=False, *args, **kwargs):
        """
        Given the fully qualified path to a directory of Bruker files, or a
        fully qualified path to a .fid or .acqus file, returns a DataRaw object
        populated with the parameters and data therein.

        For instance, all of the following are valid:
          read_raw("/home/philip/data/csi2d_03.fid")
          read_raw("/home/philip/data/csi2d_03.fid/acqus")
          read_raw("/home/philip/data/csi2d_03.fid/fid")
    
        Note. ignore_data does not work in this parser.

        """
        if not os.path.isdir(path):
            # Caller passed a directory name + file name. Strip the file name.
            path, _ = os.path.split(path)
        #else: 
            # Caller passed a directory name so I don't need to do anything

        # Check if the two files I want to read are available
        if os.path.isfile(os.path.join(path, "acqus")):
            parameters_filename = "acqus"
        else:
            raise util_exceptions.FileNotFoundError("I can't find the parameters file '%s'" % parameters_filename)

        if os.path.isfile(os.path.join(path, "fid")):
            data_filename = "fid"
        elif os.path.isfile(os.path.join(path, "ser")):
            data_filename = "ser"
        else:
            raise util_exceptions.FileNotFoundError("I can't find the ser or fid data file '%s'" % data_filename)

        # Read the params file and extract the stuff I need.
        dic, data = nmrglue_bruker.read(path, acqus_files=[parameters_filename])  #, read_prog=False)

        d = _get_parameters(dic)

        # The entire acqus file gets stuffed into the "headers"
        hdr = open(os.path.join(path, parameters_filename), "rb").read()
        hdr = hdr.decode()
        d["header"]      = hdr
        d["data_source"] = path

        # Remove the digital filtering of the first N pts in the FID that is 
        # applied by the Bruker system. Ensure that full length of the data
        # is maintained. I apply complex conjugate here so the data are aligned
        # on the Vespa x-axis properly (ie. data is 'flipped' on read).
        temp = nmrglue_bruker.remove_digital_filter(dic, data)
        data_f = np.zeros_like(data)
        data_f[0:len(temp)] = np.conj(temp)

        # Ensure the data is the right shape
        # data_f = mrs_data_raw.normalize_data_dims(data_f)
        while len(data_f.shape) < 4:
            data_f.shape = [1,] + list(data_f.shape)
        
        d["data"] = data_f

        return [mrs_data_raw.DataRaw(d),]



####################    Internal functions start here     ###############

def _get_parameters(dic):
    """
    Given the acqus file as a dict, extracts a few specific parameters 
    and returns a flat dict containing those parameters and their value.

    The returned dict is appropriate for passing to DataRaw.inflate().

    """
    # An acqus dict contains many, many keys. We're only interested in a few.
    
    nucleus = util_misc.normalize_isotope_name(dic["acqus"]["NUC1"])
    if nucleus is None:
        nucleus = ""

    params = {'sw'          : float(dic["acqus"]["SW_h"]) ,
              'frequency'   : float(dic["acqus"]["BF1"]),
              'resppm'      : 4.7,
              'echopeak'    : 0.0,
              'nucleus'     : nucleus,
              'seqte'       : 0.0,
              'seqtr'       : 1.0,
              'voxel_dimensions' : [20.0,20.0,20.0],
              'header'      : '',
              'transform'   : None,
              'data'        : None}


    return params

