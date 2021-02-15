"""
Routines for reading a Varian file pair (procpar & fid) and returning an
DataRaw object populated with the files' data.

Documentation concerning the decisions in code below is scattered across a
few different manuals. We found relevant doc in both VNMR and VNMRJ manuals.
(The latter is a successor to the former and/or a Java-based version.)
We looked at docs with the following titles:
  User Guide: Imaging
  Vnmr  Command and Parameter Reference
  VnmrJ Command and Parameter Reference
  VnmrJ User Programming

All of these docs are available at Agilent.com.

Most of the work for reading Varian files is done by modules we swiped from
the NMRGlue project (J. J. Helmus and C.P. Jaroniec, nmrglue,
http://code.google.com/p/nmrglue, The Ohio State University.)


"""
# Python modules
import os.path

# 3rd party modules
import vespa.analysis.fileio.nmrglue.varian as nmrglue_varian

# Our modules
import vespa.common.util.misc as util_misc
import vespa.common.mrs_data_raw as mrs_data_raw
import vespa.analysis.fileio.raw_reader as raw_reader 
import vespa.analysis.fileio.util_exceptions as util_exceptions




class RawReaderVarian(raw_reader.RawReader):
    def __init__(self):

        raw_reader.RawReader.__init__(self)
        self.filetype_filter = "Spectra (procpar,fid)|procpar;fid"
        self.multiple = False
        

    def read_raw(self, path, ignore_data=False, *args, **kwargs):
        """
        Given the fully qualified path to a directory of Varian files, returns
        a DataRaw object populated with the parameters and data therein.

        One can also pass the fully qualified path to a procpar or fid file. I.e.
        all of the following are valid and will result in the same output:
          read_raw("/home/philip/data/csi2d_03.fid")
          read_raw("/home/philip/data/csi2d_03.fid/procpar")
          read_raw("/home/philip/data/csi2d_03.fid/fid")
    
        - 'ignore_data' - when True, this function only reads the parameters file
              which can be much faster than reading both params & data.
        """
        if not os.path.isdir(path):
            # Caller passed a directory name + file name. Strip the file name.
            path, _ = os.path.split(path)

        # Create the name of the two files I want to read
        parameters_filename = os.path.join(path, "procpar")
        data_filename = os.path.join(path, "fid")

        if not os.path.isfile(parameters_filename):
            msg = "Parameter file not found - '%s'" % parameters_filename
            raise util_exceptions.FileNotFoundError(msg)

        if not ignore_data and not os.path.isfile(data_filename):
            msg = "Data file not found - '%s'" % data_filename
            raise util_exceptions.FileNotFoundError(msg)

        # Read the params file and extract the stuff I need.
        procpar = nmrglue_varian.read_procpar(parameters_filename)

        d = _get_parameters(procpar)

        # The entire procpar file gets stuffed into the "headers"
        hdr = open(parameters_filename, "rb").read()
        hdr = hdr.decode()
        d["header"]      = hdr
        d["data_source"] = path

        # Read data, too, if the caller wants me to do so.
        if not ignore_data:
            shape = nmrglue_varian.find_shape(procpar)

            _, data = nmrglue_varian.read_fid(data_filename, shape)

            # Ensure the data is the right shape
            while len(data.shape) < 4:
                data.shape = [1,] + list(data.shape)
            d["data"] = data

        return [mrs_data_raw.DataRaw(d),]



####################    Internal functions start here     ###############

def _get_parameters(procpar):
    """ Given a procpar file (dict) extract parameters and return a flat dict. """

    d = { }

    # isotope is given as H1, P31, etc. We prefer it in reverse order
    nucleus = util_misc.normalize_isotope_name(procpar["tn"]["values"][0])
    if nucleus is None:
        nucleus = "1H"
    resppm = 4.7 if nucleus == "1H" else 0.0


    params = {'sw'          : float(procpar["sw"]["values"][0]),
              'frequency'   : float(procpar["sfrq"]["values"][0]),
              'resppm'      : resppm,
              'echopeak'    : 0.0,
              'nucleus'     : nucleus,
              'seqte'       : float(procpar["te"]["values"][0]) * 1000.0,
              'seqtr'       : float(procpar["tr"]["values"][0]) * 1000.0,
              'voxel_dimensions' : [20.0,20.0,20.0],
              'header'      : '',
              'transform'   : None,
              'data'        : None}




    return params

