# Python modules
from __future__ import division

# 3rd party modules
import numpy as np

# Our modules
import vespa.common.constants as constants
import vespa.common.ge_read_pfile as read_pfile

from vespa.analysis.fileio.ge_pfile import RawReaderGeProbep
from vespa.common.mrs_data_raw import DataRawProbep

from vespa.common.constants import DEGREES_TO_RADIANS, RADIANS_TO_DEGREES

NUMPY_DATA_TYPE = np.float32

BYTES_PER_ELEMENT = np.zeros(1, dtype=NUMPY_DATA_TYPE).nbytes



class RawReaderGeFidcsi13c(RawReaderGeProbep):
    # This inherits from RawReaderGeProbep (q.v.). The only methods you
    # need to implement are __init__() and read_raw(). You *may* want to 
    # override or supplement some of raw_reader.RawReader's other methods.

    def __init__(self):
        RawReaderGeProbep.__init__(self)
        
        self.filetype_filter = "GE FIDCSI for 13C (*.7, *.7.00)|*.7;*.7.00"
        self.multiple = False
        
        
    def read_raw(self, filename, ignore_data=False, open_dataset=None):
        """
        Given the name of a .7 file, returns a DataRaw object
        populated with the parameters and data in the file.

        When ignore_data is True, this function only reads the parameters file
        which can be much faster than reading both params & data, depending
        on the file format.
        
        The open_dataset attribute is not used in this reader. 
        
        """
        pfile = read_pfile.Pfile(filename)
        
        if not pfile.is_svs:
            # closest thing I have to a NotSvsData exception
            raise NotImplementedError
            
        ncoil = pfile.map.get_num_coils
            
        # This function will almost always return a single object that is a 
        # instance of DataRaw or one of its subclasses (like DataRawFidsum). 
        # In unusual cases (where a single file contains multiple datasets),
        # it will return a list of DataRaw instances (or subclasses).

        d = _extract_parameters(pfile.hdr)
        
        d["header"]     = "\n".join(pfile.dump_header_strarr())
        d["data_source"] = filename

        spec_pts = pfile.hdr.rhr_rh_frame_size
        shape = (1, 1, 1, spec_pts)

        # create DataRaw objects for Suppressed and Unsuppressed
        # data if they exist. Create

        if ncoil > 1:
            phases = pfile.map.phase_first_point_deg
            phases = phases[0] - phases
            phases = np.exp(1j * phases * DEGREES_TO_RADIANS)

        data_list = []

        # At the moment, the number of Suppressed and Unsuppressed data
        # sets do not read out correctly in the ge_pfile_mapper, so we will
        # do a hack here for now ...

        data = pfile.map.avg_unsuppressed
        if ncoil > 1:
            for i in range(data.shape[3]):
                data[0,0,0,i,:] *= phases[i]
            data = np.sum(data, axis=3)
        data.shape = shape
        d["data"] = data
        data_list.append( DataRawProbep(d) )
        
        return data_list



####################    Internal functions start here     ###############

def _extract_parameters(hdr):
    """
    Given the pfile header, extracts a few specific parameters 
    and returns a flat dict containing those parameters and their value.

    The returned dict is appropriate for passing to DataRaw.inflate().
    """
    d = { }

    # A pfile contains many, many attributes. We're only interested in a few.
    
    d["sw"]         = float(hdr.rhi_user0)               # in Hz
    d["frequency"]  = float(hdr.rhr_rh_ps_mps_freq)/1e7     # in MHz
    d["seqte"]      = float(hdr.rhi_te)/1000                # in ms
    d["flip_angle"] = float(hdr.rhi_mr_flip)                # in deg
        
        
    #  Use the mps freq and field strength to determine gamma which is 
    #  characteristic of the isotop:
    gamma = ( hdr.rhr_rh_ps_mps_freq*1e-7 ) / ( hdr.rhe_magstrength/10000.0 )
   
    if abs( gamma - 42.57 ) < 0.3:
        nucleus = "1H"
    elif abs( gamma - 10.7 ) < 0.3 :
        nucleus = "13C"
    elif abs( gamma - 17.2 ) < 0.3 :
        nucleus = "31P"
    else:
        nucleus = "1H"
        
    if nucleus == "1H":
        d["resppm"] = constants.DEFAULT_PROTON_CENTER_PPM
    else:
        d["resppm"] = constants.DEFAULT_XNUCLEI_CENTER_PPM

    d["nucleus"] = nucleus

    return d


