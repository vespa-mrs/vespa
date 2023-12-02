"""
This file is a template for those who want to implement a reader for 
a file format that Vespa doesn't know about. It contains plenty of comments 
to guide you.

Focus on the class method read_raw() which mainly asks you to populate a
Python dictionary with relevant metadata. If you can do that, you can get 
Vespa to understand your file format.

Don't neglect the rest of the file which contains other useful information.

You can also look at the other files in this directory for examples of real
implementations. The files dicom_siemens.py and siemens_rda.py are examples
of formats in which the metadata and data are in the same file. For examples
of formats where the data and metadata are in separate files, have a look at
philips_spar.py, varian.py and vasf.py.

Once you have finished creating a class to read your format's header and data,
Appendix B in the Analysis User Manual can tell you how to let the program know
that it should add this format to its data import options. The Analysis User 
Manual is available from Analysis's Help menu.

"""
# Python modules

# 3rd party modules
import numpy as np

# Our modules
import vespa.analysis.fileio.raw_reader as raw_reader
import vespa.common.mrs_data_raw as mrs_data_raw
import vespa.common.util.fileio as util_fileio



class RawReaderAcme(raw_reader.RawReader):
    """
    This inherits from raw_reader.RawReader (q.v.). The only methods you
    have to implement are __init__() and read_raw(). You *may* want to
    override or supplement some of raw_reader.RawReader's other methods.

    Change the class name to something that makes sense to you. Our example
    reads the fictitious Acme format.

    """

    def __init__(self):
        """
        Set filetype_filter to the type of files you want to load.
        Set multiple to False if opening only single, not multiple, files

        """
        raw_reader.RawReader.__init__(self)
        self.filetype_filter = "Acme (*.acme)|*.acme;"
        self.multiple = False
        
        
    def read_raw(self, filename, ignore_data=False, *args, **kwargs):
        """
        Given the name of a .acme file, returns a DataRaw object
        populated with the parameters and data in the file.

        - 'ignore_data', Bool: default=False, only reads in parameters not data
             which can be much faster than reading both params & data, depending
             on the file format. It's OK to ignore the 'ignore_data' option, just
             document if it works or not.
        
        - When there are datasets already open in Vespa-Analysis, they can be
             accessed through the kwarg['open_dataset'] attribute

        This function creates a dictionary and populates it with parameters and
        the data. In the final step, the dictionary is turned into a Vespa
        DataRaw object and returned.

        This function returns a list of one or more objects that are instances
        of DataRaw or one of its subclasses, like DataRawFidsum.

        The purpose of this function is to allow you to read the parameters from
        a file, but you can hardcode them if you like. Here are the parameters
        you must populate:

        data_source: A string describing where this data came from. It's
          typically a filename, but it's free form so it could also be
          'Study database' or 'Direct from scanner'.  As with the header,
          what goes in here is entirely up to you.

        sw: sweep width in Hertz as a float, e.g. 2000.0

        frequency: frequency in MHz as a float, e.g 64.0

        resppm: float, acquisition on-resonance value in PPM

        echopeak: float [0.0 - 1.0], acquisition spectral echo peak as
          a percentage of acquistion first dimension length. A value of 0.0
          is a FID acquisition.

        nucleus: string, e.g. '1H'. This is *not* free form. Common values
          that Vespa recognizes are 1H, 13C, and 31P. The complete list is

        seqte: float, sequence echo time (TE) in milliseconds

        seqtr: float, sequence repetitions time (TR) in milliseconds

        voxel_dimensions, list, floats: SVS voxel size in mm,
          eg. [20.0,20.0,20.0]

        header: A text string containing all metadata. Vespa stores this
          but never reads it, so you can put anything you want in there. It's
          free form. Can also be and empty string, ''.

        transform, ndarray, 4x4 floats: The SVS to scanner transform for
          calculating SVS voxel overlay with MRI data sets. Set to None to
          trigger the default transform.

        data, ndarray, complex: Data must have a shape which is typically a
          4-tuple where the last element is the number of spectral points in
          the data. The other values (dimensions) in the shape should be 1 if
          you're only reading one set of data at a time.

        """
        d = { }

        shape = (1, 1, 1, 2048)

        if ignore_data:
            complex_data = np.zeros(shape, np.complex63)    # Create zero data
        else:
            complex_data = _read_data(filename)
            complex_data.shape = shape

        params = {'data_source': filename,
                  'sw': 2000.0,
                  'frequency': 127.9,
                  'resppm': 4.7,
                  'echopeak': 0.0,
                  'nucleus': '1H',
                  'seqte': 30.0,
                  'seqtr': 1500.0,
                  'voxel_dimensions': [20.0,20.0,20.0],
                  'header': 'empty',
                  'transform': None,
                  'data': complex_data}


        # This is typically the last thing that's added to the dictionary.
        d["data"] = data

        # Here's where the dictionary becomes a Vespa DataRaw object.
        raw = mrs_data_raw.DataRaw(d)

        # you have to return objects in a list, even if just one object
        return [raw,]



####################    Internal functions start here     ###############

def _read_data(filename):
    """
    Reads data from a .acme file and returns it in a numpy array. Note that
    the array is unshaped (one dimensional).

    It's up to you to write this function! This is an example for reading an
    extremely simple format. In this format, the data is written as a series
    of 32-bit floats in native byte order.

    """
    import struct

    # read data from binary file into a list of floats
    data = open(filename, "rb").read()
    element_count = len(data) // struct.calcsize('f')
    try:
        data = struct.unpack('f' * element_count, data)
    except struct.error:
        msg = "Unexpected input encountered while reading raw data"
        raise util_fileio.UnreadableDataError(msg)

    # convert list of floats into numpy array of complex floats
    data_iter = iter(data)
    data = [complex(r, i) for r, i in zip(data_iter, data_iter)]
    complex_data = np.fromiter(data, np.complex64)

    return complex_data

