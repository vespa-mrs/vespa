# Python modules
import copy


# 3rd party modules


# Our modules
from vespa.analysis.block_prep_fidsum import BlockPrepFidsum
from vespa.common.mrs_data_raw import DataRaw, DataRawFidsum
import vespa.analysis.fileio.util_exceptions as util_exceptions

# need for inline processing - no wx
try:
    from vespa.common.wx_gravy.common_dialogs import pickfile
except:
    pickfile = None

class RawReader(object):
    """
    A base class for classes that transform a specific format (e.g. Siemens
    DICOM or VASF) into a Vespa DataRaw object. It defines the minimal interface
    that subclasses must support and also provides some default implementations.

    Subclasses must override read_raw().

    Subclasses may want to override pickfile(), read_raws() and any of the
    class' "private" methods.

    """
    def __init__(self):
        self.filetype_filter = ""
        self.filenames = [ ]
        self.multiple = True
        self.raws = [ ]


    def pickfile(self, default_path=""):
        if pickfile is not None:
            self.filenames = pickfile( filetype_filter=self.filetype_filter,
                                       multiple=self.multiple,
                                       default_path=default_path)
            if not self.multiple:
                # If select one file, we get a string back. Turn it into a list.
                self.filenames = [self.filenames,] if self.filenames else []
        else:
            self.filenames = []
        return bool(self.filenames)


    def read_raw(self, filename, ignore_data, *args, **kwargs):
        # Subclasses must implement this. See template.py for instructions.
        raise NotImplementedError


    def read_raws(self, ignore_data=False, *args, **kwargs):
        """
        Calls read_raw() once for each filename in self.filenames.
        - returns either DataRaw or DataRawFidsum objects
        - returns a list with a single object or multiple objects
        - read_raw() method has to deal with the complexity of what data to
            read and how to organize it on return.

        There are five possible cases -

        Note. Return objects HAVE to be passed back in a list, regardless if
              there is only one, or more, objects to return.

        1) One filename returning ONE DataRaw object with one FID in it and data
           array shape [1,1,1,N], where N is the number of spectral points.

        2) Multiple filenames returning ONE DataRaw object with multiple FIDs in
           it and data array shape [1,1,M,N], where M is number of files. This
           case is for a 'stack of single voxel data sets' into the M dimension

        3) One filename returning ONE DataRawFidsum. The file contains multiple
           FIDs which are retuned in a DataRawFidsum object with data array shape
           [1,L,M,N] where M = number of averages, and L = 1 or 'number of coils'
           per average.

           Note. As the M dimension is used for 'coils', we can not 'stack into
           the screen'. Instead, open each file separately into a Dataset Tab.

        4) One filename containing MULTIPLE data sets each with multiple FIDs.
           At this point, read_raw() HAS to organize the collection of FIDs into
           separate objects, one for each data set. Return objects have to be
           all DataRaw or all DataRawFidsum object, you can not mix. Each object
           returned will open into its own Dataset Tab.
        
        5) (Not Implemented) Multiple filenames containing MULTIPLE data sets
           each with multiple FIDs. This is the same as 4, just more complicated.
           This method doesn't try to handle this case; it's too format-specific.
           If you want to implement this, you'll have to write your own version of
           this read_raws() method.

        If this is being called when one or more datasets is already open,
        then one of those open datasets must be passed in the open_dataset 
        param. The attributes of the raw files that this method opens must be 
        consistent with the ones that are already open. If they're not,
        this code raises MultifileAttributeMismatchError.

        Exceptions: MultifileAttributeMismatchError, UnsupportedDimensionalityError,
        IOError, SIDataError, and FileNotFoundError. Also, DataRaw.concatenate()
        can raise ValueError which is not trapped by this method.

        """
        self.raws = [ ]

        for filename in self.filenames:
            raw = self.read_raw(filename, ignore_data, *args, **kwargs)
            for item in raw: self.raws.append(item)

        nfiles = len(self.filenames)
        raw0 = self.raws[0]

        fidsum_flag = False
        if isinstance(raw0, DataRawFidsum) and not type(raw0) == DataRaw:
            fidsum_flag = True

        self._check_consistency(fidsum=fidsum_flag)
        self._check_for_si()

        if nfiles > 1:
            out = copy.deepcopy(self.raws[0])
            for raw in self.raws[1:]:
                out.concatenate(raw)
            out = [out,]
        else:
            out = self.raws

        if 'open_dataset' in list(kwargs.keys()):
            if kwargs['open_dataset'] is not None:
                self._check_compatibility(out, kwargs['open_dataset'], fidsum=fidsum_flag)

        return out


    ##################  Internal use only methods below   ##################
    #
    # These methods shouldn't be called by users of this class. However,
    # subclasses should feel free to override them as necessary.

    def _check_consistency(self, fidsum=False):
        """
        Ensures all new objects have same shape and parameters before they are
        concatenated into a single raw object

        """
        if len(self.raws) == 1 and len(self.filenames) == 1:
            pass
        elif len(self.raws) > 1:
            # Cases 2-4, all objects in the raws list must be the same type
            
            first = self.raws[0]
            if not all(type(item)==type(first) for item in self.raws):
                raise util_exceptions.MultifileTypeMismatchError

            # The attributes (data shape and sweep width) must match.
            for raw in self.raws[1:]:
                if fidsum:
                    dims_flag = (first.data_shape[-1] == raw.data_shape[-1])
                else:
                    dims_flag = (first.data_shape == raw.data_shape)
                
                if dims_flag and (abs(first.sw - raw.sw) < 0.001):
                    pass
                else:
                    raise util_exceptions.MultifileAttributeMismatchError
                    
            # Also ensure at least one dimension is 1. Datasets with all dims > 1
            # might be valid data, but Vespa doesn't support that.
            if min(first.data_shape) > 1:
                raise util_exceptions.UnsupportedDimensionalityError


    # def _check_dimensionality(self, open_dataset):
    #     """ Match new shape and sw to those of the open dataset(s) """
    #     for raw in self.raws:
    #         if not (raw.data_shape == open_dataset.raw_shape) and (raw.sw == open_dataset.sw):
    #             raise util_exceptions.OpenFileAttributeMismatchError
            
            
    def _check_compatibility(self, raws, open_dataset, fidsum=False):
        """
        Attributes and object type of new object must match those of the open dataset(s)

        """
        raw0 = self.raws[0]
        raw1 = open_dataset.blocks["prep"]

        new_is_fidsum = False
        if isinstance(raw0, DataRawFidsum) and not type(raw0) == DataRaw:
            new_is_fidsum = True

        open_is_fidsum = False
        if isinstance(raw1, BlockPrepFidsum):
            open_is_fidsum = True

        if not (new_is_fidsum == open_is_fidsum):
            raise util_exceptions.OpenFileTypeMismatchError

        for raw in raws:
            if fidsum:
                dims_flag = (raw.data_shape[-1] == open_dataset.raw_shape[-1])
            else:
                dims_flag = (raw.data_shape == open_dataset.raw_shape)

            if  dims_flag and (abs(raw.sw - open_dataset.sw) < 0.001):
                # Add the 'slight difference' equation for real world conditions. All is well!
                pass
            else:
                raise util_exceptions.OpenFileAttributeMismatchError
                   
                   
    def _check_for_si(self):
        """
        This is not a specific check for SI data, it just ensures that the
        current 4th dimension is not being used in all datasets.

        """
        for raw in self.raws:
            if self.raws[0].data_shape[0] != 1:
                raise util_exceptions.SIDataError

