# Python modules
import os

# 3rd party modules
import pydicom

# Our modules
import vespa.common.util.dir_list as dir_list

 

# MIN & MAX PIXEL are limited by what can be stuffed into a byte, because
# we turn the DICOM data into RGB data with R == G == B for each pixel.
_MAX_PIXEL = 255
_MIN_PIXEL = 0


class DicomFileInfo(object):
    """
    A convenience container for a small subset of attributes from a DICOM file.

    The advantages of this over a pydicom Dataset object are --
    1) This uses Pythonic names
    2) All attributes are guaranteed to exist (although they may be empty)
       so no need to test "if 'foo' in dataset".
    3) This offers flexible sorting. Lists of DicomFileInfo objects will sort on
       the attributes listed in the class-level attribute 'sort_attributes'.
       To sort on e.g. study_id and series_number:

           files = get_all_files(path)
           DicomFileInfo.sort_attributes = ('study_id', 'series_number')
           files.sort()

       One cannot sort on arbitrary dataset objects, only on instance attributes.

    """
    sort_attributes = ('patient_name','study_id','series_number','instance_number')

    def __init__(self, dataset, filename=""):
        """
        The dataset param must be a pydicom dataset. The filename is
        stored in this object for convenience but is otherwise unused.

        'slice_location' is calculated from ImagePositionPatient (0020 0037)
        and ImageOrientationPatient (0020 0037) if present according to the
        algorithm outlined at great length by Jolinda Smith here:
          http://www.itk.org/pipermail/insight-users/2003-September/004762.html
          http://www.cmake.org/pipermail/insight-users/2003-September/004762.html
        If those links die you can beseech Google:
          http://www.google.com/search?q=%22slightly-rotated+coronal+acquisition%22

        """
        self.dataset            = dataset
        self.filename           = filename
        self.patient_name       = dataset.PatientName       # 0010 0010
        self.patient_id         = dataset.PatientID         # 0010 0020
        self.study_instance_uid = dataset.StudyInstanceUID  # 0020 000D
        self.study_id           = ''                        # 0020 0010
        self.study_description  = ''                        # 0008 1030
        self.image_id           = 0                         # 0054 0400
        self.instance_number    = 0                         # 0020 0013
        self.series_number      = 0                         # 0020 0011
        self.series_description = ''                        # 0008 103E
        self.series_instance_uid = ''                       # 0020 000E
        self.echo_time          = 0.001                     # 0018 0081
        self.acquisition_time   = 0                         # 0008 0032

        # Here I populate the not-always-present values
        if 'ImageType' in dataset:         self.image_type          = dataset.ImageType         # 0008 0008
        if 'StudyID' in dataset:           self.study_id            = dataset.StudyID
        if 'ImageID' in dataset:           self.image_id            = dataset.ImageID
        if 'SeriesNumber' in dataset:      self.series_number       = dataset.SeriesNumber
        if 'StudyDescription' in dataset:  self.study_description   = dataset.StudyDescription
        if 'SeriesInstanceUID' in dataset: self.series_instance_uid = dataset.SeriesInstanceUID
        if 'SeriesDescription' in dataset: self.series_description  = dataset.SeriesDescription
        if 'InstanceNumber' in dataset:    self.instance_number     = dataset.InstanceNumber
        if 'EchoTime' in dataset:          self.echo_time           = dataset.EchoTime
        if 'AcquisitionTime' in dataset:   self.acquisition_time    = dataset.AcquisitionTime

        self.slice_location = 0
        if ('ImagePositionPatient'     in dataset)  and \
            ('ImageOrientationPatient' in dataset)  and \
            dataset.ImagePositionPatient            and \
            (len(dataset.ImageOrientationPatient) >= 6):

            o = dataset.ImageOrientationPatient
            slice_normals = [ (o[1] * o[5]) - (o[2] * o[4]),
                              (o[2] * o[3]) - (o[0] * o[5]),
                              (o[0] * o[4]) - (o[1] * o[3]), ]

            self.slice_location = sum([a * b for a, b in zip(slice_normals, dataset.ImagePositionPatient)])

    # bjs June 2020, Py3 no longer supports __cmp__() method, have to do all comparisons below

    def __lt__(self, other):
        """ Compares based on the class attribute sort_attributes"""
        item1, item2 = self.get_cmp_lists(other)
        return item1 < item2

    def __gt__(self, other):
        """ Compares based on the class attribute sort_attributes"""
        item1, item2 = self.get_cmp_lists(other)
        return item1 > item2

    def __eq__(self, other):
        """ Compares based on the class attribute sort_attributes"""
        item1, item2 = self.get_cmp_lists(other)
        return item1 == item2

    def __ne__(self, other):
        """ Compares based on the class attribute sort_attributes"""
        item1, item2 = self.get_cmp_lists(other)
        return item1 != item2

    def get_cmp_lists(self, other):
        item1 = [getattr(self, item) for item in self.sort_attributes]
        item2 = [getattr(other, item) for item in self.sort_attributes]
        return item1, item2


    def get_image(self):
        """Returns a wx.Image containing the image data in the dataset.
        Only a small subset of the image types described in the DICOM
        specification are supported. Specifically --
        - the dataset's PixelData attribute must be populated
        - the image must be grayscale
        - the data must be 8, 16, or 32-bit
        - the transfer syntax UID must indicate raw mode (big and little
          endian are both OK)

        If this method encounters something it doesn't understand, it raises
        a NotImplementedError.
        """
        import wx

        # I create a shortcut reference for the dataset
        ds = self.dataset

        if 'PixelData' not in ds:
            return None
        else:
            # The Transfer Syntax UID is prefixed with a noise string
            # (1.2.840.10008.1.) that's followed by the parts that matter. In
            # this case, '2' and '2.1' ==> raw little endian and '2.2' ==>
            # raw big endian.
            # ref: DICOM spec, Part 6: Data Dictionary, Annex A, Table A-1
            RAW_UIDS = ('1.2.840.10008.1.2', '1.2.840.10008.1.2.1',
                        '1.2.840.10008.1.2.2')

            # Turn the pydicom string into a list of integers
            pixels = [ord(c) for c in ds.PixelData]

            if ds.TransferSyntaxUID not in RAW_UIDS:
                raise NotImplementedError('Unsupported: TransferSyntaxUID (0002,0010) == %d' % ds.TransferSyntaxUID)
            else:
                # Data is in raw mode
                little_endian = (ds.TransferSyntaxUID != '1.2.840.10008.1.2.2')

                if ds.BitsAllocated not in (8, 16, 32):
                    raise NotImplementedError('Unsupported: BitsAllocated (0028,0100) == %d' % ds.BitsAllocated)

                # dataset.BitsAllocated is 8, 16, or 32
                bytes_per_pixel = ds.BitsAllocated / 8

                if bytes_per_pixel == 1:
                    # No need to massage the data
                    pass
                elif bytes_per_pixel == 2:
                    # Combine the byte sets into ints
                    for i in range(0, len(pixels) // 2):
                        if little_endian:
                            pixels[i] = (pixels[(i * 2) + 1] << 8) |   \
                                         pixels[ i * 2     ]
                        else:
                            # big endian
                            pixels[i] = (pixels[(i * 2)    ] << 8) |   \
                                         pixels[(i * 2) + 1]

                    # The second half of the data is now trash
                    pixels = pixels[:len(pixels) // 2]
                elif bytes_per_pixel == 4:
                    for i in range(0, len(pixels) // 4):
                        if little_endian:
                            pixels[i] = (pixels[(i * 4) + 3] << 24) |  \
                                        (pixels[(i * 4) + 2] << 16) |  \
                                        (pixels[(i * 4) + 1] <<  8) |  \
                                         pixels[(i * 4)    ]
                        else:
                            # big endian
                            pixels[i] = (pixels[(i * 4)    ] << 24) |  \
                                        (pixels[(i * 4) + 1] << 16) |  \
                                        (pixels[(i * 4) + 2] <<  8) |  \
                                         pixels[(i * 4) + 3]

                    # All but the first quarter of the data is now trash
                    pixels = pixels[:len(pixels) // 4]

            # Now pixels is a list of ints and we can forget about bits per
            # pixel and endianness.

            if ('PhotometricInterpretation' in ds) and \
               (not ds.PhotometricInterpretation.startswith('MONOCHROME')):
                # We only handle monochrome images
                raise NotImplementedError('Unsupported: PhotometricInterpretation (0028,0004) == %d' % ds.PhotometricInterpretation)
            else:
                # If Photometric Interpretation (0028,0004) isn't present,
                # I *assume* MONOCRHOME1 or MONOCRHOME2.
                if ('WindowCenter' in ds) and ('WindowWidth' in ds):
                    # Normalize pixels based on the algorithm presented in the
                    # DICOM spec (see the function for details)
                    normalize_pixels(pixels, ds.WindowCenter, ds.WindowWidth)

                if ('PhotometricInterpretation' in ds) and \
                   (ds.PhotometricInterpretation == 'MONOCHROME1'):
                    # MONOCHROME1 ==> low values are bright and high values
                    # are dark. In other words, it is the opposite of what
                    # most images libraries expect.
                    pixels = [_MAX_PIXEL - n for n in pixels]

                # Turn the pixel values into a string in preparation for
                # passing them to the image library.
                # wx.Image wants RGB data which is why I create three copies
                # of each character.
                pixels = "".join([chr(int(round(n))) * 3 for n in pixels])

                image = wx.EmptyImage(ds.Rows, ds.Columns)
                image.SetData(pixels)

            return image



def normalize_pixels(pixels, window_center, window_width):
    """
    Applies a linear conversion to the pixels based on the Window
    Center (0028,1050) and Window Width (0028,1051) values of the dataset.
    Pixels must be an iterable (e.g. list); it is altered in-place.
    The pixels returned are floating point values between _MIN_PIXEL and
    _MAX_PIXEL.

    This implements the algorithm presented in the DICOM specification
    Part 3: Information Object Definitions, Annex C.11.2.1.2.
    ftp://medical.nema.org/medical/dicom/2008/08_03pu.pdf

    Here's the pseudocode from the DICOM spec:
        if      (x <= c - 0.5 - (w-1)/2), then y = ymin
        else if (x >  c - 0.5 + (w-1)/2), then y = ymax
        else y = ((x - (c - 0.5)) / (w-1) + 0.5) * (ymax - ymin) + ymin

    The Python code below is the same as the pseudocode with some loop
    invariants calculated here to improve performance.

    """
    min_threshold = (window_center - 0.5 - ((window_width - 1) / 2))
    max_threshold = (window_center - 0.5 + ((window_width - 1) / 2))
    output_range = _MAX_PIXEL - _MIN_PIXEL

    window_width  -= 1.0
    window_center -= 0.5

    for i, value in enumerate(pixels):
        if value <= min_threshold:
            pixels[i] = _MIN_PIXEL
        elif value > max_threshold:
            pixels[i] = _MAX_PIXEL
        else:
            pixels[i] = ((((value - window_center) / window_width) + 0.5) * output_range) + _MIN_PIXEL


def is_dicom(filename):
    """
    Returns True if the file in question is a DICOM file, else False. 
    
    - a DICOM file starts with 128 reserved bytes followed by 'DICM'.
    - ref: DICOM spec, Part 10: Media Storage and File Format for Media 
           Interchange, 7.1 DICOM FILE META INFORMATION 
    """
    if os.path.isfile(filename):
        s = open(filename, 'rb').read(132)
        if isinstance(s,(bytes, bytearray)):
            try:
                s = s.decode('utf-8')
            except:
                # bugfix - was trying to read PNG and found a char utf8 did not like
                try:
                    s = s.decode('utf-16')
                except:
                    return False
        return s.endswith('DICM')
    else:
        return False


def get_all_files(path):
    """
    Gets all DICOM files in the directory indicated by the path and
    returns them as an unsorted list of DicomFileInfo instances.
    To sort the files, call the sort() method of the returned list.
    """
    return [f for f in get_files(path)]


def get_files(path, manufacturer=['SIEMENS',]):
    """
    Return DicomFileInfo instances (or a subclass thereof) for the DICOM files
    in the directory indicated by the path. This is a Python generator so they
    are returned one by one.
    
    Example, use it in a loop like so:

        for dicom_file in get_files(the_path):
            do_something(dicom_file)

    See also get_all_files().

    """
    filenames = dir_list.DirList(path, dir_list.FilterFile).fq_names

    for filename in filenames:
        if is_dicom(filename):
            dataset = pydicom.dcmread(filename)
            dataset.decode()                        # change strings to unicode
            df = DicomFileInfo(dataset, filename)
            yield df
        #else:
            # Not a DICOM file; ignore it.

