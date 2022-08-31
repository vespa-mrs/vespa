# Python modules
import os
import shutil

# 3rd party modules
import numpy as np
import pydicom as dicom


# Our modules
import vespa.common.constants  as common_constants

# List Classes then Methods ...
# - keep 'em in alphabetical order, please, unless you have a better idea.

class DicomFileItems(object):

    def __init__(self, filename):
        self.filename = filename
        self.move_path = None
        self.patient_id = None
        self.patient_name = None
        self.study_id = None
        self.study_description = None
        self.series_number = None
        self.series_description = None
        self.instance_number = None

    def __lt__(self, other):
        """
        Python docs say that, "The sort routines are guaranteed to use
        __lt__() when making comparisons between two objects. So, it is easy
        to add a standard sort order to a class by defining an __lt__() method.

        For DICOM series sort, we typically sort by the patient ID, then the
        study id, then series number attributes.

        """
        original_reverse = [self.patient_id, self.study_id, self.series_number, self.instance_number]
        other_reverse = [self.patient_id, self.study_id, self.series_number, self.instance_number]
        original_reverse.reverse()
        other_reverse.reverse()

        return original_reverse < other_reverse


# def calculate_area(dimension_size, frequencies, phases, x_left, x_right):
#     """
#     Calculates & returns the selected area.
#
#     The param dimension_size is a scalar like data.dims[0].
#     The param phases is a tuple of info.ph0 and info.ph1.
#     The params x_left and x_right are translated from event coordinates to
#     native values.
#
#     """
#     phase0 = phases[0] * common_constants.DEGREES_TO_RADIANS
#     phase1 = phases[1] * common_constants.DEGREES_TO_RADIANS * \
#         (np.arange(dimension_size, dtype=float) - (dimension_size / 2))
#     phase1 /= dimension_size
#     phase  = np.exp(complex(0, 1) * (phase0 + phase1))
#     pdata  = (frequencies * phase).real[::-1]
#
#     if x_right >= x_left:
#         area = sum(pdata[x_left:x_right + 1])
#     else:
#         area = sum(pdata[x_right:x_left + 1])
#
#     return area


def dicom_series_sort(base_path):

    # Enumerate the directory contents
    filenames = []
    for root, dirs, fnames in os.walk(base_path):
        for fname in fnames:
            # Add path to filenames
            filenames.append(os.path.join(root, fname))

    # Filter out non-files
    filenames = [filename for filename in filenames if os.path.isfile(filename)]

    headers = []
    move_paths = []

    # assume that files MAY be in one dir, or in sub-dirs under the
    # base_path that is searched (and passed in here).

    for filename in filenames:
        # go through all filenames and ensure each file is DICOM
        if is_dicom(filename):

            # create storage container for data and refresh control flags
            hdr = DicomFileItems(filename)

            data_raw = dicom.dcmread(filename, stop_before_pixels=True)

            hdr.patient_id = data_raw.PatientID
            hdr.study_id = data_raw.StudyID
            hdr.series_number = data_raw.SeriesNumber
            hdr.series_description = data_raw.SeriesDescription
            hdr.instance_number = data_raw.InstanceNumber
            hdr.study_description = data_raw.StudyDescription
            hdr.patient_name = str(data_raw.PatientName)

            if hdr.series_number < 10:
                ser_str = '000' + str(hdr.series_number)
            elif hdr.series_number < 100:
                ser_str = '00' + str(hdr.series_number)
            elif hdr.series_number < 1000:
                ser_str = '0' + str(hdr.series_number)
            else:
                ser_str = str(hdr.series_number)

            # concatenate the destination sub-directory name
            move_path = ser_str + '_'
            move_path += str(hdr.series_description).upper()
            # move_path += str(hdr.patient_id)+'_'
            # move_path += str(hdr.study_id)+'_'
            move_path = move_path.replace(" ", "_")

            # creae a list of unique sub-directory names
            if move_path not in move_paths:
                move_paths.append(move_path)

            hdr.move_path = move_path

            headers.append(hdr)

    # create all "move directories"
    #
    # Note. according to Python docs, the mkdir command
    #   only works for Windows and *nix, no Mac support
    for fdir in move_paths:
        new_path = base_path + '/' + fdir
        if not os.path.exists(new_path):
            os.mkdir(new_path)
            #print('New path - '+new_path)

    # move all files to new sub-directories
    #
    # used shutil here to be pure Python and not have to
    # worry about OS specific command names.
    for hdr in headers:
        src = hdr.filename
        path, fname = os.path.split(src)
        dst = os.path.join(base_path, hdr.move_path, fname)
        shutil.move(src, dst)
        #print('Dst file - '+dst)


def is_dicom(filename):
    """Returns True if the file in question is a DICOM file, else False. """
    # Per the DICOM specs, a DICOM file starts with 128 reserved bytes
    # followed by "DICM".
    # ref: DICOM spec, Part 10: Media Storage and File Format for Media
    # Interchange, 7.1 DICOM FILE META INFORMATION
    if os.path.isfile(filename):
        f = open(filename, "rb")
        s = f.read(132)
        f.close()
        return s.endswith(b"DICM")
    else:
        return False

    
    
    
    
    
    


