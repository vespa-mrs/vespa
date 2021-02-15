# Python modules
import struct

# 3rd party modules
try:
    import dicom
except:
    import pydicom as dicom


# Our modules
import vespa.common.util.fileio as util_fileio


def vax_to_ieee_single_float(data):
    """
    This was mainly for SPAR/SDAT parsing - added this after I put the DICOM
    code in this module as well.
    
    Converts floats in Vax format to IEEE format.

    data should be a single string of chars that have been read in from 
    a binary file. These will be processed 4 at a time into float values.
    Thus the total number of byte/chars in the string should be divisible
    by 4.
    
    Based on VAX data organization in a byte file, we need to do a bunch of 
    bitwise operations to separate out the numbers that correspond to the
    sign, the exponent and the fraction portions of this floating point
    number
    
    role :      S        EEEEEEEE      FFFFFFF      FFFFFFFF      FFFFFFFF
    bits :      1        2      9      10                               32
    bytes :     byte2           byte1               byte4         byte3    
    
    """
    f = []
    nfloats = int(len(data) / 4)
    for i in range(nfloats):
        byte2 = data[0 + i*4]
        byte1 = data[1 + i*4]
        byte4 = data[2 + i*4]
        byte3 = data[3 + i*4]
        
        # hex 0x80 = binary mask 10000000
        # hex 0x7f = binary mask 01111111

        if type(data) == type(b'1'):
            sign  =  (byte1 & 0x80) >> 7
            expon = ((byte1 & 0x7f) << 1 )  + ((byte2 & 0x80 ) >> 7 )
            fract = ((byte2 & 0x7f) << 16 ) +  (byte3 << 8 ) + byte4
        else:
            sign  =  (ord(byte1) & 0x80) >> 7
            expon = ((ord(byte1) & 0x7f) << 1 )  + ((ord(byte2) & 0x80 ) >> 7 )
            fract = ((ord(byte2) & 0x7f) << 16 ) +  (ord(byte3) << 8 ) + ord(byte4)
        
        if sign == 0:
            sign_mult =  1.0
        else:
            sign_mult = -1.0
        
        if 0 < expon:
            # note 16777216.0 == 2^24  
            val = sign_mult * (0.5 + (fract/16777216.0)) * pow(2.0, expon - 128.0)   
            f.append(val)
        elif expon == 0 and sign == 0:
            f.append(0)
        else: 
            f.append(0)
            # may want to raise an exception here ...
        
    return f




"""
#
# Shared by Sandeep Ganji by email Apr 17, 2019
# - collaborate on running Vespa-analysis inline on Philips scanner
#

Can read the file with dicom.read_file() or pydicom.dcmread()

Can use defer_size tag to not read large data initially. Note. if you try to
  access this 'deferred read' tag later, it will be read and evaluated then.
  
Can use specific_tags keyword to only retrieve certain data from entire file

Example:
        
    my_tags = [dicom.tag.Tag((0x0008, 0x0031)), 
               dicom.tag.Tag((0x2005, 0x1270)),
               dicom.tag.Tag((0x2005, 0x10c0))]
        
    dataset = dicom.read_file(filename) #, defer_size="1 KB")

    Note. At this point the data tag (0x2005,0x1270) has had its raw data read
    in but has not been evaluated. The data tag value is just a byte string.
        
        tag = dicom.tag.Tag(0x2005,0x1270)
        data_elem1 = dataset._dict[tag]         # This is tuple with raw bytes in value
        
        tmp = len(dataset[0x2005,0x1270].value) # force object to evaluate this tag
        data_elem3 = dataset._dict[tag]         # This is DataElement with evaluated numbers in value

"""

TAG_SOP_CLASS_UID      = (0x0008, 0x0016)
TAG_PHILIPS_PRIVATE_01 = (0x2005, 0x1270)
TAG_PHILIPS_PRIVATE_02 = (0x2005, 0x10c0)

# This tag holds a description of the content of a Siemens file
TAG_CONTENT_TYPE  = (0x0029, 0x1008) 

# SG_phan_20190418_Braino_Aprs28_m_8_2_raw_act.SDAT
# 4096    4    Complex (FID) - water suppressed metabolite
#
# SG_phan_20190418_Braino_Aprs28_m_8_2_raw_ref.SDAT
# 4096    4    Complex (FID) - water 




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
        return s.endswith("DICM")
    else:
        return False


def is_mrs_dicom(dataset):
    """ returns True if all criteria (Sandeep) are met """
    
    if not isinstance(dataset, dicom.dataset.Dataset):
        raise ValueError("Object passed in not a dicom Dataset.") 
    
    if not "ProtocolName" in dataset:
        return False
    
    if dataset.ProtocolName == 'ExamCard':
        return False
    
    if not (0x2005,0x10c0) in dataset:
        return False

    if dataset[0x2005,0x10c0].value != 'SPECTRO':
        return False    

    if not (0x2005, 0x1270) in dataset:
        return False

    return True    

    
    
def get_dicom_file_type(dataset):
    """ returns 'spectro','protocol', or 'other' based on Sandeeps criteria """
    
    if not isinstance(dataset, dicom.dataset.Dataset):
        raise ValueError("Object passed in not a dicom Dataset.") 
    
    if not "ProtocolName" in dataset:
        return 'other'
    
    if dataset.ProtocolName == 'ExamCard':
        return 'other'
    
    if not (0x2005,0x10c0) in dataset:
        return 'other'

    if dataset[0x2005,0x10c0].value != 'SPECTRO':
        return 'protocol'    

    return 'spectro'    



def filter_file_types(filenames):
    
    files_spectro   = []
    files_protocol  = []
    files_other     = []
    
    # only get subset to speed up dataset reads
    my_tags = [dicom.tag.Tag((0x0008, 0x0020)),     # 'Study Date' just because
               dicom.tag.Tag((0x0018, 0x1030)),     # 'Protocol Name'
               dicom.tag.Tag((0x2005, 0x1270)),
               dicom.tag.Tag((0x2005, 0x10c0))]    
    
    for filename in filenames:
        
        if not is_dicom(filename):
            continue
        
        dataset = dicom.read_file(filename, defer_size="1 KB")  #, specific_tags=my_tags)
        
        label = get_dicom_file_type(dataset)
        
        if label == 'spectro':
            files_spectro.append(filename)
        elif label == 'protocol':
            files_protocol.append(filename)
        elif label == 'other':
            files_other.append(filename)
            
    return files_spectro, files_protocol, files_other
        

def get_dicom_mrs_data(dataset):

    data = []

    if not isinstance(dataset, dicom.dataset.Dataset):
        raise ValueError("Object passed in not a dicom Dataset.") 
    
    if not is_mrs_dicom(dataset):
        raise ValueError("Dataset does not have MRS data.")
    
    tag_data = dicom.tag.Tag(0x2005,0x1270)
    
    data_elem = dataset._dict[tag_data]         # This should be a tuple with raw bytes in value

    if isinstance(data_elem, dicom.dataelem.DataElement):
        raise ValueError("RawDataElement not available in Dataset, Element has previously been accessed.")
    
    
    if isinstance(data_elem, tuple):
        # Big simplifying assumptions --
        # 0) Unpack byte stream into a series of float32 values
        # 1) Floats a series of complex numbers organized as ririri...
        #    where r = real and i = imaginary.
        # 2) Each real & imaginary number is a 4 byte float.
        # 3) Data is little endian.

        data = data_elem.value
        data = struct.unpack("<%df" % (len(data) / 4), data)

        # fyi this gave the same result as my brute force above   
        compare = dataset[0x2005,0x1270].value
        
        data = util_fileio.collapse_complexes(data)

        
    return data
        
        
        
    
def parse_philips_dicom_data(filename):

    data = []

    dataset = dicom.read_file(filename) 
        
    try:
        data = get_dicom_mrs_data(dataset)
    except ValueError as e:
        print(e)

    return data
