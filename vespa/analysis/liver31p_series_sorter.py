# Python modules
import os
import shutil

# 3rd party modules
import pydicom

# Our modules
import vespa.analysis.fileio.dicom_siemens as dicom_siemens


TAG_SOP_CLASS_UID = (0x0008, 0x0016)
# This tag holds a description of the content of a Siemens file
TAG_CONTENT_TYPE  = (0x0029, 0x1008)



def is_dicom(filename):
    """
    Returns True if the file in question is a DICOM file, else False. """
    # Per the DICOM specs, a DICOM file starts with 128 reserved bytes
    # followed by "DICM".
    # ref: DICOM spec, Part 10: Media Storage and File Format for Media 
    # Interchange, 7.1 DICOM FILE META INFORMATION 
    if os.path.isfile(filename):
        s = open(filename, "rb").read(132)
        return s.endswith("DICM")
    else:
        return False
    

def get_files(filenames):
    """
    Returns the DICOM files in the directory indicated by the path. The 
    objects returned are DicomFileInfo instances (or a subclass thereof) and 
    they're returned one by one. 
    
    This is a Python generator, which means you're expected to use it in a 
    loop like so:
    for dicom_file in get_files(the_path):
        do_something(dicom_file)

    See also get_all_files().
    """

    for filename in filenames:
        if is_dicom(filename):
            dataset = pydicom.read_file(filename)
            
            # change strings to unicode
            dataset.decode()
            
            yield filename, dataset


def file_filter_mrs(dataset):
    """
    This function is called for each file to see whether or not it 
    should be included in the tree. It returns True if the given DICOM 
    file is a Siemens spectroscopy file.
    """
    
    if TAG_SOP_CLASS_UID in dataset:
        item_uid = dataset[TAG_SOP_CLASS_UID].value.upper() 
        sop_class_uid = pydicom.uid.UID(str(item_uid))
        if sop_class_uid.name == 'MR Image Storage':
            # this is the condition for Spectroscopy Browser
            return False
        elif sop_class_uid.name == 'MR Spectroscopy Storage':
            # this is the UID value for Siemens VD files
            return True
        elif sop_class_uid.name == str(item_uid):
            # In software versions VA, VB (maybe others, but not VD),
            # Siemens uses a proprietary SOP Class UID for their spectroscopy
            # data files. I have seen 1.3.12.2.1107.5.9.1 used but am not sure
            # if it is unique. And this is not in the DICOM dictionary of UIDs
            # So, we will look for other proprietary tags that will give us 
            # solid info that these are Siemens SPEC file.
            # 
            # pydicom doesn't know about the Siemens content type tag, so I have 
            # to pass the magic numbers that represent it.
            if TAG_CONTENT_TYPE in dataset:
                content_type = dataset[TAG_CONTENT_TYPE].value.upper()
            else:
                content_type = ""
            # There are lots of values that can appear in the Siemens content 
            # type tag. The ones I care about are "SPEC NUM 4" and "Spectroscopy".
            # ref:
            #    p 132 of this PDF:
            #    http://www.medical.siemens.com/siemens/en_GLOBAL/rg_marcom_FBAs/files/brochures/DICOM/mr/syngo_MR_B17.pdf
            #    which is from here:
            #    http://www.medical.siemens.com/webapp/wcs/stores/servlet/CategoryDisplay~q_catalogId~e_-11~a_categoryId~e_16560~a_catTree~e_100003,16554,16560~a_langId~e_-11~a_storeId~e_10001.htm
            
            return ("SPEC" in content_type)
        else:
            return False    


def do_sorting(path, extn='', verbose=False):

    # this gets all files *.IMA in all subdirectories of STARTDIR
    imafiles = []
    for dirpath, dirnames, filenames in os.walk(path):
        for filename in [f for f in filenames if f.endswith(extn.upper()) or f.endswith(extn.lower())]:
            imafiles.append(os.path.join(dirpath, filename))
            if verbose: print(os.path.join(dirpath, filename))

    series_indiv  = []
    series_sum    = []

    for filename, dataset in get_files(imafiles):

        if file_filter_mrs(dataset):
            

            item_uid = dataset[TAG_SOP_CLASS_UID].value.upper()
            sop_class_uid = pydicom.uid.UID(str(item_uid))
    
            if sop_class_uid.name == 'MR Spectroscopy Storage':
    
                # First we extract a bunch of items from the dataset and put them into
                # a dict keyed by names reminiscent of VASF parameters.
                parameters = dicom_siemens._extract_parameters_from_dataset_dicom_sop(dataset)
    
            else:
                # This is how we get parameters for VA and VB system data where
                # lots of things were stored in proprietary tags
                #
                # First we extract a bunch of items from the dataset and put them into
                # a dict keyed by names reminiscent of VASF parameters.
                parameters = dicom_siemens._extract_parameters_from_dataset_siemens_proprietary(dataset)
                
            if parameters['averages'] == 16:
                series_sum.append(filename)
            elif parameters['averages'] == 1:
                series_indiv.append(filename)
        
    if series_sum:
        directory = path+'\\series_sum'
        if not os.path.exists(directory):
            os.makedirs(directory)
        for filename in series_sum:
            temp = os.path.basename(filename)
            temp = os.path.join(directory, temp)
            if not os.path.isfile(temp):
                shutil.copy2(filename, directory)

    if series_indiv:
        directory = path+'\\series_indiv'
        if not os.path.exists(directory):
            os.makedirs(directory)
        for filename in series_indiv:
            temp = os.path.basename(filename)
            temp = os.path.join(directory, temp)
            if not os.path.isfile(temp):
                shutil.copy2(filename, directory)



def main():

    # Processing of liver 31p series 'raw' directories
    STARTDIR = 'D:\\Users\\bsoher\\temp\\current\\data'

    i = 0

    # this gets all files *.IMA in all subdirectories of STARTDIR
    imafiles = []

    # all paths under STARTDIR
#     paths = [x[0] for x in os.walk(STARTDIR)]
#     extn  = '.ima'

    # limited paths under STARTDIR
    paths = ['','\\raw_fruct-005','\\raw_fruct-006','\\raw_fruct-020']
    paths = [STARTDIR+item for item in paths]
    extn  = ''                                      # NB. This does NOT exclude '.ima' files, just allows non-'.ima' files

    verbose = False

    for path in paths[1:]:

        print("processing - "+path)
        do_sorting(path)
            



        
if __name__ == '__main__':
    main()        
        