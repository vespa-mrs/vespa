# Python modules


# import packages/libraries
import GERecon
import os
import shutil
import numpy as np
#import matplotlib.pyplot as plt        # this line causes SegFault on Linux

import inline_vespa_ge as ivg

"""
Version History
--------------------------

v1 - testing in linux so we can have access to GERecon library v1.0. Calling as a direct script, no __main__

"""



# Set the input Pfile
# Data1: Phase[1] -> Slices[4] -> Echo[1] -> Channel[1]
filename =  "/home/ese_user/orchestra/orchestra-python-sdk-1.0.x86_64/testdata/P12288.7"

# Create a pfile object
pfile = GERecon.Pfile()

# Load the Pfile File
status = pfile.Load(filename)  # Additional Options: Anon-mode->'No-Anonymize', Readmode -> 'All-Available-Acquisitions'

# Create Transform object
transformer = GERecon.Transformer()

# Create Orientation object
orient = GERecon.Orientation()

# Create DICOM object
dicom = GERecon.Dicom()

# Set the Path to save the output DICOM images.
pathToSaveDicoms = "/home/ese_user/orchestra/orchestra-python-sdk-1.0.x86_64/DICOM_OUT"

# Clear the Output directory before saving the DICOM images
if os.path.isdir(pathToSaveDicoms):
    shutil.rmtree(pathToSaveDicoms, ignore_errors=True)

# Check if pfile load is successful for further processing
if status == 0:

    # Get the metadata from the loaded Pfile
    metadata = pfile.MetaData()
    print(("MetaData: ", metadata))

    acqXRes = metadata['acquiredXRes']  # spectral points
    acqYRes = metadata['acquiredYRes']  # number of FIDs
    imageXRes = metadata['imageXRes']   # 512
    imageYRes = metadata['imageYRes']   # 512
    numPhases = metadata['phases']  # 1 here
    numSlices = metadata['slices']  # 1 here
    numEchoes = metadata['echoes']  # 1 here
    numChannels = metadata['channels']  # 8 here
    slice = 0

    print(('AcqXRes is {}, AcqYRes is {}, ImageXRes is {}, ImageYRes is {}, numChannels is {}, numEchoes is {}, numSlices is {}, numPhases is {}'.format(acqXRes, acqYRes, imageXRes, imageYRes, numChannels, numEchoes, numSlices, numPhases)))

    # Pre-Initialize the kSpace and image arrays
    kSpace = np.zeros([acqXRes, acqYRes], dtype=np.complex64)
    image = np.zeros([imageXRes, imageYRes, numChannels], dtype=np.complex64)
    combinedImage = np.zeros([imageXRes, imageYRes], dtype=np.complex64)
    magnitudeImage = np.zeros([imageXRes, imageYRes], dtype=np.float64)
    orientedImage = np.zeros([imageXRes, imageYRes], dtype=np.float64)


    # Get corners and orientation for this slice location
    corners      = pfile.Corners(slice)
    orientation  = pfile.Orientation(slice)
    imageCorners = orient.OrientCorners(corners, orientation)

    if len(corners) == 0:
        print('Slice Corner Dictionary is Empty')
    else:
        print(("Slice Corners: {}".format(corners)))

    if len(orientation) == 0:
        print('Slice Orientation Dictionary is Empty')
    else:
        print(("Slice Orientation: {}".format(orientation)))

    print(("ImageCorners: ", imageCorners))

    
    dat = pfile.KSpace(0,0) # seems to give me numpy array with all [dim0, nfids, ncoils]
    dim0, nfids, ncoil = dat.shape
    nfids_wat = 2
    nfids_met = 4

    data = np.empty([ncoil,nfids,dim0], dtype=np.complex128)

    print('got here 1')

    for i in range(nfids):
        for j in range(ncoil):
            data[j,i,:] = dat[:,i,j]

    sw        = 2500.0
    freq      = 63.86
    seqte     = 35.0
    nucstr    = '1H'
    seqstr    = 'PROBE-P'
                
    fpreset       = "/home/ese_user/orchestra/orchestra-python-sdk-1.0.x86_64/testdata/preset_ge_probep_te035_v2_braino.xml"
    fname_rgb_out = "/home/ese_user/orchestra/orchestra-python-sdk-1.0.x86_64/testdata/debug_last_run_ge_P12288.bin"
    fname_out     = "/home/ese_user/orchestra/orchestra-python-sdk-1.0.x86_64/testdata/result_vespa_inline_ge_P12288.xml"                
    
    print('got here 2')

    cbuf = ivg.vespa_process(dim0, ncoil,  
                             nfids_wat, 
                             nfids_met, 
                             sw, freq, 
                             seqte, 
                             nucstr, 
                             seqstr, 
                             data, 
                             fpreset=fpreset, 
                             fname_rgb_out=fname_rgb_out, 
                             fname_out=fname_out)    

    print('got here 3')

    #plt.figure(slice + 1)
    #plt.imshow(orientedImage, cmap='gray')
    #plt.suptitle('Oriented Image: Phase[{}], Slice[{}], Echo[{}]'.format((phase + 1), (slice + 1), (echo + 1), fontsize=2))

    slice = 0
    #dicom.Write("{}/image_vespa_result_{}.dcm".format(pathToSaveDicoms, slice), orientedImage, slice, corners, orientation)

    #plt.show()


