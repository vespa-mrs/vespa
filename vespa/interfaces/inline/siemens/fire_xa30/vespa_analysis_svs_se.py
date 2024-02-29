
import ismrmrd
import os
import itertools
import logging
import numpy as np
import numpy.fft as fft
import matplotlib.pyplot as plt
import xml.dom.minidom
import base64
import ctypes
import re
import mrdhelper

# vespa imports

import io
import platform

import matplotlib as mpl
mpl.use('Agg')

import vespa.analysis.block_prep_fidsum as block_prep_fidsum
import vespa.analysis.mrs_dataset as mrs_dataset
import vespa.analysis.util_import as util_import
import vespa.analysis.figure_layouts as figure_layouts
import vespa.common.util.export as util_export
import vespa.common.constants as constants

from vespa.common.mrs_data_raw import DataRawFidsum



# Folder for debug output files
debugFolder = "/tmp/share/debug"

# RAWDATA_SCALE is based on ICE_RAWDATA_SCALE in SpecRoFtFunctor.cpp
RAWDATA_SCALE = 131072.0 * 256.0


def process(connection, config, metadata):
    logging.info("Config: \n%s", config)

    # Continuously parse incoming data parsed from MRD messages
    acqGroup = []
    imgGroup = []
    try:
        for item in connection:
            # ----------------------------------------------------------
            # Raw k-space data messages
            # ----------------------------------------------------------
            if isinstance(item, ismrmrd.Acquisition):
                # Accumulate all imaging readouts in a group
                if (not item.is_flag_set(ismrmrd.ACQ_IS_NOISE_MEASUREMENT) and
                    not item.is_flag_set(ismrmrd.ACQ_IS_PARALLEL_CALIBRATION) and
                    not item.is_flag_set(ismrmrd.ACQ_IS_PHASECORR_DATA)):
                    acqGroup.append(item)

                # When this criteria is met, run process_raw() on the accumulated
                # data, which returns images that are sent back to the client.
                if item.is_flag_set(ismrmrd.ACQ_LAST_IN_SLICE):
                    logging.info("Processing a group of k-space data")
                    image = vespa_process(acqGroup, config, metadata)
                    connection.send_image(image)
                    acqGroup = []

            # ----------------------------------------------------------
            # Image data messages
            # ----------------------------------------------------------

            # Images and waveform data are not supported in this example

            if isinstance(item, ismrmrd.Image):
                continue

            elif isinstance(item, ismrmrd.Acquisition) or isinstance(item, ismrmrd.Waveform):
                continue

            elif item is None:
                break

            else:
                logging.error("Unsupported data type %s", type(item).__name__)

        if len(imgGroup) > 0:
            logging.info("Processing a group of images (untriggered)")
            image = process_image(imgGroup, config, metadata)
            connection.send_image(image)
            imgGroup = []

    finally:
        connection.send_close()


def vespa_process(acqGroup, config, metadata):

    # orig params - need to parse from header
    
    ncol = 0
    ncha = 0
    nave = 0
    os_remove = 0
    nleft = 0
    nright = 0
    os_factor = 0
    dwell = 0
    freq = 0
    delta = 0
    seqte = 0
    nucstr = 0
    seqstr = 0
    tapstr = 0
    presetfile = 0
    data64 = 0


    if self.verbose: logging.info("currently running: vespa_process() --")

    logging.info('In method vespa_process')
    logging.info('Preset file into program is = ' + str(presetfile))

    if platform.system() == 'Windows':
        #presetfile  = "C:/bsoher/code/xmlrpc_server_vespa/preset_svs_se_water_v1.xml"
        #presetfile  = "C:/bsoher/code/xmlrpc_server_vespa/preset_svs_se_braino_v1.xml"
        presetfile   = "C:/bsoher/code/xmlrpc_server_vespa/"+presetfile
        out_filename = "C:/bsoher/code/xmlrpc_server_vespa/ice_dataset_out_8.xml"

        logging.info('seqstr = '+str(seqstr))

        if seqstr == 'svs_seF':
            out_rgb_fname = "C:/bsoher/code/xmlrpc_server_vespa/debug_last_run_seq_svs_seF.bin"
        elif seqstr == 'svs_sead':
            out_rgb_fname = "C:/bsoher/code/xmlrpc_server_vespa/debug_last_run_seq_svs_sead.bin"
        else:
            out_rgb_fname = "C:/bsoher/code/xmlrpc_server_vespa/debug_last_run_seq_svs_xxxx.bin"

    else:
        presetfile   = "/opt/med/lib/soher/"+presetfile
        out_filename = "/opt/med/lib/soher/ice_dataset_out_8.xml"

        if seqstr == 'svs_seF':
            out_rgb_fname = "/opt/med/lib/soher/debug_last_run_seq_svs_seF.bin"
        elif seqstr == 'svs_sead':
            out_rgb_fname = "/opt/med/lib/soher/debug_last_run_seq_svs_sead.bin"
        else:
            out_rgb_fname = "/opt/med/lib/soher/debug_last_run_seq_svs_xxxx.bin"

    try:
        count  = 2 * ncol * ncha * nave
        fmt = '<%d%s' % (count, 'f')

        data1 = io.StringIO(data64.data)
        data2 = data1.read(struct.calcsize(fmt))
        data3 = struct.unpack(fmt, data2)

        data_iter = iter(data3)
        data4 = [complex(r, i) for r, i in zip(data_iter, data_iter)]
        logging.info('count, len(data4) = ', count, len(data4))

        data5 = np.array(data4)
        data5.shape = nave, ncha, ncol      # acquisition order

        ncol0 = 1<<((ncol-1).bit_length()-1) # ensure power of 2

    except Exception as e:
        logging.error(str(e))

    # Set up header info ----------------------------------------

    d = { }
    d["header"] = "Siemens FIRE Vespa - vespa_process()"
    d["sw"]             = (1.0 / float(dwell))
    d["remove_os"]      = os_remove
    d["readout_os"]     = float(os_factor)
    d["sequence_type"]  = seqstr
    d["frequency"]      = float(freq)
    d["dims"]           = [1,1,1,1024]
    d["dims"][0]        = int(ncol0/d["readout_os"]) 
    d["dims"][1]        = int(nave)
    d["dims"][2]        = int(ncha)
    d["dims"][3]        = 1 
    d["seqte"]          = float(seqte)     # comes in as msec
    logging.info('seqte = ', seqte)
    d["start_point"]    = nleft
    d["nucleus"]        = nucstr

    dims        = d["dims"]
    dim0        = dims[0]
    ncoils      = dims[2]
    nfids       = dims[1]
    acqdim0     = dims[0] * int(d["readout_os"])
    remove_os   = d["remove_os"]
    start_point = d['start_point']
    end_point   = start_point + acqdim0

    scale       = RAWDATA_SCALE / float(nfids) 
    dat         = np.empty([ncoils,nfids,dim0], dtype=np.complex128)

    for i in range(nfids):
        for j in range(ncoils):

            chan = data5[i,j,start_point:end_point].copy() * scale                     # in pulseq acquisition order

            if remove_os:
                chan = np.fft.fft(chan)
                chan = np.roll(chan, int(dim0/2))
                chan = np.fft.ifft(np.roll(chan[:dim0], int(-dim0/2)))   
            else:
                chan = chan[:dim0]              
            chan = np.conjugate(chan)

            dat[j,i,:] = chan       # index coils on outside so eventually collapse to 1,1,nfid,dim0

    if d["remove_os"]: d["sw"] = d["sw"] / 2.0

    d["data"] = dat
    d["data_source"] = "Siemens FIRE Vespa - vespa_process()"
    raw = DataRawFidsum(d)

    dataset = _import_siemens_ice([raw,], open_dataset=None)
    dataset = dataset[0]

    preset = _import_preset(presetfile)

    # Update dataset with preset ------------------------------------
    dataset.apply_preset(preset, voxel=(0,0,0))

    # Process and fit data ------------------------------------------
    chain_outputs = _process_all_blocks(dataset)

    fig_call = figure_layouts.analysis_brp512
    fig = fig_call( dataset, 
                    viffpath='Analysis - Siemens FIRE Inline', 
                    vespa_version='1.1.1',
                    timestamp='',
                    fontname='Courier New',
                    minplot=0.1,
                    maxplot=4.9,
                    nobase=True,
                    extfig=None,
                    fixphase=True,
                    verbose=False, 
                    debug=False)

    buf1 = fig.canvas.tostring_rgb()


    # convert string to byte arry and write to a debug file
    logging.info("Saving degug RGB file to - " + str(out_rgb_fname))
    cbuf = np.fromstring(buf1, dtype=np.uint8)
    cbuf.tofile(out_rgb_fname)


    logging.info('buf1 = fig.canvas.tostring_rgb()  - len() = ', len(buf1))
    buf2 = xmlrpc.client.Binary(buf1)

    dataset.dataset_filename = out_filename

    util_export.export(out_filename, [dataset], None, None, False)

    logging.info('finished vespa_process()')

    return buf2   


def _process_all_blocks(dataset):
    """ for all voxels, run chain in all blocks to update """
    
    chain_outputs = {}
    
    try:
        voxel = dataset.all_voxels
        for key in list(dataset.blocks.keys()):
            if key == 'spectral':
                key = 'spectral'
                block = dataset.blocks[key]
                tmp = block.chain.run(voxel, entry='all')
                chain_outputs[key] = tmp
                
                logging.info('block._svd_outputs = ', block.get_svd_output([0,0,0]))
                
                if 'fit' in list(dataset.blocks.keys()):
                    key = 'fit'
                    block = dataset.blocks[key]
                    block.chain.run(voxel, entry='initial_only')
                    key = 'spectral'
                    block = dataset.blocks[key]
                    block.set_do_fit(True, voxel[0])
                    tmp = block.chain.run(voxel, entry='all')
                    chain_outputs[key] = tmp
            else:
                block = dataset.blocks[key]
                tmp = block.chain.run(voxel, entry='all')
                chain_outputs[key] = tmp
    except Exception as e:
        logging.error(str(e)) 
        

    return chain_outputs
    

def _import_preset(presetfile):

    try:
        msg = ""
        try:
            importer = util_import.DatasetImporter(presetfile)
        except Exception as e:
            msg = str(e)
            logging.error(str(e)) 
        
#         except IOError:
#             msg = """I can't read the preset file "%s".""" % presetfile
#         except SyntaxError:
#             msg = """The preset file "%s" isn't valid Vespa Interchange File Format.""" % presetfile

        if msg:
            logging.error(msg)
        else:
            # Time to rock and roll!
            presets = importer.go()
            preset  = presets[0]
            
            return preset

    except Exception as e:
        logging.error(str(e)) 


def _import_siemens_ice(raws, open_dataset=None):
    """
    Stolen from Analysis main.py module - trimmed for CLI usage
    
    Assumption here is that we are opening one file in the reader. If there is
    an 'open_dataset' sent in, then the current reader is for an associated 
    file. We will associate the current file with the open one at the end of
    the code.
    
    """
    datasets = [ ]

    # Convert these raw objects into fully-fledged dataset objects.
    if open_dataset:
        zero_fill_multiplier = open_dataset.zero_fill_multiplier
    else:
        zero_fill_multiplier = 0

    # Step 2
    #
    # See if any data types need special classes. We usually only
    # look for raw fidsum classes which trigger a prep fidsum block.
    block_class_specs = []
    for raw in raws:
        d = { }
        if isinstance(raw, DataRawFidsum):
            d["prep"] = block_prep_fidsum.BlockPrepFidsum
        block_class_specs.append(d)

    f = lambda raw, block_classes: mrs_dataset.dataset_from_raw(raw, block_classes, zero_fill_multiplier)
    datasets = list(map(f, raws, block_class_specs))

    if datasets:

        if open_dataset is not None:
            open_dataset.blocks['raw'].set_associated_datasets([datasets[0], ])

        return datasets[0], open_dataset
    
    else:
        return None, open_dataset


# Orig rgb.py methods -----------------------------------------------------------------------

def process_raw(group, config, metadata):
    if len(group) == 0:
        return []

    # Create folder, if necessary
    if not os.path.exists(debugFolder):
        os.makedirs(debugFolder)
        logging.debug("Created folder " + debugFolder + " for debug output files")

    # Format data into single [cha PE RO phs] array
    lin = [acquisition.idx.kspace_encode_step_1 for acquisition in group]
    phs = [acquisition.idx.phase                for acquisition in group]

    # Use the zero-padded matrix size
    data = np.zeros((group[0].data.shape[0], 
                     metadata.encoding[0].encodedSpace.matrixSize.y, 
                     metadata.encoding[0].encodedSpace.matrixSize.x, 
                     max(phs)+1), 
                    group[0].data.dtype)

    rawHead = [None]*(max(phs)+1)

    for acq, lin, phs in zip(group, lin, phs):
        if (lin < data.shape[1]) and (phs < data.shape[3]):
            # TODO: Account for asymmetric echo in a better way
            data[:,lin,-acq.data.shape[1]:,phs] = acq.data

            # center line of k-space is encoded in user[5]
            if (rawHead[phs] is None) or (np.abs(acq.getHead().idx.kspace_encode_step_1 - acq.getHead().idx.user[5]) < np.abs(rawHead[phs].idx.kspace_encode_step_1 - rawHead[phs].idx.user[5])):
                rawHead[phs] = acq.getHead()

    # Flip matrix in RO/PE to be consistent with ICE
    data = np.flip(data, (1, 2))

    logging.debug("Raw data is size %s" % (data.shape,))
    np.save(debugFolder + "/" + "raw.npy", data)

    # Fourier Transform
    data = fft.fftshift( data, axes=(1, 2))
    data = fft.ifft2(    data, axes=(1, 2))
    data = fft.ifftshift(data, axes=(1, 2))
    data *= np.prod(data.shape) # FFT scaling for consistency with ICE

    # Sum of squares coil combination
    # Data will be [PE RO phs]
    data = np.abs(data)
    data = np.square(data)
    data = np.sum(data, axis=0)
    data = np.sqrt(data)

    logging.debug("Image data is size %s" % (data.shape,))
    np.save(debugFolder + "/" + "img.npy", data)

    # Determine max value (12 or 16 bit)
    BitsStored = 12
    if (mrdhelper.get_userParameterLong_value(metadata, "BitsStored") is not None):
        BitsStored = mrdhelper.get_userParameterLong_value(metadata, "BitsStored")
    maxVal = 2**BitsStored - 1

    # Normalize and convert to int16
    data *= maxVal/data.max()
    data = np.around(data)
    data = data.astype(np.int16)

    # Remove readout oversampling
    offset = int((data.shape[1] - metadata.encoding[0].reconSpace.matrixSize.x)/2)
    data = data[:,offset:offset+metadata.encoding[0].reconSpace.matrixSize.x]

    # Remove phase oversampling
    offset = int((data.shape[0] - metadata.encoding[0].reconSpace.matrixSize.y)/2)
    data = data[offset:offset+metadata.encoding[0].reconSpace.matrixSize.y,:]

    logging.debug("Image without oversampling is size %s" % (data.shape,))
    np.save(debugFolder + "/" + "imgCrop.npy", data)

    # Format as ISMRMRD image data
    imagesOut = []
    for phs in range(data.shape[2]):
        # Create new MRD instance for the processed image
        # data has shape [PE RO phs], i.e. [y x].
        # from_array() should be called with 'transpose=False' to avoid warnings, and when called
        # with this option, can take input as: [cha z y x], [z y x], or [y x]
        tmpImg = ismrmrd.Image.from_array(data[...,phs], transpose=False)

        # Set the header information
        tmpImg.setHead(mrdhelper.update_img_header_from_raw(tmpImg.getHead(), rawHead[phs]))
        tmpImg.field_of_view = (ctypes.c_float(metadata.encoding[0].reconSpace.fieldOfView_mm.x), 
                                ctypes.c_float(metadata.encoding[0].reconSpace.fieldOfView_mm.y), 
                                ctypes.c_float(metadata.encoding[0].reconSpace.fieldOfView_mm.z))
        tmpImg.image_index = phs

        # Set ISMRMRD Meta Attributes
        tmpMeta = ismrmrd.Meta()
        tmpMeta['DataRole']               = 'Image'
        tmpMeta['ImageProcessingHistory'] = ['FIRE', 'PYTHON']
        tmpMeta['WindowCenter']           = str((maxVal+1)/2)
        tmpMeta['WindowWidth']            = str((maxVal+1))
        tmpMeta['Keep_image_geometry']    = 1

        xml = tmpMeta.serialize()
        logging.debug("Image MetaAttributes: %s", xml)
        tmpImg.attribute_string = xml
        imagesOut.append(tmpImg)

    # Call process_image() to create RGB images
    imagesOut = process_image(imagesOut, config, metadata)

    return imagesOut

def process_image(images, config, metadata):
    if len(images) == 0:
        return []

    # Create folder, if necessary
    if not os.path.exists(debugFolder):
        os.makedirs(debugFolder)
        logging.debug("Created folder " + debugFolder + " for debug output files")

    logging.debug("Processing data with %d images of type %s", len(images), ismrmrd.get_dtype_from_data_type(images[0].data_type))

    # Note: The MRD Image class stores data as [cha z y x]

    # Extract image data into a 5D array of size [img cha z y x]
    data = np.stack([img.data                              for img in images])
    head = [img.getHead()                                  for img in images]
    meta = [ismrmrd.Meta.deserialize(img.attribute_string) for img in images]

    # Reformat data to [y x z cha img], i.e. [row col] for the first two dimensions
    data = data.transpose((3, 4, 2, 1, 0))

    # Display MetaAttributes for first image
    logging.debug("MetaAttributes[0]: %s", ismrmrd.Meta.serialize(meta[0]))

    # Optional serialization of ICE MiniHeader
    if 'IceMiniHead' in meta[0]:
        logging.debug("IceMiniHead[0]: %s", base64.b64decode(meta[0]['IceMiniHead']).decode('utf-8'))

    logging.debug("Original image data is size %s" % (data.shape,))
    np.save(debugFolder + "/" + "imgOrig.npy", data)

    if data.shape[3] != 1:
        logging.error("Multi-channel data is not supported")
        return []
    
    # Normalize to (0.0, 1.0) as expected by get_cmap()
    data = data.astype(float)
    data -= data.min()
    data *= 1/data.max()

    # Apply colormap
    cmap = plt.get_cmap('jet')
    rgb = cmap(data)

    # Remove alpha channel
    # Resulting shape is [row col z rgb img]
    rgb = rgb[...,0:-1]
    rgb = rgb.transpose((0, 1, 2, 5, 4, 3))
    rgb = np.squeeze(rgb, 5)

    # MRD RGB images must be uint16 in range (0, 255)
    rgb *= 255
    data = rgb.astype(np.uint16)
    np.save(debugFolder + "/" + "imgRGB.npy", data)

    currentSeries = 0

    # Re-slice back into 2D images
    imagesOut = [None] * data.shape[-1]
    for iImg in range(data.shape[-1]):
        # Create new MRD instance for the inverted image
        # Transpose from convenience shape of [y x z cha] to MRD Image shape of [cha z y x]
        # from_array() should be called with 'transpose=False' to avoid warnings, and when called
        # with this option, can take input as: [cha z y x], [z y x], or [y x]
        imagesOut[iImg] = ismrmrd.Image.from_array(data[...,iImg].transpose((3, 2, 0, 1)), transpose=False)
        data_type = imagesOut[iImg].data_type

        # Create a copy of the original fixed header and update the data_type
        # (we changed it to int16 from all other types)
        oldHeader = head[iImg]
        oldHeader.data_type = data_type

        # Set RGB parameters
        oldHeader.image_type = 6  # To be defined as ismrmrd.IMTYPE_RGB
        oldHeader.channels   = 3  # RGB "channels".  This is set by from_array, but need to be explicit as we're copying the old header instead

        # Increment series number when flag detected (i.e. follow ICE logic for splitting series)
        if mrdhelper.get_meta_value(meta[iImg], 'IceMiniHead') is not None:
            if mrdhelper.extract_minihead_bool_param(base64.b64decode(meta[iImg]['IceMiniHead']).decode('utf-8'), 'BIsSeriesEnd') is True:
                currentSeries += 1

        imagesOut[iImg].setHead(oldHeader)

        # Create a copy of the original ISMRMRD Meta attributes and update
        tmpMeta = meta[iImg]
        tmpMeta['DataRole']                       = 'Image'
        tmpMeta['ImageProcessingHistory']         = ['PYTHON', 'RGB']
        tmpMeta['SequenceDescriptionAdditional']  = 'FIRE_RGB'
        tmpMeta['Keep_image_geometry']            = 1

         # Add image orientation directions to MetaAttributes if not already present
        if tmpMeta.get('ImageRowDir') is None:
            tmpMeta['ImageRowDir'] = ["{:.18f}".format(oldHeader.read_dir[0]), "{:.18f}".format(oldHeader.read_dir[1]), "{:.18f}".format(oldHeader.read_dir[2])]

        if tmpMeta.get('ImageColumnDir') is None:
            tmpMeta['ImageColumnDir'] = ["{:.18f}".format(oldHeader.phase_dir[0]), "{:.18f}".format(oldHeader.phase_dir[1]), "{:.18f}".format(oldHeader.phase_dir[2])]

        metaXml = tmpMeta.serialize()
        logging.debug("Image MetaAttributes: %s", xml.dom.minidom.parseString(metaXml).toprettyxml())
        logging.debug("Image data has %d elements", imagesOut[iImg].data.size)

        imagesOut[iImg].attribute_string = metaXml

    return imagesOut
