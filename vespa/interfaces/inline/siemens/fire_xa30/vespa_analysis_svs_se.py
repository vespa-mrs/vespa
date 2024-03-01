
import ismrmrd
import os
import itertools
import logging
import numpy as np
import numpy.fft as fft
import scipy as sp
import matplotlib.pyplot as plt
import xml.dom.minidom
import base64
import ctypes
import re
import mrdhelper

# vespa imports

import zlib
import xdrlib
# import io
# import platform

import matplotlib as mpl
mpl.use('wxAgg')

# import vespa.analysis.block_prep_fidsum as block_prep_fidsum
# import vespa.analysis.mrs_dataset as mrs_dataset
# import vespa.analysis.util_import as util_import
# import vespa.analysis.figure_layouts as figure_layouts
# import vespa.common.util.export as util_export
# import vespa.common.constants as constants
#
# from vespa.common.mrs_data_raw import DataRawFidsum



# Folder for debug output files
debugFolder = "/tmp/share/debug"

# RAWDATA_SCALE is based on ICE_RAWDATA_SCALE in SpecRoFtFunctor.cpp
RAWDATA_SCALE = 131072.0 * 256.0



test1 = ''


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
                    if image is not None:
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

    return None

#     # orig params - need to parse from header
#
#     ncol = 0
#     ncha = 0
#     nave = 0
#     os_remove = 0
#     nleft = 0
#     nright = 0
#     os_factor = 0
#     dwell = 0
#     freq = 0
#     delta = 0
#     seqte = 0
#     nucstr = 0
#     seqstr = 0
#     tapstr = 0
#     presetfile = 0
#     data64 = 0
#
#
#     if self.verbose: logging.info("currently running: vespa_process() --")
#
#     logging.info('In method vespa_process')
#     logging.info('Preset file into program is = ' + str(presetfile))
#
#     if platform.system() == 'Windows':
#         #presetfile  = "C:/bsoher/code/xmlrpc_server_vespa/preset_svs_se_water_v1.xml"
#         #presetfile  = "C:/bsoher/code/xmlrpc_server_vespa/preset_svs_se_braino_v1.xml"
#         presetfile   = "C:/bsoher/code/xmlrpc_server_vespa/"+presetfile
#         out_filename = "C:/bsoher/code/xmlrpc_server_vespa/ice_dataset_out_8.xml"
#
#         logging.info('seqstr = '+str(seqstr))
#
#         if seqstr == 'svs_seF':
#             out_rgb_fname = "C:/bsoher/code/xmlrpc_server_vespa/debug_last_run_seq_svs_seF.bin"
#         elif seqstr == 'svs_sead':
#             out_rgb_fname = "C:/bsoher/code/xmlrpc_server_vespa/debug_last_run_seq_svs_sead.bin"
#         else:
#             out_rgb_fname = "C:/bsoher/code/xmlrpc_server_vespa/debug_last_run_seq_svs_xxxx.bin"
#
#     else:
#         presetfile   = "/opt/med/lib/soher/"+presetfile
#         out_filename = "/opt/med/lib/soher/ice_dataset_out_8.xml"
#
#         if seqstr == 'svs_seF':
#             out_rgb_fname = "/opt/med/lib/soher/debug_last_run_seq_svs_seF.bin"
#         elif seqstr == 'svs_sead':
#             out_rgb_fname = "/opt/med/lib/soher/debug_last_run_seq_svs_sead.bin"
#         else:
#             out_rgb_fname = "/opt/med/lib/soher/debug_last_run_seq_svs_xxxx.bin"
#
#     try:
#         count  = 2 * ncol * ncha * nave
#         fmt = '<%d%s' % (count, 'f')
#
#         data1 = io.StringIO(data64.data)
#         data2 = data1.read(struct.calcsize(fmt))
#         data3 = struct.unpack(fmt, data2)
#
#         data_iter = iter(data3)
#         data4 = [complex(r, i) for r, i in zip(data_iter, data_iter)]
#         logging.info('count, len(data4) = ', count, len(data4))
#
#         data5 = np.array(data4)
#         data5.shape = nave, ncha, ncol      # acquisition order
#
#         ncol0 = 1<<((ncol-1).bit_length()-1) # ensure power of 2
#
#     except Exception as e:
#         logging.error(str(e))
#
#     # Set up header info ----------------------------------------
#
#     d = { }
#     d["header"] = "Siemens FIRE Vespa - vespa_process()"
#     d["sw"]             = (1.0 / float(dwell))
#     d["remove_os"]      = os_remove
#     d["readout_os"]     = float(os_factor)
#     d["sequence_type"]  = seqstr
#     d["frequency"]      = float(freq)
#     d["dims"]           = [1,1,1,1024]
#     d["dims"][0]        = int(ncol0/d["readout_os"])
#     d["dims"][1]        = int(nave)
#     d["dims"][2]        = int(ncha)
#     d["dims"][3]        = 1
#     d["seqte"]          = float(seqte)     # comes in as msec
#     logging.info('seqte = ', seqte)
#     d["start_point"]    = nleft
#     d["nucleus"]        = nucstr
#
#     dims        = d["dims"]
#     dim0        = dims[0]
#     ncoils      = dims[2]
#     nfids       = dims[1]
#     acqdim0     = dims[0] * int(d["readout_os"])
#     remove_os   = d["remove_os"]
#     start_point = d['start_point']
#     end_point   = start_point + acqdim0
#
#     scale       = RAWDATA_SCALE / float(nfids)
#     dat         = np.empty([ncoils,nfids,dim0], dtype=np.complex128)
#
#     for i in range(nfids):
#         for j in range(ncoils):
#
#             chan = data5[i,j,start_point:end_point].copy() * scale                     # in pulseq acquisition order
#
#             if remove_os:
#                 chan = np.fft.fft(chan)
#                 chan = np.roll(chan, int(dim0/2))
#                 chan = np.fft.ifft(np.roll(chan[:dim0], int(-dim0/2)))
#             else:
#                 chan = chan[:dim0]
#             chan = np.conjugate(chan)
#
#             dat[j,i,:] = chan       # index coils on outside so eventually collapse to 1,1,nfid,dim0
#
#     if d["remove_os"]: d["sw"] = d["sw"] / 2.0
#
#     d["data"] = dat
#     d["data_source"] = "Siemens FIRE Vespa - vespa_process()"
#     raw = DataRawFidsum(d)
#
#     dataset = _import_siemens_ice([raw,], open_dataset=None)
#     dataset = dataset[0]
#
#     preset = _import_preset(presetfile)
#
#     # Update dataset with preset ------------------------------------
#     dataset.apply_preset(preset, voxel=(0,0,0))
#
#     # Process and fit data ------------------------------------------
#     chain_outputs = _process_all_blocks(dataset)
#
#     fig_call = figure_layouts.analysis_brp512
#     fig = fig_call( dataset,
#                     viffpath='Analysis - Siemens FIRE Inline',
#                     vespa_version='1.1.1',
#                     timestamp='',
#                     fontname='Courier New',
#                     minplot=0.1,
#                     maxplot=4.9,
#                     nobase=True,
#                     extfig=None,
#                     fixphase=True,
#                     verbose=False,
#                     debug=False)
#
#     buf1 = fig.canvas.tostring_rgb()
#
#
#     # convert string to byte arry and write to a debug file
#     logging.info("Saving degug RGB file to - " + str(out_rgb_fname))
#     cbuf = np.fromstring(buf1, dtype=np.uint8)
#     cbuf.tofile(out_rgb_fname)
#
#
#     logging.info('buf1 = fig.canvas.tostring_rgb()  - len() = ', len(buf1))
#     buf2 = xmlrpc.client.Binary(buf1)
#
#     dataset.dataset_filename = out_filename
#
#     util_export.export(out_filename, [dataset], None, None, False)
#
#     logging.info('finished vespa_process()')
#
#     return buf2
#
#
# def _process_all_blocks(dataset):
#     """ for all voxels, run chain in all blocks to update """
#
#     chain_outputs = {}
#
#     try:
#         voxel = dataset.all_voxels
#         for key in list(dataset.blocks.keys()):
#             if key == 'spectral':
#                 key = 'spectral'
#                 block = dataset.blocks[key]
#                 tmp = block.chain.run(voxel, entry='all')
#                 chain_outputs[key] = tmp
#
#                 logging.info('block._svd_outputs = ', block.get_svd_output([0,0,0]))
#
#                 if 'fit' in list(dataset.blocks.keys()):
#                     key = 'fit'
#                     block = dataset.blocks[key]
#                     block.chain.run(voxel, entry='initial_only')
#                     key = 'spectral'
#                     block = dataset.blocks[key]
#                     block.set_do_fit(True, voxel[0])
#                     tmp = block.chain.run(voxel, entry='all')
#                     chain_outputs[key] = tmp
#             else:
#                 block = dataset.blocks[key]
#                 tmp = block.chain.run(voxel, entry='all')
#                 chain_outputs[key] = tmp
#     except Exception as e:
#         logging.error(str(e))
#
#
#     return chain_outputs
#
#
# def _import_preset(presetfile):
#
#     try:
#         msg = ""
#         try:
#             importer = util_import.DatasetImporter(presetfile)
#         except Exception as e:
#             msg = str(e)
#             logging.error(str(e))
#
# #         except IOError:
# #             msg = """I can't read the preset file "%s".""" % presetfile
# #         except SyntaxError:
# #             msg = """The preset file "%s" isn't valid Vespa Interchange File Format.""" % presetfile
#
#         if msg:
#             logging.error(msg)
#         else:
#             # Time to rock and roll!
#             presets = importer.go()
#             preset  = presets[0]
#
#             return preset
#
#     except Exception as e:
#         logging.error(str(e))
#
#
# def _import_siemens_ice(raws, open_dataset=None):
#     """
#     Stolen from Analysis main.py module - trimmed for CLI usage
#
#     Assumption here is that we are opening one file in the reader. If there is
#     an 'open_dataset' sent in, then the current reader is for an associated
#     file. We will associate the current file with the open one at the end of
#     the code.
#
#     """
#     datasets = [ ]
#
#     # Convert these raw objects into fully-fledged dataset objects.
#     if open_dataset:
#         zero_fill_multiplier = open_dataset.zero_fill_multiplier
#     else:
#         zero_fill_multiplier = 0
#
#     # Step 2
#     #
#     # See if any data types need special classes. We usually only
#     # look for raw fidsum classes which trigger a prep fidsum block.
#     block_class_specs = []
#     for raw in raws:
#         d = { }
#         if isinstance(raw, DataRawFidsum):
#             d["prep"] = block_prep_fidsum.BlockPrepFidsum
#         block_class_specs.append(d)
#
#     f = lambda raw, block_classes: mrs_dataset.dataset_from_raw(raw, block_classes, zero_fill_multiplier)
#     datasets = list(map(f, raws, block_class_specs))
#
#     if datasets:
#
#         if open_dataset is not None:
#             open_dataset.blocks['raw'].set_associated_datasets([datasets[0], ])
#
#         return datasets[0], open_dataset
#
#     else:
#         return None, open_dataset


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



test2 = b'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'



def get_fake_result(n256=False):
    """
    The fake data is an encoded string previously saved of an RGB file
    screenshot of a Vespa-Analysis results layout. We wrap it in and
    xmlrpclib friendly format and ship it back.

    """
    #import matplotlib.pyplot as plt

    dat = base64.b64decode(test2)
    dat = zlib.decompress(dat)
    data = np.fromstring(dat, dtype=np.uint8)
    data.shape = 512,512,3


    # fig = plt.figure(figsize=[16,16])
    # ax1 = fig.add_subplot(2, 2, 1)
    # ax1.imshow(data)
    # ax2 = fig.add_subplot(2, 2, 2)
    # ax2.imshow(dat256)
    # plt.show()

    if n256:
        data = sp.ndimage.zoom(data, [0.5, 0.5, 1.0])

    return data

if __name__ == '__main__':

    r = get_fake_result()

    bob = 11

