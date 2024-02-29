'''
This code tests Vespa processing from a miniconda2 install on either Siemens VM 
WinXP IDEA environment OR on the MRIR Linux environment. 

It does not test XMLRPC funtionality. 

It loads data from a saved Numpy array in 'data5_saved.npy' and then converts 
it into a 'mrs_data_raw_fidsum' object. Preset file is loaded and attached
to the Dataset. Dataset is then processed through an early version of 
_process_all_blocks() and then figure_layout() results are saved to PNG and 
RGB buf file/array, and dataset saved to XML.

This is the Python 2 version of this code

'''

# Python imports


import os
import sys
import zlib
import base64
import xdrlib
import xmlrpc.client
import struct
import itertools
import io
import platform

import matplotlib as mpl
mpl.use('Agg')

import numpy as np

import vespa.analysis.src.block_prep_fidsum as block_prep_fidsum
import vespa.analysis.src.mrs_dataset as mrs_dataset
import vespa.analysis.src.util_import as util_import
import vespa.analysis.src.figure_layouts as figure_layouts
import vespa.common.mrs_data_raw_fidsum as mrs_data_raw_fidsum
import vespa.common.util.export as util_export
import vespa.common.constants as constants

# RAWDATA_SCALE is based on ICE_RAWDATA_SCALE in SpecRoFtFunctor.cpp
RAWDATA_SCALE = 131072.0 * 256.0




def vespa_process():

    print("currently running: vespa_process() --")

    try:

        # debug default value
        seqstr = 'svs_seF'

        if  platform.system() == 'Windows':
            fname        = "C:/bsoher/code/xmlrpc_server_vespa/data5_saved.npy"
            presetfile   = "C:/bsoher/code/xmlrpc_server_vespa/preset_svs_se_braino_v1.xml"
            out_filename = "C:/bsoher/code/xmlrpc_server_vespa/ice_dataset_out_7.xml"
            
            if seqstr == 'svs_seF':
                out_rgb_fname = "C:/bsoher/code/xmlrpc_server_vespa/debug_last_run_seq_svs_seF.bin"
            elif seqstr == 'svs_sead':
                out_rgb_fname = "C:/bsoher/code/xmlrpc_server_vespa/debug_last_run_seq_svs_sead.bin"
            else:
                out_rgb_fname = "C:/bsoher/code/xmlrpc_server_vespa/debug_last_run_seq_svs_xxxx.bin"

        else:    
            fname        = "/opt/med/lib/soher/data5_saved.npy"
            presetfile   = "/opt/med/lib/soher/preset_svs_se_braino_v1.xml"
            out_filename = "/opt/med/lib/soher/ice_dataset_out_7.xml"

            if seqstr == 'svs_seF':
                out_rgb_fname = "/opt/med/lib/soher/debug_last_run_seq_svs_seF.bin"
            elif seqstr == 'svs_sead':
                out_rgb_fname = "/opt/med/lib/soher/debug_last_run_seq_svs_sead.bin"
            else:
                out_rgb_fname = "/opt/med/lib/soher/debug_last_run_seq_svs_xxxx.bin"

        data5 = np.load(fname)

    except Exception as e:
        print(str(e))
        return

    ncol = 4096
    ncha = 4
    nave = 64
    os_remove = 1
    os_factor = 2.0
    dwell = 0.0002
    freq = 123.251094
    delta = 0.0
    seqte = 32.0
    nucstr = '1H'
    seqstr = 'svs_seF'
    tapstr = 1
    
    

    # Set up header info ----------------------------------------

    d = { }
    d["header"] = "XmlRpc Server Ice Vespa - vespa_process()"
    d["sw"]             = (1.0 / float(dwell))
    d["remove_os"]      = os_remove
    d["readout_os"]     = float(os_factor)
    d["sequence_type"]  = seqstr
    d["frequency"]      = float(freq)
    d["dims"]           = mrs_data_raw_fidsum.DataRawFidsum.DEFAULT_DIMS
    d["dims"][0]        = int(ncol/os_factor) 
    d["dims"][1]        = int(nave)
    d["dims"][2]        = int(ncha)
    d["dims"][3]        = 1 
    d["seqte"]          = float(seqte)     # comes in as msec
    print('seqte = ', seqte)
    d["start_point"]    = 0
    d["nucleus"]        = nucstr
    if nucstr == '1H':
        d["midppm"] = constants.DEFAULT_PROTON_CENTER_PPM
    else:
        d["midppm"] = constants.DEFAULT_XNUCLEI_CENTER_PPM


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

            chan = data5[i,j,:].copy() * scale                     # in pulseq acquisition order

            if remove_os:
                chan = np.fft.fft(chan)
                chan = np.roll(chan, int(dim0/2))
                chan = np.fft.ifft(np.roll(chan[:dim0], int(-dim0/2)))                 
            chan = np.conjugate(chan)

            dat[j,i,:] = chan       # index coils on outside so eventually collapse to 1,1,nfid,dim0

    if d["remove_os"]: d["sw"] = d["sw"] / 2.0

    d["data"] = dat
    d["data_source"] = "XmlRpc Server Ice Vespa - vespa_process()"
    raw = mrs_data_raw_fidsum.DataRawFidsum(d)

    dataset = _import_siemens_ice([raw,], open_dataset=None)
    dataset = dataset[0]

    preset = _import_preset(presetfile)

    # Update dataset with preset ------------------------------------
    dataset.apply_preset(preset, voxel=(0,0,0))

    # Process and fit data ------------------------------------------
    chain_outputs = _process_all_blocks(dataset)

    fig_call = figure_layouts.analysis_brp512
    fig = fig_call( dataset, 
                    viffpath='Analysis - Siemens ICE Inline', 
                    vespa_version='0.9.x',
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
    print('buf1 = fig.canvas.tostring_rgb()  - len() = ', len(buf1))

    # convert string to byte arry and write to a debug file
    print("Saving degug RGB file to - " + str(out_rgb_fname))
    cbuf = np.fromstring(buf1, dtype=np.uint8)
    cbuf.tofile(out_rgb_fname)

    buf2 = xmlrpc.client.Binary(buf1)

    
    
    dataset.dataset_filename = out_filename

    util_export.export(out_filename, [dataset], None, None, False)

    print('finished vespa_process()')
    
    
    
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
                
                print('block._svd_outputs = ', block.get_svd_output([0,0,0]))
                
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
        print(str(e)) 
        

    return chain_outputs
    

def _import_preset(presetfile):

    try:
        msg = ""
        try:
            importer = util_import.DatasetImporter(presetfile)
        except Exception as e:
            msg = str(e)
            print(str(e)) 
        
#         except IOError:
#             msg = """I can't read the preset file "%s".""" % presetfile
#         except SyntaxError:
#             msg = """The preset file "%s" isn't valid Vespa Interchange File Format.""" % presetfile

        if msg:
            print(msg)
        else:
            # Time to rock and roll!
            presets = importer.go()
            preset  = presets[0]
            
            return preset

    except Exception as e:
        print(str(e)) 
            
    
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
        if isinstance(raw, mrs_data_raw_fidsum.DataRawFidsum):
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
    
    
#------------------------------------------------------------------------------
# Test Code 
    
def _test():
    vespa_process()    
    
    
    
if __name__ == '__main__':
    
    _test()    