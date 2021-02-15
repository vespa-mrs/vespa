# Python modules


import os
import sys
import shutil

# 3rd party modules
import numpy as np
try:
    import dicom
except:
    import pydicom as dicom

import GERecon
#import matplotlib.pyplot as plt        # this line causes SegFault on Linux

# Our modules
import vespa.analysis.util_import as util_import
import vespa.analysis.figure_layouts as figure_layouts

import vespa.common.util.export as util_export
import vespa.common.util.time_ as util_time


"""
Version History
--------------------------

v1 - testing in linux so we can have access to GERecon library v1.0. 
     Calling as a direct script, no __main__

v2 - testing in linux so we can have access to GERecon library v1.1. 
     Calling as a direct script, no __main__

v3 - always in linux, still GERecon 1.1, changed to using convention
     if __name__ == '__main__': to call this script

"""




def do_process(filename):
    """
    Data1: Phase[1] -> Slices[4] -> Echo[1] -> Channel[1]
    
    """
    msg = ''

    debug          = False
    verbose        = True

    out_set = { 'savetype'   : 'pdf',
                'minplot'    : 0.1,
                'maxplot'    : 4.9,
                'fixphase'   : True,
                'fontname'   : 'Arial',     # mayber Courier New
                'dpi'        : 128,
                'pad_inches' : 0.5
             }
    
    
    fname, fpath = os.path.split(filename) 
    
    # Create a pfile object and read header and sequence identifiers
    pfile   = GERecon.Pfile(filename)
    hdr     = pfile.Header()
    seqstr  = hdr['rdb_hdr_image']['psdname']
    seqte   = float(hdr['rdb_hdr_image']['te'])/1000.0

    fbase = '/home/ese_user/orchestra/code/python_sdk_v11/vespa_inline/'
    
    if seqstr == 'oslaser':
        # this is sLASER
        
        nfids_wat = int(hdr['rdb_hdr_rec']['user19'])   # number of reference frames
        nfids_met = int(hdr['rdb_hdr_rec']['user4'])    # number of averages (wat suppressed)
        freq      = float(hdr['rdb_hdr_rec']['rdb_hdr_ps_mps_freq'])/1e7
        sw        = float(hdr['rdb_hdr_rec']['user0'])
        nucstr    = '1H'    

        fpresets = []
        fpresets.append(fbase+"testdata/preset_ge_oslaser_metab_invivo.xml")     # metab
        fpresets.append(fbase+"testdata/preset_ge_oslaser_water.xml")            # water quant
        fpresets.append(None)    # coil combine
        fpresets.append(None)    # ecc correct
        
        fbasis_mmol = None       # mmol   fbase+'\\basis_mmol_dataset_seadMM2014_truncat2048pts_normScale100dc015.xml'
    
        fname_rgb_out  = fbase+"testdata/debug_last_run_ge_{}.bin".format(fname)
        fname_out      = fbase+"testdata/result_vespa_inline_ge_P12288.xml"                
        fname_dicom    = fbase+"DICOM_OUT/my_output_dicom_last_run.ima"
        
    elif seqstr == 'PROBE-P':
        # this is normal PROBE-P
        
        fbase = '/home/ese_user/orchestra/code/python_sdk_v11/vespa_inline/'
        
        nex       = hdr['rdb_hdr_image']['nex']
        user4     = hdr['rdb_hdr_rec']['user4']
        nfids_wat = 16 / int(nex)
        nfids_met = user4 / int(nex)
        freq      = float(hdr['rdb_hdr_rec']['rdb_hdr_ps_mps_freq'])/1e7
        sw        = float(hdr['rdb_hdr_rec']['user0'])
        nucstr    = '1H'    

        fpresets = []
        fpresets.append(fbase+"testdata/preset_ge_probep_te035_v2_braino.xml")     # metab
        fpresets.append(None)    # water quant
        fpresets.append(None)    # coil combine
        fpresets.append(None)    # ecc correct
        
        fbasis_mmol = None       # mmol   fbase+'\\basis_mmol_dataset_seadMM2014_truncat2048pts_normScale100dc015.xml'

        fname_rgb_out  = fbase+"testdata/debug_last_run_ge_{}.bin".format(fname)
        fname_out      = fbase+"testdata/provenance_last_run_vespa_inline_ge_{}.xml".format(fname)                
        fname_dicom    = fbase+"DICOM_OUT/my_output_dicom_last_run.ima"
        
    else:
        msg = 'Error: Unknown pulse sequence, returning!  seqstr = {}'.format(seqstr)

    if msg:        
        print(msg, file=sys.stderr)
        print(msg, file=sys.stdout)
        sys.exit(-1)
    
    for item in fdatasets+fpresets+[fbasis_mmol,]:
        if item is not None:
            if not os.path.isfile(item):
                msg += """\nFILE does not exist "%s".""" % item
    
    if msg:
        print(msg, file=sys.stderr)
        print(msg, file=sys.stdout)
        sys.exit(-1)
        
    
    
    # Get the metadata from the loaded Pfile
    metadata = pfile.MetaData()
    
    acqXRes     = metadata['acquiredXRes']  # spectral points
    acqYRes     = metadata['acquiredYRes']  # number of FIDs
    imageXRes   = metadata['imageXRes']     # 512
    imageYRes   = metadata['imageYRes']     # 512
    numPhases   = metadata['phases']        # 1 here
    numSlices   = metadata['slices']        # 1 here
    numEchoes   = metadata['echoes']        # 1 here
    numChannels = metadata['channels']      # 8 here
    
    print(('AcqXRes is {}, AcqYRes is {}, ImageXRes is {}, ImageYRes is {}, numChannels is {}, numEchoes is {}, numSlices is {}, numPhases is {}'.format(acqXRes, acqYRes, imageXRes, imageYRes, numChannels, numEchoes, numSlices, numPhases)))

    # Pre-Initialize the kSpace and image arrays
    orientedImage = np.zeros([imageXRes, imageYRes], dtype=np.float64)
    
    # Create Transform, Orientation and DICOM objects
    #transformer = GERecon.Transformer()
    orient      = GERecon.Orientation()
    dicom       = GERecon.Dicom()
    
    # Clear the Output directory before saving the DICOM images
    if os.path.isdir(dicom_out_path):
        shutil.rmtree(dicom_out_path, ignore_errors=True)
    
    # Get corners and orientation for this slice location
    slice = 0
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
    
    data = np.empty([ncoil,nfids,dim0], dtype=np.complex128)
    
    print('got here 1')
    
    for i in range(nfids):
        for j in range(ncoil):
            data[j,i,:] = dat[:,i,j]
    
    # Set up header info ----------------------------------------

    d = { }
    d["header"] = "Vespa Inline GE - vespa_process_ge()"
    d["sw"]             = float(sw)
    d["sequence_type"]  = seqstr
    d["frequency"]      = float(freq)
    d["dims"]           = [1,1,1,1024]
    d["dims"][0]        = int(dim0) 
    d["dims"][1]        = int(nfids_wat+nfids_met)
    d["dims"][2]        = int(ncoil)
    d["dims"][3]        = 1 
    d["seqte"]          = float(seqte)     # comes in as msec
    d["nucleus"]        = nucstr
    if nucstr == '1H':
        d["resppm"] = constants.DEFAULT_PROTON_CENTER_PPM
    else:
        print("X-nuclei detected in data parsing, returning.")
        return
 

    #-------------------------------
    # Water Dataset

    dim0        = d["dims"][0]
    ncoils      = d["dims"][2]
    nfids       = nfids_wat
    acqdim0     = d["dims"][0] 

    scale       = 1.0  #RAWDATA_SCALE / float(nfids) 
    dat         = np.empty([ncoils,nfids,dim0], dtype=np.complex128)

    for i in range(nfids):
        for j in range(ncoils):
            tmp = data[j,i,:].copy() * scale                     # in pulseq acquisition order
            # tmp = np.conjugate(chan)
            dat[j,i,:] = tmp                # index coils on outside so eventually collapse to 1,1,nfid,dim0

    d["data"] = dat
    d["data_source"] = "Vespa Inline GE - run_inline_vespa_ge() - water fids"
    d["dims"][1] = nfids_wat
    raw = DataRawFidsum(d)

    dataset = _import_ge([raw,], open_dataset=None)
    dataset_water = dataset[0]

    #-------------------------------
    # Metabs Dataset

    dim0        = d["dims"][0]
    ncoils      = d["dims"][2]
    nfids       = nfids_met
    acqdim0     = d["dims"][0] 

    scale       = 1.0  #RAWDATA_SCALE / float(nfids) 
    dat         = np.empty([ncoils,nfids,dim0], dtype=np.complex128)

    for i in range(nfids):
        for j in range(ncoils):
            tmp = data[j,i+nfids_wat,:].copy() * scale                     # in pulseq acquisition order
            # tmp = np.conjugate(chan)
            dat[j,i,:] = tmp       # index coils on outside so eventually collapse to 1,1,nfid,dim0

    d["data"] = dat
    d["data_source"] = "Vespa Inline GE - run_inline_vespa_ge() - metab fids"
    d["dims"][1] = nfids_met
    raw = DataRawFidsum(d)

    dataset = _import_ge([raw,], open_dataset=None)
    dataset_metab = dataset[0]


    data_coil = None
    data_ecc  = None        # for now

    
    #--------------------------------------------------------------------------
    # Start Processing - expand parameters
    
    presets = []
    for fname in fpresets:
        if fname is not None:
            presets.append(load_preset(fname, verbose=True, debug=debug))
        else:
            presets.append(None)
    
    preset_metab, preset_water, preset_coil, preset_ecc = presets
    
    if fbasis_mmol is not None:
        basis_mmol, msg = util_file_import.open_viff_dataset_file([fbasis_mmol,]) 
    else:
        basis_mmol = None

    
    # Load and Process - Coil Combine Dataset -----------------------

    if verbose: print("Apply Preset and Run Chain - Coil Combine")    
    if data_coil is not None and preset_coil is not None:
        try:
            msg = """applying preset - coil combine""" 
            data_coil.apply_preset(preset_coil, voxel=(0,0,0))  # update dataset object with preset blocks and chains
    
            msg = """running chain - coil combine""" 
            _process_all_blocks(data_coil)
        
        except:
            print(msg+'\n'+str(sys.exc_info()[1]), file=sys.stderr)
            sys.exit(-1)

    # Load Preset - Ecc, Water and Metab Datasets -------------------

    if verbose: print("Apply Preset - Ecc, Water and Metab Datasets")    
    try:
        #----------------------------------------------------------------------
        # Apply presets to ecc, water and metab datasets
        
        if data_ecc is not None and preset_ecc is not None:
            msg = """applying preset - ecc""" 
            data_ecc.apply_preset(preset_ecc, voxel=(0,0,0))      # chain  

        if data_water is not None and preset_water is not None:
            msg = """applying preset - water""" 
            data_water.apply_preset(preset_water, voxel=(0,0,0))  

        if data_metab is not None and preset_metab is not None:
            msg = """applying preset - metab""" 
            data_metab.apply_preset(preset_metab, voxel=(0,0,0))  

        #----------------------------------------------------------------------
        # Attach coil combine to ecc, water and metab datasets - run chain ecc

        msg = """attaching coil combine to - ecc, water and metab"""
        for dset in [data_ecc, data_water, data_metab]:
            if dset is not None and data_coil is not None:
                dset.set_associated_dataset_combine(data_coil)
        
        if data_ecc is not None:
            msg = """running chain - ecc"""
            _process_all_blocks(data_ecc)       # get combined FID for next steps
        
        #----------------------------------------------------------------------
        # Attach ecc to water and metab datasets - run chain water

        msg = """attaching ecc to - water and metab"""
        for dset in [data_water, data_metab]:
            if dset is not None and data_ecc is not None:
                dset.set_associated_dataset_ecc(data_ecc)
        
        if data_water is not None:
            msg = """running chain - water"""
            _process_all_blocks(data_water)

        #----------------------------------------------------------------------
        # Attach mmol_basis and water to metab dataset - run chain metab

        msg = """attaching mmol_basis and water to - metab"""
        for dset in [data_metab,]:
            if dset is not None and data_water is not None:
                if basis_mmol is not None:
                    dset.set_associated_dataset_mmol(basis_mmol)
                dset.set_associated_dataset_quant(data_water)

        if data_metab is not None:
            msg = """running chain - metab"""
            _process_all_blocks(data_metab)


    except:
        print('Error: '+msg+'\n'+sys.exc_info()[1].message, file=sys.stderr)
        sys.exit(-1)
    
    #--------------------------------------------------------------------------
    # Begin Output
    
    timestamp = util_time.now(util_time.DISPLAY_TIMESTAMP_FORMAT)
    
    
    # Create unique name ID for this dataset ------------------------
    
    data_metab.dataset_filename = fname_out

    # Save results to XML -----------------------------------------------------
    
    if verbose: print("""Saving dataset to XML file "%s". """ % fname_out)
    
    try:
        #bob = 10
        util_export.export(fname_out, [data_metab,], None, None, False)
    except Exception as e:
        msg = """I can't write the file "%s".""" % fname_out
        print(msg, file=sys.stderr)
        print(repr(e), file=sys.stderr)
        sys.exit(-1)


    # Save results to PDF -----------------------------------------------------

    fig_call = figure_layouts.analysis_brp_generic
    fig = fig_call( data_metab, 
                    viffpath='Analysis - Philips PRIDE Inline', 
                    vespa_version='0.9.10i',
                    timestamp=timestamp,
                    fontname=out_set['fontname'],
                    minplot=out_set['minplot'],
                    maxplot=out_set['maxplot'],
                    nobase=False,
                    extfig=None,
                    fixphase=out_set['fixphase'],
                    verbose=False, debug=False)

    buf1 = fig.canvas.tostring_rgb()

    # convert string to byte arry and write to a debug file
    cbuf = np.fromstring(buf1, dtype=np.uint8)


    if fname_rgb_out:
        print("Saving degug RGB file to - " + str(fname_rgb_out))
        cbuf.tofile(fname_rgb_out)
        fig.savefig(fname_rgb_out+'.png', dpi=out_set['dpi'], pad_inches=0.5)



    print('finished vespa_process()')

    #dicom.Write("{}/image_vespa_result_{}.dcm".format(dicom_out_path, slice), orientedImage, slice, corners, orientation)
    dicom.Write(fname_dicom, cbuf, slice, corners, orientation)

    
    print('got here 3')
    
    return 'hello GE'



def _process_all_blocks(dataset):
    """ for all voxels, run chain in all blocks to update """
    
    chain_outputs = {}
    
    voxel = dataset.all_voxels
    for key in list(dataset.blocks.keys()):
        if key == 'spectral':
            key = 'spectral'
            block = dataset.blocks[key]
            tmp = block.chain.run(voxel, entry='all')
            chain_outputs[key] = tmp
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

    return chain_outputs



def load_preset(presetfile, verbose=False, debug=False):

    # Load PRESET object ----------------------------------------------
    
    if verbose: print("""load_preset - Presetfile = "%s"."""  % presetfile) 
    if debug: return 
     
    try:
        msg = ""
        try:
            importer = util_import.DatasetImporter(presetfile)
        except IOError:
            msg = """load_preset - I can't read the preset file "%s".""" % presetfile
        except SyntaxError:
            msg = """load_preset - The preset file "%s" isn't valid Vespa Interchange File Format.""" % presetfile

        if msg:
            print(msg, file=sys.stderr)
            sys.exit(-1)
        else:
            # Time to rock and roll!
            presets = importer.go()
            preset  = presets[0]
    except:
        msg = """load_preset - Unknown exception reading Preset file "%s".""" % presetfile 
        print(msg, file=sys.stderr)
        sys.exit(-1)
    
    return preset









if __name__ == '__main__':
    
    
    if len(sys.argv) <= 1:
        # only the script name itself - no filename
        pfilename =  "/home/ese_user/orchestra/python_sdk_v10/testdata/P12288.7"
    else:
        pfilename = sys.argv[1]
        if not os.path.isfile(fname):
            msg = "Error - pfile name passed in is not a file, returning."
            print(msg, file=sys.stderr)
            print(msg, file=sys.stdout)
            sys.exit(-1)

    res = do_process(pfilename)
    
    