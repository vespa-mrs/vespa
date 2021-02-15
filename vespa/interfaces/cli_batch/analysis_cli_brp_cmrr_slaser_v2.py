# Python modules
import os
import sys
import multiprocessing

# 3rd party modules
import numpy as np

# Our modules
import vespa.analysis.util_import as util_import
import vespa.analysis.util_file_import as util_file_import
import vespa.common.util.export as util_export

# Change to True to enable the assert() statements sprinkled through the code
ASSERTIONS_ENABLED = False

# data is complex64 per Siemens documentation for Twix
NUMPY_DATA_TYPE = np.complex64
# BYTES_PER_ELEMENT expresses how many bytes each element occupies. You
# shouldn't need to change this definition.
BYTES_PER_ELEMENT = np.zeros(1, dtype=NUMPY_DATA_TYPE).nbytes

# RAWDATA_SCALE is based on ICE_RAWDATA_SCALE in SpecRoFtFunctor.cpp
RAWDATA_SCALE = 131072.0 * 256.0


DESC =  \
"""Command line interface to process MRS data in Vespa-Analysis. 
 Data filename, preset file name, data type string and CSV output 
 file name values are all required for this command to function 
 properly.
  
 Note. You may have to enclose data/preset/output strings in double 
 quotation marks for them to process properly if they have  
 spaces or other special characters embedded in them.
"""

class CliError(Exception):
    """Basic exception for errors when applying preset object"""
    def __init__(self, msg=None):
        if msg is None:
            # set default error message
            msg = 'A general cli error occured.'
        e = sys.exc_info()
        msg =  'CliErrorMessage : '+msg
        msg += '\n'
        msg += 'BaseErrorMessage: '+e[1].message
        super(CliError, self).__init__(msg)


def clean_header(header):
    """ converts all values in ICE dict into a long string"""
    return "need to write"


def analysis_cli_gulin(datasets, preset_metab,
                                 preset_coil,
                                 preset_water,
                                 preset_ecc,
                                 out_base,
                                 out_prefix,
                                 basis_mmol=None,
                                 verbose=False, debug=False, in_gui=False):
    
    
    # Sort datasets into variables ----------------------------------
    
    data_coil, data_ecc, data_water, data_ecc2, data_water2, data_metab, basis_mmol = datasets

    msg = ""
    
    # Load and Process - Coil Combine Dataset -----------------------

    if verbose: print("Apply Preset and Run Chain - Coil Combine")    
    try:
        msg = """applying preset - coil combine""" 
        data_coil.apply_preset(preset_coil, voxel=(0,0,0))  # update dataset object with preset blocks and chains

        msg = """running chain - coil combine""" 
        _process_all_blocks(data_coil)
    
    except:
        if not in_gui:
            print(msg+'\n'+str(sys.exc_info()[1]), file=sys.stderr)
            sys.exit(-1)
        else:
            raise CliError(msg)

    # Load Preset - Ecc, Water and Metab Datasets -------------------

    if verbose: print("Apply Preset - Ecc, Water and Metab Datasets")    
    try:
        # Apply presets to ecc, water and metab datasets
        
        msg = """applying preset - ecc""" 
        data_ecc.apply_preset(preset_ecc, voxel=(0,0,0))      # chain  

        msg = """applying preset - water""" 
        data_water.apply_preset(preset_water, voxel=(0,0,0))  

        msg = """applying preset - metab""" 
        data_metab.apply_preset(preset_metab, voxel=(0,0,0))  

        #----------------------------------------------------------------------
        # Attach coil combine to ecc, water and metab datasets - run chain ecc

        msg = """attaching coil combine to - ecc, water and metab"""
        for dset in [data_ecc, data_water, data_metab]:
            #dset.blocks['prep'].attach_dataset_coil_combine(data_coil)
            dset.set_associated_dataset_combine(data_coil)
        
        msg = """running chain - ecc"""
        _process_all_blocks(data_ecc)       # get combined FID for next steps
        
        #----------------------------------------------------------------------
        # Attach ecc to water and metab datasets - run chain water

        msg = """attaching ecc to - water and metab"""
        for dset in [data_water, data_metab]:
            #dset.blocks['spectral'].attach_dataset_ecc(data_ecc)
            dset.set_associated_dataset_ecc(data_ecc)
        
        msg = """running chain - water"""
        _process_all_blocks(data_water)

        #----------------------------------------------------------------------
        # Attach mmol_basis and water to metab dataset - run chain metab

        msg = """attaching mmol_basis and water to - metab"""
        for dset in [data_metab,]:
            if basis_mmol is not None:
                #dset.blocks['fit'].attach_dataset_basis_mmol(basis_mmol)
                dset.set_associated_dataset_mmol(basis_mmol)
            #dset.blocks['quant'].attach_dataset_water_quant(data_water)
            dset.set_associated_dataset_quant(data_water)

        msg = """running chain - metab"""
        _process_all_blocks(data_metab)


    except:
        if not in_gui:
            print('Error: '+msg+'\n'+sys.exc_info()[1].message, file=sys.stderr)
            sys.exit(-1)
        else:
            raise CliError(msg)
    
    
    # Create unique name ID for this dataset ------------------------
    
    outxml = out_base+out_prefix+'.xml'
    
    data_metab.dataset_filename = outxml

    # Save results to CSV file --------------------------------------

#     if verbose: print """Saving results to CSV file "%s". """ % csvfile
#     
#     fit = data_metab.blocks["fit"]
#     data_source = data_metab.blocks["raw"].get_data_source(voxel)
#     
#     val, hdr = fit.results_as_csv(voxel[0], fit.chain.fitted_lw,
#                                             fit.chain.minmaxlw[0],
#                                             fit.chain.minmaxlw[1], 
#                                             data_source, outxml)
#     nhdr = len(hdr)
#     val = ",".join(val)
#     hdr = ",".join(hdr)
#     val += "\n"
#     hdr += "\n"
#      
#     hdr_flag = True
#     if os.path.isfile(csvfile):
#         with open(csvfile, 'r+') as f:
#             data = f.readlines()
#             if len(data)>1:
#                 last = data[-1]
#                 nlast = len(last.split(','))
#                 if nlast == nhdr:
#                     hdr_flag = False
#                 
#     with open(csvfile, 'a') as f:
#         if hdr_flag:
#             f.write(hdr)
#         f.write(val)

    # Save results to XML -----------------------------------------------------
    
    if verbose: print("""Saving dataset to XML file "%s". """ % outxml)
    
    try:
        util_export.export(outxml, [data_metab,], None, None, False)
    except Exception as e:
        msg = """I can't write the file "%s".""" % outxml
        print(msg, file=sys.stderr)
        print(repr(e), file=sys.stderr)
        sys.exit(-1)
        
    return None, None
            
            
            
               
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


def analysis_kernel(param):
    
    datafname, fbase, out_base, fpreset_coil, fpreset_ecc, fpreset_water, fpreset_metab, fbasis_mmol, out_label = param
    
    debug   = False
    verbose = True

    preset_metab = load_preset(fpreset_metab, verbose=True, debug=debug) 
    preset_coil  = load_preset(fpreset_coil,  verbose=True, debug=debug)
    preset_water = load_preset(fpreset_water, verbose=True, debug=debug)
    preset_ecc   = load_preset(fpreset_ecc,   verbose=True, debug=debug)

    # Use subdir names to create unique prefix for output files
    parts = os.path.normpath(datafname).split(os.sep)
    out_prefix = out_label+parts[-3]+'_'+parts[-2]    # Ex. 'S4_V1'

#    datasets = process_twix2dataset( datafname, verbose=True, debug=debug )

    dataformat = 'siemens_twix_slaser_cmrr_ve'          
    datasets   = util_file_import.get_datasets_cli(datafname, dataformat, None) 
    
    dataset_mmol, msg = util_file_import.open_viff_dataset_file([fbasis_mmol,]) 
    for item in dataset_mmol:
        datasets.append(item)


    if verbose: 
        print("Unique Output Prefix = "+out_prefix)
        print("Unique Output Base   = "+out_base)

    if not debug:
        img0, outxml0 = analysis_cli_gulin( datasets, 
                                            preset_metab,
                                            preset_coil,
                                            preset_water,
                                            preset_ecc,
                                            out_base,
                                            out_prefix,
                                            basis_mmol=dataset_mmol,
                                            verbose=True)

        bob = 10
        bob += 1

    if verbose: 
        print("Finished - " + str(os.path.basename(datafile)))
    

    
    

if __name__ == '__main__':

    debug          = False
    verbose        = True
    single_process = False
    nprocess       = 12

    fbase    = 'D:\\Users\\bsoher\\projects\\2015_gulin_BRP\\data_sharing\\BRP_twix_v2\\'
    out_base = 'D:\\Users\\bsoher\\projects\\2015_gulin_BRP\\data_sharing\\BRP_twix_v2\\results_basis_mixed_v2\\'
    out_label = 'mixed_'

    fpreset_coil  = fbase+'preset_analysis_brp_slaser_coil_v1.xml'
    fpreset_ecc   = fbase+'preset_analysis_brp_slaser_ecc_v1.xml'
    fpreset_water = fbase+'preset_analysis_brp_slaser_water_v1.xml'
    fpreset_metab = fbase+'preset_analysis_brp_slaser_metab_mixed_v5_fixedTa_withMMol.xml'
    #fpreset_metab = fbase+'preset_analysis_brp_slaser_metab_mixed_v1.xml'
    
    fbasis_mmol   = 'D:\\Users\\bsoher\\projects\\2015_gulin_BRP\\data_sharing\\BRP_twix_v2\\basis_mmol_dataset_seadMM2014_truncat2048pts_normScale100dc015.xml'

#    out_base      = 'D:\\Users\\bsoher\\projects\\2015_gulin_BRP\\data_sharing\\BRP_twix_v2\\results_basis_indiv\\'
#    out_label     = 'indiv_'
#    fpreset_metab = fbase+'preset_analysis_brp_slaser_metab_mixed_v2_increasedInitval.xml'
#    fpreset_metab = fbase+'preset_slaser_metab_basis_mixed.xml' 
#    fpreset_metab = fbase+'preset_slaser_metab_basis_indiv.xml'
#    fpreset_coil  = fbase+'preset_slaser_combine.xml'
#    fpreset_water = fbase+'preset_slaser_water.xml'
#    fpreset_ecc   = fbase+'preset_slaser_ecc.xml'

    
    fdata = ['S1\\V1\\meas_MID00347_FID24051_CBWM_WS64.dat',
            'S1\\V2\\meas_MID00040_FID25206_CBWM_WS64_manTRA.dat',
            'S2\\V1\\meas_MID00148_FID24553_CBWM_WS64_man_wROT.dat',
            'S2\\V2\\meas_MID00520_FID26915_CBWM_WS64.dat',
            'S3\\V1\\meas_MID00329_FID25911_CBWM_WS64.dat',
            'S3\\V2\\meas_MID00613_FID28758_CBWM_WS64.dat',
            'S4\\V1\\meas_MID00257_FID29236_CBWM_WS64.dat',
            'S4\\V2\\meas_MID00146_FID30165_CBWM_WS64.dat',
            'S5\\V1\\meas_MID00062_FID30294_CBWM_WS64.dat',
            'S5\\V2\\meas_MID00259_FID30493_CBWM_WS64.dat',
            'S6\\V1\\meas_MID00393_FID39482_CBWM_WS64.dat',
            'S6\\V2\\meas_MID00148_FID39689_CBWM_WS64.dat',
            'S7\\V1\\meas_MID00090_FID22649_CBWM_WS64_Tra.dat',
            'S7\\V2\\meas_MID00257_FID22816_CBWM_WS64_traMan.dat',
            'S7\\V2_B\\meas_MID00250_FID22809_CBWM_WS64.dat',
            'S8\\V1\\meas_MID00106_FID27202_CBWM_WS64_TRA.dat',
            'S8\\V1_B\\meas_MID00105_FID27201_CBWM_WS64.dat',
            'S8\\V2\\meas_MID00536_FID28681_CBWM_WS64.dat',
            'S9\\V1\\meas_MID00218_FID38117_CBWM_WS64.dat',
            'S9\\V2\\meas_MID00541_FID41901_CBWM_WS64.dat',
             ]

    _datafiles = [fbase+item for item in fdata]
    
    datafiles = [_datafiles[8],]
#    datafiles = _datafiles 
#    datafiles = _datafiles[0:5]
    
    #----------------------------------------------------------
    # Basic file checking for existence

    msg = ''
    for datafile in datafiles:
        if not os.path.isfile(datafile):
            msg = """Main DATAFILE does not exist "%s".""" % datafile 

    for item in [fpreset_metab,fpreset_coil,fpreset_water,fpreset_ecc,fbasis_mmol]:
        if not os.path.isfile(item):
            msg += """\nPRESET FILE does not exist "%s".""" % item
    if msg:        
        print(msg, file=sys.stderr)
        print(msg, file=sys.stdout)
        sys.exit(-1)


    #----------------------------------------------------------
    # Run the processing

    if len(datafiles) == 1 or single_process:

        for datafile in datafiles:
            params = [datafile, fbase, out_base, fpreset_coil, fpreset_ecc, fpreset_water, fpreset_metab, fbasis_mmol, out_label]
            analysis_kernel(params)
    else:
        params = []
        for datafile in datafiles:
            params.append([datafile, fbase, out_base, fpreset_coil, fpreset_ecc, fpreset_water, fpreset_metab, fbasis_mmol, out_label])
            
        pool = multiprocessing.Pool(processes=nprocess)
        results = pool.map(analysis_kernel, params)



        
#if __name__ == '__main__':
#    main()        
        