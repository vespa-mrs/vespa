# Python modules
from __future__ import division
from __future__ import print_function

import argparse
import os
import sys
import platform
import imp
import xml.etree.cElementTree as ElementTree
import collections
import multiprocessing
import datetime


# 3rd party modules
import numpy as np

from matplotlib.backends.backend_pdf import PdfPages


# Our modules
import vespa.analysis.util_import as util_import
import vespa.analysis.util_file_import as util_file_import
import vespa.analysis.figure_layouts as figure_layouts

import vespa.common.util.misc as util_misc
import vespa.common.util.export as util_export
import vespa.common.util.time_ as util_time


#
#  This is for Baseline sLASER data that Dinesh sent me as a follow up test
#  after finishing the BRP_twix2 data that Joers sent me initially
#
#  More specifically, this is for reading the Siemens Twis data for just
#  the metab data, but which also has two initial FIDs that are water
#  unsuppressed that I am using for ecc and water.
#

# Change to True to enable the assert() statements sprinkled through the code
ASSERTIONS_ENABLED = False




DESC =  \
"""
 Command line interface to process MRS data in Vespa-Analysis. 
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
        msg += 'BaseErrorMessage: '+str(e)
#        msg += 'BaseErrorMessage: '+e[1].message
        super(CliError, self).__init__(msg)


def clean_header(header):
    """ converts all values in ICE dict into a long string"""
    return "need to write"


def analysis_cli(datasets,  preset_metab,
                            preset_coil,
                            preset_water,
                            preset_ecc,
                            out_base,
                            out_prefix,
                            out_set=None,
                            basis_mmol=None,
                            verbose=False, debug=False, in_gui=False):
    
    # Test for keyword values ---------------------------------------
    
    if out_set is None:
        out_set = { 'savetype'   : 'lcm_multi',
                    'minplot'    : 0.1,
                    'maxplot'    : 4.9,
                    'fixphase'   : True,
                    'fontname'   : 'Courier New',
                    'dpi'        : 300,
                    'pad_inches' : 0.5 }
    
    
    # Sort datasets into variables ----------------------------------
    
    data_coil, data_ecc, data_water, data_metab, basis_mmol = datasets

    msg = ""
    
    # Load and Process - Coil Combine Dataset -----------------------

    if data_coil is not None:
        if verbose: print(out_prefix+" - Apply Preset and Run Chain - Coil Combine")    
        try:
            msg = out_prefix+" - " + """applying preset - coil combine""" 
            data_coil.apply_preset(preset_coil, voxel=(0,0,0))  # update dataset object with preset blocks and chains
    
            msg = out_prefix+" - " + """running chain - coil combine""" 
            _process_all_blocks(data_coil)
        
        except:
            if not in_gui:
                print(msg+'\n'+str(sys.exc_info()[1]), file=sys.stderr)
                sys.exit(-1)
            else:
                raise CliError(msg)

    # Load Preset - Ecc, Water and Metab Datasets -------------------

    if verbose: print(out_prefix+"Apply Preset - Ecc, Water and Metab Datasets")    
    try:
        # Apply presets to ecc, water and metab datasets
        
        if data_ecc is not None and preset_ecc is not None:
            msg = out_prefix+" - " + """applying preset - ecc""" 
            data_ecc.apply_preset(preset_ecc, voxel=(0,0,0))      # chain  

        if data_water is not None and preset_water is not None:
            msg = out_prefix+" - " + """applying preset - water""" 
            data_water.apply_preset(preset_water, voxel=(0,0,0))  

        if data_metab is not None and preset_metab is not None:
            msg = out_prefix+" - " + """applying preset - metab""" 
            data_metab.apply_preset(preset_metab, voxel=(0,0,0))  

        #----------------------------------------------------------------------
        # Attach coil combine to ecc, water and metab datasets - run chain ecc

        if data_coil is not None:
            msg = out_prefix+" - " + """attaching coil combine to - ecc, water and metab"""
            for dset in [data_ecc, data_water, data_metab]:
                if dset is not None:
                    dset.set_associated_dataset_combine(data_coil)

        if verbose: print(out_prefix+" - " + """running chain - ecc""")
        
        if data_ecc is not None:
            msg = out_prefix+" - " + """running chain - ecc"""
            _process_all_blocks(data_ecc)       # get combined FID for next steps
        
        #----------------------------------------------------------------------
        # Attach ecc to water and metab datasets - run chain water

        if data_ecc is not None:
            msg = out_prefix+" - " + """attaching ecc to - water and metab"""
            for dset in [data_water, data_metab]:
                if dset is not None:
                        dset.set_associated_dataset_ecc(data_ecc)
        
        if verbose: print(out_prefix+" - " + """running chain - water""")
        
        if data_water is not None:
            msg = out_prefix+" - " + """running chain - water"""
            _process_all_blocks(data_water)

        #----------------------------------------------------------------------
        # Attach mmol_basis and water to metab dataset - run chain metab

        msg = out_prefix+" - " + """attaching mmol_basis and water to - metab"""
        for dset in [data_metab,]:
            if basis_mmol is not None:
                if dset is not None:
                    dset.set_associated_dataset_mmol(basis_mmol)
            dset.set_associated_dataset_quant(data_water)

        if verbose: print(out_prefix+" - " + """running chain - metab""")
        
        _process_all_blocks(data_metab)


    except:
        if not in_gui:
            print('Error: '+msg+'\n'+sys.exc_info()[1].message, file=sys.stderr)
            sys.exit(-1)
        else:
            raise CliError(msg)
    
    #--------------------------------------------------------------------------
    # Begin Output
    
    timestamp = util_time.now(util_time.DISPLAY_TIMESTAMP_FORMAT)
    
    
    # Create unique name ID for this dataset ------------------------
    
    outxml = out_base+'provenance_'+out_prefix+'.xml'
    
    data_metab.dataset_filename = outxml

    # Save provenance to XML -----------------------------------------------------
    
    if verbose: print(out_prefix+" - " + """Saving dataset to XML file "%s". """ % outxml)
    
    try:
        util_export.export(outxml, [data_metab,], None, None, False)
    except Exception as e:
        msg = """I can't write the file "%s".""" % outxml
        print(msg, file=sys.stderr)
        print(repr(e), file=sys.stderr)
        sys.exit(-1)


    # Save fitting results to PDF -----------------------------------------------------

    fig_call = figure_layouts.null_call         # default

    if out_set['savetype'] == 'lcm':
        outimg   = out_base+'plot_lcm_'+out_prefix+'.pdf'
        fig_call = figure_layouts.lcm_like
    elif out_set['savetype'] == 'lcm_multi':
        outimg   = out_base+'plots_lcm_multi_'+out_prefix+'.pdf'
        fig_call = figure_layouts.lcm_multipage_pdf

    if verbose: print(out_prefix+" - " + """Saving Results to PDF "%s". """ % outimg)
    
    try:
        figs = fig_call(data_metab, 
                        viffpath='Analysis - CLI Batch', 
                        vespa_version='0.10.0-CLI',
                        timestamp='',
                        fontname=out_set['fontname'],
                        minplot=out_set['minplot'],
                        maxplot=out_set['maxplot'],
                        nobase=False,
                        extfig=None,
                        fixphase=out_set['fixphase'],
                        verbose=False, 
                        debug=False, 
                        quantvals=True)
        
        # Create the PdfPages object to which we will save the pages:
        # The with statement endsures object closed at end of block, even if Exception
        with PdfPages(outimg) as pdf:
            for fig in figs:
                pdf.savefig(fig, 
                            dpi=out_set['dpi'], 
                            pad_inches=out_set['pad_inches'], 
                            facecolor=fig.get_facecolor(), 
                            edgecolor='none')
        
            # We can also set the file's metadata via the PdfPages object:
            today = datetime.date.today()
            d = pdf.infodict()
            d['Title']        = u'Vespa Provenance Output'
            d['Author']       = u'Brian J. Soher'
            d['Subject']      = u'Vespa results output'
            d['Keywords']     = u'PdfPages Vespa output lcm multi-page'
            d['CreationDate'] = datetime.datetime(today.year, today.month, today.day)
            d['ModDate']      = datetime.datetime.today()        
     
    except Exception as e:
        msg = """Failure to create/write file "%s".""" % outimg
        print(msg, file=sys.stderr)
        print(repr(e), file=sys.stderr)
        sys.exit(-1)


    # Save Water Quant and Fit results to CSV text file -----------------------
    
    voxel  = (0,0,0)
    outcsv = out_base+'csv_results_collated.csv'

    if verbose: print(out_prefix+" - " + """Saving Results to CSV "%s". """ % outcsv)

    try:
        raw    = data_metab.blocks["raw"]
        fit    = data_metab.blocks["fit"]
            
        val, hdr = data_metab.quant_results_as_csv(voxel, lw = fit.chain.fitted_lw, 
                                                   lwmin     = fit.chain.minmaxlw[0], 
                                                   lwmax     = fit.chain.minmaxlw[1], 
                                                   source    = raw.get_data_source(voxel),
                                                   dsetname  = data_metab.dataset_filename,
                                                   decor1    = False)
        val = "".join(val) + "\n"
        hdr = "".join(hdr) + "\n"
            
        hdr_flag = True
        if os.path.isfile(outcsv):
            with open(outcsv, 'r+') as f:
                data = f.readlines()
                if len(data)>1:
                    nlast = len(data[-1].split(','))
                    if nlast == len(hdr): hdr_flag = False

        with open(outcsv, 'a') as f:
            if hdr_flag:
                f.write(hdr)
            f.write(val)

    except Exception as e:
        msg = """Failure to create/write file "%s".""" % outcsv
        print(msg, file=sys.stderr)
        print(repr(e), file=sys.stderr)
        sys.exit(-1)

          
    return None, None
            
            
               
def _process_all_blocks(dataset):
    """ for all voxels, run chain in all blocks to update """
    
    chain_outputs = {}
    
    voxel = dataset.all_voxels
    for key in dataset.blocks.keys():
        if key == 'spectral':
            key = 'spectral'
            block = dataset.blocks[key]
            tmp = block.chain.run(voxel, entry='all')
            chain_outputs[key] = tmp
            if 'fit' in dataset.blocks.keys():
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
        pattern = "DICM"
        binary_pattern = pattern.encode()
        return s.endswith(binary_pattern)
    else:
        return False



def load_preset(presetfile, verbose=False, debug=False):

    if not presetfile:
        return None
        
    # Load PRESET object ----------------------------------------------
    
    if verbose: print("""load_preset - Presetfile = "%s"."""  % presetfile )
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
    
    try:
        
        datafname, fbase, out_base, fpreset_coil, fpreset_ecc, fpreset_water, fpreset_metab, fbasis_mmol, out_label, out_set, dformat = param
        
        debug   = False
        verbose = True
    
        # Use subdir names to create unique prefix for output files
        parts = os.path.normpath(datafname).split(os.sep)
        out_prefix = out_label+parts[-2]    # Ex. 'C009'

        if verbose:
            print('Begin - '+out_prefix)
    
        preset_coil  = load_preset(fpreset_coil,  verbose=True, debug=debug)
        preset_ecc   = load_preset(fpreset_ecc,   verbose=True, debug=debug)
        preset_water = load_preset(fpreset_water, verbose=True, debug=debug)
        preset_metab = load_preset(fpreset_metab, verbose=True, debug=debug) 
    
        datasets   = util_file_import.get_datasets_cli(datafname, dformat, None) 

        dataset_coil, dataset_water, dataset_metab = datasets
        datasets = [dataset_coil, None, dataset_water, dataset_metab]

        dataset_mmol, msg = util_file_import.open_viff_dataset_file([fbasis_mmol,]) 
        for item in dataset_mmol:
            datasets.append(item)
    
        if verbose: 
            print("Unique Output Prefix = "+out_prefix)
            print("Unique Output Base   = "+out_base)
    
        if not debug:
            img0, outxml0 = analysis_cli( datasets, preset_metab,
                                                    preset_coil,
                                                    preset_water,
                                                    preset_ecc,
                                                    out_base,
                                                    out_prefix,
                                                    out_set=out_set,
                                                    basis_mmol=dataset_mmol,
                                                    verbose=True)
        if verbose: 
            print('Finished - '+out_prefix + ", datafname - " + datafname)
        
    except Exception as e:
        if verbose:
            print('Exception - '+out_prefix)
        msg = "I am in - " + out_prefix
        raise CliError(msg)
            
    
    return (img0, out_prefix)
    

def get_time():

    now = datetime.datetime.now()
    current_time = now.strftime("%H:%M:%S")
    return current_time   
    
    

def do_main():

    print("Start Time - "+get_time()+"\n")

    debug          = False
    verbose        = True
    single_process = False
    nprocess       = 7

    out_set = { 'savetype' : 'lcm_multi',
                'minplot'  : 0.5,
                'maxplot'  : 4.2,
                'fixphase' : False,
                'fontname' : 'Courier New',
                'dpi'      : 300,
                'pad_inches' : 0.5
             }

    dformat = 'siemens_twix_svs_slaser_cmrr_vb_gulin_long'          

    fbase = "D:\\Users\\bsoher\\projects\\2015_gulin_BRP\\data_sharing\\Test_Retest_Paper_Data\\"

    out_base  = fbase + 'a_results_siemens_twix_v01\\'     # picked this so ends at top of dir hierarchy
    # out_label = 'twix_'
    out_label = 'debug_'

    fpreset_coil  = fbase + 'preset_analysis_brp_slaser_coil_v1.xml'
    fpreset_ecc   = fbase + 'preset_analysis_brp_slaser_ecc_v1.xml'
    fpreset_water = fbase + 'preset_analysis_brp_slaser_water_v1_leftShift1.xml'
    fpreset_metab = fbase + 'preset_analysis_brp_slaser_metab_indiv_v1_leftShift1.xml'
    fbasis_mmol   = fbase + 'basis_mmol_dataset_seadMM2014_truncat2048pts_normScale100dc004zf4.xml'


    fdata=[[fbase+"MJF31\\meas_MID247_sead_PCC_metab_FID116094.dat",
            fbase+"MJF31\\meas_MID248_sead_PCC_w1_FID116095.dat",
            fbase+"MJF31\\meas_MID249_sead_PCC_w3_FID116096.dat"],
            [fbase+"MJF31\\meas_MID293_sead_CBM_metab_FID116140.dat",
            fbase+"MJF31\\meas_MID294_sead_CBM_w1_FID116141.dat",
            fbase+"MJF31\\meas_MID295_sead_CBM_w3_FID116142.dat"],
            [fbase+"MJF32\\meas_MID587_sead_PCC_metab_FID117618.dat",
            fbase+"MJF32\\meas_MID588_sead_PCC_w1_FID117619.dat",
            fbase+"MJF32\\meas_MID589_sead_PCC_w3_FID117620.dat"],
            [fbase+"MJF32\\meas_MID635_sead_CBM_metab_FID117666.dat",
            fbase+"MJF32\\meas_MID636_sead_CBM_w1_FID117667.dat",
            fbase+"MJF32\\meas_MID637_sead_CBM_w3_FID117668.dat"],
            [fbase+"MJF33\\meas_MID54_sead_CBM_metab_FID118931.dat",
            fbase+"MJF33\\meas_MID55_sead_CBM_w1_FID118932.dat",
            fbase+"MJF33\\meas_MID56_sead_CBM_w3_FID118933.dat"],
            [fbase+"MJF33\\meas_MID94_sead_PCC_metab_FID118971.dat",
            fbase+"MJF33\\meas_MID95_sead_PCC_w1_FID118972.dat",
            fbase+"MJF33\\meas_MID96_sead_PCC_w3_FID118973.dat"],
            [fbase+"MJF34\\meas_MID405_sead_PCC_metab_FID120718.dat",
            fbase+"MJF34\\meas_MID406_sead_PCC_w1_FID120719.dat",
            fbase+"MJF34\\meas_MID407_sead_PCC_w3_FID120720.dat"],
            [fbase+"MJF34\\meas_MID449_sead_CBM_metab_FID120762.dat",
            fbase+"MJF34\\meas_MID450_sead_CBM_w1_FID120763.dat",
            fbase+"MJF34\\meas_MID451_sead_CBM_w3_FID120764.dat"],
            [fbase+"MKB31\\meas_MID394_sead_PCC_metab_FID117425.dat",
            fbase+"MKB31\\meas_MID395_sead_PCC_w1_FID117426.dat",
            fbase+"MKB31\\meas_MID396_sead_PCC_w3_FID117427.dat"],
            [fbase+"MKB31\\meas_MID434_sead_CBM_metab_FID117465.dat",
            fbase+"MKB31\\meas_MID435_sead_CBM_w1_FID117466.dat",
            fbase+"MKB31\\meas_MID436_sead_CBM_w3_FID117467.dat"],
            [fbase+"MKB32\\meas_MID1029_sead_PCC_metab_FID118708.dat",
            fbase+"MKB32\\meas_MID1030_sead_PCC_w1_FID118709.dat",
            fbase+"MKB32\\meas_MID1031_sead_PCC_w3_FID118710.dat"],
            [fbase+"MKB32\\meas_MID1072_sead_CBM_metab_FID118751.dat",
            fbase+"MKB32\\meas_MID1073_sead_CBM_w1_FID118752.dat",
            fbase+"MKB32\\meas_MID1074_sead_CBM_w3_FID118753.dat"],
            [fbase+"MKB33\\meas_MID593_sead_PCC_metab_FID121771.dat",
            fbase+"MKB33\\meas_MID594_sead_PCC_w1_FID121772.dat",
            fbase+"MKB33\\meas_MID595_sead_PCC_w3_FID121773.dat"],
            [fbase+"MKB33\\meas_MID639_sead_CBM_metab_FID121817.dat",
            fbase+"MKB33\\meas_MID640_sead_CBM_w1_FID121818.dat",
            fbase+"MKB33\\meas_MID641_sead_CBM_w3_FID121819.dat"],
            [fbase+"MKB34\\meas_MID38_sead_PCC_metab_FID122582.dat",
            fbase+"MKB34\\meas_MID39_sead_PCC_w1_FID122583.dat",
            fbase+"MKB34\\meas_MID40_sead_PCC_w3_FID122584.dat"],
            [fbase+"MKB34\\meas_MID81_sead_CBM_metab_FID122625.dat",
            fbase+"MKB34\\meas_MID82_sead_CBM_w1_FID122626.dat",
            fbase+"MKB34\\meas_MID83_sead_CBM_w3_FID122627.dat"],
            [fbase+"MKH31\\meas_MID138_sead_PCC_metab_FID126366.dat",
            fbase+"MKH31\\meas_MID139_sead_PCC_w1_FID126367.dat",
            fbase+"MKH31\\meas_MID140_sead_PCC_w3_FID126368.dat"],
            [fbase+"MKH31\\meas_MID180_sead_CBM_metab_FID126408.dat",
            fbase+"MKH31\\meas_MID181_sead_CBM_w1_FID126409.dat",
            fbase+"MKH31\\meas_MID182_sead_CBM_w3_FID126410.dat"],
            [fbase+"MKH32\\meas_MID1205_sead_CBM_metab_FID127432.dat",
            fbase+"MKH32\\meas_MID1206_sead_CBM_w1_FID127433.dat",
            fbase+"MKH32\\meas_MID1207_sead_CBM_w3_FID127434.dat"],
            [fbase+"MKH32\\meas_MID1246_sead_PCC2_metab_FID127473.dat",
            fbase+"MKH32\\meas_MID1247_sead_PCC2_w1_FID127474.dat",
            fbase+"MKH32\\meas_MID1248_sead_PCC2_w3_FID127475.dat"],
            [fbase+"MKH33\\meas_MID629_sead_PCC_metab_FID129204.dat",
            fbase+"MKH33\\meas_MID630_sead_PCC_w1_FID129205.dat",
            fbase+"MKH33\\meas_MID631_sead_PCC_w3_FID129206.dat"],
            [fbase+"MKH33\\meas_MID671_sead_CBM2_metab_FID129246.dat",
            fbase+"MKH33\\meas_MID672_sead_CBM2_w1_FID129247.dat",
            fbase+"MKH33\\meas_MID674_sead_CBM2_w3_FID129249.dat"],
            [fbase+"MKH34\\meas_MID166_sead_CBM_metab_FID130419.dat",
            fbase+"MKH34\\meas_MID167_sead_CBM_w1_FID130420.dat",
            fbase+"MKH34\\meas_MID168_sead_CBM_w3_FID130421.dat"],
            [fbase+"MKH34\\meas_MID206_sead_PCC_metab_FID130459.dat",
            fbase+"MKH34\\meas_MID207_sead_PCC_w1_FID130460.dat",
            fbase+"MKH34\\meas_MID208_sead_PCC_w3_FID130461.dat"],
            [fbase+"MLV31\\meas_MID128_sead_PCC_metab_FID103144.dat",
            fbase+"MLV31\\meas_MID129_sead_PCC_w1_FID103145.dat",
            fbase+"MLV31\\meas_MID130_sead_PCC_w3_FID103146.dat"],
            [fbase+"MLV31\\meas_MID83_sead_CBM_metab_FID103099.dat",
            fbase+"MLV31\\meas_MID84_sead_CBM_w1_FID103100.dat",
            fbase+"MLV31\\meas_MID85_sead_CBM_w3_FID103101.dat"],
            [fbase+"MLV32\\meas_MID341_sead_CBM_metab_FID104138.dat",
            fbase+"MLV32\\meas_MID342_sead_CBM_w1_FID104139.dat",
            fbase+"MLV32\\meas_MID343_sead_CBM_w3_FID104140.dat"],
            [fbase+"MLV32\\meas_MID382_sead_PCC_metab_FID104179.dat",
            fbase+"MLV32\\meas_MID383_sead_PCC_w1_FID104180.dat",
            fbase+"MLV32\\meas_MID384_sead_PCC_w3_FID104181.dat"],
            [fbase+"MLV33\\meas_MID592_sead_CBM_metab_FID105322.dat",
            fbase+"MLV33\\meas_MID593_sead_CBM_w1_FID105323.dat",
            fbase+"MLV33\\meas_MID594_sead_CBM_w3_FID105324.dat"],
            [fbase+"MLV33\\meas_MID632_sead_PCC_metab_FID105362.dat",
            fbase+"MLV33\\meas_MID633_sead_PCC_w1_FID105363.dat",
            fbase+"MLV33\\meas_MID634_sead_PCC_w3_FID105364.dat"],
            [fbase+"MLV34\\meas_MID674_sead_CBM_metab_FID107289.dat",
            fbase+"MLV34\\meas_MID675_sead_CBM_w1_FID107290.dat",
            fbase+"MLV34\\meas_MID676_sead_CBM_w3_FID107291.dat"],
            [fbase+"MLV34\\meas_MID720_sead_PCC_metab_FID107335.dat",
            fbase+"MLV34\\meas_MID721_sead_PCC_w1_FID107336.dat",
            fbase+"MLV34\\meas_MID722_sead_PCC_w3_FID107337.dat"],
            [fbase+"MRM31\\meas_MID1614_sead_CBM_metab_FID90952.dat",
            fbase+"MRM31\\meas_MID1615_sead_CBM_w1_FID90953.dat",
            fbase+"MRM31\\meas_MID1616_sead_CBM_w3_FID90954.dat"],
            [fbase+"MRM31\\meas_MID1657_sead_PCC_metab_FID90995.dat",
            fbase+"MRM31\\meas_MID1658_sead_PCC_w1_FID90996.dat",
            fbase+"MRM31\\meas_MID1659_sead_PCC_w3_FID90997.dat"],
            [fbase+"MRM32\\meas_MID1471_sead_CBM_metab_FID92560.dat",
            fbase+"MRM32\\meas_MID1472_sead_CBM_w1_FID92561.dat",
            fbase+"MRM32\\meas_MID1473_sead_CBM_w3_FID92562.dat"],
            [fbase+"MRM32\\meas_MID1517_sead_PCC_metab_FID92606.dat",
            fbase+"MRM32\\meas_MID1518_sead_PCC_w1_FID92607.dat",
            fbase+"MRM32\\meas_MID1519_sead_PCC_w3_FID92608.dat"],
            [fbase+"MRM33\\meas_MID3002_sead_CBM_metab_FID94091.dat",
            fbase+"MRM33\\meas_MID3003_sead_CBM_w1_FID94092.dat",
            fbase+"MRM33\\meas_MID3004_sead_CBM_w3_FID94093.dat"],
            [fbase+"MRM33\\meas_MID3045_sead_PCC_metab_FID94134.dat",
            fbase+"MRM33\\meas_MID3046_sead_PCC_w1_FID94135.dat",
            fbase+"MRM33\\meas_MID3047_sead_PCC_w3_FID94136.dat"],
            [fbase+"MRM34\\meas_MID423_sead_PCC_metab_FID95482.dat",
            fbase+"MRM34\\meas_MID424_sead_PCC_w1_FID95483.dat",
            fbase+"MRM34\\meas_MID425_sead_PCC_w3_FID95484.dat"],
            [fbase+"MRM34\\meas_MID467_sead_CBM2_metab_FID95526.dat",
            fbase+"MRM34\\meas_MID468_sead_CBM2_w1_FID95527.dat",
            fbase+"MRM34\\meas_MID469_sead_CBM2_w3_FID95528.dat"],
            [fbase+"MSH31\\meas_MID242_sead_CBM_metab_FID104041.dat",
            fbase+"MSH31\\meas_MID243_sead_CBM_w1_FID104042.dat",
            fbase+"MSH31\\meas_MID244_sead_CBM_w3_FID104043.dat"],
            [fbase+"MSH31\\meas_MID286_sead_PCC_metab_FID104085.dat",
            fbase+"MSH31\\meas_MID287_sead_PCC_w1_FID104086.dat",
            fbase+"MSH31\\meas_MID288_sead_PCC_w3_FID104087.dat"],
            [fbase+"MSH32\\meas_MID469_sead_CBM_metab_FID105199.dat",
            fbase+"MSH32\\meas_MID470_sead_CBM_w1_FID105200.dat",
            fbase+"MSH32\\meas_MID471_sead_CBM_w3_FID105201.dat"],
            [fbase+"MSH32\\meas_MID525_sead_PCC_metab_FID105255.dat",
            fbase+"MSH32\\meas_MID526_sead_PCC_w1_FID105256.dat",
            fbase+"MSH32\\meas_MID527_sead_PCC_w3_FID105257.dat"],
            [fbase+"MSH33\\meas_MID576_sead_PCC_metab_FID107191.dat",
            fbase+"MSH33\\meas_MID577_sead_PCC_w1_FID107192.dat",
            fbase+"MSH33\\meas_MID578_sead_PCC_w3_FID107193.dat"],
            [fbase+"MSH33\\meas_MID623_sead_CBM_metab_FID107238.dat",
            fbase+"MSH33\\meas_MID624_sead_CBM_w1_FID107239.dat",
            fbase+"MSH33\\meas_MID625_sead_CBM_w3_FID107240.dat"],
            [fbase+"MSH34\\meas_MID354_sead_PCC_metab_FID108499.dat",
            fbase+"MSH34\\meas_MID355_sead_PCC_w1_FID108500.dat",
            fbase+"MSH34\\meas_MID356_sead_PCC_w3_FID108501.dat"],
            [fbase+"MSH34\\meas_MID397_sead_CBM_metab_FID108542.dat",
            fbase+"MSH34\\meas_MID398_sead_CBM_w1_FID108543.dat",
            fbase+"MSH34\\meas_MID399_sead_CBM_w3_FID108544.dat"]
            ]

#    datafiles = fdata
    datafiles = fdata[0:1]
    
    #----------------------------------------------------------
    # Basic file checking for existence

    msg = ''
    for item in [fpreset_metab,None,fpreset_water,fpreset_ecc,fbasis_mmol]:
        if item is not None:
            if not os.path.isfile(item):
                msg += """\nPRESET FILE does not exist "%s".""" % item
    if msg:        
        print(msg, file=sys.stderr)
        sys.exit(-1)


    #----------------------------------------------------------
    # Run the processing

#    if False: #len(datafiles) == 1 or single_process:
    if len(datafiles) == 1 or single_process:

        for datafile in datafiles:
            params = [datafile[0], fbase, out_base, fpreset_coil, '', fpreset_water, fpreset_metab, fbasis_mmol, out_label, out_set, dformat]
            analysis_kernel(params)
    else:
        params = []
        for datafile in datafiles:
            params.append([datafile[0], fbase, out_base, fpreset_coil, '', fpreset_water, fpreset_metab, fbasis_mmol, out_label, out_set, dformat])
            
        pool = multiprocessing.Pool(processes=nprocess)
        results = pool.map(analysis_kernel, params)
    
    bob = 10
    bob += 1

    print("\nEnd Time - "+get_time())
        


if __name__ == '__main__':
    
    do_main()