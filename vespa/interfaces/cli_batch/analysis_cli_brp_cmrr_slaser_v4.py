# Python modules
import os
import sys
import multiprocessing
import datetime

# 3rd party modules
from matplotlib.backends.backend_pdf import PdfPages


# Our modules
import vespa.analysis.util_import as util_import
import vespa.analysis.util_file_import as util_file_import
import vespa.common.util.export as util_export
import vespa.common.util.time_ as util_time

import vespa.analysis.figure_layouts as figure_layouts

#
#  This is for Baseline sLASER data that Dinesh sent me as a follow up test
#  after finishing the BRP_twix2 data that Joers sent me initially
#
#  More specifically, this was for reading the Siemens DICOM individual FIDs
#  data for both metab, ecc, water.  HOWEVER, the results from this were
#  OFF by a FACTOR of 20+ so I'm reverting to Twix data in a new CLI module
#  analysis_cli_brp_cmrr_slaser_v5.py
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
    
    # test for keyword values
    
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

    if verbose: print("Apply Preset - Ecc, Water and Metab Datasets")    
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
                dset.set_associated_dataset_mmol(basis_mmol)
            dset.set_associated_dataset_quant(data_water)

        msg = out_prefix+" - " + """running chain - metab"""
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

    if verbose: print('Start -  PDF Save')

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
            d['Title']        = 'Vespa Provenance Output'
            d['Author']       = 'Brian J. Soher'
            d['Subject']      = 'Vespa results output'
            d['Keywords']     = 'PdfPages Vespa output lcm multi-page'
            d['CreationDate'] = datetime.datetime(today.year, today.month, today.day)
            d['ModDate']      = datetime.datetime.today()        
     
    except Exception as e:
        msg = """Failure to create/write file "%s".""" % outimg
        print(msg, file=sys.stderr)
        print(repr(e), file=sys.stderr)
        sys.exit(-1)


    # Save Water Quant and Fit results to CSV text file -----------------------
    
    if verbose: print('Start -  CSV Save')
    
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
        val = ",".join(val) + "\n"
        hdr = ",".join(hdr) + "\n"
            
        hdr_flag = True
        if os.path.isfile(outcsv):
            with open(outcsv, 'r+') as f:
                data = f.readlines()
                if len(data)>1:
                    nlast = len(data[-1].split(','))
                    if nlast == len(hdr.split(',')): hdr_flag = False

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
        
        datadirs, fbase, out_base, fpreset_coil, fpreset_ecc, fpreset_water, fpreset_metab, fbasis_mmol, out_label, out_set = param
        
        debug   = False
        verbose = True
    
        # Use subdir names to create unique prefix for output files
        parts = os.path.normpath(datadirs[0]).split(os.sep)
        out_prefix = out_label+parts[-2]    # Ex. 'C009'

        if verbose:
            print('Begin - '+out_prefix)
    
        preset_coil  = None  #load_preset(fpreset_coil,  verbose=True, debug=debug)
        preset_ecc   = load_preset(fpreset_ecc,   verbose=True, debug=debug)
        preset_water = load_preset(fpreset_water, verbose=True, debug=debug)
        preset_metab = load_preset(fpreset_metab, verbose=True, debug=debug) 
    
    
        # get DICOM filenames from datadirs
        
        tmpfiles = [os.path.join(datadirs[0],item) for item in os.listdir(datadirs[0])]
        fnames_metab = [item for item in tmpfiles if is_dicom(item)]
        
        tmpfiles = [os.path.join(datadirs[1],item) for item in os.listdir(datadirs[1])]
        fnames_ecc = [item for item in tmpfiles if is_dicom(item)]

        tmpfiles = [os.path.join(datadirs[2],item) for item in os.listdir(datadirs[2])]
        fnames_water = [item for item in tmpfiles if is_dicom(item)]
        
        # ecc, water, metab, ... mmol
        
        dataformat = 'import_siemens_dicom_fidsum'          
        dataset_metab = util_file_import.get_datasets_cli(fnames_metab, dataformat, None)
        dataset_coil  = None
        dataset_ecc   = util_file_import.get_datasets_cli(fnames_ecc,   dataformat, None)
        dataset_water = util_file_import.get_datasets_cli(fnames_water, dataformat, None) 
    
        if isinstance(dataset_metab, list): dataset_metab = dataset_metab[0]
        if isinstance(dataset_ecc, list):   dataset_ecc   = dataset_ecc[0]
        if isinstance(dataset_water, list): dataset_water = dataset_water[0]
    
        datasets = [dataset_coil, dataset_ecc, dataset_water, dataset_metab]
    
    
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
            print('Finished - '+out_prefix + ", datadir - " + str(datadirs[0]))
        
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
    
    

if __name__ == '__main__':

    print("Start Time - "+get_time()+"\n")

    debug          = False
    verbose        = True
    single_process = False
    nprocess       = 12

    out_set = { 'savetype' : 'lcm_multi',
                'minplot'  : 0.5,
                'maxplot'  : 4.2,
                'fixphase' : False,
                'fontname' : 'Courier New',
                'dpi'      : 300,
                'pad_inches' : 0.5
             }


    fbase = 'D:\\Users\\bsoher\\projects\\2015_gulin_BRP\\data_sharing\\BRP_twix_v3_long_SCA1_baseline_dinesh_2020\\'

    out_base  = fbase + 'a_results_indiv_v01\\'     # picked this so ends at top of dir hierarchy
    out_label = 'indiv_'

    fpreset_coil  = fbase+'preset_analysis_brp_slaser_coil_v2_zf4.xml'
    fpreset_ecc   = fbase+'preset_analysis_brp_slaser_ecc_v2_zf4_forBRP3.xml'
    fpreset_water = fbase+'preset_analysis_brp_slaser_water_v2_zf4_forBRP3.xml'

    fpreset_metab = fbase+'preset_analysis_brp_slaser_metab_indiv_v5_start4_lowInitFloor_splBase180pts_first95pct_forBRP3.xml'
    
    fbasis_mmol   = fbase + 'basis_mmol_dataset_seadMM2014_truncat2048pts_normScale100dc004zf4.xml'




    fdata = [
#              ['C001\\MR-SE018-vermis_test',             # from dinesh, wref1 or w1 is for ECC 
#               'C001\\MR-SE019-vermis_test_wref1',       # wref or lw (maybe wref3?) is for quantification
#               'C001\\MR-SE020-vermis_test_wref3'],      # by default the third directory is for metabolites
#              ['C002\\MR-SE015-vermis_64',
#               'C002\\MR-SE016-vermis_64_wref1',         # for this data, first slot is metab, second is ecc, third is quant
#               'C002\\MR-SE017-vermis_64_wref3'],        # and I have reorganized these sub-lists to be that way
#             ['C003\\MR-SE015-vermis_64',
#              'C003\\MR-SE016-vermis_64_wref1',
#              'C003\\MR-SE017-vermis_64_wref3'],
#             ['C004\\MR-SE014-vermis_test',
#              'C004\\MR-SE015-vermis_64_wref1',
#              'C004\\MR-SE016-vermis_64_wref3'],
            ['C005\\MR-SE014-vermis_64',
             'C005\\MR-SE015-vermis_64_wref1',
             'C005\\MR-SE016-vermis_64_wref3'],
            ['C006\\MR-SE013-vermis_64',
             'C006\\MR-SE014-vermis_64_wref1',
             'C006\\MR-SE015-vermis_64_wref3'],
            ['C009\\MR-SE011-vermis_64',
             'C009\\MR-SE012-vermis_64_wref1',
             'C009\\MR-SE013-vermis_64_wref3'],
            ['C010\\MR-SE012-vermis_64',
             'C010\\MR-SE013-vermis_64_wref1',
             'C010\\MR-SE010-vermis_lw'],
            ['C011\\MR-SE012-vermis_64',
             'C011\\MR-SE013-vermis_64_wref1',
             'C011\\MR-SE014-vermis_64_wref3'],
            ['C012\\MR-SE012-vermis_64',
             'C012\\MR-SE013-vermis_64_wref1',
             'C012\\MR-SE014-vermis_64_wref3'],
            ['C013\\MR-SE015-vermis_64',
             'C013\\MR-SE016-vermis_64_wref1',
             'C013\\MR-SE017-vermis_64_wref3'],
            ['C014\\MR-SE009-vermis_64',
             'C014\\MR-SE010-vermis_64_wref1',
             'C014\\MR-SE011-vermis_64_wref3'],
            ['C015\\MR-SE018-vermis_64',
             'C015\\MR-SE019-vermis_64_wref1',
             'C015\\MR-SE020-vermis_64_wref3'],
            ['C018\\MR-SE013-vermis_64',
             'C018\\MR-SE014-vermis_64_wref1',
             'C018\\MR-SE015-vermis_64_wref3'],
            ['C019\\MR-SE017-vermis_64',
             'C019\\MR-SE018-vermis_64_wref1',
             'C019\\MR-SE019-vermis_64_wref3'],
            ['C021\\MR-SE008-vermis_64',
             'C021\\MR-SE009-vermis_64_wref1',
             'C021\\MR-SE010-vermis_64_wref3'],
            ['C022\\MR-SE011-vermis_64',
             'C022\\MR-SE012-vermis_64_wref1',
             'C022\\MR-SE009-vermis_lw'],
            ['C023\\MR-SE009-vermis_64',
             'C023\\MR-SE010-vermis_64_wref1',
             'C023\\MR-SE011-vermis_64_wref3'],
            ['C024\\MR-SE012-vermis_64',
             'C024\\MR-SE013-vermis_64_wref1',
             'C024\\MR-SE010-vermis_LW'],
            ['C025\\MR-SE011-vermis_64',
             'C025\\MR-SE012-vermis_64_wref1',
             'C025\\MR-SE009-vermis_LW'],
            ['C026\\MR-SE011-vermis_64',
             'C026\\MR-SE012-vermis_64_wref1',
             'C026\\MR-SE009-vermis_LW'],
            ['C027\\MR-SE010-vermis_64',
             'C027\\MR-SE011-vermis_64_wref1',
             'C027\\MR-SE008-vermis_LW'],
            ['C028\\MR-SE011-vermis_64',
             'C028\\MR-SE012-vermis_64_wref1',
             'C028\\MR-SE009-vermis_LW'],
            ['C029\\MR-SE011-vermis_64',
             'C029\\MR-SE012-vermis_64_wref1',
             'C029\\MR-SE010-vermis_LW'],
            ['S101\\MR-SE023-vermis_64',
             'S101\\MR-SE024-vermis_64_wref1',
             'S101\\MR-SE025-vermis_64_wref3'],
            ['S103\\MR-SE017-vermis_64',
             'S103\\MR-SE018-vermis_64_wref1',
             'S103\\MR-SE019-vermis_64_wref3'],
            ['S104\\MR-SE014-vermis_64',
             'S104\\MR-SE015-vermis_64_wref1',
             'S104\\MR-SE016-vermis_64_wref3'],
            ['S105\\MR-SE014-vermis_64',
             'S105\\MR-SE015-vermis_64_wref1',
             'S105\\MR-SE016-vermis_64_wref3'],
            ['S106\\MR-SE014-vermis_64',
             'S106\\MR-SE015-vermis_64_wref1',
             'S106\\MR-SE016-vermis_64_wref3'],
            ['S107\\MR-SE015-vermis_64',
             'S107\\MR-SE016-vermis_64_wref1',
             'S107\\MR-SE017-vermis_64_wref3'],
            ['S108\\MR-SE012-vermis_64',
             'S108\\MR-SE013-vermis_64_wref1',
             'S108\\MR-SE014-vermis_64_wref3'],
            ['S109\\MR-SE012-vermis_vapor_64',
             'S109\\MR-SE013-vermis_64_wref1',
             'S109\\MR-SE014-vermis_64_wref3'],
            ['S110\\MR-SE012-vermis_64',
             'S110\\MR-SE013-vermis_64_wref1',
             'S110\\MR-SE009-vermis_lw'],
            ['S111\\MR-SE011-vermis_64',
             'S111\\MR-SE012-vermis_64_wref1',
             'S111\\MR-SE013-vermis_64_wref3'],
            ['S112\\MR-SE011-vermis_64',
             'S112\\MR-SE012-vermis_64_wref1',
             'S112\\MR-SE013-vermis_64_wref3'],
            ['S113\\MR-SE010-vermis_64',
             'S113\\MR-SE011-vermis_64_wref1',
             'S113\\MR-SE012-vermis_64_wref3'],
            ['S114\\MR-SE014-vermis_64',
             'S114\\MR-SE015-vermis_64_wref1',
             'S114\\MR-SE016-vermis_64_wref3'],
            ['S115\\MR-SE010-vermis_64',
             'S115\\MR-SE011-vermis_64_wref1',
             'S115\\MR-SE012-vermis_64_wref3'],
            ['S116\\MR-SE011-vermis_64',
             'S116\\MR-SE012-vermis_64_wref1',
             'S116\\MR-SE009-vermis_LW'],
            ['S117\\MR-SE011-vermis_64',
             'S117\\MR-SE012-vermis_64_wref1',
             'S117\\MR-SE009-vermis_LW'],
             ]
   
    datafiles = list(fdata) 
    for i in range(len(datafiles)):
        datafiles[i][0] = os.path.join(fbase, datafiles[i][0])
        datafiles[i][1] = os.path.join(fbase, datafiles[i][1])
        datafiles[i][2] = os.path.join(fbase, datafiles[i][2])
    
#    datafiles = [datafiles[7],] 
    
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

    #if False: #len(datafiles) == 1 or single_process:
    if len(datafiles) == 1 or single_process:

        for datafile in datafiles:
            params = [datafile, fbase, out_base, None, fpreset_ecc, fpreset_water, fpreset_metab, fbasis_mmol, out_label, out_set]
            analysis_kernel(params)
    else:
        params = []
        for datafile in datafiles:
            params.append([datafile, fbase, out_base, None, fpreset_ecc, fpreset_water, fpreset_metab, fbasis_mmol, out_label, out_set])
            
        pool = multiprocessing.Pool(processes=nprocess)
        results = pool.map(analysis_kernel, params)
    
    bob = 10
    bob += 1

    print("\nEnd Time - "+get_time())
        
