"""
Command line interface to batch process MRS data in Vespa-Analysis. 
 
User must provide a list of Data filenames, the preset file name, the output
directory name, and an optional prefix label that will be added to each output
filename (to help differentiate processing runs). 
  
Note. You have to enclose data/preset/output strings in double quotation marks 
for them to process properly. And it is best if they do not have spaces or odd
characters in them.

"""
# Python modules

import os
import sys
import datetime
import multiprocessing

# 3rd party modules
from matplotlib.backends.backend_pdf import PdfPages

# Our modules
import vespa.analysis.util_import as util_import
import vespa.analysis.util_file_import as util_file_import
import vespa.analysis.figure_layouts as figure_layouts

import vespa.common.util.time_ as util_time
import vespa.common.util.export as util_export


#
#  This is for WBNAA VE data that Oded sent in 2021
#
#  More specifically, this is for reading the Siemens Twix
#

FILES = [
        'D:\\Users\\bsoher\\tmp\\oded_test\\Subject_Name_01_oddData\\meas_MID00202_FID120571_fid_WBNAA.dat',
        'D:\\Users\\bsoher\\tmp\\oded_test\\Subject_Name_02\\meas_MID00202_FID120571_fid_WBNAA.dat',
        ]
PRESET   = 'D:\\Users\\bsoher\\tmp\\oded_test\\vespa_analysis_preset_vd_ve_bjs_v3.xml'
OUTDIR   = 'D:\\Users\\bsoher\\tmp\\oded_test\\results_v01\\'     
OUTLABEL = 'try1_'
SINGLECORE = True      # note - multicore processing may result in random order of CSV results





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
        super(CliError, self).__init__(msg)


def analysis_cli(datasets,  preset_metab,
                            preset_coil,
                            preset_water,
                            preset_ecc,
                            out_base,
                            out_prefix,
                            out_set=None,
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
            if data_water is not None:
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

        val, hdr = data_metab.fit_results_as_csv(voxel, lw = fit.chain.fitted_lw,
                                                   lwmin     = fit.chain.minmaxlw[0],
                                                   lwmax     = fit.chain.minmaxlw[1],
                                                   source    = raw.get_data_source(voxel),
                                                   dsetname  = data_metab.dataset_filename,
                                                   decor1    = False)

        # val, hdr = data_metab.quant_results_as_csv(voxel, lw = fit.chain.fitted_lw,
        #                                            lwmin     = fit.chain.minmaxlw[0],
        #                                            lwmax     = fit.chain.minmaxlw[1],
        #                                            source    = raw.get_data_source(voxel),
        #                                            dsetname  = data_metab.dataset_filename,
        #                                            decor1    = False)
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
        
        datafname, out_base, fpreset_coil, fpreset_ecc, fpreset_water, fpreset_metab, fbasis_mmol, out_label, out_set, dformat = param
        
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
    
        datasets = util_file_import.get_datasets_cli(datafname, dformat, None)

        dataset_metab = datasets[0]
        datasets = [None, None, None, dataset_metab]

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
    single_process = SINGLECORE
    nprocess       = 4

    out_set = { 'savetype' : 'lcm_multi',
                'minplot'  : 0.5,
                'maxplot'  : 4.2,
                'fixphase' : False,
                'fontname' : 'Courier New',
                'dpi'      : 300,
                'pad_inches' : 0.5
             }

    fdata = FILES
    out_base  = OUTDIR
    out_label = OUTLABEL

    fpreset_coil  = ''
    fpreset_ecc   = ''
    fpreset_water = ''
    fpreset_metab = PRESET
    fbasis_mmol   = ''

    datafiles = fdata
    #datafiles = fdata[0:1]  # run just first file

    dformat   = 'siemens_twix_wbnaa'

    #----------------------------------------------------------
    # Basic file checking for existence

    msg = ''
    for item in [fpreset_metab, fpreset_coil, fpreset_ecc, fpreset_water, fbasis_mmol]:
        if item is not '':
            if not os.path.isfile(item):
                msg += """\nPRESET FILE does not exist "%s".""" % item
    if msg:        
        print(msg, file=sys.stderr)
        sys.exit(-1)

    #----------------------------------------------------------
    # Run the processing

    if len(datafiles) == 1 or single_process:

        for datafile in datafiles:
            params = [datafile, out_base, fpreset_coil, fpreset_ecc, fpreset_water, fpreset_metab, fbasis_mmol, out_label, out_set, dformat]
            analysis_kernel(params)
    else:
        params = []
        for datafile in datafiles:
            params.append([datafile, out_base, fpreset_coil, fpreset_ecc, fpreset_water, fpreset_metab, fbasis_mmol, out_label, out_set, dformat])
            
        pool = multiprocessing.Pool(processes=nprocess)
        results = pool.map(analysis_kernel, params)
    
    print("\nEnd Time - "+get_time())
        


if __name__ == '__main__':
    
    do_main()