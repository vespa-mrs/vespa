# Python modules

import os
import sys
import multiprocessing
import datetime


# 3rd party modules
from matplotlib.backends.backend_pdf import PdfPages


# Our modules
import vespa.analysis.figure_layouts as figure_layouts
import vespa.analysis.util_import as util_import
import vespa.analysis.util_file_import as util_file_import
import vespa.common.util.export as util_export
import vespa.common.util.time_ as util_time




#  More specifically, this is for reading the Siemens Twis data for just
#  the metab data, but which also has two initial FIDs that are water
#  unsuppressed that I am using for ecc and water.



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



def analysis_cli(datasets, presets, 
                 out_base, out_prefix, datdir, out_set=None,
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
    
    data_water, data_metab, basis_mmol = datasets
    pset_water, pset_metab = presets

    msg = ""
    
    # Load and Process - Coil Combine and Water Dataset -----------------------

    if data_water is not None:
        if verbose: print(out_prefix+" - Apply Preset and Run Chain - Coil Combine")
        try:
            msg = out_prefix+" - " + """applying preset - coil combine"""
            data_water.apply_preset(pset_water, voxel=(0,0,0))  # update dataset object with preset blocks and chains

            msg = out_prefix+" - " + """running chain - coil combine"""
            _process_all_blocks(data_water)

        except:
            if not in_gui:
                print(msg+'\n'+str(sys.exc_info()[1]), file=sys.stderr)
                sys.exit(-1)
            else:
                raise CliError(msg)

    # Load Preset - Ecc, Water and Metab Datasets -------------------

    try:
        # Apply preset to metab dataset
      
        msg = "Apply Preset - Water and Metab Datasets"
        if verbose: print(msg)    

        if data_metab is not None and pset_metab is not None:
            msg = out_prefix+" - " + """applying preset - metab""" 
            data_metab.apply_preset(pset_metab, voxel=(0,0,0)) 

        #----------------------------------------------------------------------
        # Attach coil combine, mmol_basis and water to metab - run chain metab

        msg = out_prefix+" - " + """attaching combine, mmol_basis and water to - metab"""
        if data_water is not None:
            data_metab.set_associated_dataset_combine(data_water)
        if basis_mmol is not None:
            data_metab.set_associated_dataset_mmol(basis_mmol)
        if data_water is not None:
            data_metab.set_associated_dataset_quant(data_water)

        msg = out_prefix+" - " + """running chain - metab"""
        if verbose: print(msg)
        
        _process_all_blocks(data_metab)


    except:
        if not in_gui:
            print('Error: '+msg+'\n'+str(sys.exc_info()[1]), file=sys.stderr)
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
                        vespa_version='1.0.0-CLI',
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
    
    if verbose: print('Start -  CSV Save')
    
    voxel  = (0,0,0)
    outcsv = out_base+'csv_results_collated.csv'
    outcsv_fit = out_base + 'csv_results_collated_fit.csv'

    if verbose: print(out_prefix+" - " + """Saving Results to CSV "%s". """ % outcsv)

    try:
        raw    = data_metab.blocks["raw"]
        fit    = data_metab.blocks["fit"]

        # Save quantified fit results from Quant tab
        val, hdr = data_metab.quant_results_as_csv(voxel, lw = fit.chain.fitted_lw, 
                                                   lwmin     = fit.chain.minmaxlw[0], 
                                                   lwmax     = fit.chain.minmaxlw[1], 
                                                   source    = raw.get_data_source(voxel),
                                                   dsetname  = datdir,
                                                   #dsetname  = data_metab.dataset_filename,
                                                   decor1    = False)
        val = ",".join(val) + "\n"
        hdr = ",".join(hdr) + "\n"
            
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

        # Save fitted area results from Voigt tab
        val1, hdr1 = data_metab.fit_results_as_csv(voxel, lw=fit.chain.fitted_lw,
                                                   lwmin=fit.chain.minmaxlw[0],
                                                   lwmax=fit.chain.minmaxlw[1],
                                                   source=raw.get_data_source(voxel),
                                                   dsetname=datdir,
                                                   #dsetname=data_metab.dataset_filename,
                                                   decor1=False)
        val1 = ",".join(val1) + "\n"
        hdr1 = ",".join(hdr1) + "\n"

        hdr_flag1 = True
        if os.path.isfile(outcsv_fit):
            with open(outcsv_fit, 'r+') as f:
                data = f.readlines()
                if len(data) > 1:
                    if data[0] == hdr1: hdr_flag1 = False

        with open(outcsv_fit, 'a') as f:
            if hdr_flag1:
                f.write(hdr1)
            f.write(val1)

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
        elif key == 'prep':
            block = dataset.blocks[key]
            tmp = block.chain.run(voxel, entry='all', freq_raw=True)
            chain_outputs[key] = tmp
        else:
            block = dataset.blocks[key]
            tmp = block.chain.run(voxel, entry='all')
            chain_outputs[key] = tmp

    return chain_outputs


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
        
        fdata, fpresets, fbasis_mmol, fbase, out_base, out_label, out_set, dformat = param
        
        debug   = False
        verbose = True

        fpset_water, fpset_metab = fpresets
    
        # Use file names to create unique prefix for output files
        parts = os.path.normpath(fdata[0]).split(os.sep)
        datdir = parts[-3]
        if 'hyper' in fdata[0].lower():
            out_prefix = out_label+parts[-3]+'_deepWM_press'    # Ex. twix_2020_01_08_camrd_452_006_year4_deepWM
        else:
            out_prefix = out_label+parts[-3]+'_normGM_press'

        if verbose:
            print('Begin - '+out_prefix)
    
        pset_water = load_preset(fpset_water, verbose=True, debug=debug)
        pset_metab = load_preset(fpset_metab, verbose=True, debug=debug) 
        presets = [pset_water, pset_metab]

        data_water = util_file_import.get_datasets_cli(fdata[0], dformat[0], None)
        data_metab = util_file_import.get_datasets_cli(fdata[1], dformat[0], None)

        # only add in the datasets we want for next step
        datasets = [data_water[0], data_metab[0]]

        dset_mmol, msg = util_file_import.open_viff_dataset_file([fbasis_mmol,]) 
        for item in dset_mmol:
            datasets.append(item)
    
        if verbose: 
            print("Unique Output Prefix = "+out_prefix)
            print("Unique Output Base   = "+out_base)
    
        if not debug:
            img0, outxml0 = analysis_cli( datasets, 
                                          presets,
                                          out_base,
                                          out_prefix,
                                          datdir,
                                          out_set=out_set,
                                          verbose=True)
        else:
            img0 = None

        if verbose: 
            print('Finished - '+out_prefix + ", fdata[0] - " + str(os.path.basename(fdata[0])))
        
    except Exception as e:
        if verbose:
            print('Exception - '+out_prefix)
        msg = "I am in - " + out_prefix
        raise CliError(msg)
            
    return (img0, out_prefix)
    

def get_time():
    now = datetime.datetime.now()
    return now.strftime("%H:%M:%S")   
    

def do_main():

    print("Start Time - "+get_time()+"\n")

    debug          = False
    verbose        = True
    single_process = True
    nprocess       = 8

    out_set = { 'savetype' : 'lcm_multi',
                'minplot'  : 0.5,
                'maxplot'  : 4.2,
                'fixphase' : False,
                'fontname' : 'Courier New',
                'dpi'      : 300,
                'pad_inches' : 0.5
             }

    dformat = ['siemens_twix_svs_se',]

    fbase = 'D:\\Users\\bsoher\\projects\\2019_Priya_CNS_MRS\\data_cohort1\\'
    fpset = fbase + 'presets_press\\'

    out_base  = fbase + '_results_v05\\'
    out_label = 'presets_v2_'

    fpset_water = fpset + 'preset_analysis_press_water_v2.xml'
    fpset_metab = fpset + 'preset_analysis_press_metab_v2.xml'
    fbasis_mmol = fpset + 'basis_mmol_from_datasim_sead2014_siemens_forPRESS.xml'

    fpresets = [fpset_water, fpset_metab]


    fdata = [

            [fbase + "camrd_452_005_year6_2023_03_15\\twix\\meas_MID00382_FID18674_PRESS_TE35ms_water_hyper.dat",
             fbase + "camrd_452_005_year6_2023_03_15\\twix\\meas_MID00385_FID18677_PRESS_TE35ms_metab_hyper.dat"],
            [fbase + "camrd_452_005_year6_2023_03_15\\twix\\meas_MID00389_FID18681_PRESS_TE35ms_water_gm.dat",
             fbase + "camrd_452_005_year6_2023_03_15\\twix\\meas_MID00392_FID18684_PRESS_TE35ms_metab_gm.dat"],
            [fbase + "camrd_452_006_year7_2023_08_01\\twix\\meas_MID00186_FID13597_PRESS_TE35ms_water_hyper.dat",
             fbase + "camrd_452_006_year7_2023_08_01\\twix\\meas_MID00189_FID13600_PRESS_TE35ms_metab_hyper.dat"],
            [fbase + "camrd_452_006_year7_2023_08_01\\twix\\meas_MID00193_FID13604_PRESS_TE35ms_water_gm.dat",
             fbase + "camrd_452_006_year7_2023_08_01\\twix\\meas_MID00196_FID13607_PRESS_TE35ms_metab_gm.dat"],
            [fbase + "camrd_452_019_year4_2021_11_30\\twix\\meas_MID00308_FID13823_PRESS_TE35ms_water_hyper.dat",
             fbase + "camrd_452_019_year4_2021_11_30\\twix\\meas_MID00311_FID13826_PRESS_TE35ms_metab_hyper.dat"],
            [fbase + "camrd_452_019_year4_2021_11_30\\twix\\meas_MID00315_FID13830_PRESS_TE35ms_water_gm.dat",
             fbase + "camrd_452_019_year4_2021_11_30\\twix\\meas_MID00318_FID13833_PRESS_TE35ms_metab_gm.dat"],
            [fbase + "camrd_452_019_year5_2023_07_12\\twix\\meas_MID00206_FID11996_PRESS_TE35ms_water_hyper.dat",
             fbase + "camrd_452_019_year5_2023_07_12\\twix\\meas_MID00209_FID11999_PRESS_TE35ms_metab_hyper.dat"],
            [fbase + "camrd_452_019_year5_2023_07_12\\twix\\meas_MID00213_FID12003_PRESS_TE35ms_water_gm.dat",
             fbase + "camrd_452_019_year5_2023_07_12\\twix\\meas_MID00216_FID12006_PRESS_TE35ms_metab_gm.dat"],
            [fbase + "camrd_452_019_year6_2021_11_30\\twix\\meas_MID00308_FID13823_PRESS_TE35ms_water_hyper.dat",
             fbase + "camrd_452_019_year6_2021_11_30\\twix\\meas_MID00311_FID13826_PRESS_TE35ms_metab_hyper.dat"],
            [fbase + "camrd_452_019_year6_2021_11_30\\twix\\meas_MID00315_FID13830_PRESS_TE35ms_water_gm.dat",
             fbase + "camrd_452_019_year6_2021_11_30\\twix\\meas_MID00318_FID13833_PRESS_TE35ms_metab_gm.dat"],
            [fbase + "camrd_452_055_year1_2022_06_15\\twix\\meas_MID00036_FID30012_PRESS_TE35ms_water_hyper.dat",
             fbase + "camrd_452_055_year1_2022_06_15\\twix\\meas_MID00039_FID30015_PRESS_TE35ms_metab_hyper.dat"],
            [fbase + "camrd_452_055_year1_2022_06_15\\twix\\meas_MID00043_FID30019_PRESS_TE35ms_water_gm.dat",
             fbase + "camrd_452_055_year1_2022_06_15\\twix\\meas_MID00046_FID30022_PRESS_TE35ms_metab_gm.dat"],
            [fbase + "camrd_452_060_year1_2025_01_14\\twix\\meas_MID00037_FID13323_PRESS_TE35ms_water_hyper.dat",
             fbase + "camrd_452_060_year1_2025_01_14\\twix\\meas_MID00040_FID13326_PRESS_TE35ms_metab_hyper.dat"],
            [fbase + "camrd_452_060_year1_2025_01_14\\twix\\meas_MID00044_FID13330_PRESS_TE35ms_water_gm.dat",
             fbase + "camrd_452_060_year1_2025_01_14\\twix\\meas_MID00047_FID13333_PRESS_TE35ms_metab_gm.dat"],

    ]
    
    datafiles = fdata[0:]

    #----------------------------------------------------------
    # Basic file checking for existence

    msg = ''
    for datalist in datafiles:
        for item in datalist:
            if not os.path.isfile(item):
                msg += """DATAFILE does not exist "%s".""" % item

    for item in fpresets:
        if item:
            if not os.path.isfile(item):
                msg += """\nPRESET FILE does not exist "%s".""" % item

    for item in [fbasis_mmol,]:
        if item:
            if not os.path.isfile(item):
                msg += """\nMMol FILE does not exist "%s".""" % item
    if msg:        
        print(msg, file=sys.stderr)
        sys.exit(-1)


    #----------------------------------------------------------
    # Run the processing

    if len(datafiles) == 1 or single_process:

        for datafile in datafiles:
            params = [datafile, fpresets, fbasis_mmol, fbase, out_base, out_label, out_set, dformat]
            analysis_kernel(params)
    else:
        params = []
        for datafile in datafiles:
            params.append([datafile, fpresets, fbasis_mmol, fbase, out_base, out_label, out_set, dformat])
            
        pool = multiprocessing.Pool(processes=nprocess)
        results = pool.map(analysis_kernel, params)
    
    bob = 10
    bob += 1

    print("\nEnd Time - "+get_time())
        


if __name__ == '__main__':
    
    do_main()