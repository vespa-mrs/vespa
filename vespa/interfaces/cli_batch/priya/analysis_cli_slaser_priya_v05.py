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




#  More specifically, this is for reading the Siemens Twix data for just
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
                 out_base, out_prefix, out_set=None,
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
    pset_coil, pset_ecc, pset_water, pset_metab = presets

    msg = ""
    
    # Load and Process - Coil Combine Dataset -----------------------

    if data_coil is not None:
        if verbose: print(out_prefix+" - Apply Preset and Run Chain - Coil Combine")    
        try:
            msg = out_prefix+" - " + """applying preset - coil combine""" 
            data_coil.apply_preset(pset_coil, voxel=(0,0,0))  # update dataset object with preset blocks and chains

            msg = out_prefix+" - " + """running chain - coil combine""" 
            _process_all_blocks(data_coil)

        except:
            if not in_gui:
                print(msg+'\n'+str(sys.exc_info()), file=sys.stderr)
                sys.exit(-1)
            else:
                raise CliError(msg)

    # Load Preset - Ecc, Water and Metab Datasets -------------------

    try:
        # Apply presets to ecc, water and metab datasets
      
        msg = "Apply Preset - Ecc, Water and Metab Datasets"
        if verbose: print(msg)    

        if data_ecc is not None and pset_ecc is not None:
            msg = out_prefix+" - " + """applying preset - ecc""" 
            data_ecc.apply_preset(pset_ecc, voxel=(0,0,0))      # chain  

        if data_water is not None and pset_water is not None:
            msg = out_prefix+" - " + """applying preset - water""" 
            data_water.apply_preset(pset_water, voxel=(0,0,0))  

        if data_metab is not None and pset_metab is not None:
            msg = out_prefix+" - " + """applying preset - metab""" 
            data_metab.apply_preset(pset_metab, voxel=(0,0,0)) 

        #----------------------------------------------------------------------
        # Attach coil combine to ecc, water and metab datasets - run chain ecc

        if data_coil is not None:
            msg = out_prefix+" - " + """attaching coil combine to - ecc, water and metab"""
            for dset in [data_ecc, data_water, data_metab]:
                if dset is not None:
                    dset.set_associated_dataset_combine(data_coil)
        
        if data_ecc is not None:
            msg = out_prefix+" - " + """running chain - ecc"""
            if verbose: print(msg)
            _process_all_blocks(data_ecc)       # get combined FID for next steps
        
        #----------------------------------------------------------------------
        # Attach ecc to water and metab datasets - run chain water

        if data_ecc is not None:
            msg = out_prefix+" - " + """attaching ecc to - water and metab"""
            for dset in [data_water, data_metab]:
                if dset is not None:
                    dset.set_associated_dataset_ecc(data_ecc)
        
        if data_water is not None:
            msg = out_prefix+" - " + """running chain - water"""
            if verbose: print(msg)
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

        msg = out_prefix+" - " + """running chain - metab"""
        if verbose: print(msg)
        
        _process_all_blocks(data_metab)

    except:
        if not in_gui:
            print('Error: '+msg+'\n'+str(sys.exc_info()), file=sys.stderr)
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

    except Exception as e:
        msg = """load_preset - Exception reading Preset file "%s".\n%s""" % (presetfile, repr(e))
        print(msg, file=sys.stderr)
        sys.exit(-1)
    
    return preset



def analysis_kernel(param):
    
    try:
        
        fdata, fpresets, fbasis_mmol, fbase, out_base, out_label, out_set, dformat = param
        
        debug   = False
        verbose = True

        fpset_coil, fpset_ecc, fpset_water, fpset_metab = fpresets
    
        # Use file names to create unique prefix for output files
        parts = os.path.normpath(fdata[0]).split(os.sep)
        if '_1_slaser' in fdata[0].lower():
            out_prefix = out_label+parts[-3]+'_deepWM'    # Ex. twix_2020_01_08_camrd_452_006_year4_deepWM
        else:
            out_prefix = out_label+parts[-3]+'_normGM'

        if verbose:
            print('Begin - '+out_prefix)
    
        pset_coil  = load_preset(fpset_coil,  verbose=True, debug=debug)
        pset_ecc   = load_preset(fpset_ecc,   verbose=True, debug=debug)
        pset_water = load_preset(fpset_water, verbose=True, debug=debug)
        pset_metab = load_preset(fpset_metab, verbose=True, debug=debug) 
        presets = [pset_coil, pset_ecc, pset_water, pset_metab]

        r = util_file_import.get_datasets_cli(fdata[0], dformat, None)

        # only add in the datasets we want for next step
        data_coil, data_ecc, data_water, data_ecc2, data_water2, data_metab = r
        datasets = [data_coil, data_ecc, data_water, data_metab]

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

    dformats = ['siemens_twix_svs_slaser_cmrr_vb',
                'siemens_twix_svs_slaser_cmrr_vb',
                'siemens_twix_svs_slaser_cmrr_vb',
                'siemens_twix_svs_slaser_cmrr_vb',
                'siemens_twix_svs_slaser_cmrr_vb',
                'siemens_twix_svs_slaser_cmrr_vb',
                'siemens_twix_svs_slaser_cmrr_vb',
                'siemens_twix_svs_slaser_cmrr_vb',
                'siemens_twix_svs_slaser_cmrr_vb',
                'siemens_twix_svs_slaser_cmrr_vb',
                'siemens_twix_slaser_cmrr_ve',
                'siemens_twix_slaser_cmrr_ve',
                'siemens_twix_slaser_cmrr_ve',
                'siemens_twix_slaser_cmrr_ve',
                'siemens_twix_slaser_cmrr_ve',
                'siemens_twix_slaser_cmrr_ve',
                'siemens_twix_slaser_cmrr_ve',
                'siemens_twix_slaser_cmrr_ve',
                'siemens_twix_svs_slaser_cmrr_vb',
                'siemens_twix_svs_slaser_cmrr_vb',
                'siemens_twix_slaser_cmrr_ve',
                'siemens_twix_slaser_cmrr_ve',
                'siemens_twix_slaser_cmrr_ve',
                'siemens_twix_slaser_cmrr_ve',
                'siemens_twix_slaser_cmrr_ve',
                'siemens_twix_slaser_cmrr_ve',
                'siemens_twix_slaser_cmrr_ve',
               ]

    fbase = 'D:\\Users\\bsoher\\projects\\2019_Priya_CNS_MRS\\data\\'

    out_base  = fbase + '_results_v04\\'
    out_label = 'presets_v2_'

    fpset_coil  = fbase + 'preset_analysis_coil_v1.xml'
    fpset_ecc   = fbase + 'preset_analysis_ecc_v1.xml'
    fpset_water = fbase + 'preset_analysis_water_v1.xml'
    fpset_metab = fbase + 'preset_analysis_metab_v2.xml'
    fbasis_mmol = fbase + 'basis_mmol_from_datasim_sead2014_siemens.xml'

    fpresets = [fpset_coil, fpset_ecc, fpset_water, fpset_metab]

    fdata = [

            [fbase + "camrd_452_006_year4_2020_01_08\\twix\\meas_MID00206_FID102846_1_slaser028_avg32_wref1_resolv.dat", ],
            [fbase + "camrd_452_006_year4_2020_01_08\\twix\\meas_MID00216_FID102856_2_slaser028_avg32_wref1_resolv.dat", ],
            [fbase + "camrd_452_006_year6_2021_03_17\\twix\\meas_MID01411_FID140476_1_slaser028_avg32_wref1_resolv.dat",],
            [fbase + "camrd_452_006_year6_2021_03_17\\twix\\meas_MID01418_FID140483_2_slaser028_avg32_wref1_resolv.dat",],
            [fbase + "camrd_452_007_year4_2020_11_12\\twix\\meas_MID00355_FID130102_1_slaser028_avg32_wref1_resolv.dat",],
            [fbase + "camrd_452_007_year4_2020_11_12\\twix\\meas_MID00364_FID130111_2_slaser028_avg32_wref1_resolv.dat",],
            [fbase + "camrd_452_008_year4_2020_02_17\\twix\\meas_MID00114_FID105950_1_slaser028_avg32_wref1_resolv.dat",],
            [fbase + "camrd_452_008_year4_2020_02_17\\twix\\meas_MID00127_FID105963_2_slaser028_avg32_wref1_resolv.dat",],
            [fbase + "camrd_452_008_year5_2021_02_11\\twix\\meas_MID00032_FID137200_1_slaser028_avg32_wref1_resolv.dat",],
            [fbase + "camrd_452_008_year5_2021_02_11\\twix\\meas_MID00039_FID137207_2_slaser028_avg32_wref1_resolv.dat",],
            [fbase + "camrd_452_008_year7_2023_05_11\\twix\\meas_MID00294_FID05387_1_sLASER028_avg32_wref1_HYPER.dat", ],
            [fbase + "camrd_452_008_year7_2023_05_11\\twix\\meas_MID00315_FID05408_2_sLASER028_avg32_wref1_GM.dat", ],
            [fbase + "camrd_452_009_year3_2022_12_05\\twix\\meas_MID00071_FID10779_1_sLASER028_avg32_wref1_HYPER.dat",],
            [fbase + "camrd_452_009_year3_2022_12_05\\twix\\meas_MID00078_FID10786_2_sLASER028_avg32_wref1_GM.dat",],
            [fbase + "camrd_452_010_year4_2024_02_28\\twix\\meas_MID00367_FID09162_1_slaser_metab_avg32_hyper.dat", ],
            [fbase + "camrd_452_010_year4_2024_02_28\\twix\\meas_MID00374_FID09169_2_slaser_metab_avg32_gm.dat", ],
            [fbase + "camrd_452_020_year3_2023_08_09\\twix\\meas_MID00258_FID14445_1_sLASER028_avg32_wref1_HYPER.dat", ],
            [fbase + "camrd_452_020_year3_2023_08_09\\twix\\meas_MID00267_FID14454_2_sLASER028_avg32_wref1_GM.dat", ],
            [fbase + "camrd_452_025_year2_2020_03_02\\twix\\meas_MID00117_FID107165_1_slaser028_avg32_wref1_resolv.dat",],
            [fbase + "camrd_452_025_year2_2020_03_02\\twix\\meas_MID00126_FID107174_2_slaser028_avg32_wref1_resolv.dat",],
            [fbase + "camrd_452_025_year4_2023_08_17\\twix\\meas_MID00496_FID15352_1_dkd_slaser_avg32_wref4x2.dat", ],
            [fbase + "camrd_452_025_year4_2023_08_17\\twix\\meas_MID00517_FID15373_2_dkd_slaser_avg32_wref4x2.dat", ],
            [fbase + "camrd_452_029_year3_2022_10_10\\twix\\meas_MID00107_FID06037_1_sLASER028_avg32_wref1_HYPER.dat",],
            [fbase + "camrd_452_029_year3_2022_10_10\\twix\\meas_MID00116_FID06046_2_sLASER028_avg32_wref1_GM.dat",],
            [fbase + "camrd_452_046_year1_2023_07_27\\twix\\meas_MID00222_FID13211_1_sLASER028_avg32_wref1_HYPER.dat",],
            [fbase + "camrd_452_049_year1_2023_03_23\\twix\\meas_MID00039_FID00660_1_sLASER028_avg32_wref1_HYPER.dat",],
            [fbase + "camrd_452_049_year1_2023_03_23\\twix\\meas_MID00047_FID00668_2_sLASER028_avg32_wref1_GM.dat",],

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

        for datafile,fmt in zip(datafiles,dformats):
            params = [datafile, fpresets, fbasis_mmol, fbase, out_base, out_label, out_set, fmt]
            analysis_kernel(params)
    else:
        params = []
        for datafile, fmt in zip(datafiles, dformats):
            params.append([datafile, fpresets, fbasis_mmol, fbase, out_base, out_label, out_set, fmt])
            
        pool = multiprocessing.Pool(processes=nprocess)
        results = pool.map(analysis_kernel, params)
    
    bob = 10
    bob += 1

    print("\nEnd Time - "+get_time())
        


if __name__ == '__main__':
    
    do_main()